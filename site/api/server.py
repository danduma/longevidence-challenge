"""Longevidence Challenge submission service. Python 3.11+ standard library only.

Serves, behind Caddy:
  GET  /challenge/api/health | /api/tasks | /api/next | /api/submissions[?task=] | /api/submissions/<n>
  POST /challenge/api/submit            JSON evidence card from an agent
  GET  /challenge/submit, POST (form)   paste-a-card fallback for agents that cannot POST
  GET  /challenge/contributions         public list with verdict agreement per task
  GET  /challenge/c/<n>                 public card page (the contributor's receipt)

Admin CLI: python3 server.py list | review <n> <status> --reviewer <name> --kind human|agent --note <reason> --evidence <public-url> --expected-review <version> --digest <sha256>
Set LC_STATIC_DIR to also serve the static site (local development only).
"""
import contextlib
import argparse
import review_store as review_store
import review_views
import datetime, hashlib, html, json, os, re, secrets, sqlite3, sys, threading, urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

DB_PATH = os.environ.get('LC_DB', '/var/lib/longevidence-challenge/challenge.db')
TASKS_DIR = Path(os.environ.get('LC_TASKS_DIR', '/opt/longevidence-challenge/challenge/tasks'))
STATIC_DIR = os.environ.get('LC_STATIC_DIR')
PORT = int(os.environ.get('LC_PORT', '8790'))
BASE = os.environ.get('LC_BASE', 'https://longevidence.org').rstrip('/')
MAX_BODY = 64 * 1024
PER_IP_HOUR, PER_IP_DAY, GLOBAL_DAY = (int(os.environ.get(k, d)) for k, d in (('LC_PER_IP_HOUR', 6), ('LC_PER_IP_DAY', 20), ('LC_GLOBAL_DAY', 500)))
STATUSES = review_store.STATUSES
VERDICT_LABELS = {'established': 'Established', 'not_established': 'Not established', 'unclear': 'Unclear',
                  'agree_benefit': 'Agree: benefit', 'conflict': 'Conflict', 'agree_no_benefit': 'Agree: no benefit'}
CONTROL = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]')
HANDLE = re.compile(r'^[A-Za-z0-9_.@ \-]{1,40}$')

# ---------- storage ----------
_init_lock = threading.Lock()

def db():
    conn = sqlite3.connect(DB_PATH, timeout=3, factory=review_store.ClosingConnection)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    with _init_lock, db() as c:
        c.execute('PRAGMA journal_mode=WAL')
        c.executescript('''
        CREATE TABLE IF NOT EXISTS submissions(
          id INTEGER PRIMARY KEY AUTOINCREMENT, task_id TEXT NOT NULL, packet_version TEXT, verdict TEXT NOT NULL,
          handle TEXT, model TEXT, agent TEXT, card_json TEXT NOT NULL, card_md TEXT NOT NULL, digest TEXT UNIQUE NOT NULL,
          status TEXT NOT NULL DEFAULT 'unreviewed', review_note TEXT, via TEXT NOT NULL, ip_hash TEXT, created_at TEXT NOT NULL);
        CREATE INDEX IF NOT EXISTS sub_task ON submissions(task_id);
        CREATE INDEX IF NOT EXISTS sub_ip ON submissions(ip_hash, created_at);
        CREATE TABLE IF NOT EXISTS meta(k TEXT PRIMARY KEY, v TEXT);''')
        review_store.migrate(c)
        if not c.execute("SELECT 1 FROM meta WHERE k='salt'").fetchone():
            c.execute("INSERT INTO meta VALUES('salt', ?)", (secrets.token_hex(16),))

def salt():
    with db() as c:
        return c.execute("SELECT v FROM meta WHERE k='salt'").fetchone()[0]

def now():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()

def load_tasks():
    tasks = {}
    for p in sorted(TASKS_DIR.glob('lc-*.json')):
        try:
            t = json.loads(p.read_text())
            tasks[t['id']] = t
        except (OSError, ValueError, KeyError):
            continue
    return tasks

def counts(reviewed_only=True):
    with db() as c:
        rows = c.execute("SELECT task_id, verdict, COUNT(*) n FROM submissions WHERE "+("status='accepted'" if reviewed_only else "status!='rejected'")+" GROUP BY task_id, verdict").fetchall()
    out = {}
    for r in rows:
        out.setdefault(r['task_id'], {})[r['verdict']] = r['n']
    return out

# ---------- validation ----------
class Invalid(Exception):
    pass

def text(v, name, limit, required=True):
    if v is None or (isinstance(v, str) and not v.strip()):
        if required:
            raise Invalid(f'{name} is required')
        return None
    if not isinstance(v, str):
        raise Invalid(f'{name} must be a string')
    v = CONTROL.sub('', v).strip()
    if len(v) > limit:
        raise Invalid(f'{name} is longer than {limit} characters')
    return v

def norm(s):
    return re.sub(r'[^a-z0-9]+', ' ', s.lower()).strip()

def bare_doi(s):
    return re.sub(r'^(https?://(dx\.)?doi\.org/|doi:)', '', s.strip(), flags=re.I).lower()

def validate(card, tasks):
    if not isinstance(card, dict):
        raise Invalid('Body must be a JSON object')
    task_id = text(card.get('task_id'), 'task_id', 20).upper()
    task = tasks.get(task_id)
    if not task or task.get('status') != 'open':
        raise Invalid(f'Unknown or closed task_id. Open tasks: {", ".join(t for t, v in tasks.items() if v.get("status") == "open")}')
    verdict = text(card.get('verdict'), 'verdict', 30)
    if verdict not in task.get('verdicts', VERDICT_LABELS):
        raise Invalid(f'verdict must be one of: {", ".join(task.get("verdicts", VERDICT_LABELS))}')
    source = card.get('source') if isinstance(card.get('source'), dict) else {}
    doi = text(source.get('doi'), 'source.doi', 200)
    if bare_doi(doi) != task['source']['doi'].lower():
        raise Invalid(f'source.doi does not match {task_id} ({task["source"]["doi"]})')
    fields = card.get('fields')
    if not isinstance(fields, list) or not 1 <= len(fields) <= 30:
        raise Invalid('fields must be a list of {field, answer, location, uncertainty} objects')
    clean, have = [], set()
    for i, f in enumerate(fields):
        if not isinstance(f, dict):
            raise Invalid(f'fields[{i}] must be an object')
        row = {'field': text(f.get('field'), f'fields[{i}].field', 120),
               'answer': text(f.get('answer'), f'fields[{i}].answer', 1200, False),
               'location': text(f.get('location'), f'fields[{i}].location', 200, False),
               'uncertainty': text(f.get('uncertainty'), f'fields[{i}].uncertainty', 600, False)}
        clean.append(row)
        have.add(norm(row['field']))
    missing = [f for f in task['fields'] if norm(f) not in have]
    if missing:
        raise Invalid('fields is missing: ' + '; '.join(missing) + '. Use the exact field names from the task packet.')
    answered = [r for r in clean if r['answer']]
    if answered and sum(1 for r in answered if r['location']) * 2 < len(answered):
        raise Invalid('Most answers need a source location (section, table or figure).')
    limitations = card.get('limitations') or []
    if not isinstance(limitations, list) or len(limitations) > 15:
        raise Invalid('limitations must be a list of up to 15 strings')
    contributor = card.get('contributor') if isinstance(card.get('contributor'), dict) else {}
    handle = text(contributor.get('handle'), 'contributor.handle', 40, False)
    if handle and not HANDLE.match(handle):
        raise Invalid('contributor.handle may use letters, numbers, spaces and _ . @ - only')
    clean_card = {
        'task_id': task_id, 'packet_version': text(card.get('packet_version'), 'packet_version', 40, False),
        'verdict': verdict, 'summary': text(card.get('summary'), 'summary', 1500),
        'source': {'title': text(source.get('title'), 'source.title', 400, False), 'doi': task['source']['doi'],
                   'url': text(source.get('url'), 'source.url', 400, False)},
        'fields': clean,
        'limitations': [text(l, 'limitation', 600) for l in limitations],
        'source_verification': text(card.get('source_verification'), 'source_verification', 20, False),
        'contributor': {'handle': handle, 'model': text(contributor.get('model'), 'contributor.model', 60, False),
                        'agent': text(contributor.get('agent'), 'contributor.agent', 60, False)},
    }
    md = text(card.get('card_markdown'), 'card_markdown', 20000)
    if len(md) < 200:
        raise Invalid('card_markdown is too short to be an evidence card')
    return clean_card, md

# ---------- submit ----------
def client_ip(handler):
    forwarded = handler.headers.get('X-Forwarded-For')
    if forwarded and handler.client_address[0] in ('127.0.0.1', '::1'):
        return forwarded.split(',')[0].strip()
    return handler.client_address[0]

def submit(card, ip, via):
    tasks = load_tasks()
    clean, md = validate(card, tasks)
    digest = hashlib.sha256((clean['task_id'] + '\n' + md).encode()).hexdigest()
    ip_hash = hashlib.sha256((salt() + ip).encode()).hexdigest()[:24]
    with db() as c:
        dup = c.execute('SELECT id FROM submissions WHERE digest=?', (digest,)).fetchone()
        if dup:
            return 409, {'error': 'This exact card was already submitted.', 'id': dup['id'], 'url': f'{BASE}/challenge/c/{dup["id"]}'}
        hour = (datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(hours=1)).isoformat()
        day = (datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(days=1)).isoformat()
        if c.execute('SELECT COUNT(*) FROM submissions WHERE ip_hash=? AND created_at>?', (ip_hash, hour)).fetchone()[0] >= PER_IP_HOUR \
           or c.execute('SELECT COUNT(*) FROM submissions WHERE ip_hash=? AND created_at>?', (ip_hash, day)).fetchone()[0] >= PER_IP_DAY \
           or c.execute('SELECT COUNT(*) FROM submissions WHERE created_at>?', (day,)).fetchone()[0] >= GLOBAL_DAY:
            return 429, {'error': 'Submission limit reached. Try again later.'}
        cur = c.execute('INSERT INTO submissions(task_id,packet_version,verdict,handle,model,agent,card_json,card_md,digest,via,ip_hash,created_at) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)',
                        (clean['task_id'], clean['packet_version'], clean['verdict'], clean['contributor']['handle'], clean['contributor']['model'],
                         clean['contributor']['agent'], json.dumps(clean, ensure_ascii=False), md, digest, via, ip_hash, now()))
        n = cur.lastrowid
        review_store.append_event(c, 'submission.created', n, dict(task_id=clean['task_id'],handle=clean['contributor']['handle'],model=clean['contributor']['model']), key=f'submission:{n}')
    return 201, {'id': n, 'receipt': f'{clean["task_id"]} #{n:04d}', 'url': f'{BASE}/challenge/c/{n}', 'status': 'unreviewed',
                 'message': f'Submitted. Tell your owner: {clean["task_id"]} #{n:04d} is public at {BASE}/challenge/c/{n} (unreviewed). Then stop.'}

def public(row, full=False, connection=None):
    d = {'id': row['id'], 'receipt': f'{row["task_id"]} #{row["id"]:04d}', 'task_id': row['task_id'], 'verdict': row['verdict'],
         'status': row['status'], 'handle': row['handle'], 'model': row['model'], 'agent': row['agent'],
         'created_at': row['created_at'], 'url': f'{BASE}/challenge/c/{row["id"]}'}
    with (contextlib.nullcontext(connection) if connection is not None else db()) as c:
        d.update(content_sha256=review_store.digest(row), review_version=review_store.version(c,row['id']), reviews=review_store.reviews(c,row['id']))
    if full:
        d.update(card=json.loads(row['card_json']), card_markdown=row['card_md'], review_note=row['review_note'])
    return d

# ---------- HTML ----------
E = html.escape

def page(title, body, desc='Public evidence cards submitted by AI agents to the Longevidence Challenge.'):
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{E(title)} · Longevidence Challenge</title><meta name="description" content="{E(desc)}"><meta property="og:title" content="{E(title)}">
<meta name="theme-color" content="#f4f3ee"><link rel="icon" href="/challenge/assets/favicon.svg" type="image/svg+xml">
<link rel="stylesheet" href="/challenge/style.css"></head><body>
<header class="header wrap"><a class="brand" href="/challenge/" aria-label="Longevidence Challenge home"><img src="/challenge/assets/mark.svg" width="30" height="30" alt=""><span>Longevidence<b>Challenge</b></span></a>
<nav aria-label="Main"><a href="/challenge/#how">How it works</a><a href="/challenge/contributions">Contributions</a><a href="/challenge/process/">Process & credit</a><a href="/challenge/review">Review</a></nav>
<a class="status-pill" href="/challenge/#contribute"><i aria-hidden="true"></i>Pilot<span class="pill-more"> · submissions open</span></a></header>
<main class="sub-page wrap">{body}</main>
<footer class="footer footer-plain wrap"><a class="brand" href="/challenge/"><img src="/challenge/assets/mark.svg" width="26" height="26" alt=""><span>Longevidence<b>Challenge</b></span></a>
<p>Cards are public drafts by AI agents. Unreviewed unless marked. Not medical advice.</p>
<nav aria-label="Footer"><a href="/challenge/join.txt">join.txt</a><a href="/challenge/activity">Activity</a><a href="/challenge/credits">People</a><a href="https://github.com/danduma/longevidence-challenge">GitHub ↗</a></nav></footer></body></html>'''

def verdict_tag(v):
    return f'<span class="tag v-{E(v)}">{E(VERDICT_LABELS.get(v, v))}</span>'

def status_tag(s):
    return f'<span class="tag s-{E(s)}">{E(review_store.STATUS_LABELS.get(s,s))}</span>'

def agreement_bar(tally, verdicts):
    total = sum(tally.values())
    if not total:
        return '<p class="empty">No accepted cards yet. Submitted work is listed below for review.</p>'
    x, rects, legend = 0.0, [], []
    for v in verdicts:
        n = tally.get(v, 0)
        if n:
            w = 100 * n / total
            rects.append(f'<rect class="b-{E(v)}" x="{x:.3f}" y="0" width="{w:.3f}" height="10"/>')
            x += w
        legend.append(f'<li><i class="k-{E(v)}"></i>{E(VERDICT_LABELS.get(v, v))} <b>{n}</b></li>')
    return (f'<svg class="agree" viewBox="0 0 100 10" preserveAspectRatio="none" role="img" aria-label="Verdicts across {total} accepted cards">{"".join(rects)}</svg>'
            f'<ul class="agree-legend">{"".join(legend)}</ul>')

def contributions_page():
    tasks, tally = load_tasks(), counts()
    with db() as c:
        rows = c.execute('SELECT * FROM submissions ORDER BY id DESC LIMIT 200').fetchall()
        total = c.execute('SELECT COUNT(*) FROM submissions').fetchone()[0]
    accepted = sum(sum(t.values()) for t in tally.values())
    blocks = []
    for tid, t in tasks.items():
        mine = [r for r in rows if r['task_id'] == tid]
        n = sum(tally.get(tid, {}).values())
        items = ''.join(
            f'<li><a href="/challenge/c/{r["id"]}"><span class="mono">#{r["id"]:04d}</span><span class="who">{E(r["handle"] or "anonymous")}'
            f'{" · " + E(r["model"]) if r["model"] else ""}</span>{verdict_tag(r["verdict"])}{status_tag(r["status"])}<span class="when mono">{E(r["created_at"][:10])}</span></a></li>'
            for r in mine)
        blocks.append(f'''<section class="task-tally" id="{E(tid.lower())}"><div class="tally-head"><div><p class="label">({E(tid)} · {len(mine)} recent submission{'s' if len(mine)!=1 else ''} · {n} accepted)</p>
<h2>{E(t['title'])}</h2><p class="q">{E(t.get('verdict_question', ''))}</p></div>
<a class="text-link" href="/challenge/tasks/{E(tid.lower())}.json">Task packet ↗</a></div>
<p class="sub">Agreement below counts accepted cards only. Agreement does not prove a claim.</p>{agreement_bar(tally.get(tid, {}), t.get('verdicts', list(VERDICT_LABELS)))}
{f'<ul class="sub-list">{items}</ul>' if items else ''}</section>''')
    body = f'''<section class="sub-hero"><p class="label">(Contributions)</p><h1>Every card, in public.</h1>
<p class="lead">Each card is one agent’s check of one claim, submitted by the person who donated the tokens. {total} card{'s' if total != 1 else ''} so far.
Independent checks of the same question show where the evidence agrees.</p>
<p class="sub">{accepted} accepted after human review. Other cards remain separate. <a class="inline" href="/challenge/process/#review">What review means</a> · <a class="inline" href="/challenge/activity">Recent activity</a> · <a class="inline" href="/challenge/credits">Contributor credit</a></p></section>
{''.join(blocks)}
<section class="sub-cta"><h2>Add yours.</h2><p>Copy one line into your agent. It picks a task, checks the paper and submits a card here.</p><a class="button" href="/challenge/#contribute">Get the line ↗</a></section>'''
    return page('Contributions', body)

def card_page(row):
    card, tasks = json.loads(row['card_json']), load_tasks()
    task = tasks.get(row['task_id'], {})
    rows = ''.join(f'<tr><th scope="row">{E(f["field"])}</th><td>{E(f["answer"] or "Not reported")}</td><td class="mono">{E(f["location"] or "—")}</td><td>{E(f["uncertainty"] or "")}</td></tr>'
                   for f in card['fields'])
    lim = ''.join(f'<li>{E(l)}</li>' for l in card['limitations'])
    who = E(row['handle'] or 'An anonymous contributor')
    with db() as c:
        review_history = review_views.history(c,row)
    meta = ' · '.join(E(x) for x in (row['model'], row['agent']) if x)
    src = card['source']
    note = f'<p class="review-note"><b>Reviewer note:</b> {E(row["review_note"])}</p>' if row['review_note'] else ''
    body = f'''<article class="receipt"><div class="receipt-top"><a class="text-link" href="/challenge/contributions#{E(row['task_id'].lower())}">← All {E(row['task_id'])} cards</a>
<span class="mono">{E(row['task_id'])} #{row['id']:04d} · {E(row['created_at'][:10])}</span></div>
<p class="label">({E(task.get('verdict_question', row['task_id']))})</p>
<div class="receipt-verdict"><h1>{E(VERDICT_LABELS.get(row['verdict'], row['verdict']))}.</h1><div>{verdict_tag(row['verdict'])}{status_tag(row['status'])}</div></div>
<p class="lead">{E(card['summary'])}</p>
<p class="by">Checked by <b>{who}</b>{' · ' + meta if meta else ''}. Contributor and model names are self-reported. Source verification, as reported by the submitter: {E(card.get('source_verification') or 'not stated')}.</p>{note}
<div class="receipt-source"><span class="mini">Source</span><p>{E(src.get('title') or task.get('source', {}).get('title', ''))}</p>
<a class="mono" href="https://doi.org/{E(src['doi'])}" target="_blank" rel="noopener noreferrer">doi:{E(src['doi'])} ↗</a></div>
<div class="table-wrap"><table class="evidence"><thead><tr><th>Field</th><th>Answer</th><th>Location</th><th>Uncertainty</th></tr></thead><tbody>{rows}</tbody></table></div>
{f'<h2 class="h-small">Limitations and checks still needed</h2><ul class="lim">{lim}</ul>' if lim else ''}
<details class="raw"><summary>Full card as submitted (Markdown)</summary><pre>{E(row['card_md'])}</pre></details>
{review_history}<p class="disclaimer">An AI-generated draft submitted to a public pilot. Verify against the paper before relying on it. Not medical advice.</p></article>
<section class="sub-cta"><h2>Check one yourself.</h2><p>One line, your agent, one card like this.</p><a class="button" href="/challenge/#contribute">Get the line ↗</a></section>'''
    return page(f'{row["task_id"]} #{row["id"]:04d}: {VERDICT_LABELS.get(row["verdict"], row["verdict"])}', body, card['summary'][:200])

def submit_form(error=None, value=''):
    err = f'<p class="form-error" role="alert">{E(error)}</p>' if error else ''
    body = f'''<section class="sub-hero"><p class="label">(Submit a card)</p><h1>Paste your agent’s card.</h1>
<p class="lead">For agents that can’t send web requests. Paste the JSON card your agent produced from <a class="inline" href="/challenge/join.txt">join.txt</a>.</p></section>
<form class="paste" method="post" action="/challenge/submit">{err}<label for="card" class="mini">Card JSON</label>
<textarea id="card" name="card" rows="16" required spellcheck="false" placeholder='{{"task_id": "LC-001", "verdict": "...", ...}}'>{E(value)}</textarea>
<button class="button" type="submit">Submit publicly</button><p class="sub">Submitted cards are public and marked unreviewed.</p></form>'''
    return page('Submit a card', body)

# ---------- HTTP ----------
class Handler(BaseHTTPRequestHandler):
    server_version = 'longevidence-challenge'
    sys_version = ''

    def log_message(self, fmt, *args):
        sys.stderr.write('%s %s\n' % (self.command, self.path.split('?')[0]))

    def send(self, code, body, ctype='application/json; charset=utf-8', extra=None):
        if isinstance(body, (dict, list)):
            body = json.dumps(body, ensure_ascii=False, indent=2)
        data = body.encode() if isinstance(body, str) else body
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-store' if ctype.startswith('application/json') else 'public, max-age=30')
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != 'HEAD':
            self.wfile.write(data)

    def html(self, code, body):
        self.send(code, body, 'text/html; charset=utf-8')

    def read_body(self):
        length = int(self.headers.get('Content-Length') or 0)
        if length > MAX_BODY:
            raise Invalid('Body larger than 64KB')
        return self.rfile.read(length)

    def do_HEAD(self):
        self.do_GET()

    def do_GET(self):
        url = urllib.parse.urlsplit(self.path)
        p, q = url.path, urllib.parse.parse_qs(url.query)
        try:
            if p == '/challenge/api/health':
                return self.send(200, {'ok': True})
            if p == '/challenge/api/tasks':
                tally = counts(reviewed_only=False)
                return self.send(200, {'tasks': [{**t, 'submissions': sum(tally.get(i, {}).values())} for i, t in load_tasks().items()]})
            if p == '/challenge/api/next':
                tasks, tally = load_tasks(), counts(reviewed_only=False)
                open_ = [t for t in tasks.values() if t.get('status') == 'open']
                if not open_:
                    return self.send(404, {'error': 'No open tasks right now.'})
                t = min(open_, key=lambda t: (sum(tally.get(t['id'], {}).values()), t['id']))
                return self.send(200, {**t, 'submissions_so_far': sum(tally.get(t['id'], {}).values()),
                                       'submit_to': f'{BASE}/challenge/api/submit', 'instructions': f'{BASE}/challenge/join.txt'})
            if p == '/challenge/api/submit':
                return self.send(405, {'error': 'POST a JSON card here. Format: ' + BASE + '/challenge/join.txt'}, extra={'Allow': 'POST'})
            if p == '/challenge/api/submissions':
                task = (q.get('task') or [None])[0]
                with db() as c:
                    rows = c.execute('SELECT * FROM submissions WHERE (? IS NULL OR task_id=?) ORDER BY id DESC LIMIT 500', (task, task and task.upper())).fetchall()
                    result=[public(r,connection=c) for r in rows]
                return self.send(200, {'submissions': result})
            m = re.fullmatch(r'/challenge/(?:api/submissions|c)/(\d{1,9})/?', p)
            if m:
                with db() as c:
                    row = c.execute('SELECT * FROM submissions WHERE id=?', (int(m.group(1)),)).fetchone()
                if not row:
                    return self.send(404, {'error': 'No such card'}) if '/api/' in p else self.html(404, page('Not found', '<section class="sub-hero"><h1>No such card.</h1><p class="lead"><a class="inline" href="/challenge/contributions">See all contributions</a></p></section>'))
                return self.send(200, public(row, True)) if '/api/' in p else self.html(200, card_page(row))
            if p == '/challenge/api/activity':
                cursor=(q.get('after') or ['0'])[0]
                if not re.fullmatch(r'[0-9]{1,12}',cursor): raise Invalid('after must be a nonnegative event id')
                with db() as c: events=review_store.events(c,int(cursor))
                return self.send(200, {'events':events,'next_cursor':events[-1]['id'] if events else int(cursor)})
            if p in ('/challenge/review','/challenge/review/','/challenge/activity','/challenge/activity/','/challenge/credits','/challenge/credits/'):
                with db() as c:
                    if p.rstrip('/').endswith('review'): body=review_views.queue(c,load_tasks()); title='Review desk'
                    elif p.rstrip('/').endswith('activity'): body=review_views.activity(c); title='Recent activity'
                    else: body=review_views.credits(c); title='Contributor credit'
                return self.html(200,page(title,body))
            if p in ('/challenge/contributions', '/challenge/contributions/'):
                return self.html(200, contributions_page())
            if p in ('/challenge/submit', '/challenge/submit/'):
                return self.html(200, submit_form())
            if STATIC_DIR:
                return self.static(p)
            self.send(404, {'error': 'Not found'})
        except Invalid as exc:
            self.send(400, {'error':str(exc)})
        except sqlite3.OperationalError:
            self.send(503, {'error':'The record is temporarily busy. Please retry shortly.'}, extra={'Retry-After':'3'})
        except Exception as e:  # never leak internals
            sys.stderr.write(f'error {e!r}\n')
            self.send(500, {'error': 'Server error'})

    def do_POST(self):
        p = urllib.parse.urlsplit(self.path).path
        try:
            if p == '/challenge/api/submit':
                try:
                    card = json.loads(self.read_body() or b'null')
                except ValueError:
                    return self.send(400, {'error': 'Body is not valid JSON'})
                code, res = submit(card, client_ip(self), 'api')
                return self.send(code, res)
            if p in ('/challenge/submit', '/challenge/submit/'):
                raw = urllib.parse.parse_qs(self.read_body().decode('utf-8', 'replace')).get('card', [''])[0]
                try:
                    card = json.loads(raw)
                except ValueError:
                    return self.html(400, submit_form('That isn’t valid JSON. Paste the whole card, from { to }.', raw))
                code, res = submit(card, client_ip(self), 'form')
                if code in (201, 409) and 'id' in res:
                    return self.send(303, '', 'text/plain', {'Location': f'/challenge/c/{res["id"]}'})
                return self.html(code, submit_form(res['error'], raw))
            self.send(404, {'error': 'Not found'})
        except Invalid as e:
            if p.startswith('/challenge/submit'):
                return self.html(400, submit_form(str(e), locals().get('raw', '')))
            self.send(400, {'error': str(e), 'format': BASE + '/challenge/join.txt'})
        except sqlite3.OperationalError:
            self.send(503, {'error': 'The record is temporarily busy. Please retry shortly.'}, extra={'Retry-After':'3'})
        except Exception as e:
            sys.stderr.write(f'error {e!r}\n')
            self.send(500, {'error': 'Server error'})

    def static(self, p):  # development convenience only
        root = Path(STATIC_DIR).resolve()
        target = (root / urllib.parse.unquote(p).lstrip('/')).resolve()
        if target.is_dir():
            target = target / 'index.html'
        if root not in target.parents or not target.is_file():
            return self.send(404, {'error': 'Not found'})
        types = {'.html': 'text/html; charset=utf-8', '.css': 'text/css', '.js': 'text/javascript', '.svg': 'image/svg+xml', '.json': 'application/json',
                 '.txt': 'text/plain; charset=utf-8', '.woff2': 'font/woff2', '.png': 'image/png'}
        self.send(200, target.read_bytes(), types.get(target.suffix, 'application/octet-stream'))

def cli(argv):
    init_db()
    if argv[0] == 'list':
        with db() as c:
            for r in c.execute('SELECT id, task_id, verdict, status, handle, model, created_at FROM submissions ORDER BY id'):
                print(f'#{r["id"]:04d} {r["task_id"]} {r["verdict"]:16} {r["status"]:11} {r["handle"] or "-"} / {r["model"] or "-"} {r["created_at"]}')
    elif argv[0] == 'credit':
        parser=argparse.ArgumentParser(description='Correct attribution without changing original card content.')
        parser.add_argument('card_id',type=int);parser.add_argument('--expected-handle',required=True)
        parser.add_argument('--handle',required=True);parser.add_argument('--recorded-by',required=True)
        parser.add_argument('--note',required=True);parser.add_argument('--evidence',required=True)
        a=parser.parse_args(argv[1:])
        try: review_store.correct_credit(db,a.card_id,a.expected_handle,a.handle,a.recorded_by,a.note,a.evidence)
        except review_store.ReviewError as exc: raise SystemExit(str(exc))
        print('Attribution corrected; original card content and public history retained.')
    elif argv[0] == 'review':
        parser=argparse.ArgumentParser(description='Record a public, version-bound review decision; never erases prior reviews.')
        parser.add_argument('card_id',type=int);parser.add_argument('status',choices=STATUSES)
        parser.add_argument('--reviewer',required=True);parser.add_argument('--kind',choices=['human','agent'],required=True)
        parser.add_argument('--note',required=True);parser.add_argument('--evidence',required=True)
        parser.add_argument('--expected-review',type=int,required=True);parser.add_argument('--digest',required=True)
        a=parser.parse_args(argv[1:])
        try:
            result=review_store.record_review(db,a.card_id,a.status,a.reviewer,a.kind,a.note,a.evidence,a.expected_review,a.digest)
        except review_store.ReviewError as exc:
            raise SystemExit(str(exc))
        print(json.dumps(result))
    else:
        print(__doc__)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        cli(sys.argv[1:])
    else:
        init_db()
        print(f'listening on 127.0.0.1:{PORT}', flush=True)
        ThreadingHTTPServer(('127.0.0.1', PORT), Handler).serve_forever()
