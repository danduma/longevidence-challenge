"""Append-only public review decisions. Only the trusted operator CLI writes here."""
import hashlib
import datetime
import json
import re
import sqlite3

REPOSITORY = 'https://github.com/danduma/longevidence-challenge'
STATUSES = ('unreviewed', 'under_review', 'changes_requested', 'accepted', 'rejected')
STATUS_LABELS = {'unreviewed':'Awaiting review','under_review':'Review in progress',
                 'changes_requested':'Changes requested','accepted':'Accepted for this task','rejected':'Not accepted'}

class ClosingConnection(sqlite3.Connection):
    def __exit__(self,*args):
        try: return super().__exit__(*args)
        finally: self.close()

class ReviewError(ValueError): pass

def now():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()

def append_event(c, event_type, card_id, payload, created_at=None, key=None):
    c.execute('INSERT OR IGNORE INTO activity_events(event_key,event_type,submission_id,payload_json,created_at) VALUES(?,?,?,?,?)',
              (key,event_type,card_id,json.dumps(payload,ensure_ascii=False),created_at or now()))

def migrate(c):
    # Execute statements individually: executescript would commit the caller's transaction.
    c.execute('''CREATE TABLE IF NOT EXISTS reviews(id INTEGER PRIMARY KEY AUTOINCREMENT,
      submission_id INTEGER NOT NULL, status TEXT NOT NULL, previous_status TEXT NOT NULL,
      reviewer TEXT NOT NULL, kind TEXT NOT NULL, note TEXT NOT NULL, evidence_url TEXT NOT NULL,
      content_digest TEXT NOT NULL, created_at TEXT NOT NULL)''')
    c.execute('''CREATE TABLE IF NOT EXISTS activity_events(id INTEGER PRIMARY KEY AUTOINCREMENT,
      event_key TEXT UNIQUE, event_type TEXT NOT NULL, submission_id INTEGER NOT NULL,
      payload_json TEXT NOT NULL, created_at TEXT NOT NULL)''')
    c.execute('CREATE INDEX IF NOT EXISTS review_card ON reviews(submission_id,id)')
    for row in c.execute('SELECT id,task_id,handle,model,status,created_at FROM submissions ORDER BY id').fetchall():
        append_event(c,'submission.created',row['id'],dict(task_id=row['task_id'],handle=row['handle'],model=row['model'],
                     status_at_import=row['status']),row['created_at'],f"submission:{row['id']}")

def digest(row):
    if 'card_json' in row.keys():
        canonical=json.dumps(json.loads(row['card_json']),sort_keys=True,separators=(',',':'),ensure_ascii=False)
        return hashlib.sha256((canonical+'\n'+row['card_md']).encode()).hexdigest()
    return row['digest']

def reviews(c,card_id):
    return [dict(r) for r in c.execute('SELECT * FROM reviews WHERE submission_id=? ORDER BY id',(card_id,))]

def version(c,card_id):
    return c.execute('SELECT COALESCE(MAX(id),0) FROM reviews WHERE submission_id=?',(card_id,)).fetchone()[0]

def events(c,after=0):
    result=[]
    for row in c.execute('SELECT * FROM activity_events WHERE id>? ORDER BY id LIMIT 200',(after,)):
        r=dict(row);r['payload']=json.loads(r.pop('payload_json'));r.pop('event_key');result.append(r)
    return result

def record_review(connect,card_id,status,reviewer,kind,note,evidence_url,expected_review_id,expected_digest):
    if status not in STATUSES: raise ReviewError('Choose a documented review status.')
    if kind not in ('human','agent'): raise ReviewError('Review kind must be human or agent.')
    if status=='accepted' and kind!='human': raise ReviewError('An agent check cannot confer human acceptance.')
    if not isinstance(reviewer,str) or not re.fullmatch(r'[A-Za-z0-9_.@ -]{1,80}',reviewer.strip()):
        raise ReviewError('Name the reviewer using a public handle.')
    if not isinstance(note,str) or not 20<=len(note.strip())<=4000:
        raise ReviewError('Explain what was checked, the outcome and its limits (20–4000 characters).')
    prefix=re.escape(REPOSITORY)
    if not isinstance(evidence_url,str) or not re.fullmatch(prefix+r'/(issues|pull)/[1-9][0-9]*(#[A-Za-z0-9_-]+)?',evidence_url):
        raise ReviewError('Link the public review discussion or pull request in the project repository.')
    if not isinstance(expected_review_id,int) or expected_review_id<0:
        raise ReviewError('Supply the review_version read from the public card API.')
    with connect() as c:
        c.execute('BEGIN IMMEDIATE')
        row=c.execute('SELECT * FROM submissions WHERE id=?',(card_id,)).fetchone()
        if row is None: raise ReviewError('No such card.')
        if digest(row)!=expected_digest or version(c,card_id)!=expected_review_id:
            raise ReviewError('The content or review changed. Read the current card before deciding again.')
        if status=='accepted' and row['handle'] and row['handle'].casefold().lstrip('@')==reviewer.strip().casefold().lstrip('@'):
            raise ReviewError('Acceptance needs a reviewer other than the named contributor.')
        stamp=now()
        cur=c.execute('''INSERT INTO reviews(submission_id,status,previous_status,reviewer,kind,note,evidence_url,content_digest,created_at)
          VALUES(?,?,?,?,?,?,?,?,?)''',(card_id,status,row['status'],reviewer.strip(),kind,note.strip(),evidence_url,digest(row),stamp))
        review_id=cur.lastrowid
        c.execute('UPDATE submissions SET status=?,review_note=? WHERE id=?',(status,note.strip(),card_id))
        append_event(c,'review.recorded',card_id,dict(review_id=review_id,previous_status=row['status'],status=status,
                     reviewer=reviewer.strip(),kind=kind,note=note.strip(),evidence_url=evidence_url,task_id=row['task_id']),stamp,f'review:{review_id}')
    return dict(id=review_id,submission_id=card_id,status=status)

def correct_credit(connect,card_id,expected_handle,new_handle,recorded_by,note,evidence_url):
    if not re.fullmatch(re.escape(REPOSITORY)+r'/issues/[1-9][0-9]*(#[A-Za-z0-9_-]+)?',evidence_url):
        raise ReviewError('Link a public attribution request in the project repository.')
    for value in [new_handle,recorded_by]:
        if not re.fullmatch(r'[A-Za-z0-9_.@ -]{1,80}',value): raise ReviewError('Use a public handle.')
    if not 20<=len(note)<=4000: raise ReviewError('Explain the attribution correction.')
    with connect() as c:
        c.execute('BEGIN IMMEDIATE')
        row=c.execute('SELECT * FROM submissions WHERE id=?',(card_id,)).fetchone()
        if row is None or (row['handle'] or '')!=expected_handle:
            raise ReviewError('Attribution changed; inspect the current record.')
        key=f'credit:{card_id}:{evidence_url}'
        if c.execute('SELECT id FROM activity_events WHERE event_key=?',(key,)).fetchone():
            raise ReviewError('This attribution request has already been recorded.')
        c.execute('UPDATE submissions SET handle=? WHERE id=?',(new_handle,card_id))
        append_event(c,'credit.corrected',card_id,dict(task_id=row['task_id'],previous_handle=row['handle'],
                     handle=new_handle,recorded_by=recorded_by,note=note,evidence_url=evidence_url),key=key)
