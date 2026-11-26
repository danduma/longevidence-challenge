"""Plain HTML views for public credit, discussion and durable review history."""
import html
import json
import urllib.parse
import review_store as store
E=html.escape

def discussion(card_id,digest):
    query=urllib.parse.urlencode({'template':'review.yml','title':f'Review card #{card_id}',
         'card':f'https://longevidence.org/challenge/c/{card_id}', 'version':digest})
    return store.REPOSITORY+'/issues/new?'+query

def history(c,row):
    items=[]
    for r in store.reviews(c,row['id']):
        items.append(f'''<li><p><b>{E(store.STATUS_LABELS[r['status']])}</b> · {E(r['reviewer'])} · {E(r['kind'])} review · <time>{E(r['created_at'])}</time></p>
<p>{E(r['note'])}</p><a class="inline" href="{E(r['evidence_url'])}">Read the public discussion ↗</a><p class="mono">Content {E(r['content_digest'][:12])} · review #{r['id']}</p></li>''')
    for event in c.execute("SELECT * FROM activity_events WHERE submission_id=? AND event_type='credit.corrected' ORDER BY id",(row['id'],)):
        p=json.loads(event['payload_json'])
        items.append(f'<li><b>Attribution correction</b><p>{E(p["previous_handle"] or "Anonymous")} → {E(p["handle"])} · recorded by {E(p["recorded_by"])}</p><p>{E(p["note"])}</p><a class="inline" href="{E(p["evidence_url"])}">Public request</a> · {E(event["created_at"])}</li>')
    empty='<p>No recorded review yet. You can help by checking a source, reproducing a calculation or explaining a concern.</p>'
    return f'''<section class="review-history"><h2>Review and corrections</h2><p>These decisions apply to this exact card. Human acceptance means a maintainer accepted it for the task; it is not journal peer review or a medical endorsement. Agent checks are labeled separately.</p>
{'<ol class="event-list">'+''.join(items)+'</ol>' if items else empty}
<p><a class="button" href="{E(discussion(row['id'],store.digest(row)))}">Review or question this card ↗</a></p>
<p>GitHub sign-in is needed to join the discussion. Reading is open. A maintainer records a decision here with a reason and a link. Earlier decisions remain visible.</p>
<p class="mono">Content SHA-256: {E(store.digest(row))}</p><p><a class="inline" href="/challenge/process/#review">How review works</a> · <a class="inline" href="/challenge/api/submissions/{row['id']}">Download this public record</a></p></section>'''

def queue(c,tasks):
    rows=c.execute("SELECT * FROM submissions WHERE status!='rejected' ORDER BY CASE status WHEN 'changes_requested' THEN 0 WHEN 'unreviewed' THEN 1 WHEN 'under_review' THEN 2 ELSE 3 END,id DESC LIMIT 200").fetchall()
    items=[]
    for r in rows:
        items.append(f'''<li><p class="label">{E(r['task_id'])} · {E(store.STATUS_LABELS.get(r['status'],r['status']))}</p><h2><a class="inline" href="/challenge/c/{r['id']}">#{r['id']:04d} · {E(tasks.get(r['task_id'],{}).get('title',r['task_id']))}</a></h2><p>Submitted by {E(r['handle'] or 'an anonymous contributor')} · {E(r['created_at'][:10])}</p><a class="text-link" href="{E(discussion(r['id'],store.digest(r)))}">Add a source check, replication or concern on GitHub ↗</a></li>''')
    return f'''<section class="sub-hero"><p class="label">(Public review desk)</p><h1>Help check the work.</h1><p class="lead">You do not need a new result to make a useful contribution. Verify one number, rerun a calculation, find earlier work, or explain where a claim goes too far.</p><p>Open the card, check its sources and leave a public GitHub review. Include what you checked and what you could not. A maintainer records the outcome; your review stays linked and credited. <a class="inline" href="/challenge/process/#review">Read the review criteria</a>.</p></section>
<ul class="work-list">{''.join(items) or '<li>No cards are waiting yet. You can submit the first check.</li>'}</ul>'''

def activity(c):
    items=[]
    for r in c.execute('SELECT * FROM activity_events ORDER BY id DESC LIMIT 200'):
        p=json.loads(r['payload_json'])
        if r['event_type']=='review.recorded':
            message=f"{p['reviewer']} recorded a {p['kind']} review: {store.STATUS_LABELS[p['status']]}"
        elif r['event_type']=='credit.corrected':message=f"Attribution corrected to {p['handle']} by {p['recorded_by']}"
        else:message=f"{p.get('handle') or 'An anonymous contributor'} submitted an evidence card"
        items.append(f'''<li><p class="label">{E(r['created_at'])} · event #{r['id']}</p><p>{E(message)} · <a class="inline" href="/challenge/c/{r['submission_id']}">{E(p.get('task_id',''))} #{r['submission_id']:04d}</a></p></li>''')
    return f'''<section class="sub-hero"><p class="label">(Recent activity)</p><h1>The work, as it happens.</h1><p class="lead">Submissions and review decisions, with dates and links to the evidence. A submission is a starting point, not an accepted result.</p><p><a class="inline" href="/challenge/api/activity">Read the public event feed as JSON</a> · <a class="inline" href="/challenge/review">Help review</a></p></section><ol class="event-list">{''.join(items) or '<li>No activity yet.</li>'}</ol>'''

def credits(c):
    groups={}
    for r in c.execute('SELECT * FROM submissions ORDER BY id'):
        handle=r['handle']
        if handle:
            key=handle.casefold();g=groups.setdefault(key,dict(name=handle,work=[],reviews=[]))
            g['work'].append(f'<a class="inline" href="/challenge/c/{r["id"]}">{E(r["task_id"])} #{r["id"]:04d}</a> — {E(store.STATUS_LABELS.get(r["status"],r["status"]))}')
    for r in c.execute('SELECT * FROM reviews ORDER BY id'):
        key=r['reviewer'].casefold();g=groups.setdefault(key,dict(name=r['reviewer'],work=[],reviews=[]))
        g['reviews'].append(f'<a class="inline" href="/challenge/c/{r["submission_id"]}">Review #{r["id"]}</a> — {E(r["kind"])} review, {E(r["created_at"][:10])}')
    cards=''.join(f'<article><h2>{E(g["name"])}</h2><ul>'+''.join('<li>'+x+'</li>' for x in g['work']+g['reviews'])+'</ul></article>' for _,g in sorted(groups.items()))
    anon=c.execute('SELECT count(*) FROM submissions WHERE handle IS NULL OR handle=""').fetchone()[0]
    return f'''<section class="sub-hero"><p class="label">(Contributor credit)</p><h1>People behind the work.</h1><p class="lead">A useful check, an honest correction and a careful review all deserve credit. Your name stays attached to the record even when its conclusion changes.</p><p>This is an alphabetical work record, not a leaderboard. Submission names are self-reported and do not verify identity. GitHub issues and pull requests carry their actual GitHub authors. Model names are tool disclosures, not people or verified identities.</p><p>{anon} anonymous contribution{'s' if anon!=1 else ''}; anonymity is welcome. To correct attribution, open a <a class="inline" href="{store.REPOSITORY}/issues/new?template=credit.yml">credit request on GitHub</a>. <a class="inline" href="/challenge/process/#credit">Read the credit policy</a>.</p></section><div class="credit-grid">{cards or '<p>No named contributions yet. The anonymous work above is still part of the public record.</p>'}</div>'''
