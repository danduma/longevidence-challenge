import json
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from pathlib import Path
from http.server import ThreadingHTTPServer
import server
import review_store

class PublicWorkflow(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        server.DB_PATH=str(Path(self.tmp.name)/'cards.db')
        server.TASKS_DIR=Path(__file__).resolve().parents[1]/'dist/challenge/tasks'
        server.init_db()
        self.http=ThreadingHTTPServer(('127.0.0.1',0),server.Handler)
        self.thread=threading.Thread(target=self.http.serve_forever,daemon=True);self.thread.start()
        self.base=f'http://127.0.0.1:{self.http.server_port}'
        self.addCleanup(self.stop)
        task=server.load_tasks()['LC-001']
        self.card=dict(task_id='LC-001',verdict='not_established',source=dict(doi=task['source']['doi']),
          summary='A test of a specific source-based conclusion.',fields=[dict(field=f,answer='Unknown',location='Methods') for f in task['fields']],
          contributor=dict(handle='alice',model='test-agent'),card_markdown='Source-linked evidence draft. '*15)
    def stop(self):
        self.http.shutdown();self.http.server_close();self.thread.join()
    def request(self,path,data=None):
        raw=None if data is None else json.dumps(data).encode()
        req=urllib.request.Request(self.base+path,data=raw,headers={'Content-Type':'application/json'})
        try:
            with urllib.request.urlopen(req) as r:return r.status,r.read().decode()
        except urllib.error.HTTPError as e:
            with e: return e.code,e.read().decode()
    def test_submit_review_and_readback_use_real_public_event_contract(self):
        status,body=self.request('/challenge/api/submit',self.card);self.assertEqual(status,201)
        card=json.loads(self.request('/challenge/api/submissions/1')[1])
        self.assertEqual(card['status'],'unreviewed');self.assertEqual(card['review_version'],0)
        with server.db() as c:self.assertEqual(server.counts(),{})
        result=review_store.record_review(server.db,1,'accepted','bob','human','I checked the primary source and its limits.',review_store.REPOSITORY+'/issues/1',0,card['content_sha256'])
        card=json.loads(self.request('/challenge/api/submissions/1')[1])
        self.assertEqual(card['status'],'accepted');self.assertEqual(card['review_version'],result['id'])
        self.assertEqual(card['reviews'][0]['reviewer'],'bob')
        self.assertEqual(server.counts()['LC-001']['not_established'],1)
        feed=json.loads(self.request('/challenge/api/activity')[1]);self.assertEqual([r['event_type'] for r in feed['events']],['submission.created','review.recorded'])
        self.assertNotIn('ip_hash',json.dumps(feed));self.assertNotIn('salt',json.dumps(card))
        html=self.request('/challenge/c/1')[1];self.assertIn('bob',html);self.assertIn('Review and corrections',html)
    def test_duplicate_and_invalid_cards_do_not_create_extra_events(self):
        self.assertEqual(self.request('/challenge/api/submit',self.card)[0],201)
        self.assertEqual(self.request('/challenge/api/submit',self.card)[0],409)
        self.assertEqual(self.request('/challenge/api/submit',{})[0],400)
        self.assertEqual(len(json.loads(self.request('/challenge/api/activity')[1])['events']),1)
    def test_public_review_writes_are_not_authorized(self):
        self.assertEqual(self.request('/challenge/api/review',dict(status='accepted'))[0],404)
    def test_empty_pages_and_bad_cursor_have_clear_responses(self):
        for path in ['/challenge/review','/challenge/activity','/challenge/credits']:
            code,body=self.request(path);self.assertEqual(code,200);self.assertIn('<h1>',body)
        self.assertEqual(self.request('/challenge/api/activity?after=bad')[0],400)
    def test_public_total_does_not_confuse_recent_list_limit_with_all_contributions(self):
        self.request('/challenge/api/submit',self.card)
        with server.db() as c:
            original=c.execute('SELECT * FROM submissions WHERE id=1').fetchone()
            for i in range(2,202):
                c.execute("INSERT INTO submissions(task_id,verdict,card_json,card_md,digest,via,created_at) VALUES(?,?,?,?,?,?,?)",('LC-001','not_established',original['card_json'],original['card_md'],str(i),'test',original['created_at']))
        code,html=self.request('/challenge/contributions')
        self.assertEqual(code,200)
        self.assertTrue('201 cards so far' in html, 'Public total must include cards beyond the 200-item recent list.')

    def test_database_context_closes_connection_after_read(self):
        c=server.db()
        with c: c.execute('SELECT 1')
        import sqlite3
        with self.assertRaises(sqlite3.ProgrammingError): c.execute('SELECT 1')

    def test_wal_read_paths_still_work_while_another_writer_holds_lock(self):
        self.request('/challenge/api/submit',self.card)
        with server.db() as c:
            c.execute('BEGIN IMMEDIATE')
            for path in ['/challenge/api/health','/challenge/api/tasks','/challenge/api/submissions/1','/challenge/api/activity']:
                self.assertEqual(self.request(path)[0],200)
            c.rollback()

    def test_full_content_fingerprint_includes_structured_answers(self):
        self.request('/challenge/api/submit',self.card)
        with server.db() as c:
            row=dict(c.execute('SELECT * FROM submissions').fetchone())
        original=review_store.digest(row)
        data=json.loads(row['card_json']);data['fields'][0]['answer']='Changed'
        row['card_json']=json.dumps(data)
        self.assertNotEqual(original,review_store.digest(row))

if __name__=='__main__':unittest.main()
