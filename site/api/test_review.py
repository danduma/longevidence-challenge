import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
import review_store as store


class ReviewWorkflow(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = str(Path(self.tmp.name) / 'cards.db')
        with self.connect() as c:
            c.execute('CREATE TABLE submissions(id INTEGER PRIMARY KEY,task_id TEXT,verdict TEXT,handle TEXT,model TEXT,status TEXT,review_note TEXT,digest TEXT,created_at TEXT)')
            c.execute("INSERT INTO submissions VALUES(1,'LC-001','not_established','alice','test-model','unreviewed',NULL,'abc','2026-11-24T12:00:00+00:00')")
            store.migrate(c)
        self.addCleanup(self.tmp.cleanup)

    def connect(self):
        c = sqlite3.connect(self.path, timeout=.1)
        c.row_factory = sqlite3.Row
        return c

    def record(self, status='accepted', **kw):
        args = dict(connect=self.connect, card_id=1, status=status, reviewer='bob', kind='human',
                    note='Checked source table and sample counts; acceptance applies only to this card.',
                    evidence_url='https://github.com/danduma/longevidence-challenge/issues/1',
                    expected_review_id=0, expected_digest='abc')
        args.update(kw)
        return store.record_review(**args)

    def test_migration_keeps_existing_card_and_is_idempotent(self):
        with self.connect() as c:
            store.migrate(c)
            self.assertEqual(c.execute('SELECT status FROM submissions').fetchone()[0], 'unreviewed')
            self.assertEqual(c.execute('SELECT count(*) FROM activity_events').fetchone()[0], 1)

    def test_decision_preserves_named_history_and_reopens_without_erasing(self):
        first = self.record()
        self.record(status='changes_requested', expected_review_id=first['id'], reviewer='carol')
        with self.connect() as c:
            history = store.reviews(c, 1)
            self.assertEqual([x['status'] for x in history], ['accepted', 'changes_requested'])
            self.assertEqual([x['reviewer'] for x in history], ['bob', 'carol'])
            self.assertEqual(c.execute('SELECT status FROM submissions').fetchone()[0], 'changes_requested')
            self.assertEqual(len(store.events(c)), 3)

    def test_stale_decision_does_not_overwrite_newer_review(self):
        self.record()
        with self.assertRaises(store.ReviewError): self.record(status='rejected')
        with self.connect() as c: self.assertEqual(len(store.reviews(c, 1)), 1)

    def test_wrong_content_version_cannot_be_reviewed(self):
        with self.assertRaises(store.ReviewError): self.record(expected_digest='changed')

    def test_rationale_reviewer_and_public_evidence_are_required(self):
        for kw in [dict(note=''),dict(reviewer=''),dict(evidence_url=''),dict(evidence_url='javascript:alert(1)'),dict(evidence_url='https://github.com/other/private/issues/1')]:
            with self.subTest(kw=kw), self.assertRaises(store.ReviewError): self.record(**kw)

    def test_agent_review_cannot_award_human_acceptance(self):
        with self.assertRaises(store.ReviewError): self.record(kind='agent')
        self.record(status='under_review',kind='agent')

    def test_named_author_cannot_accept_own_card(self):
        with self.assertRaises(store.ReviewError): self.record(reviewer='ALICE')

    def test_unknown_card_and_bad_status_do_not_create_events(self):
        for kw in [dict(card_id=999),dict(status='verified_forever')]:
            with self.assertRaises(store.ReviewError): self.record(**kw)
        with self.connect() as c: self.assertEqual(len(store.events(c)),1)

    def test_public_event_contains_no_private_submission_data(self):
        self.record()
        with self.connect() as c:
            body=json.dumps(store.events(c))
            self.assertNotIn('ip_hash',body)
            self.assertNotIn('salt',body)
            self.assertIn('review.recorded',body)

    def test_receipt_credit_is_retained_after_rejection(self):
        self.record(status='rejected')
        with self.connect() as c:
            row=c.execute('SELECT * FROM submissions').fetchone()
            self.assertEqual(row['handle'],'alice')
            self.assertEqual(store.reviews(c,1)[0]['reviewer'],'bob')

    def test_credit_correction_is_versioned_and_keeps_original_attribution_in_history(self):
        store.correct_credit(self.connect,1,'alice','alice-fixed','maintainer','Confirmed public attribution correction.',store.REPOSITORY+'/issues/2')
        with self.connect() as c:
            self.assertEqual(c.execute('SELECT handle FROM submissions').fetchone()[0],'alice-fixed')
            change=store.events(c)[-1]
            self.assertEqual(change['event_type'],'credit.corrected')
            self.assertEqual(change['payload']['previous_handle'],'alice')
        with self.assertRaises(store.ReviewError):
            store.correct_credit(self.connect,1,'alice','other','maintainer','Outdated request cannot replace credit.',store.REPOSITORY+'/issues/3')

    def test_event_cursor_replays_only_new_records(self):
        with self.connect() as c: cursor=store.events(c)[0]['id']
        self.record()
        with self.connect() as c:
            result=store.events(c,after=cursor)
            self.assertEqual(len(result),1)
            self.assertEqual(result[0]['event_type'],'review.recorded')

if __name__=='__main__': unittest.main()
