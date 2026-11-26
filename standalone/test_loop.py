import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from loop import Store, Runner, RunnerLock, execute, process_birth
from protocol import validate_task, validate_gap, validate_review, ValidationError, evidence_manifest, REVIEW_CHECKS


def task(key="test", question="Does the repair change predictions?"):
    return dict(id=key, title=key, question=question,
                increment="Compare repaired and original features, not a repeat of identity verification.",
                why_useful="Decides whether the correction changes benchmark interpretation.",
                context=["https://example.org/primary"], depends_on=[], watches=[])


class Contracts(unittest.TestCase):
    def test_reject_empty_or_duplicate_question_identity(self):
        with self.assertRaises(ValidationError):
            validate_task({"id": "empty"})

    def test_data_gate_needs_discriminating_evidence_and_design(self):
        with self.assertRaises(ValidationError):
            validate_gap({"missing_data": "more experiments"})

    def test_evidence_paths_cannot_escape_job(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(ValidationError):
                evidence_manifest(Path(folder), ["../private.txt"])

    def test_review_cannot_accept_with_a_failed_scientific_gate(self):
        review = dict(verdict="accept", reason="Mostly works", strongest_objection="Duplicate claim",
                      checks={k: True for k in REVIEW_CHECKS}, approved_followups=[])
        review["checks"]["distinct"] = False
        with self.assertRaises(ValidationError):
            validate_review(review, {"kind": "contribution"})

    def test_review_cannot_invent_followups(self):
        with self.assertRaises(ValidationError):
            validate_review(dict(verdict="reject", reason="Known", strongest_objection="Known",
                                 approved_followups=["invented"]), {"kind": "not_useful"})

    def test_evidence_directories_expand_to_individual_hashed_files(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "sources").mkdir()
            (root / "sources" / "data.csv").write_text("1,2,3")
            manifest = evidence_manifest(root, ["sources/"])
            self.assertEqual(set(manifest), {"sources/data.csv"})


class Persistence(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.store = Store(self.root)

    def tearDown(self):
        self.store.close()
        self.temp.cleanup()

    def test_duplicate_questions_and_ids_do_not_inflate_queue(self):
        self.assertTrue(self.store.add(task()))
        self.assertFalse(self.store.add(task()))
        self.assertFalse(self.store.add(task("other", " DOES the repair CHANGE predictions? ")))
        self.assertEqual(self.store.status()["tasks"], {"pending": 1})
        self.assertEqual(self.store.status()["verified_contributions"], 0)

    def test_claim_and_transition_are_durable_named_events(self):
        self.store.add(task())
        claimed = self.store.claim()
        self.assertEqual(claimed["id"], "test")
        self.assertIsNone(self.store.claim())
        self.store.finish("test", "exhausted", {"reason": "Already established in the source."})
        self.store.close()
        self.store = Store(self.root)
        self.assertEqual(self.store.status()["tasks"], {"exhausted": 1})
        self.assertEqual([e["name"] for e in self.store.events()],
                         ["task.queued", "task.started", "task.exhausted"])

    def test_unfinished_task_recovered_without_false_acceptance(self):
        self.store.add(task())
        self.store.claim()
        self.store.recover()
        self.assertEqual(self.store.status()["tasks"], {"pending": 1})
        self.assertEqual(self.store.events()[-1]["name"], "task.recovered")
        self.assertEqual(self.store.status()["verified_contributions"], 0)

    def test_dependencies_wait_for_accepted_results(self):
        self.store.add(task())
        dependent = task("dependent", "Does the downstream model change?")
        dependent["depends_on"] = ["test"]
        self.store.add(dependent)
        self.assertEqual(self.store.claim()["id"], "test")
        self.assertIsNone(self.store.claim())

    def test_invalid_transition_does_not_create_fake_event(self):
        self.store.add(task())
        with self.assertRaises(ValueError):
            self.store.finish("test", "accepted", {})
        self.assertEqual(self.store.status()["verified_contributions"], 0)

    def test_single_runner_lock(self):
        with RunnerLock(self.root):
            with self.assertRaises(RuntimeError):
                with RunnerLock(self.root):
                    pass

    def test_export_excludes_private_worker_logs(self):
        self.store.add(task())
        (self.root / "private.log").write_text("a private transcript")
        self.store.export()
        exported = (self.root / "public" / "ledger.json").read_text()
        self.assertNotIn("a private transcript", exported)
        self.assertTrue((self.root / "public" / "README.md").exists())

    def test_data_gap_does_not_count_and_reopens_explicitly(self):
        self.store.add(task())
        self.store.claim()
        self.store.finish("test", "data_gated", {"reason": "Need matched assay measurements"})
        self.assertEqual(self.store.status()["verified_contributions"], 0)
        self.store.reopen("test", "Matched data released")
        self.assertEqual(self.store.claim()["id"], "test")
        self.assertTrue(any(e["name"] == "task.reopened" for e in self.store.events()))

    def test_operator_interrupts_do_not_exhaust_technical_retry_budget(self):
        self.store.add(task())
        for _ in range(4):
            claimed = self.store.claim()
            self.store.retry("test", claimed["attempt"], KeyboardInterrupt())
        self.assertEqual(self.store.status()["tasks"], {"pending": 1})
        self.assertEqual(self.store.tasks()[0]["failures"], 0)

    def test_three_actual_technical_failures_need_attention_not_data(self):
        self.store.add(task())
        for _ in range(3):
            claimed = self.store.claim()
            self.store.retry("test", claimed["attempt"], RuntimeError("provider unavailable"))
            with self.store.db:
                self.store.db.execute("UPDATE tasks SET next_at=0")
        self.assertEqual(self.store.status()["tasks"], {"needs_attention": 1})
        self.assertEqual(self.store.status()["verified_contributions"], 0)

    def test_same_changed_evidence_reopens_all_affected_questions(self):
        evidence = self.root / "measurement.csv"
        evidence.write_text("old")
        for key in ["first", "second"]:
            spec = task(key, f"Question {key}")
            spec["watches"] = [str(evidence)]
            self.store.add(spec)
            self.store.claim()
            self.store.finish(key, "data_gated", {"reason": "Need updated measurements"})
        runner = Runner(self.store, self.root)
        runner.refresh()
        evidence.write_text("new")
        runner.refresh()
        self.assertEqual(self.store.status()["tasks"], {"pending": 2})

    def test_recovery_stops_matching_orphan_before_requeuing(self):
        self.store.add(task())
        self.store.claim()
        job = self.root / "jobs" / "orphan"
        job.mkdir(parents=True)
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"], start_new_session=True)
        try:
            (job / "produce.process.json").write_text(json.dumps(dict(pid=child.pid, birth=process_birth(child.pid))))
            self.store.recover()
            child.wait(timeout=5)
            self.assertIsNotNone(child.returncode)
            self.assertEqual(self.store.status()["tasks"], {"pending": 1})
            self.assertTrue(any(e["name"] == "worker.orphan_stopped" for e in self.store.events()))
        finally:
            if child.poll() is None:
                child.kill()
                child.wait()


class Processes(unittest.TestCase):
    def test_command_timeout_stops_its_process(self):
        with tempfile.TemporaryDirectory() as folder:
            log = Path(folder) / "worker.log"
            with self.assertRaises(subprocess.TimeoutExpired):
                execute([sys.executable, "-c", "import time; time.sleep(60)"], Path(folder), .05, log)
            marker = json.loads(log.with_suffix(".process.json").read_text())
            self.assertIn("finished", marker)
            self.assertIsNotNone(marker["returncode"])


class Verification(unittest.TestCase):
    def test_failed_reproduction_cannot_be_accepted(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            store = Store(root)
            runner = Runner(store, context=root)
            job = root / "jobs" / "trial"
            job.mkdir(parents=True)
            (job / "check.py").write_text("raise SystemExit(1)\n")
            with self.assertRaises(ValidationError):
                runner.check(job, {"verification": [["python3", "check.py"]]})
            self.assertEqual(store.status()["verified_contributions"], 0)
            store.close()


class EpisodeIntegration(unittest.TestCase):
    """Offline protocol fixtures test infrastructure, never real research output."""
    def run_fixture(self, mode):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        driver = root / "fixture_agent.py"
        driver.write_text('''import json,sys
from pathlib import Path
output=Path(sys.argv[1]); mode=sys.argv[2]
if output.name == "review.json":
 if mode == "tamper": Path("report.md").write_text("altered after verification")
 verdict="reject" if mode == "reject" else "accept"
 answer=dict(verdict=verdict,reason="Offline infrastructure fixture only",strongest_objection="Synthetic fixture has no scientific value",checks={k:True for k in ["distinct","source_supported","reproducible","useful","scope_honest","prior_art_checked"]},approved_followups=[])
else:
 Path("source.txt").write_text("Offline fixture source; not actual research")
 Path("report.md").write_text("A synthetic infrastructure test, not a scientific contribution")
 Path("data.json").write_text("[1,2,3]")
 Path("results.json").write_text('{"mean":2.0}')
 Path("verify.py").write_text("import json,statistics\\nfrom pathlib import Path\\nassert statistics.mean(json.loads(Path('data.json').read_text())) == json.loads(Path('results.json').read_text())['mean']\\n")
 answer=dict(kind="contribution",title="Synthetic fixture",summary="Infrastructure test only",prior_art="Known arithmetic",limitations="No scientific value; tests protocol only",sources=[dict(url="https://example.org/fixture",path="source.txt",read_scope="synthetic fixture",role="test only")],artifacts=["report.md","data.json","results.json","verify.py"],verification=[["python3","verify.py"]],claim="Fixture mean equals two",increment="fixture",decision_changed="fixture",verification_scope="fixture",followups=[])
output.write_text(json.dumps(answer))
''')
        store = Store(root / "run")
        self.addCleanup(store.close)
        store.add(task())
        runner = Runner(store, root, [sys.executable, str(driver), "{output}", mode])
        runner.episode(store.claim())
        return store

    def test_real_subprocess_protocol_accepts_only_after_checks_and_review(self):
        store = self.run_fixture("accept")
        self.assertEqual(store.status()["verified_contributions"], 1)
        self.assertTrue(store.audit()["ok"])
        names = [e["name"] for e in store.events()]
        self.assertLess(names.index("agent.finished"), names.index("contribution.verified"))

    def test_review_rejection_never_counts(self):
        store = self.run_fixture("reject")
        self.assertEqual(store.status()["verified_contributions"], 0)
        self.assertEqual(store.status()["tasks"], {"pending": 1})
        self.assertTrue(any(e["name"] == "candidate.rejected" for e in store.events()))

    def test_reviewer_cannot_mutate_evidence_and_still_accept(self):
        store = self.run_fixture("tamper")
        self.assertEqual(store.status()["verified_contributions"], 0)
        self.assertTrue(any(e["name"] == "error.surfaced" for e in store.events()))

    def test_zero_verification_commands_is_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            store = Store(Path(folder))
            runner = Runner(store, context=Path(folder))
            with self.assertRaises(ValidationError):
                runner.check(Path(folder), {"verification": []})
            store.close()

    def test_accepted_evidence_tampering_is_detected(self):
        store = self.run_fixture("accept")
        record = store.tasks()[0]["result"]
        (store.root / record["job"] / "report.md").write_text("changed after acceptance")
        self.assertFalse(store.audit()["ok"])


if __name__ == "__main__":
    unittest.main()
