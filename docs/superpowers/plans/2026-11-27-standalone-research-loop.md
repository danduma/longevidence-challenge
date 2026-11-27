# Standalone research loop implementation plan

**Goal:** Build and run a minimal research agent without OmniHarness, pursuing 100 distinct verified contributions without lowering the evidence standard.
**Architecture:** A Python standard-library CLI owns a SQLite task/event ledger. Disposable external agent processes investigate and separately review each candidate. Reproducible checks and immutable evidence manifests gate acceptance. Public JSON/Markdown exports expose work, review, data gaps and failures.
**Tech Stack:** Python 3, SQLite, subprocess, configurable agent command; direct Codex CLI for this installation.
**Current Milestone:** Working standalone runner, failure/recovery tests, real research execution, and an honest verified contribution ledger.
**Spec:** This document and `standalone/README.md`.

No OmniHarness code, configuration, credentials bridge, goal API or process is required. Existing public research is read-only context. Existing contributions are not recounted. One investigation is one contribution; parameter rows, source downloads, routine reproductions, generated tasks and gap cards do not inflate the target. Automated verification and agent review do not establish human acceptance, novelty or experimental validation.

## Files and ownership

- `standalone/loop.py`: CLI, single runner lock, SQLite state/events, process lifecycle, verification and export.
- `standalone/protocol.py`: validated task/result/review contracts and evidence hashing.
- `standalone/test_loop.py`: persistence, recovery, verifier, deduplication and export tests.
- `standalone/seeds.json`: distinct, source-grounded initial investigations, not 100 invented results.
- `standalone/README.md`: setup, operations, counting, review and data-gap explanations.
- `standalone/.gitignore`: operational DBs, worker transcripts and temporary artifacts excluded from source control.
- Runtime directory outside tracked sources: jobs, evidence, SQLite ledger, exported public records, runner PID/log.

## Task 1: Persistence and evidence contracts

- [x] Write and run failing tests for durable transitions, duplicates, failed checks, tampering and crash recovery.
- [x] Implement contracts, transactionally recorded named events and single-runner ownership.

**Verification:** `python3 -m unittest discover -s standalone -v`: initial import failure confirmed missing implementation; the expanded 31-test suite passes (27 November). It covers real subprocess fixtures, orphan cleanup, shared-source reopening, reviewer rejection, post-review tampering, timeouts and interruption/failure distinction.

## Task 2: Real agent execution and independent review

- [x] Implement configurable producer/reviewer processes, bounded episodes, preserved logs, timeouts and retry states.
- [x] Implement fresh reproducibility checks, hash-pinned evidence, separate reviewer verdicts, honest target counts and safe restart. **Resume repair verified:** immutable numerical receipts and append-only investigation assessments are separate; the live target counter is 7, not 30.
- [x] Implement validated follow-up queue expansion, explicit data gates, re-open and watch operation.

**Verification:** focused integration tests plus a real direct-provider episode; no mocked success is counted as research.

## Task 3: Research frontier and public records

- [x] Seed genuine unresolved investigations from the controlling reassessment; document counting and review rules.
- [x] Export contributions, gaps, experiments, producer candidates, source provenance and event history as readable Markdown and JSON.

**Verification:** inspect real exports, link them to job evidence and verify export does not include operational transcripts or credentials.

## Task 4: Execute toward 100 verified contributions

- [x] Run the standalone coordinator with target 100 and record real outcomes.
- [x] Audit accepted records for shared investigations and superseded conclusions; enforce one count per investigation without discarding evidence.
- [ ] Reach 100 distinct contributions passing evidence, reproducibility and separate review gates.
- [ ] Verify final ledger and update this checklist with evidence and any remaining gates.

**Verification:** fresh `status`, `audit`, test suite and inspection of accepted records. If the target is not reached, leave its checkbox unchecked and state the observed count and blockers; running jobs are not completed contributions.

## Execution evidence (27 November, ongoing)

- Runtime: `research/standalone-run-2026-11-27` in the operational workspace, outside the public checkout. Start command: `python3 standalone/loop.py --root /absolute/path/to/run run --target 100 --watch`. Normal `~/.codex` login; parent harness credential/config variables removed.
- Six source-grounded initial questions; a seventh follow-up was approved by the separate reviewer.
- Real protease candidate passed fresh computational verification, then received **revise**, not acceptance: its reviewer identified another feasible sequence-informed comparator. The candidate does **not** count. Its evidence and review remain in the runtime ledger.
- The original first job exposed a shell-network restriction; explicit workspace network access resolved it. A later artifact-directory contract failure was reproduced in a failing test, fixed by hashing individual contained files, and safely retried. Neither failure became a scientific data gate.
- Initial publication snapshot: **0 / 100**. On automatic session resume, the detached coordinator remains alive (PID 93735) and has 30 accepted computational receipts, two data gates, one task needing attention and one running task. All 30 receipt manifests pass a fresh integrity audit. Distinct investigation count is under an independent cross-record audit; 30 accepted receipts must not be presented as 30 distinct contributions. The 100-contribution and final-ledger gates remain unchecked.
- Verification: 25 standalone tests pass; 20 existing API/review tests pass; compilation and diff checks pass. A real earlier protease verifier passed, and a deliberately corrupted result copy failed. Real public export includes the revision history and excludes operational prompts/credential identifiers.
- No claim of human acceptance, biological validation, 100 completed findings, website deployment or automatic public outreach.
- First source publication: commit `6e9e495`; [CI run 45913545083](https://github.com/danduma/longevidence-challenge/actions/runs/45913545083) passed the standalone suite, API suite and existing thermal benchmark.
- Live-run correction: the reviewer rejected the Rubisco data gate because public checkpoint regeneration and a known target-row discrepancy remain untested. A regression test first failed when rejection closed the question; it now passes with rejection returning the question to pending work. The task was explicitly reopened with those required actions. Three unsuccessful reviews become `needs_attention`, never scientific closure. The 25-test suite passes after the correction; the updated coordinator was restarted. This lesson is recorded in the research-audit skill and its case study.

## Resumption verification — 27 November

- Three independent domain audits group the first 30 accepted receipts into seven investigations; pRR1 lineage is retained as supporting source validation without another count. Original receipts and all cached evidence remain unchanged. The assessments, audit reports, independent numerical-check results and two data-gap/experiment briefs are published in `research/standalone-2026-11-27/`.
- `standalone/investigations.py` adds append-only counting assessments. Acceptance and assignment are transactional. Existing receipts remain uncounted until assessed. New controls/refinements reuse their investigation identity; independently justified questions can add a count. Five focused regression scenarios cover shared-family counting, genuinely independent questions, missing legacy assessments, duplicate question identities and rollback of an invalid assignment.
- Fresh verification: 31 standalone tests and 20 existing API tests pass; compilation, skill validation and diff checks pass. Fresh independent checks reproduce selected raw-workbook, prediction and stress-score results. Full inference, nested fitting and all mechanical solves were not rerun by the resumption audit.
- The active cryo job completed its fresh verifier successfully, then its reviewer rejected independent-contribution status while retaining it as useful evidence within the existing thermomechanical investigation. It remains open for scoped validation under the corrected protocol, not another count.
- A consistent operational SQLite backup was saved before importing assessments. The import assigned all 30 receipts, yielding **7 / 100** counted investigations, **30** accepted receipts, **0** unassessed receipts and **0** human acceptances. A fresh full receipt hash audit passes.
- The protease inventory was reopened with the independent audit's precise scope correction and ordinary review gates intact. Updated standalone coordinator PID **98736** is running with target **100**, independent of OmniHarness. The 100-contribution and final-ledger gates remain unchecked.

- Subsequent real-provider gate check: the reopened PERRC release-eligibility audit passed its fresh verifier and separate review with an explicit independent investigation identity. Live state advanced to **8 / 100** across **31** accepted receipts, with **0** unassessed receipts and **0** human acceptances. Its declared source scripts, report and review are included in the public snapshot; the earlier unaccepted attempts are not counted again. The cryo continuation remains within its original investigation.
