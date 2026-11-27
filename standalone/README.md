# Standalone research loop

Run a persistent research queue without OmniHarness. Python 3.11+ and its standard library own all coordination and storage. A configurable command supplies inference. The default is the ordinary Codex CLI and local provider login; no OmniHarness installation, process, API, configuration or credential bridge is used.

The runner investigates a question, runs its reproducibility checks, asks a separate agent process to challenge the evidence, records the result and picks the next question. It can propose additional questions, subject to a separate task-selection review. When useful executable work runs out, it waits rather than inventing contributions.

## Start

From the repository root:

```sh
python3 -m unittest discover -s standalone -v
python3 standalone/loop.py --root /absolute/path/to/research-run init
python3 standalone/loop.py --root /absolute/path/to/research-run run --target 100 --watch
```

The default worker needs `codex` on PATH and a valid ordinary login in `~/.codex`. It uses `--ignore-user-config` so a parent application's configuration is not required. No provider SDK or Python package installation is needed by the coordinator. A particular investigation may need scientific packages in its own job directory; it must document those dependencies.

To try one episode, pass `--max-episodes 1`. To drain the current frontier and return when no work is executable, omit `--watch`. The target is a stop condition, not a promise that 100 useful findings exist.

`--watch` keeps the process alive while it waits. Its default polling interval is one hour, configurable with `--poll`. Explicitly configured watched files/HTTPS endpoints can reopen data-gated tasks when their content changes. It does not discover unpublished laboratory data or automatically poll every paper on the internet. An operator can reopen a task when new evidence arrives.

Use a process manager appropriate for your machine to start the same command after reboot. The filesystem lock prevents two coordinators from running this queue at once. Restarting resumes unfinished tasks from their durable evidence and previous review, using a new job directory. Each worker/check has a configurable timeout (`--timeout`, default 1800 seconds); three failed attempts require attention instead of retrying indefinitely.

## Inspect and operate

```sh
python3 standalone/loop.py --root /absolute/path/to/research-run status
python3 standalone/loop.py --root /absolute/path/to/research-run events
python3 standalone/loop.py --root /absolute/path/to/research-run audit
python3 standalone/loop.py --root /absolute/path/to/research-run export
python3 standalone/loop.py --root /absolute/path/to/research-run stop
python3 standalone/loop.py --root /absolute/path/to/research-run reopen TASK_ID --reason 'New matched measurements released at SOURCE'
python3 standalone/loop.py --root /absolute/path/to/research-run assess --file investigation-assessments.json
```

The live console reports each completed episode. `research.sqlite` is the state owner. Named events record task selection, worker starts/finishes, verification, rejection, recovery, errors, queue expansion, waiting and target completion. The event sequence provides a resume cursor for another consumer. Short SQLite transactions and WAL separate live status reads from long worker processes; no database transaction remains open during inference.

The runtime directory contains:

- `jobs/`: one directory per attempt, with evidence, code, results, submission, verification logs and review.
- `scouts/`: frontier-expansion proposals and their separate reviews.
- `research.sqlite`: durable tasks, contributions and events.
- `public/ledger.json` and `public/README.md`: readable status, evidence references, gaps and review history, without raw worker transcripts.

Exports are an inspectable local read model. They do not deploy a website, send messages, open upstream issues or imply that a partner has agreed to help. Before publishing third-party snapshots, check their licenses; downloading evidence does not give us permission to relicense it. The original source URL, retrieval/read scope and local snapshot hash remain part of the evidence record.

## What counts toward 100

One distinct, useful investigation counts once after all of these gates:

1. The question and proposed increment are explicit and differ from previously published work, including our own C01–C10.
2. Original source evidence is retained with provenance and a clear account of what was actually read.
3. A fresh executable check recomputes the decisive result and fails if that result is wrong.
4. A separate reviewer checks distinctness, source support, reproducibility, usefulness, honest scope and prior art. Every check must pass.
5. Evidence hashes are unchanged through review and are frozen in the acceptance record.

This is **computational verification plus agent review**. It does not prove worldwide novelty, establish biological validation, count as independent experimental replication or confer human acceptance. The human-accepted count stays zero; this minimal runner has no mechanism that awards it.

A well-supported correction or useful negative result can count. Parameter rows, folds, downloads, routine reproductions, generated tasks, generic literature summaries, rejected candidates, infrastructure tests and data-gap cards cannot count. Exact duplicate questions and claims are blocked by SQLite uniqueness constraints; semantic duplication is an explicit reviewer responsibility. There is no perfect automatic novelty detector.

`verified_records` counts accepted computational receipts. `verified_contributions` counts distinct investigations with explicit counting assessments; this is the counter used for the target. A reviewer must assign a stable investigation identity and explain whether a receipt answers an independent question, validates an existing investigation, supersedes earlier evidence or is excluded from the target. Controls, solver repairs, convergence checks and scoring/tuning variants normally remain in their parent investigation. Their receipts and reviews remain visible.

Older receipts are preserved on upgrade and remain uncounted until assessed. Stop the runner before `assess` or manual `reopen`, then restart it with the same run command. Assessment files contain a list of objects with `task_id`, investigation `id`, `question`, `relation` and `rationale`. Corrections append new assessments rather than overwriting old ones. The [first cross-record audit](../research/standalone-2026-11-27/README.md) groups 30 receipts into seven investigations; numerical acceptance alone had overcounted the work.

The audit command checks hashes of accepted evidence. A hash audit establishes integrity, not scientific correctness. Read the scripts, review the assumptions and rerun the recorded commands when reviewing a finding independently.

## Data gates, experiments and suitable collaborators

A data gate is a separate outcome, never a contribution-count shortcut. Its public record must specify:

- the question, existing evidence and why it does not settle the alternatives;
- exactly which measurements, sample identities, metadata or access are missing;
- outcomes that support or refute the claim, and outcomes that remain inconclusive;
- a minimum discriminating design, independent experimental unit, controls and the basis for precision or sample size;
- computational/search work still possible, and the event that would reopen the question;
- potential producers, documented capabilities, strategic fit, source URLs and confidence. An empty candidate list is better than an invented match.

An inference about a company's interests is labelled as such. It is not a partnership, an outreach record or proof the data already exist. No contact or wet-lab execution follows automatically.

Technical failures, quota limits and unavailable compute appear as pending retries or `needs_attention`, not scientific data gates. Dependencies remain explicit. Waiting means the present frontier has no executable task; it is not a claim that all science is solved.

A reviewer can reject a submission while the underlying question remains open. Both rejection and revision return the task for further work with the objection attached. After three unsuccessful reviews it requires attention; it is not relabelled as resolved or data-gated. A question closes only through an accepted contribution, an accepted evidence gap or an accepted finding that further work is not useful. Operator interruptions do not consume the technical-failure budget.

## Other inference providers

The worker interface is a command, not a proprietary coordinator integration. It receives its prompt on stdin and must write the final JSON object to `{output}`. It runs with the job directory as cwd. Configure another agent CLI, a local-model tool agent, or your own API adapter:

```sh
python3 standalone/loop.py --root /absolute/path/to/research-run run \
  --command-json '["/absolute/path/to/your-agent", "--output", "{output}"]' \
  --reviewer-command-json '["/absolute/path/to/other-agent", "--output", "{output}"]' \
  --target 100 --watch
```

The adapter needs actual file/research tools to perform investigations; a text-completion endpoint alone is not a research agent. `protocol.py` contains the contracts supplied to both agents. A different reviewer provider is supported but not required. Separate processes using the same model do not imply independent scientific judgment.

## Execution boundaries

Workers are instructed to write only inside their job directory and use public research as read-only context. The default producer uses the CLI's workspace-write sandbox and the reviewer its read-only sandbox. Fresh verification scripts execute as ordinary local processes, **not in a security sandbox**. Only run this coordinator with trusted agents and research sources; external text must be treated as evidence, never instructions. Arbitrary adapter commands also have the permissions of the invoking user.

No automatic Git pushes, external outreach, destructive actions or laboratory experiments are part of this program. Public source and evidence review remain normal repository workflows.
