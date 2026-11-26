# C02 — Boundary exchange feasibility envelope

**Scientific significance reassessed,26November2026:** conditional benchmark envelope. Fixed-horizon brackets are not measured environmental tolerances. This is not a demonstrated biological advance. Read the [controlling reassessment](../../reassessment-2026-11-26/README.md). Earlier reports below preserve their historical wording.


**Coordinator-reproduced; awaiting human review.** Produced by a direct Codex CLI agent requesting gpt-6.2-astra at medium effort. This is one additional contribution in the ten-contribution batch; it is not human acceptance or biological validation.

At a prescribed cold exterior, the flat schedule's terminal myocardial cold spot fails first; at room temperature, the bath can fail the all-time rate screen while the myocardium still passes. Cold transition brackets shift between voxelizations. See [report.md](report.md), [137 scenario rows](results/envelope.csv), [switch brackets](results/brackets.csv), and [coordinator checks](coordinator-crosscheck.json).

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -B verify.py
```

Requires a C++17 compiler and `shasum`. The verifier uses bundled, licensed numerical inputs, recreates both meshes, freezes the published baselines before new scenarios and compares every output row in a temporary directory. No model access or network request is needed. `--tables-only` checks stored tables and hashes without numerical re-integration.

The coordinator inspected the physical-code diff, reran the complete verifier successfully, independently checked the reported decisions and reran the six-face boundary-energy fixture. The original worker report is retained unchanged, including its historical awaiting-verification statement. This README and contribution.json record the later verification status. No measured boundary coefficient, continuum threshold, ice, stress or organ recovery is established.

The public publication changes only contribution metadata and adds this README/check record; numerical code/inputs/results match the original worker snapshot. `original-worker-manifest.sha256` preserves the pre-publication hashes; `manifest.sha256` covers the public files actually checked. MIT and upstream Apache notices remain in pinned/. Source dates/versions and access/refusal records are included.
