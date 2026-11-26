# C01 — Equal-dose schedule family

**Coordinator-reproduced; human review pending.** The coordinator inspected the physical solver diff and reran all 66 integrations, scalar quadrature checks and report generation. Independent equal-dose, energy, finite-family and narrower-window checks passed.

Run `PYTHON=python3 sh reproduce.sh` after installing `requirements.txt`; a C++17 compiler is required. All runtime inputs are bundled with their licenses. No network or model call is needed. `python3 -B verify.py` checks stored results without reintegration.

The observed wider-window null and narrower-window tradeoff are limited to the declared finite family, grids and constitutive assumptions. They establish no ice, mechanical or recovery outcome. [Report](report.md) · [Coordinator checks](coordinator-checks.json).

Publication changes replace workstation interpreter paths with python3 and refresh the source-provenance hashes accordingly. Original pre-publication hashes are retained in original-worker-provenance.json. The report retains its historical awaiting-verification status; this README and publication.json record the later coordinator check.
