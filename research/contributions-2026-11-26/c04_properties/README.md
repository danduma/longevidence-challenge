# C04 — Constitutive sensitivity

**Coordinator-reproduced; awaiting human review.** The coordinator inspected the source/capacity/enthalpy changes and diagnostics, reran all40 full-field cases and scalar checks, and reproduced every original result file byte for byte. Separate energy, analytic onset and bracket checks passed. Requested tool: gpt-6.2-astra, medium effort.

Install requirements.txt and run `python3 -B verify.py`. Requires a C++17 compiler. Licensed runtime inputs and source attribution are bundled; no network or model access is needed. [Report](report.md) · [Coordinator checks](coordinator-checks.json) · [Sources](sources.md).

The source/capacity multipliers are hypothetical sensitivity probes, not measured uncertainty or probabilities. The nominal conditional pass has a small margin, and the paired post-glass ranking survived only the finite sampled grid. No ice, stress or recovery was evaluated. Original worker hashes are retained; only interpreter paths are made portable.
