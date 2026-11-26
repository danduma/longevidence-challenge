# C10 — Held-out reporter calibration

**Coordinator-reproduced; awaiting human review.** A direct CLI agent produced the grouped calibration and normalization audit. The coordinator inspected split/normalization logic, regenerated all four output tables and independently reconstructed the 13 linear fits using centered closed-form OLS. Signed negative estimates remain included.

Install requirements.txt and run `python3 -B verify.py`. All licensed inputs are bundled; no network or model access is required. See [report](report.md), [source and data-license attribution](sources.md), and [coordinator checks](coordinator-checks.json). Original requested tool: gpt-6.2-astra, medium effort.

The finite held-group error advantage does not establish absolute turnover, whole-proteome performance or reliable predictions at unseen folded sites. Normalization and shared assay context limit generalization. The report/contribution retain the historical agent-only status; publication.json records coordinator reproduction. Only the interpreter path is made portable; input/protocol hashes and scientific outputs are unchanged.
