# C09 — Depth and candidate stability

**Coordinator-reproduced; awaiting human review.** The coordinator inspected filtering, population definitions, selection and shared-motif bootstrap logic, reran the complete analysis and reproduced all twelve frozen scientific output hashes. An independent direct-input calculation reproduced the depth500 attrition and replacement counts.

Install requirements.txt and run `C09_PYTHON=python3 sh verify.sh`. All consumed processed-data snapshots are bundled with CC BY attribution in [sources/README.md](sources/README.md); no network or model is needed. [Report](report.md) · [Coordinator checks](coordinator-checks.json). Requested tool: gpt-6.1-sol, medium effort.

This filters existing corrected scores; it does not re-estimate read noise or rerun barcode corrections at new depths. The unchanged ordering on fixed keys is structural. The measured selection attrition and reference-panel sensitivities are the contribution. No new biological candidate validation or optimal cutoff is claimed. The report retains its historical agent-only status; publication.json records later coordinator reproduction.
