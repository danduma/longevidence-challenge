# Rubisco control benchmark

Read report.md for scientific scope, new results, prior art and unresolved anchor accesses. Run `./run.sh` for a local, network-free numerical reproduction with the specified shared Python.

- `protocol.md`: frozen historical-control analysis (SHA256 in verification).
- `protocol_2025_calibration.md`: extension frozen before overlap values were scored.
- `notes/amendment_01.md`, `notes/primary_reuse_audit.md`: explicit post hoc choices.
- `derived/measurement_ledger.csv`: 623 endpoint rows, including 12 from the 2026 table, with origin, conditions, error kind, identity warning and validation status.
- `derived/control_ledger.csv`: 62 target/comparator ratios; imported comparators are explicitly distinguished from contemporaneous controls.
- `derived/primary_source_ledger.csv`: corrected primary identifiers and remaining unknowns.
- `derived/assay_frontier_benchmark.csv`: 145 same-row/table trait pairs and within-origin/condition mean dominance; no certified challenge-success labels.
- `derived/library2025_measurement_ledger.csv`: all 469 screened records; numeric, below-threshold and unreported states are distinct.
- `results/control_rectangles.csv`, `source_pair_metrics.csv`, `leave_one_source_out.csv`, `control_sensitivity.csv`: portability benchmark and sensitivities.
- `results/library2025_calibration_pairs.csv`: all 11 author-defined overlaps with original temperature strings, exact identifiers/sequences/references and correction residuals.
- `results/reused2023_2026_audit.csv`: all 12 shared-cell rounding/SD-to-SE checks.
- `raw/shared_inputs/`: exact local input snapshots. Existing acquisition scripts retrieve public sources only and are not called by run.sh.

Units: S is conventional dissolved-concentration CO2/O2 specificity (mol/mol); vC is carboxylation turnover (s^-1 per catalytic site). Source-reported gas-ratio orientation ambiguity is retained. The 2025 corrected coupled-assay rates are at 30°C and cannot enter the 25°C conventional frontier. Missing numeric cells mean unknown. No error column is presumed SD solely from its historical `_SD` name. No rate is inferred from bacterial growth or affinity.

This is a local research artifact. Original papers, figures, supplementary data and code retain their own licenses/rights; no blanket license is assigned to third-party snapshots. A future redistribution should separate factual derived tables from original copyrighted content and respect each source's terms.
