# Protease contribution audit — 27 November 2026

The two accepted protease records support **one conditional benchmark investigation**, not two independent scientific contributions. The second record is an important robustness check for the first: replacing same-study reporter-derived rules with external qualitative context changes the conclusion about when restricted fitting is necessary. Neither demonstrates transfer beyond the frozen ProtRec panel. The stalled inventory has a useful, independently checkable PERRC coverage result that can be resubmitted as a bounded technical contribution; it does not establish the broader data gap.

This audit used the research-audit skill, its case study, the original 26 November controlling reassessment, all three current task outcomes/reviews, accepted reports/code/results/source excerpts, and the latest inventory report/code/workbook. SQLite was opened with `mode=ro`. No runtime, job or public-source files were changed. The independent check imports no producer module and writes only this directory.

## What survives and how to count it

| Task | Assessment | Counting recommendation |
|---|---|---|
| protease-unseen-context-baseline | Supported conditional computational finding; cannot substantiate completely unseen reporter context | 1 investigation in `protease-held-position-calibration` |
| protease-external-position-context-sensitivity | Supported sensitivity analysis with independently sourced qualitative context; necessary control for the same investigation | 0 additional contributions; attach to the family above |
| protease-reporter-calibration-independent-benchmark-inventory | Useful release-specific eligibility audit; unresolved submission classification and incomplete global inventory | 0 accepted now; potentially 1 distinct **technical audit** after an ordinary new review, if useful technical eligibility audits meet the agreed target criteria |

There is no universal claim that a useful null or data audit cannot count. The decision here is that the two accepted calibration rows share the same question, data, outcome, folds and deployment decision. Altering context provenance and fitting weights resolves an objection to that one question. Counting a follow-up task as another advance would split necessary controls. The PERRC audit has a different experimental object and decision (whether this external release can support validation), so it can form a second family, but is not an external validation result or a biological advance. The quota must not determine its classification.

## Decisive independent checks

`python3 -B independent_check.py` independently parses the common-key CSV, reconstructs the seven held-position groups and directly fits pooled means and OLS using ordinary sums. It retains all 780 noncanonical pairs and 170 negative estimates. With the same hard context applied in both arms:

| Fitting population | Context mean RMSE | Linear reporter RMSE |
|---|---:|---:|
| All training rows | 0.19982924 | 0.21229410 |
| Rows left open by rules | 0.17943022 | 0.16062843 |

This reproduces the decisive reversal by separate code. The saved producer verification covers nested selection and isotonic fits; this audit does not pretend to have independently refitted those additional models. The narrower surviving interpretation is that fitting and evaluating the same deployment subpopulation matters for this hard-rule benchmark.

A separate diagnostic takes the external-context saved outer scores and replaces P5 prediction by 1 in **all** arms. P5 RMSE is independently recomputed as 0.47969464 from the CSV. The matched fixed-regime mean errors become 0.21804447 (context), 0.19825583 (linear), 0.20135565 (isotonic). Thus the reported reviewer objection is real: a simple P5 comparator beats the original fixed reporter scores, but applying it fairly to all arms leaves a smaller reporter increment. This is exploratory arithmetic, not a new benchmark contribution or a prespecified model-selection result. Other external folds were taken from saved scores, not refitted in this audit.

The independent workbook parser finds exactly `TEVEp` and `TEVESNp`, each on `ENLYFQS` and `ENLYFES`. Both noncanonical comparisons vary P1. There is one possible nonempty position group and zero folds with any remaining noncanonical training position. Missing paired reporter joins or peptide-fluorescence calibration cannot rescue that coverage. The 72-trace/70-archived-slope reconstruction is covered by the existing separate verifier; this check deliberately targets the decisive positional impossibility rather than repeating every trace fit.

## Provenance, leakage and significance

The hard rules use the source paper's reporter profiles at held positions. The outer biochemical target is held out, and test-pair AUC is explicitly permitted as a predictor, so this is not direct held-target leakage. It is a broader information-set problem: those folds cannot establish generalization to positions absent from all prior reporter profiling. Calling them simply “genuinely unseen positions” conceals this condition. The manuscript reports it; future headlines and task descriptions must carry the qualification.

The second analysis removes that particular same-panel context source, but not all analyst adaptivity. Its numbers 0.8/0.2/0.2/0.1 are analyst anchors, not independently measured rates. Context origin and rule form change together. Nested training-only hyperparameter choice protects against direct outer-target fitting; it does not retrospectively make a sequence of designs chosen after seeing this panel a prospective experiment. Shared denominators, six related enzymes, seven overlapping training sets, and frozen upstream AUC correction remain limitations. No practical gain margin, external validation, adoption or human review has been established.

A fresh primary-source inspection supports qualitative priors and directly challenges exact hard rules. The [2010 structural comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC3005794/) describes solvent-facing P5 and strong P4 dependence; it does not provide the second analysis's four numerical anchors. The [2024 recombinant-protein study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11129917/) documents exceptions and context effects, including P1 substitutions, opposing universal hard-zero claims. On 27 November the browser exposed full primary text for those two sources; relevant results/discussion were inspected. A fresh [PERRC PMC request](https://pmc.ncbi.nlm.nih.gov/articles/PMC12301959/) returned a browser interstitial, so the inventory's cached XML and publisher workbook were used; this is recorded as an access failure, not negative evidence. The cached specialist ResearchGate excerpt concerns tag accessibility and does not establish novelty or benchmark validation. No exhaustive novelty search is claimed.

## Smallest concrete unblock

The third inventory review already states that the verified PERRC finding can stand as a bounded computational contribution while the broader data gate remains open. All six ordinary contribution checks were true; the reviewer rejected `data_gap`, not the workbook result. The smallest next change is therefore a **new contribution submission with corrected scope**, reusing the verified evidence, not weakening the data gate or restarting a broad search.

Exact resumed-task instructions are in `resume_instructions.md`. Required revisions are `kind: contribution`, an exact release-specific claim and decision, the contribution-schema fields, and a clear statement that no independent held-position validation was run. Preserve unresolved Yi supplement access and endpoint equivalence as limitations, not resolved gates. Preserve all prior artifacts. Rerun the existing verifier before submission. Let the separate reviewer decide acceptance under unchanged checks; do not auto-convert a rejected gate into acceptance.

If the broader inventory question is pursued separately, achievable remaining work includes a bounded alternative primary-host retrieval of Yi's supplement and a precise enzyme–substrate/construct join for PERRC reporter and FRET releases. An endpoint-specific P1 descriptive replay can be done now, with uncertainty and concentration dependence; it cannot become multi-position validation. A quantitative fluorescence/product bridge must be demonstrated, not assumed from canonical normalization. These steps are unnecessary to establish the narrower release-specific positional ceiling and must not be counted as independent discoveries merely for completing the inventory.
