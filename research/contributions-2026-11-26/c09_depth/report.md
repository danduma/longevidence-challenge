# C09 — Corrected protease depth and candidate stability

**Executed 26 November 2026.** Higher final-pair depth requirements substantially change which candidates are available and selected, while leaving scores and ordering on common keys unchanged. Biochemical agreement does not show a supported directional improvement. This is a frozen five-cutoff sensitivity contribution on corrected data, not a new protease, scoring method, independent assay or prospective success estimate. **Agent-produced; awaiting coordinator verification. No new work here is human expert reviewed.**

## Question and scope

[protocol.md](protocol.md) fixed cutoffs 20/50/100/200/500, selection rules and endpoints before C09 scoring. Earlier depth100 and correction results were already known. Inputs are pinned copies of the final corrected pair tables from the accepted reconstruction. The [accepted synthesis](../../FINDINGS.md), corrected reconstruction report, original raw-data review and earlier timecourse null were read as evidence. Their correction/recurrence/kinetics calculations are prior work and are not counted again.

Minimum depth means the minimum of pooled N+F across five time points. The verifier checks that definition and pooled count sums for every input pair. It uses the existing corrected AUC rather than averaging barcode AUCs. This analysis filters **final pairs conditional on the published barcode-depth20 corrections**. It does not rerun barcode correction at higher cutoffs, downsample reads or estimate sequencing uncertainty. AUC measures recorder activity, not catalytic specificity.

Primary candidates are script05 parent/single-mutant enzyme × parent/single-mutant intact standard-amino-acid substrate pairs. All author labels, including terminal235 substitutions, are preserved: **355,242 pairs, 4,440 enzymes, 134 motifs** at depth20. The familiar Fig.2 count 355,097/4,437 additionally excludes three terminal235 labels and their 145 pairs. An incorrect initial assertion borrowing that figure count stopped before scoring; the protocol records its repair without changing the frozen population. Off-panel motifs are excluded from candidate ranking and separately enumerated in [off_panel_coverage.csv](results/off_panel_coverage.csv): the parent/single-enzyme off-panel tail falls from 15,527 pairs/26 motifs at20 to 1,688/25 at500.

## Candidate ordering and coverage

Each substrate's eligible enzymes are ranked directly by corrected AUC. Calls use the inclusive upper decile within that substrate, with linear quantiles and all ties retained. A second endpoint keeps all maximum-AUC enzyme ties per substrate. These are deterministic retrospective selection rules for sensitivity testing, not a claim to reproduce the authors' hit nomination algorithm. Global-decile results are also supplied.

| Minimum depth | Eligible pairs | Enzymes | Within-motif decile calls | Depth20 decile calls still selected | Depth20 maximum calls still selected |
|---:|---:|---:|---:|---:|---:|
| 20 | 355,242 | 4,440 | 35,810 | 35,810 | 135 |
| 50 | 321,457 | 4,439 | 32,454 | 28,328 | 81 |
| 100 | 274,198 | 4,439 | 27,696 | 20,379 | 58 |
| 200 | 206,009 | 4,438 | 20,847 | 12,284 | 39 |
| 500 | 95,963 | 4,432 | 9,748 | 3,946 | 15 |

All 134 motifs remain represented at every cutoff, but **only wt has a complete 134-motif profile**, even at20. Retaining almost all enzyme labels therefore does not establish adequate coverage for specificity profiles or off-target evaluation. See [enzyme_coverage.csv](results/enzyme_coverage.csv).

At500, 27.01% of baseline pairs remain. Of 35,810 baseline within-motif decile calls, **31,864 lose eligibility, none are demoted while eligible, and 3,946 remain selected** (11.02%). The new set includes **5,802 newly selected keys**, with baseline/new Jaccard 0.09483. This is turnover caused by the available population; it does not establish that lost low-depth calls are false positives. For the maximum-AUC endpoint, 120/135 original calls disappear and 121 replacement calls enter (ties give 136 at500).

The 15 original maximum calls surviving every cutoff are: D127W×EHLYFQS, D127W×ENLYFQC, N177F×ENLYFQS, Y11A×ENLYRQS, F217V×ENLYFQK, H46R×ENLYYQS, L55R×ENLYFYS, M124H×ENLYFSS, M82L×ENLYEQS, M82P×ENLYFFS, N177Y×ENLYVQS, Q145M×ENSYFQS, V216Q×ENLYAQS, V219C×ENLRFQS and W143Q×ENLYHQS. Their pooled minimum depths and scores are in [top_candidates_stable.csv](results/top_candidates_stable.csv); the complete baseline top-call ledger documents every disappearance in [top_candidates20.csv](results/top_candidates20.csv).

On the **95,963 fixed common depth500 keys**, all five cutoff views have maximum AUC difference zero, Spearman rho=1, zero within-motif rank changes, the same 9,748 upper-decile calls and the same 136 maxima when selection is recomputed on common keys. This structural null is expected because the operation filters existing scores. Absolute rank numbers in the changing populations can change as competitors disappear; no pairwise score inversion occurs. It is not biological evidence that low-depth AUC estimates are accurate. [Candidate summary](results/candidate_summary.csv), [fixed-key checks](results/fixed_common.csv), [full eligibility/rank/selection ledger](results/candidate_trajectories.csv.gz).

## Matched biochemical comparisons

The original relative initial-rate file has 840 rows. Collapsing repeated canonical representations gives 804 enzyme×motif keys; these representations are not independent replicates. Signed noisy rate estimates are retained. Primary parent provenance follows script03: replace control parent0 with DMS wt. The other five panels use corrected controls. These six panels are pre-existing reference enzymes, **not an independent biochemical validation of the newly ranked DMS candidate calls above**. The published assay agreement itself is prior art. [Huber et al. (2025), data and code availability](https://www.nature.com/articles/s41467-025-60622-7).

| Enzyme, primary branch | Matched n:20/50/100/200/500 | rho20 → rho500 | Recorder/rate decile overlap20 → 500 | Baseline recorder calls lost at500 |
|---|---|---|---|---:|
| 0 (DMS wt) | 134/134/134/134/134 | .8116 → .8116 | 10/14 → 10/14 | 0 |
| NX1 | 133/132/126/111/68 | .7067 → .7523 | 6/14 → 3/7 | 4 |
| NX2 | 133/132/131/117/59 | .8347 → .8137 | 8/14 → 4/6 | 11 |
| VI | 134/134/132/127/114 | .7575 → .7960 | 8/14 → 8/12 | 0 |
| VI* | 134/133/133/124/94 | .7044 → .6640 | 10/14 → 7/10 | 5 |
| VII | 134/134/130/127/105 | .6648 → .6527 | 9/14 → 7/11 | 3 |

Overlap denominators are the recorder call counts, not sample sizes or prospective precision. Both tails are redefined on the eligible population. Frozen depth20 biochemical-high flags are also retained in the identity ledger. Primary total matched coverage falls **802 → 799 → 786 → 740 → 574**. NX1×ENLYLQS and NX2×ANLYFQS are already missing at20. No absent rates or recorder values are imputed.

DMS-wt's entire comparison is stable. VI loses no baseline selected keys to depth, but two become unselected because its upper-decile threshold changes with population size. NX1 has four depth losses and three such demotions; NX2 has eleven depth losses and no demotions. VI* loses five and VII loses three. At500, the surviving and still-selected baseline recorder calls number 14/7/3/12/9/11 for 0/NX1/NX2/VI/VI*/VII. Thus unchanged AUC order does not guarantee unchanged percentile membership.

**In the primary branch, eleven calls that were both recorder-selected and biochemical-high at20 lose depth eligibility at500:** NX1×EDLYFQS/EELYFQS; NX2×EDLYFQS/EHLYFQS/EPLYFQS/EQLYFQS/ETLYFQS; VI*×EHLYFQS/EPLYFQS/ETLYFQS; VII×ETLYFQS. These calls' biochemical-high labels are frozen at20; their loss is coverage attrition, not a failed biochemical measurement. [Lost-call identities](results/biochemical_high_calls_lost500.csv), [all baseline selected identities](results/biochemical_selected20.csv), [complete matched ledger](results/biochemical_trajectories.csv).

Across six panels, mean correlation change500−20 is **+0.001752**, with paired shared-motif bootstrap 95% percentile interval **[−0.038832, +0.034172]**. Two panels improve, three worsen and the DMS-wt parent is unchanged. The 2,000 draws resample 134 motif identities jointly across all panels, preserving shared-motif multiplicity and each panel's eligibility. The interval measures conditional sensitivity to observed motifs at fixed measurements; it omits biochemical/cultivation replicate error and independent-enzyme sampling. It supports neither a directional improvement nor equivalence. On each panel's fixed500 keys, all five views give exactly the same correlations and overlap.

Parent provenance matters. With the manually cloned parent, matched coverage is **132/132/128/121/70**, rho changes **.7111 → .8119**, overlap **9/14 → 4/7**, four baseline selected keys lose eligibility and three are demoted. Its missing-at20 parent keys are ENLVFQS and RNLYFQS. The other five panels are unchanged. The six-panel mean change in this branch is +.018550, interval **[−.027129, +.056565]**. Neither branch supports directional improvement. The parent branch difference must not be attributed to the depth cutoff alone. [All panel/cutoff results](results/biochemical_summary.csv), [dependence-aware uncertainty](results/biochemical_uncertainty.csv).

## Contribution, limits and verification

C09 adds a frozen multi-cutoff **coverage and candidate-identity audit** of final corrected data and a matched, provenance-preserving biochemical sensitivity. It does not rebrand the known depth100 controls, all-depth recurrence, published corrections, previous kinetics null or source reconstruction as new outcomes. The ordering null is algebraic; the informative measured result is the extent and identity of selection attrition, plus its differing effect on selected biochemical comparisons. No cutoff was tuned on biochemical targets, and500 is not recommended as optimal.

Candidates share enzymes, motifs, barcodes/library and experimental context. Counts are a census of this snapshot, not independent replicates. The biochemical references are selected related enzymes on short peptides; endogenous substrate accessibility, cellular cleavage, proteome-wide off-target effects and prospective accuracy remain unsupported. The test cannot separate read-count precision from coverage selection, rerun all barcode-stage decisions at higher cutoffs, or supply missing biochemical replicate uncertainty.

Run **`./verify.sh`** in this directory with the shared Python. It verifies pinned input SHA256, pooled depth/sum, unique keys, exact intended-population subset, fixed-key score/call invariance, the previously known depth100 biochemical values as a sanity check, and complete selection-count reconciliation. All result tables and fixed-seed bootstrap intervals are regenerated offline without models or network access. Environment: Python3.14.10, NumPy2.5.5, pandas3.1.1, SciPy1.18.1. [Verification record](results/verification.json), [input provenance hashes](sources/input_manifest.json), [output hashes](results/output_hashes.json), [public source/download and redistribution instructions](sources/README.md), [source failures and refusal log](refusals.md).
