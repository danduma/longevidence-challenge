# Rubisco Table 3 rank transfer under genus and sequence exclusion

Removing training rows related to each fixed Table 3 test fold exposed a limit of the repaired embedding model. For specificity (S), pooled ESM2 Spearman fell from 0.720 with exact-sequence/likely-species groups to 0.527 with genus exclusion and to 0.129 with joint genus and 95% global identity component exclusion. For the CO₂ Michaelis constant (KC), it went from 0.503 to 0.390 and 0.239, with a rebound to 0.466 under the intervening 99% rule. The strictest holdout retained enough rows to fit every fold, but only 26–37 training rows remained per fold. Across the three new arms, trimer ridge beat ESM2 on pooled rank and log10 RMSE for both endpoints. Under the 95% rule, identity KNN achieved the strongest KC rank, while trimer ridge retained slightly lower KC RMSE.

All comparisons use the same 183 positive S and 197 positive KC rows, including historical duplicate weights. Values in model columns are pooled Spearman / log10 RMSE.

| Endpoint | Exclusion | Outer train rows by fold | ESM2 | Trimer ridge | Identity KNN | Mean RMSE |
|---|---|---|---|---|---|---:|
| S | Exact/species | 146, 146, 146, 147, 147 | .720 / .128 | .734 / .122 | .787 / .110 | .173 |
| S | Genus | 95, 90, 88, 114, 98 | .527 / .134 | .678 / .118 | .658 / .121 | .175 |
| S | Genus + 99% | 80, 79, 74, 107, 89 | .512 / .150 | .664 / .120 | .534 / .116 | .178 |
| S | Genus + 95% | 33, 33, 31, 28, 31 | .129 / .219 | .433 / .123 | .398 / .130 | .212 |
| KC | Exact/species | 157, 157, 158, 158, 158 | .503 / .269 | .615 / .218 | .492 / .286 | .322 |
| KC | Genus | 108, 89, 101, 113, 109 | .390 / .339 | .495 / .233 | .252 / .301 | .333 |
| KC | Genus + 99% | 99, 80, 93, 95, 101 | .466 / .269 | .575 / .227 | .396 / .255 | .336 |
| KC | Genus + 95% | 31, 37, 27, 28, 26 | .239 / .401 | .336 / .265 | .467 / .268 | .444 |

Rank aggregation does not explain the strictest arm's model ordering. Mean within-fold Spearman under that rule was .217 ESM2, .506 trimer, and .356 identity KNN for S, and .207, .399, and .475 for KC. The minimum inner training count was 17 for each endpoint. The training mean has tied scores within each fold.

## What the training boundary removed

The earlier exact/species split left 105 S and 108 KC test rows with a literal genus represented in training. The genus purge removed all such exposure at both CV levels, but the maximum test-to-training identity remained 99.784% on both endpoints. The joint 99% and 95% rules reduced that maximum to 98.952% and 94.802%. Median nearest-training identities for S were 98.90%, 97.27%, 96.54%, and 91.89% across the four arms. For KC they were 98.45%, 96.54%, 95.90%, and 92.07%. The largest component across 214 benchmark input rows grew from 7 to 19, 20 and finally 166. At 95%, the model comparisons describe transfer with sharply reduced training coverage. The rank changes combine reduced sample size, changed training composition and separation. A test of their separate effects would hold training counts fixed and compare transitive exclusion with a direct-pair exclusion.

## Prior evidence and scope

The preceding grouped replay isolated exact sequences and likely conspecific rows and tested the same trimer model. The public reassessment and ledger already credit the CGA009 pairing repair and recovery of the historical ESM2 path. The deposit's tuner supports grouped CV and uses training-only preprocessing. The Table 3 run used row-wise folds. The prior replay recovered the released predictions with ESM2-t12 residue means. This job uses those repaired t12 features. This analysis extends the grouped replay using fixed test assignments and nested relatedness purges rather than claiming a new general validation method. [Young et al.](https://www.mdpi.com/3042-6723/1/2/7) · [Version 3 deposit](https://doi.org/10.5281/zenodo.19546367)

Iqbal et al. used genus holdouts for plant Form IB kinetics and compared against Levenshtein KNN under leave-one-out evaluation. De Pins et al. tested sequence-identity restrictions for carboxylation rate. Petti and Eddy described single-linkage bridging and alternatives that remove bridge sequences. The repository issue/PR snapshot includes a revision for models and balancing. The FAIR and Biohub ESM practitioner threads examined pooling and downstream representations. None of those inspected sources supplies the same test rows, paired metrics and training counts for this Table 3 comparison. [Iqbal et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) · [De Pins et al.](https://doi.org/10.1073/pnas.2501433122) · [Petti and Eddy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8929697/) · [Revision PR](https://github.com/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/pull/3) · [Biohub ESM discussion](https://github.com/Biohub/esm/discussions/243)

Literal first-token genera can miss synonyms and group ambiguous names. Global identity includes terminal gaps, so short fragments can have lower apparent identities. A 95% ceiling still permits very similar training candidates and does not establish deep phylogenetic transfer. Pooling historical assays leaves conditions and provenance heterogeneous. Results are descriptive for one partition on a historical panel. The KNN and trimer wins concern benchmark ranking, with no external validation on new assays. Trimer ridge remains the more consistent rank comparator across the two traits, and identity KNN is necessary for a fair KC comparison under the strictest tested boundary.

## Reproduction

With the pinned dependencies in requirements.txt, run `sh verify.sh`. Run `python3 -B analysis.py` with single-thread BLAS to regenerate results.json. The verifier reads the local release and feature snapshots, recomputes identity from alignment coordinates, reconstructs components with a sparse-graph implementation, independently constructs trimer and KNN fits, and checks the nested choices, row coverage and metrics. ESM2 ridge/PCA fitting reuses the prior replay estimator. Embedding inference itself begins from cached feature arrays. The archive preserves the original code MIT license, CC BY 4.0 data license and third-party notices. Historical metadata, UniProt record and FAIR ESM features retain their source attributions in submission.json and the cached license and provenance files.
