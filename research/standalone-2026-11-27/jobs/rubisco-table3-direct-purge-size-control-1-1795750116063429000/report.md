# Rubisco Table 3 rank loss survives direct 95% exclusion

Enforcing the same pairwise boundary without transitive purging kept more training rows in all ten endpoint-fold combinations. ESM2 specificity rank recovered modestly, from 0.129 to 0.207, while KC rank fell from 0.239 to 0.062. Both were well below the exact/species reference. At the component arm's row counts, the most consistent composition effect was on the KC identity-neighbor baseline. Its component fit ranked above all 64 thinner direct-pool fits and also above the full direct fit. Equal row counts therefore leave an important difference between the training pools.

## Results

All rows are pooled out-of-fold Spearman / log10 RMSE. For matched draws the entries are medians calculated separately for each metric across identical sets of 64 draws.

| Endpoint | Training rule | Rows by outer fold | ESM2 ridge | Trimer ridge | Identity KNN | Mean RMSE |
|---|---|---|---|---|---|---:|
| Specificity | Exact/species reference | 146, 146, 146, 147, 147 | .720 / .128 | .734 / .122 | .787 / .110 | .173 |
| Specificity | Components | 33, 33, 31, 28, 31 | .129 / .219 | .433 / .123 | .398 / .130 | .212 |
| Specificity | Direct | 56, 52, 38, 61, 59 | .207 / .195 | .482 / .124 | .270 / .137 | .190 |
| Specificity | Matched direct draws | 33, 33, 31, 28, 31 | .187 / .205 | .401 / .142 | .277 / .144 | .193 |
| KC | Exact/species reference | 157, 157, 158, 158, 158 | .503 / .269 | .615 / .218 | .492 / .286 | .322 |
| KC | Components | 31, 37, 27, 28, 26 | .239 / .401 | .336 / .265 | .467 / .268 | .445 |
| KC | Direct | 62, 48, 42, 43, 49 | .062 / .368 | .379 / .257 | .232 / .298 | .381 |
| KC | Matched direct draws | 31, 37, 27, 28, 26 | .174 / .367 | .298 / .274 | .216 / .320 | .385 |

Full direct fits recovered 7–33 rows per fold. Maximum cross-fold test-to-training identity stayed below 95% for both endpoints. Median nearest-training identity increased from 91.89% to 92.99% for specificity and from 92.07% to 93.20% for KC. A neighbor model can therefore see closer allowed pairs when bridge-mediated exclusions are removed.

Matched-draw ESM2 pooled-rank 5th–95th percentiles were .057–.345 for specificity and .048–.357 for KC. The component results fall within both ranges. The best matched-draw ESM2 ranks were .392 and .382, respectively, below the exact/species reference for each endpoint. The specificity gain looks smaller after ranking within folds. Component and full direct mean within-fold ESM2 ranks were .217 and .256, and the matched-draw median was .212. For KC the corresponding ranks were .207, .121 and .192. Direct KC ESM2 RMSE improved even as rank fell, so the direction depends on the evaluation metric.

KC identity KNN reached .467 pooled rank with component purging, compared with .232 with all eligible direct rows and a matched-draw 5th–95th percentile interval of .011–.342. Its component RMSE was also lower than all matched draws. Mean within-fold rank favored the component fit over 63 of the 64 draws. Trimer ridge led both pooled rank and RMSE for both full direct fits. These results preserve the conclusion that the embedding model needs strong sequence baselines. Choose the purge to match the evaluation setting and report its training coverage. The component fit alone cannot estimate what happens under the pairwise rule.

## Design and scope

The 183 positive specificity and 197 positive KC rows follow the archived Table 3 label order, with the earlier CGA009 repair and exact/species grouped test assignments. The exact/species and component arms reproduce the previous audit. Component memberships are computed on the full prior 214-row sequence graph with exact/species, literal first-token genus and global identity edges at or above 95%. The direct arm removes a candidate when any held row has the same genus or reaches 95% identity. Exact/species groups remain disjoint by design. Given the fixed test rows, the direct rule retains every eligible candidate. Global alignment uses the first optimal alignment, the prior additive scoring scheme and all alignment columns including gaps.

Each draw samples uniformly without replacement from the direct pool, separately by fold, at the corresponding component arm's outer row count. All four training rules use four inner exact/species group folds, and each arm applies its own exclusion policy against inner validation rows. Every draw is retuned. Scaling, PCA and model fitting use only retained fit rows. Grids and baselines match the preceding audit. The inner training counts vary, reaching as low as eight in thinned samples. The matched contrast captures the consequences of changing eligible composition and the resulting retuning at fixed outer counts. The full-versus-thinned contrast probes reducing training size within the direct policy. Differences from the exact/species reference still mix training size and boundary changes.

The historical panel uses row weights and heterogeneous assays. Genus tokens and global identities with terminal gaps are operational proxies. Results concern one fixed partition. Draw quantiles describe subsampling variability on these folds and provide no confidence interval for future taxa or assays. The contrast does not locate a causal clade or show that an allowed closer neighbor is inherently misleading. A follow-up can freeze neighbor count across arms and trace per-row changes to selected neighbors to distinguish label transfer from retuning.

## Prior work and provenance

Petti and Eddy already describe single-linkage bridges and construct alternative pairwise-separated benchmark sets. GraphPart moves or removes records to create threshold-compliant partitions and explicitly checks pairwise identities. The MMseqs2 practitioner thread #398 reports a short-fragment clustering issue and underscores the importance of checking pairwise separation directly. It concerns a different alignment pipeline. The new result here is the fixed-test, matched-row-count benchmark contrast and its endpoint- and model-dependent outcome. The ledger credits the grouped replay and component audit separately. The archived tuner and S1 input were checked in the cached Version 3 release. [Petti and Eddy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8929697/) · [GraphPart](https://pmc.ncbi.nlm.nih.gov/articles/PMC10578201/) · [MMseqs2 discussion](https://github.com/soedinglab/MMseqs2/issues/398) · [Young et al.](https://www.mdpi.com/3042-6723/1/2/7) · [Version 3](https://doi.org/10.5281/zenodo.19546367) · [Flamholz et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6686151/)

The release archive retains its MIT code license, CC BY 4.0 original data and figure license, and third-party notices. The benchmark measurements come from Flamholz et al., the repair sequence from UniProt, and cached residue-mean features from the prior replay of FAIR ESM2-t12. The cache contains features and checkpoint digests. Verification starts from those arrays. No model weights are bundled.

## Reproduction

Run `sh verify.sh` with the packages in `requirements.txt` installed. The script uses a local `.venv` when present or accepts a Python interpreter through `PYTHON`. It checks source hashes, recomputes identities from alignment coordinates and components by sparse graph traversal, reconstructs folds and draws, refits the models and asserts predictions and aggregate contrasts. The verifier reuses the archived-input parser and the ESM2 ridge/PCA estimator. To regenerate `results.json`, run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B analysis.py` with the same environment.
