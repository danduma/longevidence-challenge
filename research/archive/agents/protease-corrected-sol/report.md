# A corrected-population audit of ProtRec specificity evaluation

26 November 2026. **Reconstruction succeeded. The authors’ final single-mutant population is reproduced exactly, including every deposited AUC. The earlier additive-model extrapolation failure becomes substantially weaker after correction and eligibility matching; residual negative-control recurrence persists but does not establish a transferable background correction.** This is a reproducible assay/evaluation contribution, not a new protease mechanism or a solution to the programmable-protease challenge.

The decision-relevant result is that the released mapped input is unsuitable as a substitute for the authors’ final evaluation population. In the original 10,018-row substrate holdout, a single author-excluded backbone associated with ENLCFRG accounts for 369 rows and **54.0% of the ridge model’s excess squared error over the training-mean baseline**. This attribution uses the original fixed fit; it is a post hoc diagnostic, explicitly distinguished from the locked matched-population refits. On corrected common keys, the RMSE difference is +0.00201 dimensionless AUC, with a motif-bootstrap interval spanning zero. It is therefore unwarranted to carry forward the original robust-failure claim as though it had been demonstrated on final data. This does not prove successful extrapolation.

**What was reconstructed and verified.** The primary source is [Huber et al. (2025)](https://doi.org/10.1038/s41467-025-60622-7), with ProtRec pinned to `8b9a1865108fa44c3e0a24ac20a59754222847ca`. The complete Git tree was inspected before acquisition: unmapped DMS was 35,178,112 bytes compressed / 188,797,883 bytes for its data member, and controls were 3,057,164 bytes. The required manual suspicious-barcode list was deposited; no flags were guessed. The source workbook was acquired after a HEAD check (35,182,821 bytes). Repository blobs, input hashes, response metadata and exact local files are retained in sources (archival input; see source-downloads.json and the source URLs). The release/Zenodo metadata were checked; the needed inputs were present, so the missing-data fallback was not invoked. No raw-read download was necessary.

[R was unavailable](access-failures.md). [reconstruct.py](reconstruct.py) ports the `min20DNAlevel` branch actually read by scripts 03 and 05, including the following author-defined transformations:

- Compute barcode-level fractions from F/(F+N), minimum depth across five samples, and trapezoidal time-average over 0–4 h. The AUC is dimensionless, not a catalytic rate.
- Exclude the first 11 non-stop high-median substrate barcodes under the script’s rules, and the deposited manual list minus its nine explicitly confirmed-correct barcodes, plus the separately named barcode. This yields 62 unique backbones. The high-median rank-11 boundary is untied (0.11195 versus 0.11130).
- Require each retained barcode pair to have at least 20 reads at every time point; pool counts, rather than averaging barcode AUCs. Reject entire eligible duplicate groups when max/min AUC ≥3 and at least one is above the relevant noise cutoff.
- For active groups with ≥3 eligible barcode pairs, remove the first maximum only if pooled AUC / pooled AUC without it ≥2, excluding zero denominators. Reaggregate these groups with the script’s rounding stages even when no outlier is removed. Exclude intact-substrate values at or above the specified positive control.
- Fit the published pooled normal noise model to the designated controls. Its mean is 0.0621537313432836, population SD 0.00581247573214537 and 95th-percentile cutoff 0.0717144031328703. This is the authors’ model, including their normality assumption. Published noise blanking is retained separately as `dnAUC`; corrected AUC is the primary endpoint here.

All 413,779 deposited mapped DMS keys and all 3,482 mapped control keys match upstream pooled **counts exactly**. Their stored AUCs agree to ≤0.0001. The final reconstruction has stronger validation:

| Primary deposited check | Result |
|---|---:|
| Fig. S8e single-mutant keys, including stops | 358,058/358,058, no missing or extra keys |
| Every S8e AUC | Exact equality, maximum error 0 |
| Fig. 2 intact-substrate population | 355,097 pairs, 4,437 enzymes, 134 substrates |
| Fig. S7a negative-control distributions | All 402 values exact, separately by construct |
| Fig. S7c, manually cloned parent branch | All 836 nonmissing values reproduced; max numerical error 5.2×10⁻¹⁵ |

The first port differed from S8e by one AUC rounding unit on 1,296 rows. A recorded numerical audit identified R’s `fround`, minute scaling in reaggregation, and extended-precision row summation as the cause. The repaired arithmetic reproduces every deposited value without substituting the source table as an output. [Source checks](results/source_validation.json), [rounding audit](rounding_audit.py), and [protocol amendments](protocol.md) document this post hoc implementation repair.

There is an actual provenance fork: current script03 replaces manually cloned parent `0` with DMS `wt`, but deposited Fig. S7c uses the manually cloned parent. With the code’s replacement, 55 nonmissing source values differ; with the manual branch, all values and missingness match. Both are retained. The primary biochemical comparison follows the locked code-replacement branch; a separate manual-parent sensitivity follows the figure. This is not evidence of a biological discrepancy.

The correction-stage tables retain 400,912 DMS pairs and 3,448 control pairs. They include off-panel motifs; they must not be called the final intended single-substrate analysis. Script05’s intended population restricts both axes to parent/single substitutions. The 25 corrected multi-mutant holdout motifs below are outside that population.

Every upstream pair has a key and an explicit sequential disposition:

| Upstream key disposition | DMS | Controls |
|---|---:|---:|
| All barcodes excluded as author-flagged backbones | 9,931 | 21 |
| No remaining barcode meets depth20 | 121,314 | 244 |
| Discordant duplicate | 731 | 2 |
| Intact substrate at/above positive control | 47 | 0 |
| Retained | 400,912 | 3,448 |
| Total upstream keys | 532,935 | 3,715 |

The barcode ledgers plus separate backbone-row ledgers cover all 1,725,713 DMS and 24,103 control input rows. They retain original row identifiers, both barcode keys, enzyme/substrate keys, relevant depth/AUC and exclusion flags. Final tables and all-key ledgers are in [results](results/); [exclusion_waterfall.csv](results/exclusion_waterfall.csv) and [reconstruction.json](results/reconstruction.json) summarize them.

**Locked holdouts and thresholds.** The original review’s raw scores were known before this job. [protocol.md](protocol.md) was written before reconstructed outcomes were computed. Eligible model rows have depth ≥100, intact standard-amino-acid seven-residue substrates, and parent/single-substitution enzyme labels. Common keys require eligibility in both raw and corrected data. The original seed and raw eligible order determine the random mask before intersection. Each matched model is refitted, using additive enzyme identity and substrate-residue indicators and unit ridge regularization, including the constant after target centering. No model or threshold was optimized on the results.

| Evaluation population | Test rows / motifs | Ridge RMSE | Mean RMSE | Paired RMSE difference, motif-bootstrap 95% interval |
|---|---:|---:|---:|---:|
| Original raw multi-mutant holdout | 10,018 / 27 | .029063 | .021497 | +.007566 [.001572, .012657] |
| Common keys, raw values/refit | 9,614 / 25 | .018722 | .016324 | +.002398 [−.000365, .005271] |
| Common keys, corrected values/refit | 9,614 / 25 | .018397 | .016384 | +.002013 [−.000724, .004939] |

The last two use exactly the same train/test keys. The corrected random-pair fit has RMSE .012989 versus mean .019611. Thus numerical blocked degradation remains, but the evidence for a positive blocked excess is no longer robust to the specified motif resampling. The intervals use 3,000 whole-motif draws, fixed fits and complete enzyme rows per sampled motif. They quantify sensitivity to these 25/27 observed motifs, not population coverage, training variation or biological replication. A new model has not been validated.

The post hoc fixed-fit decomposition assigns 54.0% of original excess squared error to ENLCFRG alone. All 369 of its original eligible rows lose every barcode to the authors’ backbone exclusions. Across all 404 rows failing corrected eligibility, the fraction is 54.7%; the remaining exclusions include nine discordant duplicates, three above-positive rows, 22 below corrected depth100, and one with no eligible barcode. This isolates a concrete technical dependency without claiming that all remaining error is technical or that deletion proves generalization. [Per-motif results](results/blocked_motif_impact.csv) and [fixed-fit decomposition](results/blocked_survival_decomposition.csv) preserve the calculation.

For the locked empirical 99th-percentile reference, raw all-control cutoff .136788 becomes .102604 after correction; eligibility-matched cutoffs are .096580 raw and .078334 corrected. These are empirical recorder thresholds, **not measured biochemical false-positive rates**.

| Common-key blocked endpoint | AP | Positive prevalence |
|---|---:|---:|
| Raw, all negative rows as reference | .11878 | .00135 |
| Raw, eligibility-matched reference | .15060 | .01467 |
| Corrected, all negative rows as reference | .14352 | .00988 |
| Corrected, eligibility-matched reference | .37594 | .09538 |

Both corrected references give AP above their own prevalence. Different cutoffs define different endpoints, so AP changes are not a model improvement estimate. Correcting the earlier AP interpretation is not a new biological discovery. As a descriptive candidate-impact check, inclusive upper-decile selection on 282,646 common candidates changes 1,061 previously selected candidates and adds 1,110; ties make the two set sizes slightly different. This demonstrates a concrete selection consequence, not prospective precision. [All model/threshold metrics](results/models.csv), common keys (archival input; see source-downloads.json and the source URLs), [candidate impact](results/candidate_decile_impact.json).

**How much background recurrence remains?** The locked leave-one-negative-construct-out test averages the other two constructs, predicts the held-out construct’s ranks, and uses inclusive upper deciles on 138 common complete depth100 motifs.

| Held-out negative construct | Raw ρ | Corrected ρ [motif-bootstrap interval] | Raw overlap | Corrected overlap / flagged |
|---|---:|---:|---:|---:|
| (1–13)-mCherry | .396 | .256 [.087, .410] | 7/14 | 3/14 |
| I C151A 30*31* | .422 | .226 [.056, .384] | 8/14 | 5/15 |
| No TEVp | .453 | .210 [.054, .354] | 6/14 | 2/14 |

The held-out high-tail base rate is 14/138. All three still meet the prespecified positive-direction and above-base-rate criteria. The association therefore **persists but attenuates**. Paired 3,000-motif intervals for corrected-minus-raw ρ are respectively [−.258, −.023], [−.323, −.081], [−.385, −.108]. Tail estimates are much less stable: corrected precision bootstrap intervals include zero for mCherry and No TEVp. The three rotations share motifs, constructs, library and cultivation; they are not three independent screens. WT/single-only and all-depth sensitivities retain positive point estimates and enrichment; the all-depth truncated-control ρ interval includes zero. No depth threshold was lowered to rescue cross-screen validation.

For corrected MBP, only 52 depth-qualified common-with-the-three motifs remain: ρ=.147 and one of six calls overlaps six held-out high motifs. Screen3 still supplies only three and one complete depth100 motifs by replicate, below the locked ten-motif minimum. These checks **do not support transport across constructs or screens**. Residual within-screen recurrence cannot distinguish substrate effects from shared barcode/library effects and does not justify routine motif-wise subtraction. [Transfer results](results/transfer.csv), [paired uncertainty](results/transfer_paired_uncertainty.csv), [MBP](results/MBP.csv), [cross-screen coverage](results/cross_screen.csv).

**Orthogonal biochemical evidence.** The six pre-existing relative initial-rate panels are dimensionless rates normalized to canonical substrate, with canonical representations collapsed to one enzyme×motif key. The primary comparison replaces control parent `0` with DMS parent `wt`, as the current code specifies, and requires depth100 in both versions. No multi-mutant holdout motif has a biochemical label here.

| Enzyme | Common motifs | Raw ρ with relative rate | Corrected ρ | Raw → corrected biochemical top-decile hits |
|---|---:|---:|---:|---:|
| 0 | 134 | .7263 | .8116 | 9 → 10 of 14 |
| NX1 | 126 | .7350 | .7426 | 6 → 6 of 13 |
| NX2 | 131 | .8191 | .8383 | 8 → 8 of 14 |
| VI | 132 | .7428 | .7722 | 8 → 8 of 14 |
| VI* | 133 | .6987 | .7068 | 10 → 10 of 14 |
| VII | 130 | .6227 | .6701 | 7 → 8 of 13 |

All six rank correlations improve, with mean change +.0328. A joint motif bootstrap that resamples each motif across all enzyme panels gives [.0065, .0680]. This preserves shared-motif dependence, but omits unprovided cultivation and biochemical replicate uncertainty; the six selected enzymes are not six independent assay-development experiments. Using manually cloned parent instead gives .7050 → .7160 for that parent on 128 motifs and 8 → 9 top-tail hits; the other five panels are unchanged. Thus orthogonal evidence supports modest technical ranking improvements after the **published correction**, not a novel background correction or prospective enzyme-selection improvement. Published blanking and negative-mean subtraction are supplied as secondary analyses; their effects vary across panels and do not establish a new transferable estimator. [Biochemical results](results/biochemical.csv), [joint uncertainty](results/biochemical_uncertainty.json), [manual-parent sensitivity](results/biochemical_manual_parent_sensitivity.csv).

**Contribution and remaining gap.** The bounded additional finding is a source-validated population-dependence of a concrete public-data extrapolation diagnosis, together with quantified attenuation of negative-control transfer and an orthogonal check of candidate-ranking reliability. It changes how this release should be used: benchmark on the corrected keyed population, state whether the intended single-mutant panel or an off-panel tail is being evaluated, and do not use strong raw-control recurrence to infer a broadly transferable subtraction. The corrected tables, exact arithmetic port, source reconciliation and complete exclusion ledgers make those decisions executable and falsifiable.

The original authors already published the recorder artifacts, correction rules, control blanking and biochemical agreement. The preceding independent review already supplied the raw holdout and transfer hypotheses. Our AP correction, port, and additional evaluation quantities do not constitute a novel biological discovery. The prior-art check is bounded to the sources and searches recorded in [prior-art.md](prior-art.md); it does not establish first-ever novelty. The candidate-reliability evidence is retrospective and panel-specific. A prospective advantage, independent cross-screen correction, and calibrated biochemical cleavage probabilities remain unestablished.

The full Millennium problem requires blinded prospective design and cellular, catalytic and off-target performance on endogenous targets. The peptide-only partial milestone is also prospective. Neither is tested here; no new enzyme was designed or validated.

**Reproduction and verification.** Run `./run.sh` from this job directory, using the supplied shared Python by default. Dependencies and executable path are recorded in requirements.txt (archival input; see source-downloads.json and the source URLs) and environment.json (archival input; see source-downloads.json and the source URLs). `acquire.py` verifies pinned Git blobs if re-acquisition is needed; `acquire_source_data.py` checks sizes and fetches primary supplements. All scientific writes remain inside this directory. No shared input was edited, no subagent/model session was launched, and no external contacts, posts or experiments were performed. Input hashes (archival input; see source-downloads.json and the source URLs), output hashes (archival input; see source-downloads.json and the source URLs), final assertions (archival input; see source-downloads.json and the source URLs), [access failures](access-failures.md), and [refusal log](refusals.md) are separate artifacts.
