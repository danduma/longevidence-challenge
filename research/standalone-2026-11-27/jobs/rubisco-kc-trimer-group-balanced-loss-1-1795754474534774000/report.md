# The pooled KC lead survives balancing the trimer ridge loss

On the 197 Table 3 KC rows and the same direct eligibility pools, giving each eligible exact-sequence/species group one unit of training weight slightly lowers trimer ridge's pooled Spearman and leaves its lead over both capped KNN variants intact. A fixed-alpha check shows a small error reduction from weighting. Retuning by inner rank reverses that error reduction, with the selected alpha changing in one outer fold.

| Predictor | Pooled Spearman | Pooled log10 KC RMSE |
|---|---:|---:|
| Unweighted trimer ridge, nested alpha | 0.378783 | 0.257000 |
| Balanced trimer ridge, nested alpha | 0.364963 | 0.257579 |
| Balanced trimer ridge, original alpha | 0.373630 | 0.255990 |
| Nearest eligible row per group KNN, nested k | 0.239145 | 0.279895 |
| Eligible group median KNN, nested k | 0.242032 | 0.279651 |

Against median-group KNN, which has the strongest pooled results of the two capped comparators, the pooled rank lead changes from 0.136752 to 0.122932 and the RMSE advantage changes from 0.022651 to 0.022072 after retuning. Balanced ridge beats median-group KNN by within-fold rank in two of five folds. Keep the ordinary trimer ridge as the primary pooled-rank comparator for this panel. The opposite directions of the small RMSE changes at fixed and retuned alpha give little basis for choosing training weights on prediction error here.

## Paired fit

The input order and labels come from the deposited cross-species KC benchmark. We inherit the audited CGA009 sequence replacement and the predecessor's exact-sequence/species graph, global alignment identity and outer grouping. The seven duplicate FASTA headers have identical records. Each outer fit drops rows sharing genus or at least 95% identity with any outer test row. The outer training counts are 62, 48, 42, 43 and 49 rows from 47, 44, 41, 40 and 36 groups. Split each eligible pool into four group validation folds and apply the same purge again for each inner fit. Group sizes are recounted from the eligible rows separately in every fit. The cached group graph was originally built before filtering to KC-positive rows.

Use the dot product of unit-L2 overlapping amino-acid trimer count vectors. For fit rows, let W be diagonal with weight 1 per row or 1 divided by its current eligible group size. Center training targets at their W-weighted mean m. The dual solves (K + alpha W⁻¹)c = y − m and predicts m + K(test, train)c. The centering convention keeps the unweighted fit identical to the predecessor. The objective uses a weighted **sum** of squared residuals plus a ridge penalty. Fixed numerical alpha therefore keeps the explicit penalty coefficient fixed while total loss weight falls. Tune over 0.001, 0.01, 0.1, 1, 10 and 100 by pooled inner Spearman, taking smaller alpha on ties. Both validation and test scoring give each row equal weight. Predictions are rounded to ten decimals before scoring, as in the predecessor. Unweighted alpha choices are 0.01, 0.001, 0.001, 0.1 and 0.01. Balancing changes only the first to 0.001.

Both KNN comparators use the same eligible rows and splits. Each group takes the identity of its closest member, with archive order breaking representative ties. Its label is either the closest member's log10 label or the group median of the eligible log10 labels. Groups tied at the cutoff all vote. Pooled inner Spearman tunes k over 1, 3, 5 and 10 for each rule, taking smaller k on ties.

## Prior evidence and scope

The ledger's closest predecessor changed KNN votes and reused a frozen unweighted trimer ridge. Young et al. define the original Table 3 target and tuning reference. Their revision tooling balances target quantile bins in separate benchmark experiments. Iqbal et al. used genus holdouts and edit-distance neighbors for plant Rubisco. Redundancy-weighted protein training losses already exist in work on secondary structure. The BioStars thread discusses the dependence of cluster-representative choice on downstream use, and the MMseqs2 maintainer thread discusses cluster membership and pairwise identity. The increment here is the paired test of whether donor groups can explain the stronger comparator's lead through its own training loss. [Young et al.](https://www.mdpi.com/3042-6723/1/2/7) · [Version 3](https://doi.org/10.5281/zenodo.19546367) · [Revision PR](https://github.com/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/pull/3) · [Iqbal et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) · [Redundancy-weighted training](https://academic.oup.com/bioinformatics/article/36/12/3733/5809527) · [BioStars](https://test.biostars.org/p/11987/) · [MMseqs2](https://github.com/soedinglab/MMseqs2/issues/398)

This is a sensitivity result for one historical panel, fixed split, exact/species grouping and row-weighted scoring. No sampling uncertainty was estimated. The pooled and fold-level rank results differ, and equal treatment of held-out groups could give a different interpretation. The analysis does not establish prospective transfer or biochemical effects.

## Reproduce

Run `python3 -B analysis.py` to regenerate `results.json`, then `python3 -B verify.py` using standard-library Python 3.10+. The verifier rebuilds the folds and eligibility from cached groups and identities, checks the archived target order and values, recounts trimer features, uses a separately implemented rescaled Gaussian solve, recomputes both KNN variants, checks every saved prediction and metric, and matches the archived unweighted predecessor vector. It inherits the global alignments and group IDs. The source directory preserves the deposit notices and scoped snapshots.
