# Equal-group inner validation preserves the main KC rank–error trade-off

Equal-group rank tuning changed 8 of 35 fold-by-model settings on the inherited 197-row, 155-group panel. Ordinary trimer ridge in both eligible pools and component-purged identity KNN kept every setting and prediction. Component KNN still led pooled rank among the nested arms, with Spearman 0.46924, while direct ridge had lower log10 KC error. Changes were concentrated in full-direct alternatives. Group-balanced ridge improved slightly on both metrics. The direct median-group KNN lost 0.04326 in rank and gained 0.01212 in RMSE.

| Eligible pool and fitted model | Choice changes out of 5 | Outer rank, row → group tuned | Outer RMSE, row → group tuned |
|---|---:|---:|---:|
| Direct, ordinary trimer ridge | 0 | 0.39977 → 0.39977 | 0.24630 → 0.24630 |
| Direct, group-balanced trimer ridge | 2 | 0.39204 → 0.40298 | 0.24640 → 0.24451 |
| Direct, identity KNN with row votes | 1 | 0.24042 → 0.25785 | 0.29473 → 0.29371 |
| Direct, nearest group member KNN | 2 | 0.31690 → 0.30447 | 0.27326 → 0.28176 |
| Direct, median group label KNN | 3 | 0.31985 → 0.27659 | 0.27314 → 0.28526 |
| Component, ordinary trimer ridge | 0 | 0.36737 → 0.36737 | 0.25117 → 0.25117 |
| Component, identity KNN with row votes | 0 | 0.46924 → 0.46924 | 0.26776 → 0.26776 |

Both outer columns use equal group weight in both arms. On the common direct pool, ordinary ridge’s rank lead over the best nested capped KNN widened from 0.07992 to 0.09530. The preferred capped variant within this panel also switched from median labels to nearest-member labels on both metrics. Balanced ridge moved just ahead of ordinary ridge among nested direct models. The difference was 0.00321 in rank and 0.00179 in log10 RMSE under aligned tuning.

## Strong comparators

Retrospective fixed-grid references keep the tuning result in perspective. Component tie-inclusive KNN at fixed k=3 reached 0.23139 RMSE and 0.39422 rank. The earlier strict-count k=1 component comparator reached 0.48602 rank and 0.27712 RMSE. Direct balanced ridge at fixed alpha 0.01 reached 0.40611 rank and 0.24409 RMSE, both slightly better than its group-rank-tuned version. These settings are descriptive because the same outer outcomes identify their strengths. The frozen row-tuned embedding ridge references scored below ordinary trimer ridge on both metrics within their respective pools. The predecessor's size-matched draws addressed training-pool size, a separate comparison whose settings and sampled pools were outside this replay.

The paired comparison therefore supports keeping the rank–error distinction when reporting the nested models. Metric alignment alone provided no gain for the leading rank model, and targeting equal-group rank gave no assurance of lower error. If error is the intended endpoint, inner RMSE selection is a separate prospective rule worth comparing with the fixed k=3 diagnostic.

## Design and source checks

The five outer tests and four inner validation splits per outer pool come from the earlier benchmark audits. Direct training excludes any same-genus or at-least-95%-identity held pair. Component purging excludes entire connected components of the exclusion graph. Direct outer pool sizes are 62, 48, 42, 43 and 49. Component sizes are 31, 37, 27, 28 and 26. The equal-weight scoring units are the inherited exact-sequence/species groups, distinct from the wider component-purge blocks. The cache contains the inherited mapping correction, identities, and splits.

Each setting produces pooled inner out-of-fold predictions. Row tuning maximizes ordinary Spearman. Group tuning uses inverse validation group size in both marginal midranks and their Pearson correlation across the pooled inner predictions. Equal group RMSE takes the square root of the mean of within-group mean squared errors. Ridge uses cosine-normalized trimer counts with alpha in {0.001, 0.01, 0.1, 1, 10, 100}. KNN uses k in {1, 3, 5, 10}, including cutoff ties. Training weighting, donor rules, precision and ascending-grid tie handling are inherited. Ordinary ridge weights each training row equally, while its balanced variant gives equal weight to each eligible training group. The original and new selections share the exact same candidate predictions.

Young et al. supply the historical Table 3 benchmark and Version 3 tuner. In the deposited code, `inner_oof_spearman` uses unweighted pooled ranks even with a grouped splitter. Table 3 reports 0.605 mean Spearman for 197 KC rows. That published metric applies to its own modeling and split design. Iqbal et al. already use Rubisco genus holdouts and edit-distance KNN. Tu et al. describe the equal-cluster total rank estimator. The MMseqs2 maintainer discussion illustrates why sequence clustering details matter and supplies no evidence about the effect of validation weighting. The cached upstream issue/PR index covers release and benchmark changes. The closest precursor in the coordinator ledger applied equal-group scoring to frozen outer predictions. Earlier ledger entries weighted donor votes or the ridge fit loss. This audit varies the inner selection objective while keeping outer scoring identical across the pair. [Young et al.](https://www.mdpi.com/3042-6723/1/2/7) · [Version 3](https://doi.org/10.5281/zenodo.19546367) · [Iqbal et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) · [Tu et al.](https://doi.org/10.1002/sim.10326) · [MMseqs2 discussion](https://github.com/soedinglab/MMseqs2/issues/398) · [Upstream issues and PR index](https://api.github.com/repos/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/issues?state=all)

## Scope and reproduction

The result concerns one historical panel and one set of folds. Pooled ranks compare predictions across folds. Sampling uncertainty, different splits, equal-group tuning of the embedding models or size-matched draws, measurement heterogeneity and external transfer were outside scope. This is a benchmark audit without new biochemical measurements.

Run `python3 -B analysis.py` to regenerate the outputs, then `python3 -B verify.py`. The verifier independently refits via Gaussian elimination and a separate KNN routine, selects settings via ordinary tied ranks of integer-replicated observations, checks deposited targets and inherited eligibility, and matches the original cached row-tuned predictions. It asserts every saved inner score, choice, paired prediction and outer metric. The source parser and trimer Gram construction are shared. [Input provenance](sources/provenance.md)
