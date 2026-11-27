# Inner error tuning recovers part of the fixed k=3 KC advantage

On the inherited 197-row, 155-group KC panel, component-purged identity KNN tuned by inner equal-group RMSE improved pooled outer log10 KC RMSE from **0.26776 to 0.25580**. That recovers 32.9% of the RMSE gap between rank tuning and the retrospective fixed k=3 reference (0.23139). Pooled equal-group Spearman fell from 0.46924 to 0.36832. Fixed k=3 scored 0.39422 on rank and remained ahead of error-tuned KNN on both endpoints. Keep k=3 as a retrospective diagnostic.

| Pool and model | Outer RMSE, rank tuned | Outer RMSE, error tuned | Outer rank, rank tuned | Outer rank, error tuned |
|---|---:|---:|---:|---:|
| Direct trimer ridge | 0.24630 | 0.24705 | 0.39977 | 0.41850 |
| Direct balanced trimer ridge | 0.24451 | 0.24604 | 0.40298 | 0.41468 |
| Direct row-voting KNN | 0.29371 | 0.29373 | 0.25785 | 0.25566 |
| Direct nearest group member KNN | 0.28176 | 0.27902 | 0.30447 | 0.28832 |
| Direct median group label KNN | 0.28526 | 0.27899 | 0.27659 | 0.29610 |
| Component trimer ridge | 0.25117 | 0.25110 | 0.36737 | 0.37069 |
| Component row-voting KNN | 0.26776 | 0.25580 | 0.46924 | 0.36832 |

Component KNN under error tuning had higher RMSE than even same-pool error-tuned trimer ridge (0.25110), with slightly lower rank (0.36832 versus 0.37069). Across both pools the strongest nested RMSE comparator remained rank-tuned direct balanced trimer ridge (0.24451), when the direct pool is eligible for the application. That is a descriptive comparison across model families with different eligible training information. The retrospective strict-count component k=1 rank benchmark remains 0.48602. Frozen row-tuned embedding ridge references scored 0.35514 and 0.39515 RMSE in the direct and component pools, respectively. All fixed settings in the inherited grids are included in the result record.

| Outer fold | Component KNN rank-selected k | Error-selected k | Outer RMSE, rank | Outer RMSE, error | Outer RMSE, fixed 3 |
|---|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 0.28076 | 0.28076 | 0.23569 |
| 2 | 1 | 5 | 0.19912 | 0.15948 | 0.15454 |
| 3 | 1 | 3 | 0.23812 | 0.19901 | 0.19901 |
| 4 | 1 | 1 | 0.30264 | 0.30264 | 0.23865 |
| 5 | 3 | 3 | 0.30284 | 0.30284 | 0.30284 |

The two changed choices both improved outer error. Most of the remaining gap came from folds 1 and 4, where inner error preferred k=1 and outer error favored k=3. This is a comparison on shared held-out rows. The pooled rank change cannot be recovered by averaging fold rank changes.

## Design and evidence

Five fixed outer folds each have four inherited inner validation splits. Direct eligibility purges same-genus and at-least-95%-identity pairs. Component eligibility purges the full connected components. Exact-sequence/species groups are the scoring units. Within each outer training pool the new rule minimizes the mean across validation groups of the within-group mean squared error on the pooled inner out-of-fold predictions. Taking its square root gives the same selection. The comparison rule maximizes weighted midrank correlation on the same predictions. The grids are k ∈ {1, 3, 5, 10} with cutoff ties and alpha ∈ {0.001, 0.01, 0.1, 1, 10, 100} for cosine trimer ridge. Ascending grid order breaks selection ties. The chosen settings predict outer test rows from the eligible outer training pool. Both arms use pooled equal-group outer RMSE and rank. The saved fold metrics also include a fixed comparator for every model (k=3 or alpha 0.01). [Cache provenance](sources/provenance.md)

Young et al.'s Version 3 tuner separates grouped splitting from the selection metric, which is pooled unweighted inner Spearman. The paper's Table 3 benchmark reports a mean-fold rank on its own splits and model. Iqbal et al. previously used genus holdouts, absolute prediction error and an edit-distance neighbor comparator. Tu et al. describe population rank estimation with equal cluster weights. The MMseqs2 specialist issue and maintainer comments concern short-fragment clustering behavior and illustrate the dependence of clustering on the rule. They do not test KC selection. The cached upstream issue and PR listing supplied release and benchmark context. The immediate predecessor changed the weighting of the inner rank objective and already established the retrospective k=3 error reference. The ledger also includes frozen outer rescoring, donor weighting, and training-loss weighting. The new result isolates the effect of switching the selection endpoint to squared error on the same grid. [Published benchmark](https://www.mdpi.com/3042-6723/1/2/7) · [Release and tuner](https://doi.org/10.5281/zenodo.19546367) · [Iqbal et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) · [Tu et al.](https://doi.org/10.1002/sim.10326) · [MMseqs2 discussion](https://github.com/soedinglab/MMseqs2/issues/398) · [Upstream issue and PR index](https://api.github.com/repos/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/issues?state=all)

## Reproduction and scope

Run `python3 -B analysis.py` and then `python3 -B verify.py` with Python 3.10 or later and no external packages. The verifier checks released target rows, cached eligibility and fold coverage. It independently refits the candidate grids using Gaussian elimination and a separate KNN implementation, computes groupwise mean squared errors and integer-replicated tied ranks, and compares choices, saved predictions and metrics. The input parser and trimer kernel are shared. The verifier also asserts the decisive metric orderings. This is one historical panel with inherited identities, mapping and splits. Model-family comparisons are descriptive, sampling uncertainty and independent transfer are untested, and the frozen embedding references retain their original tuning.
