# One vote per donor group leaves trimer ridge ahead for direct-pool KC

On the same 197 held-out Table 3 rows, capping the identity KNN at one vote per eligible exact-sequence/species group lowers its direct-pool log-scale RMSE from 0.2980 to 0.2799. Pooled Spearman rises from 0.2317 to 0.2391 after tuning. Trimer ridge remains ahead on both endpoints at 0.3788 and 0.2570, respectively.

| Direct-pool predictor | Pooled Spearman | Log10 KC RMSE |
|---|---:|---:|
| Row KNN, nested k | 0.2317 | 0.2980 |
| One nearest donor per group, row-selected k | 0.2222 | 0.2816 |
| One nearest donor per group, nested k | 0.2391 | 0.2799 |
| One eligible group median per vote, nested k | 0.2420 | 0.2797 |
| Trimer ridge, frozen | **0.3788** | **0.2570** |
| Embedding ridge, frozen | 0.0616 | 0.3681 |

The rank improvement in the nested result depends on retuning. Both group methods choose five instead of three neighbors in the first outer fold and keep all other choices. With the original row-selected k values, both group methods lower RMSE and lower pooled Spearman. The strongest direct group-weighted result across all tested fixed and nested settings still trails trimer ridge on both endpoints. Across individual folds, nearest-group KNN beats trimer ridge in three of five by within-fold rank. The pooled lead therefore describes this panel as scored.

## What the group cap changes

The eligible direct training pools have 62, 48, 42, 43 and 49 rows, representing 47, 44, 41, 40 and 36 groups. The cap changes 85 nested predictions. Both group methods leave the nested component-purge predictions identical to row KNN, with pooled Spearman 0.4672 and RMSE 0.2680. Thus the component-minus-direct rank gap moves from 0.2355 to 0.2281 under the nearest-group rule, and to 0.2252 under the median-label rule. Multiple votes per exact/species group explain little of that observed contrast under either rule. The component comparison is a diagnostic under another eligibility rule. The headline model comparison uses a common direct pool.

Each direct-pool fit excludes training rows sharing a genus or reaching 95% global identity with any held row. Component purging also excludes rows connected transitively through the predecessor's genus/95%/exact-species graph. For each query, eligible rows are ordered by descending identity and ascending archived row index. A group takes the identity and label of its closest row. The median variant uses the same group identity and the median log10 label of eligible group members. Distinct groups tied at the kth identity all contribute one vote. The same cutoff-tie rule applies to the row comparator. Predictions are rounded to ten decimal places to match the baseline. This procedure can fill neighbor slots with additional groups that repeated rows previously occupied.

The five outer test folds and four inner validation folds per outer training pool come from the predecessor's unshuffled GroupKFold design. In each inner fit, eligibility is determined again against the validation fold before votes are collapsed. The inner validation rows and all outer test rows keep their original row weights. For each outer fit and method, pooled inner-validation Spearman selects k from 1, 3, 5 and 10, favoring smaller k on ties.

## Evidence and limits

Young et al.'s published Table 3 and deposited tuner define the original cross-species task. Their reported 0.605 is a mean fold rank score under the paper's original design. Here the metric is pooled rank under stricter eligibility. Iqbal et al. previously used genus holdouts and a Levenshtein KNN baseline for plant Rubisco. The fetched Young repository issues/PR listing mentions data balancing in revision tooling but yielded no equivalent KC donor-vote test. BioStars contributors discuss choosing protein cluster representatives and the dependence of that choice on the downstream task. The MMseqs2 maintainer discussion illustrates how cluster membership and pairwise identity can differ. The closest internal predecessor attributed the fixed-small-k gap to newly eligible donor labels and influential test groups. This test changes the number of votes per eligible donor group while preserving the folds, labels and trimer predictions. [Young et al.](https://www.mdpi.com/3042-6723/1/2/7) · [Version 3](https://doi.org/10.5281/zenodo.19546367) · [Iqbal et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) · [Repository history](https://api.github.com/repos/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/issues?state=all&per_page=100) · [BioStars](https://test.biostars.org/p/11987/) · [MMseqs2](https://github.com/soedinglab/MMseqs2/issues/398)

This is one historical panel and partition. Exact/species grouping captures a narrow form of redundancy. The nearest representative can depend on archive order when identities tie, and the median sensitivity changes the label aggregation rule. Held-out groups with multiple rows retain their influence on pooled rank. The error reduction establishes a change in retrospective prediction error on this panel. Its prospective importance and any biochemical cause of the remaining gap are unknown.

## Reproduction

Run `python3 analysis.py` and `python3 verify.py` with Python 3.10+ and the bundled cache. Both use the standard library. The verifier rebuilds the graph and fold allocations from cached identities and group IDs, checks row order and observed values against the deposited CSV, independently selects neighbors and k, and checks every reported metric and prediction vector. It also checks the predecessor row-KNN vectors and the comparison directions. The alignment, group construction from raw sequences and ridge training are inherited. The verifier scores the frozen ridge vectors without refitting them.

The input cache is derived from the coordinator's KC neighbor-transfer and direct-purge size-control analyses of the Version 3 release at commit 92ce749. It retains only the KC rows, identities, group IDs, model vectors and fold data needed here. The alignment used match 2, mismatch −1, gap open −5, extension −1 and included gap columns in the identity denominator. The inherited targets are Flamholz species medians and the sequence metadata includes the audited CGA009 correction. The release code/data licenses and third-party notices are bundled. [Flamholz et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6686151/)
