# Equal held-out group weight preserves the KC rank–error trade-off

On 197 frozen out-of-fold predictions from 155 exact-sequence/species groups,
ordinary trimer ridge beats the strongest nested KNN trained on the full direct
pool. Across eligible training policies, component-purged nested KNN leads
pooled Spearman under both weighting rules. Ordinary ridge retains lower log10
RMSE in that nested comparison. A retrospective fixed three-neighbor component
fit achieves lower RMSE than either nested fit.

| Frozen fit | Row Spearman | Group Spearman | Row RMSE | Group RMSE |
|---|---:|---:|---:|---:|
| Full direct ordinary trimer ridge | 0.378783 | 0.399771 | 0.257000 | 0.246301 |
| Full direct balanced trimer, fixed alpha | 0.373630 | 0.400310 | 0.255990 | 0.244650 |
| Full direct balanced trimer, retuned | 0.364963 | 0.392039 | 0.257579 | 0.246404 |
| Full direct nested median-group KNN | 0.242032 | 0.319849 | 0.279651 | 0.273136 |
| Component nested identity KNN | 0.467217 | 0.469239 | 0.267971 | 0.267762 |
| Component nested trimer ridge | 0.335713 | 0.367372 | 0.265015 | 0.251172 |
| Component fixed k=1, strict count | 0.473330 | 0.486019 | 0.278805 | 0.277118 |
| Component fixed k=3, inclusive ties | 0.451547 | 0.394221 | 0.235781 | 0.231388 |

Under the common full direct pool, equal test-group weights shrink ordinary
ridge's rank advantage over nested median-group KNN from 0.13675 to 0.07992.
Its RMSE advantage grows from 0.02265 to 0.02684. The best retrospective
full-direct fixed-k group rank is 0.336015. Balancing the ridge loss at the
original alpha choices reverses the two ridge fits' pooled rank order by
0.000539 under equal-group scoring. Treat these ridge ranks as practically
unresolved on this panel.

Component KNN's rank advantage over full-direct ordinary ridge shrinks from
0.08843 to 0.06947 and persists. Its equal-group within-fold Spearman averages
0.47029 versus 0.42691 for ridge. Among the 64 previously sampled size-matched
direct pools, the highest group rank across their frozen nested models is
0.40649 and the lowest RMSE is 0.24244. Those extrema are retrospective across
draws. Within the fixed component grid, one-neighbor strict count gives the
highest rank and three neighbors with cutoff ties gives the lowest RMSE.
Neither fixed setting was selected by an inner validation rule for this audit.
For reporting the historical benchmark, pair rank with error and state the
training-pool policy beside the model. Under the frozen nested policies the
preferred fit depends on the endpoint used for scoring.

## Evaluation and evidence

Each group contributes unit mass, divided evenly among its observed rows.
The weighted midrank is the mass below a value plus half the tied mass.
Spearman is weighted Pearson correlation of the two weighted midranks.
RMSE is the square root of the mean of group-level mean squared errors on
log10 KC. The group IDs join rows linked by exact sequence or the same species
label, including transitive links. Row scores give every measured row unit
weight. This application uses the established total rank estimator.
[Tu et al.](https://doi.org/10.1002/sim.10326)

The same five held-out folds are used throughout. Direct purging excludes
training candidates that share a genus or reach 95% global identity with a
test row. The component rule removes transitive graph neighbors as well. Its
fold training counts are 31, 37, 27, 28 and 26 against 62, 48, 42, 43 and 49
under direct purging. The revision adds the component arm and matched direct
draws after the reviewer identified the omitted rank comparator. The ledger
already contains dominance comparisons at uniform row weight, donor weighting,
training-loss weighting and single test-group omissions. The new comparison
gives every held-out group equal mass while retaining the observed labels.

Young et al.'s Table 3 reports mean-fold rank. The deposited tuner computes
pooled rank and errors over log-transformed predictions. Iqbal et al. provide
earlier Rubisco genus holdouts and edit-distance KNN. Petti and Eddy describe
transitive bridges in protein benchmark splits. The MMseqs2 specialist thread
concerns short-fragment clustering behavior in a different pipeline. It
illustrates why a stated identity threshold alone is insufficient evidence of
group membership. [Young et al.](https://www.mdpi.com/3042-6723/1/2/7) ·
[Version 3](https://doi.org/10.5281/zenodo.19546367) ·
[Iqbal et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) ·
[Petti and Eddy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8929697/) ·
[MMseqs2 discussion](https://github.com/soedinglab/MMseqs2/issues/398)

This is a descriptive comparison for one panel, grouping rule and fold
assignment. Pooled rank includes cross-fold comparisons. Nested choices still
reflect their original row-weighted tuning, and the fixed-grid extrema use
outer outcomes. Sampling uncertainty, performance on new panels and biological
effects remain unmeasured.

Run `python3 -B analysis.py` to regenerate the score table and
`python3 -B verify.py` to verify it. The verifier independently integer-replicates
rows and computes ordinary tied ranks. It averages within-group mean squared
errors independently, checks the archived species order and targets, fold
partition, frozen cache hashes, eligible pool subset relationships, and prior
pooled row scores. Group construction, alignments and model training are
inherited. See [cache provenance](sources/frozen-provenance.md).
