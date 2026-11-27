# Registered before model fitting

Preserve the historical raw species key, observed positive log10 median targets,
deposited row order and weights, and prior Q6N0W9 feature replacement. Rebuild
the benchmark table from the archived observations, historical metadata and FASTA
and assert exact species and numerical target agreement with the deposit.
On the complete joined row table, substitute the repaired sequence and connect
rows using exact sequence or matching case-folded first two species tokens
(including numeric second tokens). Take transitive components, then apply each
endpoint mask. Compare the deposited shuffled outer 5/inner 4 KFold at seed 0
to outer 5/inner 4 GroupKFold. Assert component isolation in both levels and
count random-fold exposure. Choices, tie order, train-only standardization,
randomized PCA seed 0, and pooled inner Spearman match the release.

Compare paired training-only log means and a supervised alignment-free baseline
using cosine similarity of raw overlapping amino acid trimer counts. Tune kernel
ridge regularization over the released alpha grid extended by 0.001 and 0.01
using identical inner folds and pooled Spearman. The kernel involves no
corpus-fitted transform. The mean is an RMSE control and ranks no rows within a
fold. Evaluate pooled Spearman, RMSE and mean within-fold Spearman for ranking
models. Also examine paired random and shuffled group split sensitivity for
seeds 1 through 4, retaining PCA seed 0. Treat repeated splits as descriptive,
with the deterministic grouped split as primary. Confirm random-fold replay
against archive before replacing the focal input. Refit using cached
embeddings, as checkpoint regeneration is outside numerical verification scope.

Verification adjustment after the first fit: 10-decimal rounding of trimer
predictions reduces differences in ties at floating-point precision. The
verifier checks the saved ranks exactly and bounds independent rank differences. The analysis was rerun with this rule at both CV levels.
