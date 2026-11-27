# Which transferred KC labels drive the Table 3 rank gap?

On the 197 historical benchmark rows, the component-purged pool beats the full direct-purge pool at exactly one or three neighbors per prediction. At one neighbor, six held-out *Limonium* rows in one exact-sequence/species group account for the largest shift. Replacing only their direct-pool predictions with their component-pool predictions recovers 0.0701 of the 0.2028 pooled Spearman gap. The gap remains 0.1214 when those six test rows are omitted from scoring.

## Comparison on fixed test rows

Every arm uses the same five held-out folds. Direct exclusion removes any outer training candidate that shares a genus or has global identity at least 95% with any test row. Component exclusion also removes candidates connected transitively by the prior graph. The transferred label is the training neighbor's log10 KC, or the arithmetic mean of the selected training labels at larger k. Exact counts use descending global identity and ascending archived row index for ties. The tie-inclusive sensitivity keeps every training neighbor at the cutoff.

| Exact k | Component | Full direct | Size-matched direct median | Size-matched direct maximum |
|---:|---:|---:|---:|---:|
| 1 | 0.4733 | 0.2706 | 0.2273 | 0.3503 |
| 3 | 0.4525 | 0.2483 | 0.2048 | 0.4251 |
| 5 | 0.3224 | 0.2606 | 0.1887 | 0.3637 |
| 10 | 0.0174 | 0.2318 | 0.0661 | 0.2557 |

All values are pooled Spearman. The 64 cached direct subsamples match the component pool size within each fold. The prior nested-tuned tie-inclusive KNN gives 0.4672 versus 0.2317. Under a common tie-inclusive k of one, the scores are 0.4733 versus 0.3051. At three they are 0.4515 versus 0.2298. The tracked direct-pool trimer ridge reaches 0.3788 on these test rows, above all direct-pool KNN settings in the table. Its prediction vector is carried over from the prior nested fit.

## Which rows move the score?

At exactly one neighbor, 115 of 197 rows switch to a newly eligible training donor and receive a different label. The component fit transfers approximately 9.5 µM from *Carpobrotus edulis* to six *Limonium* rows whose observed KC values range from 8.0 to 9.4 µM. The full direct fit transfers 25.2 µM from *Macrotyloma uniflorum*. For *L. antonii‐llorensii*, the new donor increases global identity from 89.43% to 91.41% while moving the prediction away from the observed 8.7 µM.

Swapping predictions for the three most influential groups jointly recovers 0.1102 (54%) of the gap. Those groups contain the six *Limonium* rows, *Puccinellia lemmonii* and *Aphanizomenon flos aquae*. Swapping the top five groups jointly recovers 0.1417 (70%). Groups were ordered by the pooled Spearman improvement when their predictions were swapped individually. Since Spearman reranks the entire vector after a swap, these effects need not add. All single-group scoring omissions leave a positive gap. The substitutions and omissions are retrospective diagnostics with frozen training fits.

## Evidence and scope

Young et al.'s Table 3 and Version 3 tuner establish the released benchmark. The tuner fits ridge on embeddings and deposits log-scale prediction tables. The prior direct-purge audit already established the fixed-small-k persistence against matched draws. The surviving increment is the row-level label transfer and group concentration analysis, with the full direct-pool and trimer comparisons anchoring its interpretation. Iqbal et al. previously combined plant Rubisco genus holdouts with an edit-distance KNN baseline in their broader modeling study. Petti and Eddy explain the transitive bridge problem in protein benchmark splitting. The MMseqs2 maintainer thread concerns clustering short fragments, and illustrates a distinct reason to audit pairwise boundaries directly. [Young et al.](https://www.mdpi.com/3042-6723/1/2/7) · [Version 3](https://doi.org/10.5281/zenodo.19546367) · [Iqbal et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) · [Petti and Eddy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8929697/) · [MMseqs2 discussion](https://github.com/soedinglab/MMseqs2/issues/398)

The nearest admissible donor can make a worse label proxy even as its sequence identity increases. Retuning alone therefore cannot explain this fixed-small-k pattern. The reversal at k=10 and concentration in a few test groups limit the interpretation to this historical panel and fold assignment. The substitutions do not identify biochemical causes or prospective performance.

Run `python3 analysis.py` to regenerate the results and `python3 verify.py` to recompute the decisive comparisons independently. Both use the Python standard library. Row IDs are zero-based indices into the archived KC prediction CSV. The cached rows and exact/species groups were reconstructed with the prior benchmark input parser, which uses species medians from Flamholz et al. and the previously audited CGA009 pairing. The cached identities were rebuilt with the prior global alignment scoring (match 2, mismatch −1, gap open −5, extension −1, denominator including gap columns). Fold memberships, eligible training pools, tuned choices and trimer predictions came from the direct-purge predecessor. Verification checks the archived target order and values and recomputes predictions from the cached identities and memberships. It does not repeat the sequence alignment, trimer training or nesting. The deposited code and data notices are included. [Flamholz et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6686151/)
