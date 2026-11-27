# A matched replay of the Rubisco Table 3 input repair

Reconstructing the archived species keys and regenerating final-layer, residue-mean ESM2-t12 embeddings reproduced every Table 3 out-of-fold prediction within 0.000015 log10 units. It also recovered all ten PCA-dimension and ridge-penalty choices across the two reported endpoints. This is numerical reproduction of the archived benchmark under the released fit functions and default fold settings. [Paper and Table 3](https://www.mdpi.com/3042-6723/1/2/7) · [Version 3 archive](https://doi.org/10.5281/zenodo.19546367)

## Paired Table 3 replay

Only the focal feature vector changes. All target masks, log10 labels and fold memberships stay fixed, and scaling, PCA, tuning and ridge fitting are rerun. Focal predictions and RMSE are in log10 units. The table reports the regenerated baseline for a paired comparison.

| Endpoint | Focal observed | Focal prediction before → after | Pooled Spearman before → after | Pooled RMSE before → after |
|---|---:|---:|---:|---:|
| Specificity S | 1.07918 | 1.45601 → 1.50268 | 0.753349 → 0.753313 | 0.134949 → 0.135728 |
| CO₂ affinity KC | 2.01703 | 1.58121 → 1.54609 | 0.605059 → 0.597812 | 0.239530 → 0.243188 |

The focal absolute errors increase by 0.04667 for S and 0.03512 for KC. The updated KC rank correlation rounds to 0.598, down from 0.605, and S still rounds to 0.753. Both retain positive rank agreement and S stays ahead of KC. The KC rerun selects ridge alpha 100 in place of 10 in one outer fold, with the same PCA dimension. The largest KC change across all predictions is 0.13080 log10 units and the median change is 0.00328.

## Effect with archived hyperparameters

Freezing the archived PCA dimension and ridge penalty per fold isolates the replacement from the tuning decision. Preprocessing and ridge coefficients are still refit using the repaired input.

| Endpoint | Original Spearman | Repaired with archived settings | Repaired with retuning |
|---|---:|---:|---:|
| S | 0.753349 | 0.753313 | 0.753313 |
| KC | 0.605059 | 0.598916 | 0.597812 |

The focal predictions agree exactly across the two corrected fits. The KC rank correlation falls by 0.006144 with archived settings and a further 0.001103 after retuning. KC RMSE moves from 0.239530 to 0.241482 with archived settings and to 0.243188 after retuning. Retuning changes an individual prediction by as much as 0.064855 log10 units. The largest change relative to the original predictions is 0.130799 under either corrected fit. The KC decline and both focal error increases therefore persist under archived settings.

This is a replacement submission for the earlier attempt whose verification command was rejected for its format. The published reassessment ended at the pairing repair. The earlier attempt already measured the matched baseline and retuned effect. This submission adds the archived-parameter comparison and supplies a shell verifier.

## Archived benchmark controls

The fixed-ridge output is a separate archived control. The t12 replay also meets the prediction tolerance for all six endpoints in that branch. Its paired S and KC Spearman changes are 0.64188 → 0.62321 and 0.59775 → 0.59167. The nested turnover endpoint vC exceeds the acceptance cutoff in the primary replay at 0.00002755 maximum error, despite recovering the archived PCA dimensions and penalties.

The 1,280-dimensional t33 candidate produces maximum Table 3 baseline errors of 0.965 for S and 1.130 for KC. Its selected hyperparameters also disagree with the archive. Substituting the archived hyperparameters still leaves maximum errors of 1.077 and 0.977. This comparison points to final-layer t12 for reproducing the archived Table 3 branch with the restored inputs. It establishes numerical agreement on those outputs within the measured residuals. The historical run command and exact serialized cache remain unverified.

## Historical inputs

The 26 November release audit had already reproduced the target rows by choosing the last raw spelling for each trimmed species name. The repository history supplies the actual input file. Commit `71534b19c5019a715790ade1ab4f161ee49c226c` preserves thirteen trailing-space species fields that Version 3 trimmed. The FASTA, kinetics input and prediction CSVs have identical Git blob IDs across the two revisions. The released scripts use the final metadata entry for a repeated label and retain every corresponding FASTA record. Combining the earlier metadata with the Version 3 FASTA and observations reproduces the archived target names, order and log-median values in each branch. The literal Version 3 join yields 177 S and 191 KC target rows. Restoring the historical join yields 183 and 197, matching the archived prediction tables. [Earlier metadata](https://github.com/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/blob/71534b19c5019a715790ade1ab4f161ee49c226c/supplementary/doubling_time_package_v2/manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv) · [Version 3 normalization commit](https://github.com/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/commit/92ce749692703e8560776f28d18de7fe5aa42694)

## Replay

The Version 3 README names a t12 mean cache for fixed ridge, and the paper's methods name t33. The release's supplementary FASTA extractor calls FAIR ESM, runs one sequence per batch, mean-pools the last layer over residues, excludes both boundary tokens and stores float32 vectors. Identical sequences share one extraction here, with the original FASTA IDs restored before fitting. The script also embeds the Q6N0W9 input identified by C07. [Version 3 archive](https://doi.org/10.5281/zenodo.19546367) · [Paper and Table 3](https://www.mdpi.com/3042-6723/1/2/7) · [UniProt](https://rest.uniprot.org/uniprotkb/Q6N0W9.json)

The replay uses Version 3's fixed ridge and nested-CV fit functions with shuffled five-fold outer CV and seed 0. The nested analysis also searches the released grid with four inner folds, Spearman selection and train-only scaling and PCA. It separately tests the impact of forcing the archived per-fold hyperparameters. The acceptance rule uses full-precision predictions, a maximum absolute error below 0.00002 log10 units and identical nested PCA dimensions and penalties. For any accepted replay, the paired comparison changes only the Q6N0W9 input vector and refits with the same splits, masks and kinetic targets. [Version 3 archive](https://doi.org/10.5281/zenodo.19546367)

FAIR's specialist discussion on sequence pooling raises the general mean-versus-BOS question and supplies no additional information about the Young et al. run. [FAIR ESM discussion #591](https://github.com/facebookresearch/esm/discussions/591)

## Reproduction

From this directory, run `bash verify.sh` with the packages in `requirements-verify.txt` installed in the active Python environment. This job also supports a local Python 3.14 package directory at `runtime/verify-deps`. The verifier checks deposit and Git blob hashes, independently rebuilds both metadata joins, and refits the estimators from cached NumPy features. It checks both feature candidates against the archive and verifies the retuned and archived-parameter replacement comparisons.

For full regeneration, install `requirements.txt` into a local Python environment. Run `python3 -B analysis.py embed 12` and `python3 -B analysis.py embed 33`, each followed by `python3 -B analysis.py replay` with the same layer, then `python3 -B analysis.py summarize`. Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` and `VECLIB_MAXIMUM_THREADS=1` for the replay. Embedding regeneration loads the public FAIR weights into the job's runtime cache when needed. In the full environment, `python3 -B verify.py --verify-embedding` also checks both focal vectors for each checkpoint. The standard shell verifier uses the generated arrays in `sources/`.

The downloaded deposit retains the original MIT, CC BY 4.0 and third-party notices. The earlier metadata and FASTA remain attributed to the benchmark authors and their third-party sources. The replacement record is attributed to the UniProt Consortium under CC BY 4.0. The archive, accession evidence and derived vectors are cached in this job for verification.

To pass a repaired input file to the released scripts, run `python3 -B analysis.py export 12`. It writes the locally constructed `ids` and `emb` NumPy files to `runtime/`. Use `sources/71534b19c-flamholz_rubisco_wt_meta.csv` for the species join. The exported cache keeps all 214 FASTA IDs in order, and the two files differ in one embedding row. [Earlier metadata](https://github.com/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/blob/71534b19c5019a715790ade1ab4f161ee49c226c/supplementary/doubling_time_package_v2/manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv)

CPU checkpoint inference used the original FAIR model files. Both file sizes and multipart ETags were verified, and `sources/checkpoint-digests.json` records SHA-256 digests. The original FAIR ESM library license is cached separately. [FAIR ESM license](https://github.com/facebookresearch/esm/blob/main/LICENSE)
