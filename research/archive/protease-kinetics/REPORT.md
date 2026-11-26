# A locked recorder-timecourse score did not improve assay transfer

26 November 2026. **The prespecified primary hypothesis failed.** A negative-control-adjusted cumulative-hazard slope did not improve rank agreement with relative biochemical initial rates over the existing recorder AUC on the selected validation panel. No new calibration method or biological improvement is claimed, and the best-looking secondary score is not being substituted as the primary outcome.

## Question and scope

The recorder has five measurements of flipped/unflipped DNA counts. AUC summarizes them, but it might mix substrate-dependent background and timecourse shape. Before computing any new score, [PROTOCOL.md](PROTOCOL.md) fixed a target-free transformation and falsification criterion. The parent reporter was replaced by the corrected DMS parent, following the original source comparison. Standard seven-residue motifs, depth ≥100 at each time, complete measurements of the three named negative constructs and matched biochemical rates were required. Signed noisy biochemical slope estimates were retained.

The primary score fits a fixed through-origin slope to `H(t)-H(0)`, where `H=-log(1-f)` and `f=(F+0.5)/(F+N+1)`, then subtracts the mean slope of the three matched negative constructs. Sampling indices 0–4 are used; this is an empirical recorder score, **not kcat or a validated cleavage-rate constant**. There are no fitted target parameters or tuned thresholds.

## Results

The eligible comparison contains **773 enzyme–motif pairs**, **131 motifs**, and **six pre-existing selected enzyme panels**. The mean within-enzyme Spearman difference, primary score minus raw AUC, is **−0.000385**. Three panels improve. A bootstrap using the same resampled motif list across all six panels gives **[−0.01582, +0.01725]** for the mean difference. That conditional sensitivity interval crosses zero and the prespecified five-of-six direction criterion fails.

| Enzyme label | Raw AUC ρ | Primary adjusted-hazard ρ |
|---|---:|---:|
| 0 | 0.8146 | 0.8011 |
| NX1 | 0.7395 | 0.7359 |
| NX2 | 0.8363 | 0.8306 |
| VI | 0.7815 | 0.7816 |
| VI* | 0.7117 | 0.7157 |
| VII | 0.6754 | 0.6918 |

All prespecified diagnostic scores, upper-decile retrieval counts and signed matched observations are saved in [metrics.csv](results/metrics.csv) and [matched_scores.csv](results/matched_scores.csv). An unadjusted hazard score looks somewhat better on several panels, but that observation is secondary and has not received independent validation. No alternative was chosen after seeing these outcomes.

## Interpretation

This rules out the claimed improvement under the locked panel criterion; it does not establish equivalence, prove that all timecourse methods fail, or validate recorder-to-cleavage mechanisms. The related enzymes and shared motifs are selected assay panels, not a probability sample or independent biological replicates. The bootstrap does not include experimental rate uncertainty or training/dataset uncertainty. The source reconstruction's detailed agreement with the authors' final data is being reviewed separately; exact input snapshots and their original/derived hashes are preserved in inputs/manifest.json (archival input; see source-downloads.json and the source URLs).

The measurements concern local peptide preference in engineered reporter/FRET contexts, not accessibility of folded endogenous sites or proteome-wide off-target cleavage. The source paper already discusses binding background and selection of the AUC interval. This failed alternative score is a retained exploratory result, not satisfaction of the goal's novelty gate. [Primary paper](https://doi.org/10.1038/s41467-025-60622-7).

## Reproduce

From the workspace root:

```sh
research/millennium/.venv/bin/python -m unittest discover -s research/millennium/protease-kinetics -p 'test_*.py'
research/millennium/.venv/bin/python research/millennium/protease-kinetics/timecourse.py
```

Three focused tests check constant-fraction cancellation, recovery of a known cumulative-hazard slope and invalid inputs. The script reuses pinned local inputs after their initial acquisition, writes all six scores/panels and 3,000 fixed-seed shared-motif bootstrap draws, and asserts finite results. [Numerical summary](results/summary.json).
