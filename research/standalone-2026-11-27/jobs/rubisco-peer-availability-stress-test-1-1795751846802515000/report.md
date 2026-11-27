# A single contemporaneous peer erases nearly all of the fixed midpoint’s gain

The earlier label-holdback audit found that a fixed 50:50 blend of the direct *R. rubrum* (Rr) ratio and the Rr-inclusive peer fit beat both endpoints. That result used every available peer on the held assay label. I kept the same purified repeat-cell targets, retained every eligible observation on other labels and every Rr reading, and replayed each possible choice of exactly one surviving candidate peer on the held label.

| Available information and evaluation support | Equal-candidate log-rate RMSE |
|---|---:|
| Direct Rr ratio, 375 matched cells | 0.265287 |
| Reconstructed published normalization, same cells | 0.267195 |
| Unadjusted candidate history, same cells | 0.300014 |
| Full-peer fixed midpoint, same cells | 0.257556 |
| One-peer fixed midpoint, connected masks on the same cells | 0.265101 |
| One-peer Rr-inclusive fit, connected masks on the same cells | 0.277857 |
| One-peer fixed midpoint with direct fallback, all masks on the same cells | 0.263828 |
| Direct Rr ratio, all 385 ratio-eligible cells | 0.263417 |
| One-peer fixed midpoint with direct fallback, all 385 eligible cells | 0.262196 |

The connected one-peer midpoint retained about 2.4% of the full midpoint’s absolute RMSE improvement over the direct ratio. The policy that falls back to the direct ratio had a small descriptive gain, with somewhat different averaging because disconnected masks contributed direct-ratio losses. The direct ratio remains the practical benchmark for sparse coverage. The larger full-peer result gives little basis for expecting a comparable gain when only one peer survives.

## Connectivity and scoring

A cell pools same-candidate readings on one recorded assay label by averaging natural logs. Eligibility uses the predecessor’s candidate/accession and count reconciliation, excludes ancestral constructs from its pool, infers purified preparation from Form II and II/III, and requires at least one other candidate reading on an Rr-covered label. The target is a measured raw log rate. Both reference baselines use the same surviving candidate repeats. The joint endpoint fits additive candidate and label effects on the target candidate’s connected component of the **candidate-only** graph, then adds Rr observations for component labels. This preserves the predecessor’s support rule. Prediction uses the candidate and test-label effects. The direct ratio uses the test label’s Rr log reading plus the mean candidate-minus-Rr log reading on the remaining labels.

All 375 original cells had a connected singleton option, as expected from the last edge of a path from the candidate to the test label. Which peer survives matters. Of 4,239 singleton masks, 215 disconnected the target. At least one mask disconnected 127 of 375 cells. Choosing an original cell uniformly and then choosing its surviving peer uniformly gives a 7.26% disconnection probability. The pooled mask fraction is 5.07% because cells have different numbers of peers. With zero peers, all 375 lost peer-graph connectivity. Of the wider set of 385 cells from 111 candidates across 65 labels, ten already lacked peer-graph connectivity before masking. The direct ratio covered all 385. An alternative graph fit that permits paths through Rr would still be identifiable at zero peers. That would relax the predecessor’s component rule and use a different endpoint.

For each cell and policy I averaged squared errors over masks before averaging cells within candidates, then averaged candidates and took the square root. The connected-only comparison conditions on mask success within each cell. The fallback comparison averages over every singleton choice and uses the direct ratio for disconnected masks or when no peer exists. Losses share observations and assay labels. These exact finite-data averages leave uncertainty about real missingness, new batches, underlying kinetics and rankings.

## Evidence and reproducibility

The original methods document contemporaneous Rr normalization. The pinned authors’ notebook fits a separate common literature scaling factor on aggregated rates. I checked those sources, the released measurement workbook and candidate table, the repository root and bounded release history, and the preceding holdback artifacts and ledger. Iñiguez et al. studied method-to-method disagreement across labs. Plantae’s specialist coverage discussed the earlier activity survey, and the ResearchGate discussion covered activity readouts. These are contextual precedents. The closest quantitative predecessor is the 375-cell full-peer blend, which the independent replay checks row by row. The new result estimates error and coverage under test-label peer removal with candidate histories and controls fixed. [de Pins et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/) · [Pinned calibration notebook](https://gitlab.com/milo-lab-public/rubisco_is_slow/-/blob/e1c9ed945864a0744e81b9cfdbfeb609ccf67907/1_regression_measured_lit.ipynb) · [Iñiguez et al.](https://doi.org/10.1093/jxb/erab356) · [Plantae](https://plantae.org/highly-active-rubiscos-discovered-by-systematic-interrogation-of-natural-sequence-diversity-embo-j/)

Run `python3 -B analysis.py` to regenerate the output and `python3 -B verify.py` to check it through an independent workbook parse, union-find and a different least-squares reduction. The two programs use the Python standard library. Cached inputs and previous replay outputs retain their original attribution and terms.
