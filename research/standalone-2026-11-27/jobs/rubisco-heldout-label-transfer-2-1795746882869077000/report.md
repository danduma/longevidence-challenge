# The value of adding peer data to Rubisco reference adjustment differs by preparation

Across 535 held-out candidate–assay cells, the peer-only fit had higher log-rate RMSE than both reference baselines in both preparation groups. The Rr-inclusive fit, which estimates shared label effects jointly from peers and the *Rhodospirillum rubrum* standard, had the lowest RMSE and MAE among the tested predictors for the 375 purified cells. On the 160 lysate cells, the published-reference reconstruction had the lowest RMSE and MAE. This revises the first attempt's general preference for reference normalization after the review identified the omitted Rr-inclusive comparator. [de Pins et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/).

## Matched prediction errors

Each entry gives cell-weighted RMSE / cell-weighted MAE / equal-candidate RMSE in natural-log units. Equal-candidate RMSE averages the candidate-level mean squared errors before taking the square root.

| Predictor | All (535 cells, 172 candidates) | Lysate (160, 63) | Purified (375, 109) |
| --- | ---: | ---: | ---: |
| Candidate only | 0.2738 / 0.2011 / 0.2827 | 0.2409 / 0.1740 / 0.2500 | 0.2867 / 0.2127 / 0.3000 |
| Published reference | 0.2464 / 0.1735 / 0.2601 | 0.2376 / 0.1705 / 0.2474 | 0.2500 / 0.1748 / 0.2672 |
| Direct Rr ratio | 0.2464 / 0.1738 / 0.2602 | 0.2406 / 0.1726 / 0.2510 | 0.2488 / 0.1744 / 0.2653 |
| Peer only | 0.2870 / 0.1827 / 0.3181 | 0.3172 / 0.2079 / 0.3314 | 0.2731 / 0.1720 / 0.3101 |
| Half Rr, half peer | 0.2491 / 0.1647 / 0.2728 | 0.2567 / 0.1727 / 0.2711 | 0.2459 / 0.1612 / 0.2738 |
| Rr-inclusive fit | 0.2440 / 0.1624 / 0.2673 | 0.2573 / 0.1716 / 0.2752 | 0.2381 / 0.1584 / 0.2625 |

The Rr-inclusive fit reduced purified RMSE by 0.0107 relative to the direct ratio and by 0.0028 with equal candidate weights. It increased lysate RMSE by 0.0167 and equal-candidate RMSE by 0.0242. The pooled RMSE advantage of 0.0024 reversed under equal candidate weights. As an exploratory sensitivity check, we omitted each label from scoring in turn, keeping fitted predictions fixed. The pooled advantage reversed on omission of CA88 or CA93. The purified gain survived all 57 such omissions, and the lysate penalty survived all 49. These are descriptive sensitivity results with shared training data.

On the 522 cells supporting both pooled and same-preparation peer fits, restricting nonreference training candidates to the target preparation changed Rr-inclusive pooled RMSE from 0.2457 to 0.2459. Purified RMSE changed from 0.2389 to 0.2393, and lysate RMSE from 0.2615 to 0.2614. This exploratory restriction retained Rr readings at the labels represented in each training component.

## Design

Supplementary Data 3 supplies form, accessions and inclusion status. The analysis includes extant candidates with released numeric corrected rates and reserves Rr for the control. Each cell is the arithmetic mean of the logs of all its readings. Of 650 eligible nonreference cells, 644 had a same-label Rr observation. The whole target cell is held out. Scoring requires a path from the target candidate to the target label in the remaining nonreference graph on control-supported labels. All six predictors use the same held cells. The additive fits use the remaining connected component with equal cell weights. The Rr-inclusive fit adds one mean-log Rr observation per label in that component, including the held label, and estimates an Rr intercept just as it does for the other enzymes. Rr measurements and peer measurements are available contemporaneously in this prediction setting. [Supplementary Data 2](https://pmc-oa-opendata.s3.amazonaws.com/PMC12663927.1/pnas.2501433122.sd02.xlsx), [pinned Supplementary Data 3](https://gitlab.com/milo-lab-public/rubisco_is_slow/-/raw/e1c9ed945864a0744e81b9cfdbfeb609ccf67907/outputs/Supplementary_Data_3.csv).

The candidate-only predictor averages the candidate's remaining log raw rates. The direct ratio baseline averages corresponding log ratios to same-label Rr, then adds the target label's log Rr raw rate. The published-reference reconstruction averages the candidate's deposited normalized log rates in training and adds the target label's difference between log raw and log normalized Rr. The rounding in the released workbook makes this reconstruction approximate. No predictor uses either released rate of the held candidate cell. The fixed midpoint averages the direct ratio and peer-only predictions on the log scale. The common downstream multiplier cancels when prediction and target are expressed on the same log scale.

Preparation is inferred from forms using the paper's methods. The tested groups can differ in candidate mix, protocol and graph support. The outcome is a released reading on a repeat candidate with contemporaneous label data. The comparison leaves true rates, performance on new candidates or batches without available peers, and leader-order accuracy unresolved. Fold scores are dependent and the differences have no confidence coverage. For this replay, keep the published-reference reconstruction as the lysate squared-error benchmark and carry the Rr-inclusive fit forward for purified repeat-cell prediction. The pooled mean alone gives a fragile basis for choosing one adjustment across both groups.

## Closest precedent and evidence

The [preceding overlap audit](sources/prior-overlap.md) fit the models with and without Rr to examine in-sample leader ordering. The transfer comparison withholds entire cells and measures predictive error, and this revision extends the review's pooled Rr-inclusive spot check to independently verified cell predictions, matched preparation-restricted fits, and label scoring sensitivity. C06 addresses calibration assumptions and C08 addresses form-level hit fractions. The original study uses Rr normalization, and the pinned calibration notebook fits a common literature multiplier to aggregated rates. The cached root listing, 12 recent commit headings and empty public issue listing provide bounded coverage. [Pinned calibration notebook](https://gitlab.com/milo-lab-public/rubisco_is_slow/-/raw/e1c9ed945864a0744e81b9cfdbfeb609ccf67907/1_regression_measured_lit.ipynb).

The earlier assay paper documents use of Rr as a control. Iñiguez et al. examine measurement comparability across methods. Plantae discusses the preceding survey, and the specialist ResearchGate thread discusses assay readout choices. The inspected material did not supply the matched held-cell comparison. [Earlier assay study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11251275/), [Iñiguez et al.](https://academic.oup.com/jxb/article/72/22/7846/6331685), [Plantae](https://plantae.org/highly-active-rubiscos-discovered-by-systematic-interrogation-of-natural-sequence-diversity-embo-j/), [ResearchGate](https://www.researchgate.net/post/Is-there-any-protocol-to-measure-activity-of-RuBisCO-without-involving-radioactivity-Non-radioactive-assay-method).

## Reproduce

Run `python3 -B analysis.py` to regenerate the results and `python3 -B verify.py` to recompute and check them using the standard library. The verifier independently parses the workbook and metadata, builds connected components with union-find, solves label effects after profiling out candidate intercepts, and checks every fold prediction and aggregate including the label omission diagnostics. The analysis instead traverses the graph and solves candidate effects after profiling out label intercepts.

The de Pins et al. article is CC BY 4.0. The supplement and pinned GitLab files are cached as released and retain upstream attribution and terms. The community excerpts retain attribution. The caches used here were carried forward from the first attempt. No human review or biological validation is claimed.
