# A fixed blend beats label-tuned shrinkage in the purified Rubisco replay

On 375 purified candidate–label cells from 109 candidates and 57 recorded labels, fitting a blend on separate labels reduced equal-candidate log-rate RMSE from 0.26529 for the direct *Rhodospirillum rubrum* (Rr) ratio and 0.26253 for the Rr-inclusive peer fit to 0.25947. A fixed 50:50 average of those same two predictors scored 0.25756. This replay gives a reason to carry the midpoint forward as the simpler benchmark for predicting a missing reading when contemporaneous peers and reference readings are available.

| Predictor on the same 375 cells | Equal-candidate RMSE |
| --- | ---: |
| Reconstructed published normalization | 0.26719 |
| Direct Rr ratio | 0.26529 |
| Rr-inclusive peer fit | 0.26253 |
| Label-tuned blend | 0.25947 |
| Fixed 50:50 blend | **0.25756** |
| Retrospective best single blend weight (diagnostic) | 0.25748 |

The retrospective optimum assigned 0.555 weight to the inclusive predictor. It was fitted on the evaluation outcomes and serves only as a diagnostic. The midpoint came within 0.00008 RMSE of that optimum. The tuned weights ranged from 0.498 to 0.804 across the five blocks.

## Design and scope

The target is the mean natural log of the released raw readings for one candidate at one recorded label. The eligible extant candidates have a released corrected rate in Supplementary Data 3. Forms II and II/III define the purified group through the paper's sample-preparation methods. The fit can use both preparations as peers, as in the preceding audit. The comparison is limited to cells with a contemporaneous Rr control and a connected path through other nonreference cells after deletion of the entire target cell. Every scored predictor uses identical evaluation support. [Paper and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/), [Supplementary Data 2](https://pmc-oa-opendata.s3.amazonaws.com/PMC12663927.1/pnas.2501433122.sd02.xlsx).

The direct predictor averages the target candidate's remaining log ratios to Rr and adds the test-label log Rr reading. The inclusive predictor fits additive candidate and label effects by equal-cell least squares, adding Rr as another candidate with one mean log reading per label. The reconstruction uses released normalized training readings and the target label's raw-minus-normalized Rr difference. Release rounding makes this approximate.

Labels from the control-supported eligible pool were sorted by SHA-256 of their UTF-8 names and assigned round-robin to five blocks. For each block, all its candidate and control readings were purged before constructing the inner leave-one-cell-out tuning predictions. A single coefficient between zero and one minimized the mean over candidates of their within-candidate mean squared log errors. The evaluation-stage endpoints used the full control-supported graph minus the target cell, including available peers and Rr at the evaluation label. Thus a label was held back from selecting its coefficient. Its available peer and reference readings still informed the prediction. The final RMSE averages squared errors first within each candidate and then across candidates before taking the square root.

## Closest work and interpretation

The previous transfer audit compared the two fixed endpoints and a midpoint between the direct ratio and a **peer-only** fit. The present midpoint uses the **Rr-inclusive** endpoint and gives a stronger benchmark. The ledger's C06 and C08 cover assay-calibration assumptions and library-level hit fractions. The pinned upstream notebook fits a common literature multiplier on aggregated rates. Its root inventory, recent commit titles and cached public issue listing provide a bounded code check. The earlier form I work used Rr as an internal reference, and Iñiguez et al. studied cross-method corrections. Plantae's survey commentary and the specialist ResearchGate assay discussion address survey results and readout choice. These checks did not identify an equivalent comparison. [Pinned calibration code](https://gitlab.com/milo-lab-public/rubisco_is_slow/-/raw/e1c9ed945864a0744e81b9cfdbfeb609ccf67907/1_regression_measured_lit.ipynb), [earlier assay paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11251275/), [Iñiguez et al.](https://academic.oup.com/jxb/article/72/22/7846/6331685), [Plantae](https://plantae.org/highly-active-rubiscos-discovered-by-systematic-interrogation-of-natural-sequence-diversity-embo-j/), [ResearchGate](https://www.researchgate.net/post/Is-there-any-protocol-to-measure-activity-of-RuBisCO-without-involving-radioactivity-Non-radioactive-assay-method).

This is an exploratory replay on one workbook and one deterministic partition. Candidates recur between tuning and evaluation, folds share measurements, and recorded labels have no validated mapping to independent future batches. The observed differences have no confidence coverage. The outcome is a noisy repeat reading. True turnover, performance without contemporaneous peers, and candidate ranking remain unresolved.

Run `python3 -B analysis.py` to regenerate the results and `python3 -B verify.py` to recompute and assert the support, coefficients, predictions and scores from the original caches using only the standard library. The scripts adapt the previous transfer audit's workbook readers and independent least-squares solvers. The de Pins et al. article is CC BY 4.0. Its supplemental data and pinned repository snapshots retain upstream attribution and terms. Specialist excerpts retain attribution.
