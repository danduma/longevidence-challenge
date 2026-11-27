# Assay overlap connects the two extant Rubisco leaders, but their adjusted order depends on the reference

RBC2_36 and RBC8_39 have no assay label in common. They can still be compared under a model with a common effect for each label. With the study’s R. rubrum reference removed, one shortest path is:

```text
RBC2_36 → CA106 → RBC2_45 → CA96 → Hn → CA188 → RBC8_39
```

The contrast is estimable after every single-label omission among the 127 labels in the nonreference graph. Of the 8,001 distinct two-label omissions, only the pair CA187 and CA188 disconnects the candidates. Those two labels are RBC8_39’s entire observed assay support. The result holds in the reference-inclusive graph as well, among 131 labels and 8,515 pairs. Thus no single intermediate label is essential, even with the reference removed. This establishes support for an indirect comparison under the common-label-effect model. [Published Supplementary Data 2](https://pmc-oa-opendata.s3.amazonaws.com/PMC12663927.1/pnas.2501433122.sd02.xlsx), [Bailey and Sajjad, incomplete-block estimability](https://link.springer.com/article/10.1007/s42519-021-00215-x).

## What happens to the ordering

Every signed contrast is RBC2_36 minus RBC8_39. Rate differences are on the workbook scale before the common literature multiplier. Log differences use natural logs.

| Comparison | Full contrast | Range across single-label omissions |
|---|---:|---:|
| Published reference-normalized arithmetic means | +0.0714 s⁻¹ | −0.4286 to +0.5714 s⁻¹ |
| Raw rates divided by the same-label Rr observation, geometric means | +0.00621 log units | −0.03859 to +0.05102 |
| Candidate and label effects on log raw rates, excluding Rr | −0.09328 log units | −0.19005 to −0.01194 |
| Candidate and label effects on log raw rates, including Rr | +0.05272 log units | +0.01427 to +0.11136 |

The peer-only log fit estimates an RBC2_36/RBC8_39 ratio of 0.911. Including Rr raises the fitted ratio to 1.054. Both signs survive every single-label omission within their own model. Among omissions from the peer graph, CA276 gives the most negative log contrast even though neither leader was measured there. Omitting CA188 moves the estimate closest to zero. The same-label Rr ratio relies only on the control measurements paired with the two candidates. It puts the full ratio at 1.006 and reverses when CA187 is removed. This ratio is a fresh baseline formed from rounded raw observations. The deposited reference-normalized means are the baseline for comparison with the authors’ processing. [Published Supplementary Data 2](https://pmc-oa-opendata.s3.amazonaws.com/PMC12663927.1/pnas.2501433122.sd02.xlsx), [de Pins et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/).

The two log fits on the already normalized column also disagree across reference inclusion, giving −0.12096 without Rr and +0.04500 with it. They have the same sign pattern across every single-label omission as the log raw fits. On the original raw scale, the peer fit gives −0.8819 s⁻¹ and turns slightly positive after omitting CA188. Including Rr gives −0.0021 s⁻¹ on the raw scale. The additive fits on normalized rates give −0.8478 without Rr and +0.5706 with it. Additional fits are exploratory checks of scale and adjustment choices. The model dependence means the overlap audit gives no basis for promoting either candidate on the strength of a single adjusted rank.

## How the comparison was made

The unit in the calculation is a released Supplementary Data 2 reading. The graph collapses repeated candidate–label pairs to one edge. Least squares gives each reading equal weight. Supplementary Data 3 determines eligibility and reconciles released counts and accessions. The nonreference graph contains 660 readings across 258 extant candidates and 127 labels. It excludes the two reconstructed ancestors, five low-flag records and the Rr standard. The leaders’ component contains 237 candidates and 105 labels. The reference-inclusive graph adds 130 Rr readings and connects all 259 extant candidates. Assay labels are used as recorded. [Supplementary Data 3 at the pinned release](https://gitlab.com/milo-lab-public/rubisco_is_slow/-/raw/e1c9ed945864a0744e81b9cfdbfeb609ccf67907/outputs/Supplementary_Data_3.csv).

The fitted mean is a candidate intercept plus a label intercept. Along a chain, differences within each label cancel the label intercept and sum to the terminal contrast. With disconnected components, each component admits its own arbitrary shift of candidate and label effects. The fit uses all readings in the leaders’ component, including multiple paths. The log model treats label effects as common multiplicative shifts on raw rates. Including Rr adds both paths and data that influence the estimated offsets. The paper reports distinct preparation methods for the two forms, and label effects can interact with candidate or preparation. The release also gives rounded readings and limited preparation provenance. The deletion ranges describe sensitivity to recorded labels and carry no confidence coverage. [de Pins et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/).

## Relation to previous work

The earlier selection replay already established the small reference-normalized mean gap, lack of direct overlap and instability of the extant point leader to assay-label omissions. This analysis reuses that mean baseline and adds a check of indirect estimability, exhaustive one- and two-label connectivity, and the sign reversal when reference observations enter the fitted log model. C06 and the assay × temperature audit concern transport of calibration across protocols. C08 concerns form-level hit fractions. The released calibration notebook starts from already aggregated rates and applies a common multiplier derived from the literature comparison. Its calculation cannot change the signs of the same-scale contrasts here. The repository tree and cached commit and issue listings were inspected for a candidate-by-label reanalysis. [Pinned calibration notebook](https://gitlab.com/milo-lab-public/rubisco_is_slow/-/raw/e1c9ed945864a0744e81b9cfdbfeb609ccf67907/1_regression_measured_lit.ipynb).

The later Kehl survey discusses turnover across under-sampled clades. Specialist Plantae commentary on Davidi et al. highlights an earlier screening leader, while the indexed ResearchGate assay thread discusses readout choices. The Rapamycin.news thread is about food protein. The inspected material did not supply an adjusted comparison of these two candidates. Search coverage is limited. [Kehl et al.](https://www.nature.com/articles/s41467-026-73982-5), [Plantae](https://plantae.org/highly-active-rubiscos-discovered-by-systematic-interrogation-of-natural-sequence-diversity-embo-j/), [ResearchGate](https://www.researchgate.net/post/Is-there-any-protocol-to-measure-activity-of-RuBisCO-without-involving-radioactivity-Non-radioactive-assay-method), [Rapamycin.news](https://www.rapamycin.news/t/the-rise-of-rubisco-a-perfect-plant-protein/25110).

## Reproduce

Run `python3 analysis.py` to regenerate the results and `python3 verify.py` to check them with Python’s standard library. The verifier uses union-find for connectivity and profiles out candidates to solve for assay offsets by Cholesky factorization. The analysis uses breadth-first traversal and profiles out assay offsets to solve for candidate effects by conjugate gradients. The two calculations agree on each deposited estimate and deletion outcome. Both use the same workbook reader, adapted from the earlier selection replay.

The de Pins et al. article is CC BY 4.0 and the supplement is cached as released. GitLab files retain their upstream status. The Bailey and Sajjad and Kehl snapshots retain their original notices. Community excerpts or snapshots are attributed in the submission.
