# C06 — Rubisco temperature-calibration evidence audit

**Executed 26 November 2026. Agent-produced; awaiting coordinator verification. No human expert review.**

C06 finds **no verified temperature bridge for any of the eleven 2025 calibration enzymes in the audited finite inputs**. The full historical workbook contains no eligible within-origin in-vitro turnover temperature pair. Public primary tables do supply measurable temperature responses: C06 extracted 20 crop series from the original Hermida-Carrera tables, comprising 60 measured turnover estimates and 60 correlated contrasts. These plant series constrain their own source-defined temperature responses. Transferring them to the entirely microbial calibration cohort requires an additional biological assumption.

This is an evidence-coverage result, **not proof that relevant microbial measurements do not exist publicly**, and not evidence that the authors' Q10=2.2 is wrong. The known eleven-overlap rank-1 calculation is baseline context and is not counted as C06's contribution.

## Question, protocol and extraction

The question and eligibility rules were saved in [protocol.md](protocol.md) before scoring new outcomes. The audit reads all 634 rows and all 34 columns of the original workbook's sole `kinetic_data` sheet, including secondary sources, mutants, excluded rows and missing endpoints. Every cell is reconciled with the prior unfiltered CSV; no temperature is filled or inferred. All eleven overlaps are re-extracted from the original deposited 2025 library, not taken on trust from the earlier analysis.

Candidate enumeration first uses literal species names to avoid overlooking recorded temperature differences. This is deliberately permissive for inventory, not evidence of exact construct identity. Eligibility then checks origin, endpoint availability, recorded expression/isoform/variant/pH, and in-vivo flags. Missing identity fields remain unknown. No cross-species alias, PCC strain equivalence or isoform equivalence is guessed. The [row ledger](README.md), [all candidate pairs](results/workbook_candidate_pairs.csv) and [67-source inventory](results/source_inventory.csv) expose both included and excluded information.

The public expansion parses the original [25°C Table I](https://pmc.ncbi.nlm.nih.gov/articles/PMC4972260/) and the authors' deposited [15/35°C Table S1](https://researchonline.ljmu.ac.uk/id/eprint/3863/3/Supplemental%20Material.pdf) directly. The methods identify the native crop populations/cultivars, 14CO2 fixation, Michaelis–Menten Vmax divided by 14C-CABP active sites, pH 8 adjusted at each temperature, and nonoxygenic measurements; oxygenated assays also determine Kcair. Their fitted saturation Vmax is the source's measured kcat estimate, not a temperature-model prediction. Measurements at different temperatures belong to the same source-defined cultivar/native-enzyme population. They do **not** establish a homogeneous molecular isoform, identical individual molecules, matched replicate samples or equivalence to any recombinant microbial construct. Two crop entries lack a named cultivar; that missingness is retained, and no across-study identity link is made.

## Finite workbook coverage

| Audit quantity | Executed result |
|---|---:|
| Original rows / source labels | 634 / 67 |
| Rows with missing temperature | 3 |
| Rows with a reported temperature other than 25°C | 23 |
| Literal-species cross-temperature candidate pairs | 46 |
| Candidates within one source | 14 |
| Within-source candidates with turnover at both temperatures | 1 |
| Eligible within-source in-vitro turnover pairs | **0** |

The 14 within-source candidates consist of six soybean Wayne comparisons with missing turnover, six Trifolium comparisons with missing turnover, one Thalassiosira comparison with missing turnover, and one Fragilariopsis comparison at −1/20°C explicitly labeled **in-vivo estimate** in the workbook. C06 does not convert the latter into an in-vitro enzyme Q10. These counts describe dependent row combinations, not biological replicates.

The other 32 candidates cross assay origins: 30 Rhodospirillum rubrum combinations, one Prochlorococcus MIT9313 combination and one Thermococcus combination. Five Rubrum cross-origin combinations have both turnover values; they compare 25°C literature entries with a 30°C Satagopan entry. Different origins, incomplete construct correspondence and condition mismatches prevent interpreting those divisions as measured temperature responses. The ledger leaves their ratio and Q10 cells blank. The non-calibration Thermococcus 25/83°C comparison is likewise not a matched temperature series.

## What the retrieved primary series measure

For each source-defined crop population, the code calculates only

`R = kcat(T_high) / kcat(T_low)` and `Q10_interval = R ** (10 / (T_high − T_low))`.

| Observed interval | Source-defined series | Minimum effective Q10 | Median | Maximum |
|---|---:|---:|---:|---:|
| 15→25°C | 20 | 1.3333 | 1.8819 | 4.4000 |
| 25→35°C | 20 | 1.3600 | 1.5663 | 2.2500 |
| 15→35°C | 20 | 1.4268 | 1.7579 | 2.5690 |

For example, maize Carella has measured estimates 2.4, 4.1 and 5.7 s⁻¹ at 15, 25 and 35°C: its observed 25→35°C response is 5.7/4.1 = **1.39024**. These are calculations from the published measurements, not new biological findings. Their value to C06 is establishing what temperature evidence can actually be extracted and linked to the calibration.

All 60 source estimates, reported SEs and source locations are in [primary_temperature_points.csv](results/primary_temperature_points.csv); the 60 dependent contrasts are in [measured_temperature_pairs.csv](results/measured_temperature_pairs.csv). There are only **20 population series in one study**, not 60 independent observations of Q10. Printed zero SEs remain rounded reported zeros. Exact endpoint n and covariance are unavailable, so no confidence intervals or significance test are manufactured. The sample summaries above do not estimate a universal Rubisco Q10. In particular, no 30°C value is measured in this source and none is interpolated, even though 30°C lies inside its measured range.

## Calibration coverage, condition mismatch and rank

The eleven deposited overlaps are all marked prokaryotic: Rubrum, Methanococcoides, Thermosynechococcus, Acidithiobacillus, Prochlorococcus, Anaerolineae, Synechococcus, Calditrichaeota, Halothiobacillus, Promineofilum and Riftia endosymbiont constructs. **None of the 20 extracted crop species even matches an overlap species label**, so none can establish exact accession/LSU/SSU correspondence. Verified same-enzyme series coverage is **0/11**, and exact measured 25→30°C pairs are **0** in the retrieved primary crop tables. This is finite audited coverage, not a global absence estimate.

The [eleven-row coverage ledger](results/calibration_overlap_coverage.csv) also preserves assay and condition mismatch. The [2025 primary methods](https://doi.org/10.1073/pnas.2501433122) use a 30°C NADH coupled assay at 4% CO2/0.4% O2 and inhibitor titration for active-site concentration. The crop measurements use radiometric saturation kinetics and CABP binding. Sharing the CABP concept is insufficient to establish interchangeable active-site estimates. Exact pH compatibility with these recombinant constructs is not established. Source origin, enzyme identity, gas regime, saturation estimation and temperature interval therefore differ or remain unknown.

The executable [rank tests](results/rank_tests.csv) separate three claims:

1. **Admit only verified calibration bridges:** no new admissible row is available; the two-parameter design remains rank **1**. C06's addition is the executed coverage/exclusion basis for this result, not rediscovery of the baseline rank.
2. **Force one Q10 across the plant series and microbial calibration:** an actual temperature contrast adds `[0, ΔT/10]`, and the two-parameter design becomes rank **2**. This is a conditional algebraic rank gain; the equality of temperature responses is not independently demonstrated.
3. **Allow plant Q10 to differ from microbial Q10:** with columns `(assay alpha, microbial beta, plant beta)`, the enlarged design has rank **2 of 3**. The explicit null vector `[0.5, 1, 0]` survives. Thus the plant data can constrain their own temperature coefficient while microbial assay bias and temperature response remain inseparable.

No fitted plant coefficient is substituted into the 2025 correction and no microbial rate is temperature-normalized by C06. A fully heterogeneous model with separate crop coefficients would not repair the microbial null direction either.

## What remains unsupported

The exhaustive absence statement applies to **all supplied workbook rows**, and the positive extraction statement applies to **all 20 rows of the successfully retrieved primary crop table**. The audit does not exhaust every latent table in all 67 historical sources or all public literature.

The original [Orr 2016 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5047088/) explicitly reports 75 plant species at 20/25/30°C. Its Table S1 was inaccessible at the logged endpoints, so C06 does not invent its values, count them as extracted series or treat the access failure as missing biological evidence. Those plant series could improve direct interval coverage but would still require a cross-enzyme transfer assumption for the microbial overlaps. Thermosynechococcus and Rubrum temperature-related papers were identified as potentially relevant leads; original full numerical tables were not retrieved here, so their ability to supply a compatible turnover bridge remains **unresolved**, not disproved. Specificity-response data alone cannot identify turnover Q10 without further measurements or a kinetic model.

The plant temperature response literature, matched-temperature comparisons and elementary rank arguments are prior art. C06's bounded contribution is the full workbook exclusion ledger, primary-table extraction, exact-overlap coverage audit and explicit demonstration that apparent rank gain from plant data depends on a transfer assumption. It identifies the information still needed to validate calibration transport, without proposing a better enzyme or claiming a biological challenge is solved.

## Verification and provenance

Run **`./run.sh`** in this directory. Main verifier: **[verify.py](verify.py)**, executable with `python3`. It uses pinned local snapshots only, with bytecode writes disabled and no network/model calls. Fresh verification checks all input hashes, all workbook cells, byte-identical re-extraction, six independently specified primary-table sentinels, log-ratio arithmetic and telescoping contrasts for every crop series, the in-vivo exclusion and the rank/null-vector distinction. See [verification.json](results/verification.json), [run.log](coordinator-checks.json) and [environment.json](results/environment.json).

Public acquisition URLs, versions, retrieval dates, redistribution restrictions and SHA-256 locations are in [sources.md](sources.md). Raw third-party material is retained for local verification only; no public redistribution permission is inferred from download access. [refusals.md](refusals.md) distinguishes source failures from model refusals. No human expert has reviewed this new contribution.
