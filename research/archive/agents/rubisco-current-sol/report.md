# Rubisco measurement/control benchmark and calibration audit

**Research date: 26 November 2026. Local computational contribution; no experimental challenge solved.**

The deliverable is an executable, endpoint-specific measurement/control resource with two current, decision-relevant findings. **The 2025 library's eleven calibration overlaps cannot independently identify assay bias and temperature response. The 2026 cyanobacterial table also contains an imported turnover uncertainty that cannot be transported from its cited source as a simple SD-to-SE conversion.** These findings change which uncertainty and control links are defensible. They do not identify a better enzyme or invalidate the papers' biological results.

The release contains **623 historical/current S or turnover endpoint records**, a separate **469-record 30°C screen ledger**, **145 same-original-row or same-table trait pairs**, **62 explicitly typed comparator ratios**, all historical source/control rectangles, and the full numerical checks. The 145 pairs include one mutation-name/WT-flag conflict, explicitly ineligible for WT use. Empty numeric cells mean unknown, never zero. No missing uncertainties or temperature corrections were supplied by this analysis.

## Main current result: the 2025 calibration is conditional on temperature

I inspected the primary [2025 tree-of-life paper](https://doi.org/10.1073/pnas.2501433122), its downloaded source tables, and its actual calibration notebook at [GitLab commit e1c9ed945864a0744e81b9cfdbfeb609ccf67907](https://gitlab.com/milo-lab-public/rubisco_is_slow/-/blob/e1c9ed945864a0744e81b9cfdbfeb609ccf67907/1_regression_measured_lit.ipynb). The notebook fits a fixed unit slope in log space and bootstraps the common intercept. It multiplies both rates and their SD columns by the fitted factor. The retrieved notebook's Git blob, `bafbc0633ccf8ec9f9608181711a00a9182bc110`, matches the scout's source inventory.

A protocol extension was saved before inspecting these overlap values. All **11/11 original literature temperature fields contain exclusively 25°C**, including every entry in the slash-separated multi-study fields. The coupled measurements are at 30°C. The authors standardized the literature with Q10=2.2 before fitting; the actual deposited factor is **2.1292109151**, not exactly the rounded 2.1.

Let `L_i` denote the original 25°C literature rate, `A_i` the 30°C coupled rate, `q` the common Q10 and `c` the assay multiplier. Their calibration equation is

`log(L_i) − log(A_i) = log(c) − 0.5 log(q) + residual_i`.

Every temperature covariate is 0.5. The two-column design has **rank 1**: only `c/sqrt(q)` is identified, equal to **1.4355137061** here. An arbitrary positive q can be offset by changing c, leaving the calibration residuals unchanged. This is a release-specific identifiability result established from the deposited temperature fields, not an assertion that the stipulated Q10 is wrong. Additional temperature information outside this table can supply the missing constraint. I did not convert any measured rate to another temperature.

The overlap observations also distinguish mean calibration uncertainty from individual-enzyme portability:

| Quantity | Result |
|---|---:|
| Effective per-overlap factor, authors' standardized literature / raw coupled rate | 1.1902–4.7987 |
| Median residual factor after deposited correction | 1.4122 |
| Largest residual factor | 2.2538 |
| RMS natural-log error, uncorrected / corrected | 0.8744 / 0.4397 |
| Leave-one-overlap-out fitted correction RMS log error | 0.4837 |
| Effective factors outside the paper's common-factor CI 1.7–2.8 | 8 / 11 |

The eight observations outside that interval **do not invalidate a confidence interval for a common mean factor**. They show why it cannot be used as an enzyme-level prediction interval. The overlapping strains/holoenzymes were defined by the authors; exact biological equivalence of every slash-separated literature entry has not been independently certified. Residuals can combine identity, assay, literature and temperature effects; they do not estimate a pure biochemical noise distribution. The cross-validation is conditional on this small, reused calibration set.

**Decision:** use the library as a separately labeled 30°C, corrected coupled-assay rate cohort, with a shared calibration group and residual portability information. Its common-factor CI is insufficient to assign uncertainty to a candidate's absolute rate, and its calibration does not independently establish a 25°C maize comparison. The 469 screened records include 261 numeric rates, 26 explicit below-0.5 flags and 182 unreported/unsuccessful states; none supplies measured S. [Calibration pairs](results/library2025_calibration_pairs.csv), [identifiability calculation](results/calibration_identifiability.json), figure (archival input; see source-downloads.json and the source URLs), [screen ledger](derived/library2025_measurement_ledger.csv).

## A primary-source uncertainty audit of the 2026 cohort

The [2026 primary paper](https://doi.org/10.1093/plphys/kiag259) contains three newly characterized cyanobacterial kinetic pairs, two pairs explicitly imported from 2023, and a newly measured **wheat** control. The supplied request described spinach; neither Table 1 nor the inspected methods supports that identification. The modeling paragraph explicitly says the wheat kinetics were measured in the present work. The [2023 primary supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs11120-023-01008-y/MediaObjects/11120_2023_1008_MOESM1_ESM.pdf) labels its errors SD; the 2026 caption labels errors SE.

| Table entity | Identity retained | S | kcat,c / s⁻¹ | Kinetic origin |
|---|---|---:|---:|---|
| C. aponinum | PCC10605 | 58.6 ± 2.1 | 7.8 ± 0.1 | 2026 |
| C. thermalis | KOMAREK1964/111; CCALA048 in 2023 | 66.0 ± 1 | 9.1 ± 0.2 | 2023 reuse |
| G. violaceus | PCC7421 | 74.6 ± 2.2 | 11.8 ± 0.1 | 2026 |
| O. acuminata | PCC6304 | 60.4 ± 1.9 | 9.9 ± 0.1 | 2026 |
| Synechococcus sp. | PCC6301; UTEX625/CCALA188 in 2023 | 48.2 ± 2.6 | 8.7 ± 0.6 | 2023 reuse |
| T. aestivum | cultivar unknown | 96.5 ± 0.9 | 3.1 ± 0.1 | 2026 control |

The deterministic reuse audit was specified **post hoc**, after reading both tables. It checks all six common kinetic columns for both reused strains, preserving strings and units. For C. thermalis turnover the old SD and new SE are both printed as **0.2**. Allowing the full rounding ranges [0.15,0.25], `(SD/SE)^2` is at most **25/9 = 2.7778**, below every reported replicate count (old n=3–6; new n=3–5). Therefore the two error entries cannot be an unchanged-data SD-to-SE conversion, even under conservative rounding.

Across the twelve audited cells, two error pairs are incompatible without any unit conversion; a third efficiency pair is incompatible conditional on the 1,000-fold unit scale implied by its quotient definition (the 2023 SI header incompletely prints s⁻¹). Two imported Kc21%O2 means also change beyond rounding. A changed fit, sample subset or other reprocessing could explain these differences. The benchmark records the issue as **unresolved statistical provenance**, without substituting a guessed SE or deleting the mean. The particularly clear turnover field is excluded from uncertainty propagation. The [January 2026 correction to the 2023 paper](https://doi.org/10.1007/s11120-025-01185-y) addresses a photosynthesis figure and text, not this kinetic SD-to-SE transport. [All twelve checks](results/reused2023_2026_audit.csv), [exact table typography](derived/cyano2026_table1_cells.json), primary SI transcription (archival input; see source-downloads.json and the source URLs).

Endpoint conditions are retained separately: carboxylation pH 8.1, specificity pH 8.3, 25°C; kcat from N2 rates divided by 14C-CABP active-site concentration; specificity from 3H-RuBP product partition with 99.95% O2/0.05% CO2. The primary text calls 0.038 a CO2/O2 solubility ratio; its orientation requires care because that magnitude corresponds to the reciprocal solubility ordering. This textual ambiguity is retained, not used to numerically “correct” S. Per-trait exact n and covariance remain unknown.

## Historical anchors: what was actually rechecked

| Anchor | Correct primary source and numerical check | Conditions/error provenance and remaining gap |
|---|---|---|
| Galdieria partita | [Uemura 1997, 10.1006/bbrc.1997.6497](https://pubmed.ncbi.nlm.nih.gov/9144578/). Primary abstract verifies rounded S=238. Workbook 238.1 ± 8.8 and kcat 1.60 retained. | Corrected the misdirected 1996 method link. Publisher/API access did **not** expose the original numerical table. Tokara identity, detailed pH/gases, active-site quantification and error kind remain unverified. The API's HTTP200 metadata response is not full text. |
| Griffithsia monilis | [Whitney 2001 primary PDF](https://ripe.illinois.edu/sites/ripe.illinois.edu/files/2018-06/Whitney_et_al-2001-The_Plant_Journal.pdf), Figure 1a visually checked: S 167 ± 3, turnover 2.6 ± 0.1. | Sydney Harbour field collection; S at pH 8.3 on purified enzyme; turnover uses CABP sites. Figure 1a error kind/n are not defined in the inspected material. The tobacco comparator is explicitly imported from Whitney 1999, not a fresh shared control. Assay pH for turnover is unresolved. |
| Zea mays | [Kubien 2008, 10.1093/jxb/erm283](https://doi.org/10.1093/jxb/erm283). Primary text verifies S 74.9 and spinach 79.8; methods and Figure 1 caption checked. | All measurements 25°C; specificity 500 ppm CO2 in O2; turnover 14CO2/CABP sites. Figure 1 caption says SE, with specificity n = 3–5 repeated assays and turnover n = 5–7 leaf samples. Exact figure entry 4.05 ± 0.04 remains compilation-derived because the original numerical figure could not be fully transcribed. Cultivar unspecified; activation/extraction pH 8 is not evidence that S was assayed at pH 8. |

These access limits are material. **This is not a completed primary numerical certification of all three anchors.** In particular, I do not promote the historical Galdieria and maize numbers to universal acceptance cutoffs. The original paper links, retained values, validation status and endpoint conditions are explicit in [the anchor/current audit](derived/anchor_and_current_audit.csv). The inherited `_SD` column is not evidence that every underlying publication reported SD.

## What the source/control graph establishes

The predeclared historical diagnostic uses the original workbook's `use=True`, primary, mutant=False and exactly 25°C fields, positive measured endpoints, and whitespace normalization only. It retains cultivar strings. It gives 337 S records and 274 turnover records before aggregation. Multiple records are geometrically averaged only within species × source × recorded pH, with member IDs preserved. Recorded pH can itself combine endpoint conditions, so this remains a broad bridge graph, not a certified control graph.

A two-species/two-source rectangle tests a source-only multiplicative correction. Its discrepancy is `exp(abs(log(a1/b1) − log(a2/b2)))`; unity would mean that one control exactly ports the other enzyme's ratio. Rectangles sharing enzymes/studies are dependent, so no significance test treats them as independent experiments.

| Endpoint | Enzymes / source-condition nodes | Components / independent cycle rank | Rectangles / source pairs | Median / maximum discrepancy factor | Above10% / above 25% |
|---|---:|---:|---:|---:|---:|
| S | 261 / 45 | 12 / 27 | 49 / 22 | 1.057 / 1.389 | 10 / 2 |
| kcat,c | 237 / 39 | 18 / 13 | 30 / 7 | 1.118 / 1.795 | 17 / 5 |

Equal-source-pair-weighted RMS log portability error decreases from **0.1613 to 0.1074** for S and **0.3318 to 0.2569** for turnover after one-control correction. This quantifies conditional utility, not universal adequacy. The turnover squared-error benefit reverses when either Hermida-Carrera 2016 or Orr 2016 is left out.

Post hoc checks matter. Excluding measurements already annotated “suspect” in the workbook removes the largest turnover rectangle: maximum factor falls to 1.336 and RMS log error falls from 0.3365 to 0.1411 with correction. Restricting to native-expression edges gives the same turnover result. Thus the strongest claim that scalar correction generally fails is **not supported**. Residual differences still exceed 25% in four remaining turnover rectangles. The two 2016 crop studies also have cultivar/identity differences: the primary Hermida-Carrera methods identify wheat Cajeme, soybean Ransom and spinach Butterfly. Their literal species-name bridges are not verified common enzyme preparations. [Full sensitivity output](results/control_sensitivity.csv), [rectangles](results/control_rectangles.csv), [leave-one-source-out results](results/leave_one_source_out.csv).

Source identity changes the apparent anchor connectivity. In the weak paper-label S graph, Galdieria, Griffithsia and maize are connected. Reassigning the Whitney 2001 tobacco comparator to its stated 1999 origin disconnects Griffithsia from both other anchors. Reassigning the Uemura 1997 Porphyridium turnover to the Read 1994 origin recorded in the workbook also prevents it being called a fresh paired measurement. Appending the modern table with imported rows assigned to 2023 creates no redundant common-control cycles. These are statements about the finite audited graph; missing publications or documented aliases could reconnect it. [Origin graph](results/origin_graph_summary.csv).

## Reusable benchmark, challenge relevance and limitations

The [frontier benchmark](derived/assay_frontier_benchmark.csv) preserves both endpoints and their separate pH/origin fields. Dominance is calculated **within origin/condition cohorts only**. Singleton frontiers are uninformative; plus/minus values are not converted to confidence bounds. The 2025 library never enters the 25°C S–turnover frontier. The 2026 imported records are not new independent assays.

At the contextual historical diagnostic S≥238.1 and kcat≥4.05, the 145 table/row pairs contain **zero positives**; one pair has a mutation-name/WT-flag conflict. This is not a collection of certified experimental challenge failures. The [official criterion](https://millenniumproblems.bio/) requires paired assays with G. partita and maize. No appended cohort supplies both named controls, and lack of positive labels prevents validating joint-success precision or extreme-tail calibration. Independent release reviews had already established sequence/target defects; this work does not rebrand those findings or another held-out split as a discovery.

The graph uses species labels, not a complete holoenzyme sequence/isoform audit. It cannot separate genetic, culture, purification and measurement effects. Neither retained CABP normalization nor shared taxonomy guarantees comparability. No sequence model was trained and no candidate was called improved.

## Bounded prior-art and novelty assessment

Known work already establishes weak/simple tradeoff interpretations ([Flamholz 2019](https://doi.org/10.1021/acs.biochem.9b00237)), phylogenetic structure ([Bouvier 2021](https://doi.org/10.1093/molbev/msab079)), and assay-method variability ([Iñiguez 2021](https://doi.org/10.1093/jxb/erab356)). The 2025 library already uses a common control and fitted scale correction, compares sequence methods against mean/simple controls, and examines identity-dependent prediction. Its calibration notebook explicitly computes a mean-factor bootstrap. None of those methods is new here. The [April 2025 correction](https://doi.org/10.1038/s44318-025-00419-y) separately corrects a dilution error and inaccurate active-site measurements; it is not evidence that every cohort should receive a twofold adjustment.

Current developments are also distinct: [2026 Gallionella fitness/engineering](https://pmc.ncbi.nlm.nih.gov/articles/PMC13544237/) includes growth selection and selected biochemical validation; [2026 chloroplast editing](https://doi.org/10.1038/s41467-026-73783-w) includes improved measured plant kinetics, but its in-vivo specificity and plant outcomes are separate assay layers; the [September 2026 design preprint](https://doi.org/10.64898/2026.09.04.749267) explicitly distinguishes endpoint apparent specificity from conventional S. These studies do not provide the named-control labels needed to close this benchmark's joint threshold. Their results are not dismissed by this analysis.

The defensible contribution is **the versioned measurement/control resource and its specific calibration/uncertainty audit**, including (i) the all-eleven temperature confounding diagnosis, (ii) the exhaustive 2023→2026 statistical-provenance check and quarantined turnover uncertainty, and (iii) origin-aware rather than paper-label control connectivity. I did not find those exact checks in the inspected papers, supplements, downloaded code or supplied reviews. That bounded comparison is not a claim of universal first publication. The rank argument is elementary; its contribution is identifying its concrete applicability to this current release. The resource is useful incremental computational work, not a new enzyme mechanism, frontier escape or completed Millennium challenge.

## Next falsification tests

1. **Uncertainty provenance:** obtain a publicly deposited per-replicate 2023/2026 kinetic table or explicit reprocessing record. For C. thermalis, a documented different fit/sample set could explain the unchanged printed 0.2; exact original SD, n and resulting SE can resolve it. Until then the field remains flagged. No experiment or contact was initiated.
2. **Calibration identifiability:** locate an independent same-enzyme same-temperature direct/coupled measurement pair, or enzyme-specific temperature-response measurements. A bridge at a different temperature offset raises design rank and can test the shared correction separately from Q10. Validate on a different enzyme/preparation rather than recycling its fitted control. Prospective assay work is outside this computational session.
3. **Anchor portability:** a public study measuring G. partita, maize, wheat/spinach and Griffithsia under documented common conditions could test the currently disconnected comparison. A second independently measured control would make source-only correction testable rather than simply assumed.

## Reproduction and verification

Run `./run.sh` in this directory. It uses the specified shared scientific Python and local snapshots; numerical reproduction makes no network requests. Optional acquisition scripts are separate. The local `vendor/` contains pypdf 7.1.1 extracted from the supplied wheel, used only for primary PDF extraction. Scientific package versions and protocol hashes are in environment.json (archival input; see source-downloads.json and the source URLs). Shared-input snapshots (archival input; see source-downloads.json and the source URLs) preserve all seventeen read-only inputs/context files byte-for-byte; source URLs, HTTP status, sizes and digests are in `raw/retrievals.jsonl`.

Fresh independent checks re-enumerate every rectangle using a separate pivot-based calculation, reproduce the deposited correction and SD scaling, verify the temperature rank, conservative rounding bound, cohort reuse, notebook Git blob and shared-input hashes. The final verification output and execution log are [verification.json](results/verification.json) and reproduction.log (archival input; see source-downloads.json and the source URLs). These are computational/provenance checks, not independent biological replication. [Refusals](refusals.md): none observed. [Access failures](access_failures.md) are recorded separately, including the unresolved primary anchor accesses. No subagents, other model sessions, external contacts/posts/deployments or wet-lab/human/animal work were used.
