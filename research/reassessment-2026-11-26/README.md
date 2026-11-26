# What survived the skeptical audit

26 November 2026. **This assessment supersedes earlier implications that ten reproduced analyses establish ten meaningful scientific advances.** The original calculations remain available; their scientific importance is a separate question. All checks here are agent checks. No external human acceptance, upstream adoption or biological validation is claimed.

The strongest surviving result is **C07: one specific sequence–measurement pairing correction**. **C08 supplies a modest conditional inference about a selected library.** The most consequential negative finding is **C10's omitted comparator**. These are the findings we put forward for review, with the remaining work retained as technical audit material.

## C07: a small, actionable data correction

The released cross-species benchmark maps *Rhodopseudomonas palustris* to Q07N61, a 460-residue BisA53 protein. The cited kinetic study identifies the assayed organism as CGA009. The authors' structure deposition 4LF1 identifies CGA009 and links Q6N0W9; its native 461-residue core matches the deposition exactly after removal of the documented tag.

- [Released metadata, pinned row182](https://github.com/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/blob/ade87f378aa25c24a24c19915722aa04be92936e/supplementary/doubling_time_package_v2/manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv#L182).
- [Original assay, Satagopan et al. (2014)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4118107/?pdf=1), DOI 10.1074/jbc.M114.578625.
- [Author deposition 4LF1](https://www.rcsb.org/structure/4LF1), including CGA009 and Q6N0W9 reference.
- [Machine-readable patch and evidence](../contributions-2026-11-26/c07_mapping/README.md).

This supports repairing one benchmark pairing, not a discovery of the correct sequence itself: the accession and structure were already public. The two nonfragment quarantines in C07 concern unresolved provenance, not proven protein inequality. We have not regenerated the replacement embedding, demonstrated changed model performance, or obtained upstream acceptance. The public script checks accession/deposition consistency; the paper's assay identity remains a separately linked primary-source check.

## C08: a bounded library-level inference

In the selected 469-record screen, Form II's lower bound on the fraction reaching10 s⁻¹ exceeds Form I's upper bound even when unresolved outcomes are assigned adversarially. The new audit enumerates every relevant exact threshold breakpoint and intervening interval for a **single common calibration multiplier from 1.7 to2.8**, the interval reported by the [original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/?pdf=1).

| Cutoff | Factor breakpoints/interior values | Smallest Form II lower minus FormI upper | Ordering survives? |
|---|---:|---:|---|
| ≥ 10 s⁻¹ |43|+1.50593 percentage points|Yes, under the stated assumptions|
| ≥ 15 s⁻¹ |73|−11.09551 percentage points|No|

This is a finite-library sensitivity result. It assumes comparable latent endpoints for unresolved entries and a common calibration factor. It does not account for individual-enzyme measurement error, establish a 95% confidence statement for the ordering, or estimate natural-population frequencies. The expanded469-record collection must not replace the paper's 221-representative-variant survey. No improved enzyme is identified.

## C10: the omitted baseline weakens the headline

On exactly the same 780 noncanonical pairs and six leave-enzyme-out folds, predict each held enzyme's substrate activity using the mean activity of that same substrate in the other five enzymes. Every test motif has training measurements; no held-enzyme labels enter fitting. Signed negative estimates remain included.

| Predictor | Mean-fold RMSE, lower is better |
|---|---:|
| Original pooled training mean |0.430206|
| C10 AUC-linear calibration |0.308142|
| Training-only motif mean |**0.121256**|
| Training-only motif median |0.119981|

The motif mean beats the AUC-linear model in all six folds. The original arithmetic is reproducible; the practical inference from beating the weak pooled mean was overstated. The comparator uses motif identity and applies to this known-motif held-enzyme design. It does not refute the separate unseen-position test or an explicitly restricted application where only AUC is available. These additional comparisons are exploratory audit diagnostics, not a newly prespecified model-selection experiment. The underlying recorder/biochemical association is prior work: [Huber et al. (2025)](https://www.nature.com/articles/s41467-025-60622-7).

## Reclassification of all ten artifacts

| ID | Current assessment |
|---|---|
|C01|Finite schedule/objective audit. No stress recalculation or demonstrated preservation improvement.|
|C02|Fixed-horizon numerical boundary brackets. Late rate minima and the source's temperature-based stopping rule limit practical interpretation.|
|C03|Restricted, overlapping-patch null; spatial discretization dominates the small placement contrasts.|
|C04|Numerical fragility warning, substantially predictable from the initial source/capacity ratio.|
|C05|Compartment/crossing ledger and bounded null; no changed feasibility decision or ranking.|
|C06|Source/coverage audit. Temperature assumptions and assay-calibration limitations were already disclosed.|
|C07|One supported pairing correction; downstream value and upstream acceptance untested.|
|C08|Modest conditional library inference, with the above additional sensitivity check.|
|C09|Potentially reusable coverage ledger. General depth-quality concerns are prior art; fixed-score rank invariance is structural.|
|C10|Secondary calibration analysis; omitted comparator weakens held-enzyme value. Position holdout remains a narrower result.|

The thermal studies do not improve the upstream mechanical objective by demonstrating a higher thermal floor. [Vitrify](https://github.com/Marc-Dvci/Vitrify) deliberately throttles early heating to reduce stress. Slow initial warming below the glass transition was already studied in [2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5651684/), and the [2022 heart model](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/?pdf=1) explicitly distinguishes temperature stages and ends nanowarming when the coldest heart point exceeds its endpoint. A rate-screen failure late in a fixed 70 s run is not itself evidence of devitrification. Higher assumed loading also changes the upstream design space. These are possible software-review findings, not demonstrated new biological advances.

## Other sources and longevity forums

Three domain-specific skeptical reviews searched primary papers/SI, code and issue histories, later work, and field/community archives. Direct reinspection and new calculations supported the major conclusions. No claim of exhaustive worldwide novelty is made.

| Source | What it contributes to this assessment |
|---|---|
|[CryoNet,14April2002](https://www.cryonet.org/cgi-bin/dsp.cgi?msg=18925)|Mike Darwin discusses heterogeneous warming and regionally increasing magnetic coupling to lagging regions: concrete prior discussion of the broad spatial idea.|
|[Biostasis,2017](https://www.biostasis.com/new-warming-breakthrough-for-cryopreserved-organs/)|Context on nanowarming history, temperature dependence and glass-transition stress.|
|[LongeCity,January2022](https://www.longecity.org/forum/topic/114184-fight-aging-newsletter-january-17th-2022/)|Existing circulation of heart/kidney nanowarming work. A newsletter mirror with no replies, not independent expert consensus.|
|[Fight Aging!,2017 and comments](https://www.fightaging.org/archives/2017/03/a-method-of-rapidly-warming-vitrified-organs-with-minimal-damage/)|Practical limitations, including the need to load particles before preservation.|
|[Rapamycin.news cryonics discussion](https://www.rapamycin.news/t/how-come-people-here-are-not-very-interested-in-cryonics-brain-preservation/5309)|Community debate and expectations; opinions are not scientific validation.|
|[Rapamycin.news Rubisco discussion,June2026](https://www.rapamycin.news/t/the-rise-of-rubisco-a-perfect-plant-protein/25110)|Food-protein commercialization, not an exact calibration/missingness predecessor.|

Some LongeCity/PMC/publisher routes returned403, interstitials or redirects. Indexed excerpts were distinguished from complete reads; selected PMC primary text was separately fetched in memory. Forum gaps do not prove novelty. No private community was accessed, nobody was contacted, and no forum post or upstream issue was published by the audit.

The original de Pins paper already discloses Q10=2.2, the 11-enzyme calibration and its uncertainty. [Iñiguez et al. (2021)](https://doi.org/10.1093/jxb/erab356) addresses methodological disagreement, and [rbcmodel](https://cran.r-project.org/web/packages/rbcmodel/index.html) supplies preexisting temperature/kinetics infrastructure. ProtRec/MLDEEP already accounts for read quality; training-stage threshold tuning is not identical to C09's final-pair filtering, so exact duplication of its ledger was not established. [Original MLDEEP data handling](https://github.com/BorgwardtLab/MLDEEP/blob/main/src/data.py).

## Reproduce this reassessment

From this directory, install `requirements.txt`, then run:

```sh
python3 -B verify.py --download-inputs
```

Two original Rubisco CSVs are downloaded from the hash-pinned owner URLs in the earlier C08 source manifest. They remain local/ignored and are not relicensed or mirrored. Subsequent runs use `python3 -B verify.py`. To use copies acquired separately, set `--rubisco-inputs /path/to/input-directory`. The code checks input hashes and recomputes all three numerical checks against [results.json](results.json). The two other checks use licensed snapshots already in the public batch. Numerical agreement does not authenticate the papers' interpretation or establish scientific importance. `--write` explicitly regenerates the result record.

Read [LESSONS.md](LESSONS.md), [NEXT.md](NEXT.md) and the reusable [research-audit skill](../../skills/auditing-research-contributions/SKILL.md). Original reports and numerical outputs are preserved so the correction remains traceable.
