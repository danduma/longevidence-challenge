**Independent review of the computational cryopreservation contribution in `../cryo-robust-astra` — 26 November 2026.**

**Decision: a precise, useful, reproducible new computational contribution is supported, with a major correction to the assay-derived claim and a narrow novelty statement.** I found no fatal implementation error in the principal thermal comparisons. Independent calculations reproduce the equal-mean spatial effects, a schedule-ranking reversal between the two rate screens, and the conditional high-loading clamped-SAR pass. The evidence does **not** establish that a measured triplicate necessarily contains a failing heart: the mathematical bound is correct, but the source linkage required for that empirical inference is insufficiently documented. Remove that empirical inference from the accepted contribution.

**Supported novelty sentence:** Within the pinned Vitrify analytic rat-heart benchmark, this work supplies a reproducible comparison showing that mean-preserving myocardial iron redistributions can defeat the released grid-margin design, and that replacing its ramp with an equal-integrated-duty constant schedule can reverse thermal-screen rankings and yield a narrow clamped-SAR pass at a changed, high uniform loading.

This is bounded novelty relative to the primary sources and public repository history searched through 26 November 2026. It does not establish universal priority. The result is more than a source audit or reproduction of published findings: the particular perturbations, schedule comparison and counterexample were not found in the inspected prior work. It is a thermal-model contribution, without a demonstrated improvement in ice suppression, fracture resistance or organ recovery.

**Prioritized findings.**

1. **Major correction — retain the triplicate result only as a conditional lemma.** Gao reports 1.47 ± 0.36 mg Fe/mL in the text; Fig.1e labels the muscle group n=3 but plots mg Fe/mg dry weight. The methods describe drying and digestion, and identify figure error bars as SD. I inspected fresh article XML, the original figure and all supplement text. None supplies the three volumetric observations or an explicit conversion establishing that their sample SD is 0.36 on the common volumetric scale. It is plausible that the summaries refer to the same hearts, but this does not document how the moments were transformed. The report itself states this condition, yet its bold assertion that the “reported triplicate already excludes” universal feasibility exceeds the verified evidence. Delete that empirical assertion unless the volumetric triplet or its conversion is provided. The deterministic sensitivity inputs and the solver comparisons survive this correction. [Primary article](https://doi.org/10.1002/admt.202100873), archived supplement (archival input; see source-downloads.json and the source URLs).
2. **Principal numerical findings survive an independent implementation.** My own finite-volume source, face-flux and aggregation code, using midpoint RK2 and an exact 70 s horizon, reproduces the main signs and classifications. Reusing the pinned geometry and six design coefficients does not make the integration independent of those modeling assumptions; it does make this a separate calculation of their consequences. Numerical tables below are fresh reviewer results.
3. **The clamped pass is conditional and fragile.** At 0.30 mm the flat schedule has only 0.302 °C/min clearance. The loading, spatial and boundary tests confirm that a bulk mean of 1.83 mg/mL cannot certify it. Changing myocardial loading also changes the input of the upstream infeasibility statement. This counterexample does not contradict that statement at its fixed 1.47 mg/mL input. [Pinned parameterization](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/coupler/workflow.py).
4. **Broad conceptual novelty is unavailable.** Regional iron effects, local thermal constraints and staged warming precede this work. An additional 2026 author proposal explicitly discusses weakest-region margins and testing spatial compensation at equal total iron. It contains no equivalent numerical Vitrify result, but should be cited to narrow any broader placement/mean-versus-local claim. [Skow, sections 6 and 9](https://truthinresearch.pw/index.php/TIR/article/download/83/66/269).
5. **Reproduction is strong; physical validation of the new cases is absent.** Baseline agreement establishes fidelity to the pinned calculations. It cannot validate synthetic anatomy, SAR extrapolation, a 50 °C/min instantaneous ice screen, or viability. Those limits are compatible with accepting the stated computational progress.

**Independent numerical reproduction.**

I wrote [independent.cpp](independent.cpp) and [run_review.py](run_review.py). The latter imports only the pinned NumPy geometry module. It separately constructs radial shell weights, iron assignments, distance-ranked myocardial splits and active-cell connectivity. The compiled solver separately evaluates the primary heat equation, material correlations, source conversion, harmonic face conductivities, stage-wise rate minima, local crossing times and heat/enthalpy balances. It does not call the contribution's solver, setup functions or assertions. Primary constitutive expressions were checked against JR23 Table 2. [Joshi & Rabin, equations and Table 2](https://doi.org/10.1371/journal.pone.0290063).

Each internal face contributes equal and opposite power density; the exterior is adiabatic unless explicitly perturbed. Duty is interpolated at physical times 0, 35 and 70 s. The RK2 midpoint stage and step endpoints are scored, including the final RHS at exactly 70 s. Local first upward crossings are linearly interpolated, and incomplete crossings are counted. All 22 main full-field runs had positive minimum rates and every specimen voxel completed both crossings, so first crossings are unambiguous here. CPA bath, chambers and muscle are screened; the conducting wall stores energy but is excluded from specimen extrema and has zero iron.

The reviewed preregistration's SHA-256 matches its saved digest (`ffcd44855d01889a59d6467da73d8be458fee85739379345926e4696cea683b4`). It specifies the main designs, perturbation signs, two rate screens and equal-duty alternative. Its local dated file has no external timestamp, so the digest does not independently prove chronology. The assay lemma, loading brackets and extended fine-grid checks are labeled follow-ups rather than preregistered discoveries. The scenarios are deterministic sensitivities, without a sampling distribution or a statistical robustness guarantee.

The main rows below use Δt=0.010 s at 0.50/0.75 mm and 0.004 s at 0.30 mm. Rates are °C/min; temperatures are °C. Complete results, locations/times, energy, iron mass and geometry volumes are in [summary.csv](results/summary.csv) and [all.json](results/all.json).

| 0.50 mm, extrapolated SAR | All-time minimum | Local −100 to −40 minimum | Cold endpoint | All-time / window thermal decision |
|---|---:|---:|---:|---|
| Ramp, uniform 1.47 | 52.31482 | 74.15701 | −29.47073 | Pass / pass |
| Ramp, deeper half low, mean 1.47 | 39.51116 | 74.07088 | −29.75153 | Fail / pass |
| Ramp, deeper half high, mean 1.47 | 39.51262 | 74.23620 | −29.19405 | Fail / pass |
| Flat, uniform 1.47 | 52.68140 | 71.61565 | −29.15117 | Pass / pass |
| Ramp, uniform 1.83 | 55.20368 | 75.27620 | −27.81579 | Pass / pass |
| Flat, uniform 1.83 | 52.77857 | 72.69343 | −27.52106 | Pass / pass |

Both spatial maps preserve myocardial mean 1.47 and total iron to floating-point precision; each uses 1.11/1.83 mg/mL in complementary halves, with any median voxel unchanged. The depleted ramp minima occur at t=0.030 and 0 s, around −150 °C. This independently confirms failure of the released all-time screen without failure of its terminal constraints or the secondary rate screen. The uniform ramp versus flat nominal comparison reverses ranking: flat wins all-time by 0.36659, while ramp wins in-window by 2.54137. This reversal is **not universal**: at high loading with extrapolated SAR, the ramp wins both rate comparisons. The more limited wording matters.

| Uniform 1.83, clamped SAR | Grid (mm) | All-time minimum | Local-window minimum | Final coldest / hottest | Decision |
|---|---:|---:|---:|---:|---|
| Ramp | 0.75 | 44.43576 | 74.63029 | −34.75194 / −12.90371 | Fail |
| Flat | 0.75 | 50.32540 | 72.00984 | −34.38394 / −13.06235 | Pass |
| Ramp | 0.50 | 44.43576 | 75.76681 | −32.99651 / −12.84565 | Fail |
| Flat | 0.50 | 50.31942 | 73.05326 | −32.66986 / −12.99979 | Pass |
| Ramp | 0.30 | 44.43576 | 75.63253 | −33.43840 / −13.08047 | Fail |
| Flat | 0.30 | 50.30153 | 72.93807 | −33.09177 / −13.23953 | Pass |

At 0.30 mm, the flat minimum occurs at 0.182 s and −149.847 °C; the ramp minimum is at onset. Slowest local transit rates are 95.02683 versus 97.51724 for flat versus ramp. The window criterion prefers the ramp while the all-time criterion prefers flat. I independently tested three spacings, rather than claiming to rerun every reported refined mesh. The reviewed tables additionally show failure of the flat case on the original 1.4 mm grid.

| 0.30 mm, clamped ramp, mean muscle Fe 1.47 | All-time minimum | Window minimum | Cold endpoint | Terminal lower threshold |
|---|---:|---:|---:|---|
| Deeper half low | 26.95284 | 74.36609 | −35.66553 | Fail |
| Deeper half high | 26.95284 | 74.65730 | −34.72652 | Pass |

The placement effect on the cold endpoint is 0.93902 °C. Both maps still fail the all-time rate test. Only the endpoint/window classification separates them. The contribution's 0.25 mm pair, which I inspected but did not independently rerun, has the same temperature-difference sign. Changing voxelization changes anatomy and tie-ranked interfaces, so these results do not establish a continuum error bound or a mesh-independent hard pass for the favorable placement.

**Units, energy, crossings and numerical boundaries.**

At an initially uniform −150 °C, internal conduction is zero. With Fe in mg/mL and SAR in W/g Fe,

`q = 1000 × Fe × SAR × duty` W/m³, and `initial rate = 60000 × Fe × SAR × duty / C` °C/min.

I obtain `C(−150)=1,442,177.75 J/(m³ K)`. The factor 1000 is correct: 1 mg/mL = 1000 g/m³. Under clamped SAR, necessary onset loadings are 2.05915250 for the ramp, 1.81776134 for flat, and 1.73923993 even at full duty. Thus the nominal fixed-input clamped impossibility already follows at onset, without relying on optimizer exhaustion. Extrapolated necessary loadings are 1.40461439 and 1.23995369 for ramp and flat. These are necessary conditions for this all-time screen, not physical ice-formation thresholds. [Independent arithmetic and checks](results/analytic.json).

Continuous commanded duty is exactly equal: flat=0.9568032350959227, and both integrals are 66.97622645671458 duty-seconds. It is **not equal absorbed energy**. My 0.30 mm clamped runs deposit 687.09256 J (ramp) and 687.18116 J (flat), a 0.08860 J difference caused by the temperature-dependent source. At 0.50 mm the original convention runs 5766 steps through 70.00181772 s. Its rectangular sums are 66.97754694 versus 66.97796566 duty-seconds, differing by 0.00041872. The equality claim is therefore exact for the intended continuous schedules and approximate for the original discrete implementation. The independent exact-horizon comparison removes this small discrepancy and preserves the result.

The maximum full-field enthalpy imbalance across my 22 runs was 1.98×10⁻⁷ of integrated net input (0.0000198%). A separate scalar quadrature check gave an exact endpoint of −101.80463142 °C versus −101.80481494 in the solver, an absolute error of 0.00018352 °C. Two identical connected cells remained identical; zero source preserved −150 exactly. These checks test conversion, nonlinear capacity integration, symmetry and conservation. They do not validate the empirical property fits.

Halving Δt from 0.010 to 0.005 s at 0.50 mm in both high-loading clamped schedules changed the flat minimum by 0.000176 °C/min, the largest window minimum by 0.002724, and the largest endpoint by 0.000014 °C. The stage-sampled window floor can slightly overestimate the continuous local-window minimum. I therefore also wrote [crossing_review.py](crossing_review.py), which evaluates the local RHS at interpolated −100/−40 crossings using interpolated neighbor temperatures. For the four central schedule cases, the crossing correction was below 0.005 °C/min and preserved the ranking and decisions. The statistic remains a numerical instantaneous derivative, distinct from a measured finite-difference slope or a window transit average. [Instrumented results](results/crossing.json).

Two new stress tests show the pass's practical limits. At 0.50 mm, flat/clamped with mean 1.83 and complementary ±0.010 mg/mL loading gives 50.04498 °C/min, whereas ±0.015 gives 49.90775. An exterior heat-transfer coefficient of 10 W/(m² K) to a −150 °C environment changes the uniform flat case to 46.56141 °C/min and a −37.19273 °C cold endpoint, with 71.44434 J lost. At 5 W/(m² K), the cold endpoint is already −35.00360; this marginal threshold crossing should not be overinterpreted. The same coefficient to 22 °C instead warms the cold endpoint to −29.05210. These are explicit boundary perturbations, not estimates of the actual experimental environment. Adiabatic feasibility is a specified model result, not a guarantee under arbitrary boundary exchange.

The fresh scalar fixtures intentionally do not reach −100 °C. The code records missing crossings and a large sentinel for their absent window statistic. Such a sentinel alone must never be interpreted as passing a hazard-window test. In the main cases, the positive rates, completed crossings and endpoint tests avoid that ambiguity.

**Triplicate arithmetic and interpretation.**

For three values with mean μ and sample SD s, let d=μ−min(x). The centered values sum to zero and are each at least −d. Their maximum squared sum is 6d², attained at (−d,−d,2d). Consequently `min(x) ≤ μ−s/√3`. I verified the equality triplet and 100,000 independently generated centered/scaled triplets. At μ=1.47, s=0.36, the sharp bound is 1.26215390 mg/mL; adverse hundredth-place rounding gives 1.27004065 and an extrapolated-ramp onset rate of 45.20958430 °C/min.

Those numbers are correct **if** these are the moments of the relevant volumetric triplet. Dry-weight observations obey `volumetric Fe_i = dry-normalized Fe_i × dry_mass_i / tissue_volume_i`. Sample-dependent conversion factors need not preserve the sample mean, SD, or ordering. The figure's n cannot by itself certify the volume-normalized moments. I found no explicit documentation resolving this in the primary paper or supplement. This is an evidence gap, not proof that the authors' volumetric summary is wrong. It justifies retaining the bound as conditional mathematics and dropping the empirical “at least one measured heart must fail” inference. Neither the ±SD sensitivity maps nor a pass/fail screen estimates a population failure probability.

**Original-code reproduction versus physical validation.**

I compiled the contribution's `thermal.cpp` into this job's `thermal.so` and used [reproduce_original.py](reproduce_original.py) to execute its unmodified numerical functions with read-only absolute upstream paths. All 70 committed baseline metrics reproduced with maximum absolute discrepancy 9.0097×10⁻¹². Four original main comparisons reproduced separately. Baseline results (archival input; see source-downloads.json and the source URLs), [original comparisons](results/original_comparison.json), log (archival input; see source-downloads.json and the source URLs).

For those four cases, the independent exact-horizon RK2 implementation differed from the original by at most about 0.004 °C/min in the all-time floor and 0.012 °C in endpoints, without changing decisions. The contribution's own second integrator shares more setup and constitutive code with its first solver than mine; that is useful verification but not a fully independent model.

I did not run the JAX/Tesseract/FEniCSx stack: JAX is absent from the shared interpreter. These fresh results reproduce the C++ transcription against pinned outputs and independently solve the primary thermal equations; they do not claim execution of the original container stack or validation of its mechanical adjoints.

JR23 Table 5 averages rates to the instant of maximum stress. The pinned thermal benchmark instead stops at a mean-temperature glass crossing. Its extrapolated I/II/III means, freshly reproduced here, are 336.23/212.34/197.07 against published 340/198/213. The II/III ordering is reversed; the upstream boolean checks I against the other cases rather than all pairwise ordering. The contribution correctly flags this issue. Agreement of scale under a different averaging interval and synthetic geometry is not experimental validation. [JR23 Table 5](https://doi.org/10.1371/journal.pone.0290063).

**Prior work, scientific value and limits.**

The inspected public Vitrify history contains 19 commits through the pinned August 2026 head, one branch and no all-state issues or PRs. Its history already includes SAR transfer/reoptimization, horizon checks, raised-floor design and grid studies. Those are not new here. The fixed myocardial input is 1.47 mg/mL; the new perturbation and high-loading schedule tables extend that benchmark rather than overturning its fixed-input result. Exact URLs, source versions, history diffs and search queries are preserved in [sources.md](sources.md).

Prior cardiac studies already connect compartment loading to local heating/stress; the 2021 kidney study uses regional iron estimates and compares modeled temperatures with probes. [Cardiac thermal study](https://doi.org/10.1115/1.4053105), [kidney study](https://doi.org/10.1002/advs.202101691). The October 2026 pig-kidney v1 preprint already uses staged RF warming and different major-vessel/surrounding concentrations. [Primary preprint](https://doi.org/10.64898/2026.10.07.756990). None of these inspected sources reports the precise frozen-design, mean-preserving myocardial perturbations plus equal-integrated-duty ranking/pass comparison verified here. The novelty is an incremental numerical benchmark result, not discovery of heterogeneity, slow-then-fast warming, or spatial control in general.

Its useful scientific role is to separate numerical grid margin, delivery sensitivity and the temperature domain of a rate criterion. It supplies concrete regression cases and an analytic necessary-start-rate check for later optimization. A grid-safe design at one fixed loading is insufficient to establish delivery robustness, and criterion choice can change which schedule looks preferable. These are useful new calculations even after dropping the assay inference.

The geometry is an analytic matched-volume phantom, not measured coronary anatomy; chamber fraction and myocardial placement are assumed. Muscle and CPA share properties; concentration is prescribed per fixed voxel, and SAR is spatially uniform apart from temperature and loading. Neither SAR reading is newly measured. The window floor borrows a concern interval and uses a fixed 50 °C/min threshold, without nucleation/growth history, latent heat, CPA transport, feedback or injury kinetics. Kidney calorimetry shows that tissue and CPA dilution affect CWR substantially; it cannot calibrate this myocardium. [Primary tissue study, Table 1](https://doi.org/10.1016/j.cryobiol.2024.104977). Mechanics matters especially near storage/glass temperatures, exactly where flat duty adds power; no mechanical benefit of the alternative is established. A terminal gradient also cannot substitute for a stress calculation.

No requirement for a biological breakthrough or a whole-body solution was imposed in reaching the positive computational decision. The accepted result remains conditional on the named equations, geometry, loading maps, time horizon and thermal criteria.

**Exact executed numerical commands, from this job directory.**

```sh
c++ -O3 -std=c++17 independent.cpp -o independent
PYTHONDONTWRITEBYTECODE=1 ../../.venv/bin/python run_review.py > results-run.log 2>&1
c++ -O3 -std=c++17 -shared -fPIC ../cryo-robust-astra/thermal.cpp -o thermal.so
PYTHONDONTWRITEBYTECODE=1 ../../.venv/bin/python reproduce_original.py > original-run.log 2>&1
PYTHONDONTWRITEBYTECODE=1 ../../.venv/bin/python analytic_review.py > analytic-run.log 2>&1
PYTHONDONTWRITEBYTECODE=1 ../../.venv/bin/python crossing_review.py > crossing-run.log 2>&1
PYTHONDONTWRITEBYTECODE=1 ../../.venv/bin/python verify_review.py > verification.log 2>&1
```

Source retrieval used `../../.venv/bin/python fetch_sources.py > fetch.log 2>&1` and `../../.venv/bin/python fetch_extras.py > fetch-extra.log 2>&1`, plus the exact S3 supplement request recorded in `sources/supplement-retrieval.json`. PDF text extraction used `swiftc -module-cache-path .swiftcache extract_pdf.swift -o extract_pdf` and `./extract_pdf sources/heart-supplement.pdf > sources/heart-supplement.txt`. Runtime: Python 3.14.10, NumPy 2.5.5, macOS arm64, Apple Clang 21; environment (archival input; see source-downloads.json and the source URLs). The first reviewer reproduction attempt hit a NumPy-integer JSON serialization error after its baseline calculation; the serialization was fixed and the complete command rerun successfully.

Final verification.log (archival input; see source-downloads.json and the source URLs) checks the saved classifications, equal duty and mean iron, conservation, crossing corrections, original baseline agreement, unchanged read-only input hashes and document links. Code and result hashes are retained in manifest.sha256 (archival input; see source-downloads.json and the source URLs).

All generated code, binaries, downloads, results and review documents were written inside this job directory. Other jobs were read-only. Session metadata records `gpt-6.1-sol` at medium effort; no subagents or additional model calls were used. No person was contacted, and no posting, deployment or laboratory experiment occurred. [Refusals](refusals.md): none. [Access failures and fallbacks](access-failures.md) are separate from the scientific conclusions.
