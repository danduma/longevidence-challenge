# C05 — tissue versus bath constraint attribution

C05 found an informative null for domain-dependent feasibility: **requiring the surrounding bath as well as myocardium changed no thermal decision or paired schedule ranking in the predefined comparison.** Attribution still matters: myocardium set every whole-specimen cold endpoint; surrounding CPA set every hot endpoint and most restricted-window rate minima. Excluding the bath changes reported margins and locations, but does not rescue a failing myocardial screen here.

Status: agent-produced; awaiting coordinator verification. Executed 2026-11-26; no human expert review. Computational thermal results only.

The new contribution is the systematic compartment ledger, full-cell failure/crossing accounting, and the tested domain-comparison null for frozen schedules. Reproducing the previously accepted ramp/flat result is a prerequisite, not claimed as new. Equal-iron spatial testing, staged heating, and temperature-window concepts are prior art. The rejected measured-heart triplicate inference is absent.

The frozen [protocol](protocol.md), recorded before new scoring, defines the exact question, comparisons and falsification. “Original” means the released **margin ramp**, not the earlier non-margin optimum. Flat retains that ramp's surrounding iron field and uses its time-average commanded duty, 0.9568032351. Each delivers 66.97622646 duty-seconds over 70 s; absorbed energy is not forced equal. Myocardium is uniformly nominal 1.47 or high 1.83 mg Fe/mL; chamber CPA is 10. Both low-temperature SAR readings are tested. No schedule or loading was tuned to these results.

The conservative finite-volume midpoint RK2 solver integrates wall, surrounding CPA, myocardium, and chamber CPA in every run. Primary screening requires local instantaneous rate ≥50 °C/min throughout 0–70 s and every terminal temperature between −35 and 0 °C. Secondary diagnostics restrict the rate screen to local −100 to −40 °C or −123 to −35 °C, require complete crossings of the respective interval, and keep the same endpoint constraints. Stage states, midpoint states and the final RHS contribute to instantaneous minima. The complete trajectories, not mean temperature, determine crossing times by linear interpolation.

C05's primary comparison is label 2 (myocardium) versus labels 1+2+3 (whole specimen). The retained diagnostic domains are wall 0, bath 1, chamber 3, and heart contents 2+3. Upstream `specimen_mask` explicitly excludes wall, while `active_mask` includes it. This distinction comes from the [pinned Vitrify geometry](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/src/vitrify/geometry.py), and the thermal implementation uses the specimen mask for reported minima. The zero-iron wall starts at zero warming rate and would fail an all-time ≥50 rule in every case if included. It remains in heat transport and energy accounting; excluding it from specimen screening is an explicit source definition, not removal of an adverse simulation domain.

[Joshi & Rabin 2023](https://doi.org/10.1371/journal.pone.0290063) models CPA-filled chambers and a heart immersed in surrounding CPA, with different iron concentrations by domain. It assigns CPA properties to tissue because tissue-specific properties were unavailable, and discusses differing tissue/pure-CPA critical rates. Therefore applying the same numerical threshold to both domains is a declared benchmark assumption, not proof of their biological equivalence. The present geometry is Vitrify's synthetic matched-dimension/volume phantom, not the article's MRI anatomy. [Sources, permissions and downloads](sources.md).

Discovery used all eight combinations at 0.75 mm, dt≤0.01 s. The predefined 0.5-mm subset used both schedules at high/clamped and nominal/extrapolated settings, at dt≤0.01 and ≤0.005 s: 16 deterministic integrations total, not 16 independent experiments. Every integration has six domain summaries. Spatial cells are quadrature volumes, not biological replicates. There are four discovery schedule pairs and four additional pairs from spatial/time checks, not eight independent observations.

Discovery results below give myocardial / whole-specimen minima in °C/min. Cold is myocardial (also the specimen minimum); hot is whole-specimen surrounding CPA. All-time and window decisions are identical under both domain definitions. The −123 to −35 °C decisions also coincide.

| Loading / SAR | Schedule | All-time, muscle / whole | −100 to −40, muscle / whole | Cold / hot °C | All-time / window |
|---|---|---:|---:|---:|---|
| 1.47 / extrapolated | original | 52.320 / 52.320 | 73.174 / 73.052 | -31.197 / -9.147 | pass / pass |
| 1.47 / extrapolated | flat | 54.181 / 52.767 | 70.708 / 70.603 | -30.826 / -9.309 | pass / pass |
| 1.47 / clamped | original | 35.694 / 35.694 | 73.603 / 73.513 | -36.517 / -13.219 | fail / fail |
| 1.47 / clamped | flat | 40.429 / 40.429 | 71.012 / 70.936 | -36.109 / -13.381 | fail / fail |
| 1.83 / extrapolated | original | 56.669 / 55.296 | 74.324 / 74.161 | -29.520 / -8.845 | pass / pass |
| 1.83 / extrapolated | flat | 54.163 / 52.861 | 71.815 / 71.680 | -29.186 / -9.003 | pass / pass |
| 1.83 / clamped | original | 44.436 / 44.436 | 74.760 / 74.630 | -34.752 / -12.904 | fail / pass |
| 1.83 / clamped | flat | 50.325 / 50.325 | 72.118 / 72.010 | -34.384 / -13.062 | pass / pass |

No decision differs in any of the three screens over all 16 integrations (also true for adding only surrounding CPA to heart contents 2+3; see `results/bath_addition.csv`); no paired rate-ranking sign differs over all eight numerical pairs. Across discovery, whole-specimen all-time minima lie in myocardium in five runs and surrounding CPA in three; all restricted-window minima lie in surrounding CPA. Every cold endpoint lies in myocardium and every hot endpoint lies in surrounding CPA. Chamber CPA sets none of these whole-specimen extrema, although its own endpoints/transits are retained. These are descriptions of this finite design set, not universal rules.

At nominal/extrapolated/flat, discovery's myocardial all-time minimum is 54.180882 versus whole-specimen 52.767221 °C/min, the largest domain difference (1.413660). At 0.5 mm, dt=.005 s, those minima are 53.211295 and 52.681423, respectively. Rankings still coincide: flat wins the all-time comparison at nominal/extrapolated; original wins both restricted-temperature comparisons. High/extrapolated retains the prior exception in which original also wins all-time.

The refined high/clamped case shows the domains separately (0.5 mm, dt=.005 s). Rates and local transit ranges are °C/min; endpoints are °C. Transit rates are interval width divided by each cell's crossing duration, not the instantaneous minimum.

| Schedule / domain | All-time min | −100 to −40 min | −123 to −35 min | Terminal cold / hot | −100 to −40 transit min–max |
|---|---:|---:|---:|---:|---:|
| original / myocardium | 44.4358 | 75.7727 | 71.5096 | -32.9965 / -14.2074 | 97.9409–118.2883 |
| original / surrounding_CPA | 59.0141 | 75.7641 | 71.4794 | -30.2772 / -12.8457 | 97.9205–120.0527 |
| original / chamber_CPA | 59.5727 | 76.9010 | 72.5065 | -28.2062 / -14.4681 | 99.7218–117.5120 |
| original / whole_specimen | 44.4358 | 75.7641 | 71.4794 | -32.9965 / -12.8457 | 97.9205–120.0527 |
| flat / myocardium | 50.3196 | 73.0537 | 68.8615 | -32.6699 / -14.3535 | 95.4410–115.9327 |
| flat / surrounding_CPA | 56.4162 | 73.0533 | 68.8377 | -30.0780 / -12.9998 | 95.3522–117.7733 |
| flat / chamber_CPA | 56.9560 | 74.0903 | 69.7576 | -28.0597 / -14.6224 | 97.0420–115.1161 |
| flat / whole_specimen | 50.3196 | 73.0533 | 68.8377 | -32.6699 / -12.9998 | 95.3522–117.7733 |

Locations are not unique when symmetry or a uniform start ties cells. The detailed table records representative voxel indices, physical voxel-center coordinates, extremum time/temperature, tie counts and rate-tie bounding indices. Coordinates use the phantom array origin, in mm; they are not anatomical coordinates. For refined high/clamped/original:

- The 44.435757 °C/min all-time minimum occurs at −150 °C, t=0, throughout all 3,700 myocardial cells. One representative is (4.75, 8.25, 5.75) mm; it is not a uniquely identified cold spot.
- The −100 to −40 minimum is 75.764081 °C/min in surrounding CPA at (11.75, 8.75, 1.25) mm, t=61.3525 s, T=−40.002311 °C. The broader-window minimum is 71.479443 there at t=65.43 s, T=−35.000914 °C.
- The terminal minimum is −32.996510 °C in myocardium at (11.75, 8.75, 2.75) mm; the maximum is −12.845650 °C in surrounding CPA at (3.25, 8.75, 6.25) mm, both at 70 s.

Exact bottleneck labels are not uniformly mesh-stable. In nominal/extrapolated/flat at 0.5 mm, dt=.005 s, the −100 to −40 minimum is myocardial (71.613167), just below surrounding CPA (71.614939). At 0.75 mm it is in surrounding CPA. Thus C05 supports the observed decision concordance but does not assert that the bath always sets the restricted-rate minimum. The refined high/clamped flat bath–muscle gap is also tiny (0.000430 °C/min); fine distinctions in tied or nearly tied locations should not be interpreted as anatomical discoveries.

Failure volumes and crossing accounting are essential to the null. At 0.75 mm myocardium is 0.461109375 mL and specimen 2.19965625 mL. Both nominal/clamped schedules have all myocardial cells below 50 °C/min at some point, so their all-time union failure volume is exactly 0.461109375 mL. Only myocardial cells fail the cold endpoint: 0.013921875 mL (33 cells) for original and 0.010968750 mL (26 cells) for flat. These same cells have not reached −35 °C, and are explicitly incomplete in the −123 to −35 interval. They are not omitted to produce a favorable transit statistic. All have reached −40 °C, so their narrow-window endpoint failure remains separate from rate clearance.

High/clamped/original also fails the all-time rate over all myocardium (0.461109375 mL at discovery; 0.462500 mL at 0.5 mm), but clears the endpoint and both restricted-rate screens. High/clamped/flat has zero failed specimen volume for all three screens. All other discovery combinations have zero failed specimen volume. Surrounding and chamber CPA have zero all-time/window/endpoint failure volumes in this set. No specimen cell exceeds 0 °C at the terminal time. Wall failures are separately tabulated, not counted as specimen failures.

Every nonwall cell in every run completes the −100 to −40 crossing. The only incomplete −123 to −35 specimen crossings are the nominal/clamped discovery cells just listed. Wall crossing accounting is also retained. Complete-only transit ranges are labeled alongside missing counts and volumes in the table; they must not be read as full-domain traversal when an exit is missing.

Verification was executed with the shared scientific Python using pinned local inputs, no model calls or network dependency:

```sh
./verify.sh
# Portable alternative:
PYTHON=/absolute/path/to/python ./verify.sh
```

The verifier compiles the local C++17 solver, checks protocol/input hashes, runs analytic checks, reproduces the central baseline before new scoring, executes all model comparisons, and regenerates domain/ranking ledgers. The central high/clamped pair matches the pinned independent RK2 rates/endpoints to maximum absolute error 3.65×10⁻¹³ (tolerance 0.02). This is a code baseline reproduction, not independent experimental validation.

The largest relative enthalpy/source-energy imbalance across full-model runs is 1.230×10⁻⁷, including wall energy. The scalar nonlinear-heating quadrature differs by 0.0000635 s at dt=.001 s; an independent −123 °C crossing calculation differs by 2.49×10⁻¹⁰ s. The unit check validates q=1000 C_Fe SAR u and the conversion to °C/min. Zero-source uniform storage is preserved and symmetric heated cells remain equal. All observed rates are nonnegative; crossing order and complete/missing counts reconcile.

Halving the full-model time step at 0.5 mm changes any reported instantaneous minimum by at most 0.002872 °C/min and either endpoint by at most 0.00001481 °C. None of the predefined subset's domain decisions or rankings changes with mesh or time step. Endpoint temperatures do change materially with voxelization (e.g. high/clamped/original cold −34.752 at 0.75 mm versus −32.997 at 0.5 mm); C05 does not claim continuum convergence. Nominal/clamped was discovery-only, as predefined, so its cold-failure volumes have no new 0.5-mm confirmation here.

The full evidence is in [domains.csv](results/domains.csv), [paired comparisons](results/paired.csv), [domain decisions](results/domain_comparison.csv), [domain rankings](results/domain_rankings.csv), [constraint attribution](results/attribution.csv), [time refinement](results/time_refinement.csv), [baseline checks](results/baseline_check.csv), and [verification](results/verification.json). Per-cell CSVs and exact solver input snapshots remain local under `results/`; these permit independent re-aggregation of any mask. Six overlapping domains do not constitute independent replicates. Versions and source hashes accompany the tables in [provenance.json](provenance.json) and `manifest.sha256`.

C05 does not establish that it is appropriate to ignore surrounding CPA biologically. The bath still conducts heat into tissue and can matter mechanically; neither mechanism is tested by changing an output mask. Source-defined whole-specimen screening remains the primary benchmark, with myocardium-only as an explicitly isolated thermal diagnostic. The uniform loading, phantom geometry, SAR extrapolation/clamp, CPA property substitution, adiabatic exterior and numerical criteria are model assumptions. Co-registered tissue loading/CPA/thermometry and outcome measurements are unavailable. No ice fraction, stress, fracture, toxicity, survival, organ recovery or probability of failure follows. No human expert has reviewed this new work. [Source-access failures and development exceptions](refusals.md) are recorded separately from scientific negative results.
