# C05 protocol: tissue versus bath constraint attribution

Frozen before new scoring, 2026-11-26. Local timestamp/hash only; no external registration or human review.

Question: In the pinned Vitrify synthetic heart phantom, does screening myocardium alone rather than all CPA+tissue change the thermal decision or the ranking of the released margin ramp and its equal-integrated-duty flat counterpart? Which domain, location, time and temperature sets each constraint?

Primary comparison: original means the released **margin ramp** (not the earlier non-margin optimum); flat retains its three iron-shell coefficients and replaces its three RF ordinates by (u0+2u1+u2)/4. Uniform myocardium is nominal 1.47 or high 1.83 mg Fe/mL; chamber is 10; wall is zero; surrounding CPA uses the pinned radial basis. Cross both schedules, both loadings and both SAR readings (linear extrapolation below -80 C, or clamp to -80 C). Duration 70 s, initial -150 C, adiabatic outside, nonlinear published properties, exact-horizon midpoint RK2, conservative harmonic face fluxes. No parameter fitting or schedule tuning.

Before scoring: reproduce the central 0.5-mm high/clamped original and flat baseline against the portable independent solver's pinned core, at dt<=0.01 s (rate tolerance 0.02 C/min, endpoint tolerance 0.02 C). This is a gate, not the new contribution.

Discovery: all 8 combinations at 0.75 mm, dt<=0.01 s. Predetermined decision subset: high/clamped and nominal/extrapolated, both schedules, at 0.5 mm, dt<=0.01 and <=0.005 s. These include the previously reported conditional pass and ranking reversal; selection precedes new compartment results. Report all outcomes without promoting a favorable case.

Domains: myocardium label 2; chamber CPA label 3; surrounding CPA label 1; heart contents 2+3; whole specimen 1+2+3; wall label 0 diagnostic only. Geometry explicitly excludes wall from specimen. Every run still integrates every active cell. Neither bath nor wall is removed physically.

Primary screen: instantaneous minimum >=50 C/min throughout [0,70] and all terminal temperatures in [-35,0] C. Secondary diagnostic screens restrict the local rate temperature to [-100,-40] and [-123,-35] C with identical endpoint constraints; require all selected cells to complete the respective upward crossing pair. Record missing entry/exit counts explicitly and do not impute a transit for incomplete cells. These are thermal screens, not biological safety tests. Rate minima use both full and midpoint stage states plus final RHS. Crossings use linear interpolation of full step states; no mean-temperature substitute.

For every domain report instantaneous minima with label, coordinates, time and temperature; per-cell transit min/max durations and rates with missing counts; terminal cold/hot locations; volumes/fractions failing rate, cold, hot, crossing or union criteria. Voxel counts are quadrature cells, never independent replicates. Compare paired ramp-flat signs separately for each rate definition and domain; report decision dependence on including bath explicitly.

Falsification: domain dependence is unsupported if all corresponding decisions and rate-ranking signs coincide; attribution to bath is false if myocardial cells set the same bottleneck (ties reported). Any apparent decision switch lost under the predefined space/time checks is only mesh-dependent. No claim of universal novelty, ice, stress, fracture, survival, recovery, or measured-heart failure. The rejected triplicate inference is excluded.

Verification: source-unit checks (mg/mL to g/m3 factor 1000; 60 s/min), zero-source uniform-state preservation, one-cell nonlinear heating quadrature, face conservation via global enthalpy/source balance, monotone crossing order, complete-cell accounting and predeclared time/space comparisons. Record versions, SHA-256 input/code/output hashes, primary-source definitions and retrieval attempts. Stronger unavailable experimental data are a limitation, not a model refusal.
