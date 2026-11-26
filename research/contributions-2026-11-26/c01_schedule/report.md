# C01 — Equal-dose three-ordinate schedule family

**Executed result: the wider post-glass criterion does not force a sampled all-time/post-glass tradeoff, but the narrower criteria do.** The prespecified middle-deficit schedule is the only nondominated member of the nine-schedule discovery family for whole-specimen all-time versus observed local −123 to 0 °C minima, using the frozen 0.02 °C/min comparison tolerance, in each of four fixed loading/SAR scenarios. Under terminal constraints this remains true in three scenarios; nominal/clamped at 0.75 mm has no eligible schedule. With the narrower −123 to −35 °C or −100 to −40 °C criterion, four schedules remain nondominated in each terminal-eligible discovery scenario. These are finite deterministic model comparisons, not independent replicates or a continuous optimality certificate. No new work here is human expert reviewed.

The incremental C01 contribution is the executed equal-command family, conditional nondominance and localization of its objective tradeoffs. The accepted ramp/flat comparison, spatial testing, rate-window principle and staged heating are prior art. The rejected measured-heart triplicate inference is not used. Read the frozen [protocol](protocol.md), [primary sources and rights](sources.md), [provenance hashes](provenance.json), and [correction/access ledger](refusals.md).

The fixed shell and synthetic geometry come from Vitrify commit `1d498ee5c583e0112a762e60892d576ec2e27f57`. Myocardium is uniformly 1.47 or 1.83 mg Fe/mL; chambers are 10 mg/mL. The exterior is adiabatic, initial temperature −150 °C, duration exactly 70 s. Both clamped and extrapolated low-temperature SAR readings are evaluated. These are deterministic loading assumptions, without empirical population implications. The portable independent finite-volume midpoint RK2 solver was adapted locally; wall energy is included but wall cells are excluded from specimen rate and temperature decisions.

All schedules are linear between 0, 35 and 70 s, with duty in [0,1]. Their exact continuous integral is 66.9762264567146 duty-seconds. The duty-deficit simplex has d0+2d1+d2=0.1727870596163092. Its three vertices, three edge midpoints and barycenter were fixed before scoring, alongside original and flat controls. A deficit name means *less* power at that ordinate.

| Schedule | u(0) | u(35) | u(70) |
|---|---:|---:|---:|
| original | 0.844638719 | 0.991287111 | 1.000000000 |
| flat | 0.956803235 | 0.956803235 | 0.956803235 |
| deficit_start | 0.827212940 | 1.000000000 | 1.000000000 |
| deficit_mid | 1.000000000 | 0.913606470 | 1.000000000 |
| deficit_end | 1.000000000 | 1.000000000 | 0.827212940 |
| edge_start_mid | 0.913606470 | 0.956803235 | 1.000000000 |
| edge_start_end | 0.913606470 | 1.000000000 | 0.913606470 |
| edge_mid_end | 1.000000000 | 0.956803235 | 0.913606470 |
| barycenter | 0.942404313 | 0.971202157 | 0.942404313 |

The central known control was reproduced first with the retained unmodified independent solver at 0.5 mm: ramp 44.43576 and flat 50.31942 °C/min, with complete −100/−40 crossings and conservation within the protocol tolerance. This is validation, not new contribution. New scoring comprised 36 runs at 0.75 mm/dt=0.02 s, 20 checks at 0.5 mm/dt=0.01 s (controls and all three vertices in all scenarios), then 10 checks at 0.5 mm/dt=0.005 s (the same five schedules, high loading, both SAR readings).

At high loading/clamped SAR, the time-refined whole-specimen outcomes are below. Floors are °C/min; endpoints are °C. `terminal` tests every specimen cell against [−35,0] °C at 70 s. An all-time pass additionally requires ≥50 °C/min.

| Schedule | All-time | −123 to 0 | −123 to −35 | −100 to −40 | Final cold / hot | Terminal / all-time pass |
|---|---:|---:|---:|---:|---|---|
| original | 44.43576 | 59.01411 | 71.47867 | 75.76209 | -32.997 / -12.846 | True / False |
| flat | 50.31960 | 56.41623 | 68.83558 | 73.05122 | -32.670 / -13.000 | True / True |
| deficit_start | 43.51900 | 58.96249 | 71.56606 | 75.91616 | -33.022 / -12.835 | True / False |
| deficit_mid | 52.57649 | 59.47610 | 70.64143 | 74.28940 | -32.775 / -12.939 | True / True |
| deficit_end | 47.69477 | 47.69477 | 62.70649 | 68.02605 | -32.114 / -13.290 | True / False |

The middle-deficit schedule improves the high/clamped all-time floor by 8.14073 °C/min relative to the released ramp and the wider post-glass floor by 0.46199. It nevertheless lowers the two narrower floors by 0.83724 and 1.47269 °C/min. This retains the informative secondary tradeoff; it does not support the stronger hypothesis that a tradeoff is unavoidable for every post-glass support. The middle-deficit schedule was a prespecified vertex, not a schedule fitted after scoring. Its relative success is conditional and is not promoted to a robust protocol.

Across fixed scenarios, the 0.5-mm checks give the following wider-primary comparisons (original → middle deficit). Use dt=0.005 s for high-loading rows and 0.01 s for nominal rows. All selected schedules pass terminal bounds at this mesh, but every nominal/clamped schedule fails the all-time 50 °C/min floor.

| Fe mg/mL | SAR | Whole all-time | Whole −123 to 0 | Myocardial all-time | Myocardial −123 to 0 |
|---:|---|---|---|---|---|
| 1.47 | clamped | 35.6943 → 42.2399 | 58.9193 → 59.3880 | 35.6943 → 42.2399 | 59.5631 → 59.9760 |
| 1.47 | extrapolated | 52.3148 → 55.5610 | 55.1085 → 55.5610 | 52.3148 → 56.0674 | 55.6649 → 56.0674 |
| 1.83 | clamped | 44.4358 → 52.5765 | 59.0141 → 59.4761 | 44.4358 → 52.5765 | 59.5914 → 60.0042 |
| 1.83 | extrapolated | 55.2037 → 55.6495 | 55.2037 → 55.6495 | 55.6966 → 56.0985 | 55.6966 → 56.0985 |

In high/clamped middle deficit, the all-time minimum is myocardial at t=0.1125 s, −149.9014 °C, coordinate (11.75, 5.75, 4.25) mm. The limiting cell reaches −123 °C only at 13.2164 s. Thus its all-time objective rewards storage-temperature power. Its myocardial wider-post-glass floor is 60.00425, while the whole-specimen floor is 59.47610 at 70 s in surrounding CPA, −13.1381 °C, coordinate (19.75, 8.75, 6.25) mm. The narrower floors occur at the upper temperature crossing: −35 °C at 67.3610 s, and −40 °C at 63.2193 s, in myocardium. Endpoints differ by domain and are all retained in the tables.

At 0.75 mm/high-clamped, the raw all-time maximum instead belongs to edge_mid_end (52.59128), whose minimum occurs near −149.92 °C; middle deficit gives 52.58544. This 0.00585 difference is below the predefined 0.02 tolerance, while edge_mid_end loses 5.85204 °C/min in the wider post-glass floor. Maximizing the raw all-time scalar can select an essentially tied storage-temperature result with poorer later performance. That edge midpoint was not among the predefined 0.5-mm subset, so its tiny advantage is not a refined ranking. At nominal/clamped, the raw all-time maximum belongs to deficit_end, also near storage, but all discovery endpoints fail. In extrapolated scenarios the primary winner is limited at 70 s in warm CPA instead. [Objective locations and choices](results/objective_choices.csv) distinguish these cases.

The narrower conditional discovery frontier comprises original, deficit_start, edge_start_mid and deficit_mid for each of the three terminal-eligible scenarios. At 0.5 mm the selected five-schedule subset retains original, deficit_start and deficit_mid on those secondary frontiers, including nominal/clamped once terminal constraints pass. These subset frontiers do not exclude dominance by untested intermediate schedules. [Primary nondominance](results/nondominance.csv) and [secondary nondominance](results/secondary_nondominance.csv) include both domains, all phases and terminal eligibility.

Equal command does not mean equal energy: at 0.75 mm the absorbed-source ranges across the family are 681.97501–682.13516 J (nominal/clamped), 712.86669–712.93166 J (nominal/extrapolated), 687.01250–687.16194 J (high/clamped), and 717.80150–717.86036 J (high/extrapolated). At time-refined high/clamped, original, flat and middle deficit absorb 686.82728, 686.91784 and 686.95187 J. Temperature-dependent SAR makes those unequal despite identical integrated command; equal iron is fixed only within each mesh/loading comparison.

Numerical checks: maximum full-field relative enthalpy/source mismatch is 5.15691×10⁻⁷. Independent scalar quadrature agrees within 0.000237 °C, and completed-crossing transit rates within 1.27×10⁻⁷ °C/min; connected uniform cells agree exactly and zero source preserves −150 °C. Actual constitutive C(−150)=1,442,177.75 J/m³/K is used for the onset unit identity q=1000 CFe SAR u and conversion from K/s to °C/min. Window minima include linearly interpolated rate at local support boundaries rather than dropping crossing-edge values. All 66 runs complete −100 and −40 °C crossings; every run intentionally leaves some 0 °C crossings incomplete, since the endpoints are below zero. The wider floor is on the observed 70-s trajectory inside the temperature support, not a complete −123-to-0 °C transit.

Time halving changes any checked rate by at most 0.000260751 °C/min and any endpoint by 0.000018802 °C, with no terminal flips. Changing 0.75 to 0.5 mm changes checked rates by up to 1.21793 °C/min and endpoints by up to 1.79139 °C across both domains; all five nominal/clamped terminal classifications flip in both domains. Whole-specimen nominal/clamped cold endpoints move from approximately −36.55…−35.46 to −34.76…−33.77 °C. No continuum certificate follows from two meshes. [Refinement table](results/refinement.csv) records every difference. Two check-development corrections, a noncrossing scalar fixture and rounded analytic metadata, are explicitly recorded in refusals.md; neither changed a scored schedule or trajectory.

Full offline verification: `PYTHON=python3 ./reproduce.sh`. It uses pinned local inputs, a C++17 compiler and Python NumPy/SciPy, with no model/network calls. The initial computation completed all 66 trajectories; after correcting analytic metadata, the verifier and scalar checks pass and a fresh selected trajectory confirms the central result. See [verification](results/verification.json), [checks](results/numerical_checks.json), [summary](results/summary.csv), and [minimum locations](results/minima.csv). Source download instructions and exact versions accompany the snapshots.

Limitations: synthetic voxelized phantom; uniform assumed Fe; adiabatic boundary; two uncertain SAR extensions; no empirical local temperature validation for these schedules; no constitutive uncertainty sweep, anatomical validation, mechanical re-evaluation, ice kinetics, stress, fracture, toxicity or recovery calculation. Nondominance is conditional on selected supports, horizon, family, mesh and 0.02 tolerance. The new contribution is an informative finite benchmark, including a primary null about an unavoidable tradeoff, not a biological success or a claim of general schedule novelty.

**C01 status: agent-produced; awaiting coordinator verification.**
