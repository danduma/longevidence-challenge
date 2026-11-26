# C02 · Boundary exchange feasibility envelope

At a prescribed −150 °C exterior, the constant equal-duty schedule first lost feasibility because its coldest myocardial voxel finished below −35 °C. At +22 °C, the first failure occurred in the surrounding cryoprotectant solution (CPA). Its instantaneous warming rate fell below 50 °C/min at the end of the run, even while the myocardium passed all three criteria. The cold transition moved substantially between meshes. The room-temperature transition fell in the same sampled bracket on both.

## Where the constraints switched

Each interval below gives a **passing lower coefficient and a failing upper coefficient** for the indicated constraint, in W/(m² K). “Cold” and “hot” refer to terminal whole-specimen limits of −35 and 0 °C. “Rate” refers to the 50 °C/min instantaneous floor across the whole 70 s trajectory. Additional switches show when a second failure appeared as the coefficient increased. The released ramp failed the rate constraint in every scenario, including the adiabatic baseline. Its released heating command at onset gives only 44.43576 °C/min in myocardium, before any external heat can propagate there.

| Ambient °C | Mesh mm | Constant equal duty, first failure | Additional constant-duty failure | Released ramp, additional terminal failure |
|---:|---:|---|---|---|
| −150 | 0.75 | Cold, 0–2 | Rate, 6.25–7.5 | Cold, 0–2 |
| −150 | 0.50 | Cold, 3.5–5 | Rate, 6.25–7.5 | Cold, 3.5–5 |
| −80 | 0.75 | Rate, 7.5–8.75 | Cold, 28.75–30 | Cold, 27.5–28.75 |
| −20 | 0.75 | Rate, 10–11.25 | Hot, 52.5–53.75 | Hot, 47.5–48.75 |
| +22 | 0.75 | Rate, 12.5–13.75 | Hot, 15–16.25 | Hot, 15–16.25 |
| +22 | 0.50 | Rate, 12.5–13.75 | Hot, 15–16.25 | Hot, 15–16.25 |

The −150 °C flat baseline had 0.616 °C of terminal cold clearance on the 0.75 mm mesh and 2.330 °C on the 0.50 mm mesh. At h=5 the fine mesh missed the threshold by just 0.00360 °C. The disjoint cold brackets show why a transferable h limit cannot be inferred even though both meshes identify the same limiting tissue and constraint. Exposed staircase areas were similar, 1559.25 and 1558.50 mm² respectively. Voxelization also changes the internal geometry and resolves the wall differently.

At +22 °C and h=13.75 on the fine mesh, the specimen rate minimum occurred in the CPA at 70 s, at a local temperature of −2.44 °C. The myocardial minimum was 50.319 °C/min. By h=16.25 the warmest CPA reached +0.166 °C while the warmest myocardium was still −4.200 °C. Increasing the net heat supplied can therefore improve the cold endpoint and still violate the **instantaneous** rate floor. The heat source depends on temperature and the rate includes local conductive flux. The model does not identify a biological event at this rate crossing.

## Compartment endpoints and energy

Selected 0.50 mm constant-duty results at dt=0.01 s. Heat in joules is integrated over the bag, including the conducting wall. Positive boundary energy denotes net heat entering. Temperatures are terminal extrema in °C and rates are minima in °C/min.

| Ambient °C | h | Specimen rate | Specimen cold / hot | Myocardial rate | Myocardial cold / hot | Source J | Boundary J |
|---:|---:|---:|---|---:|---|---:|---:|
| −150 | 0 | 50.319 | −32.67 / −13.00 | 50.319 | −32.67 / −14.35 | 686.92 | 0 |
| −150 | 3.5 | 50.319 | −34.32 / −15.35 | 50.319 | −34.32 / −16.05 | 692.79 | −27.04 |
| −150 | 5 | 50.319 | −35.004 / −16.24 | 50.319 | −35.004 / −16.74 | 695.18 | −37.92 |
| +22 | 12.5 | 50.319 | −24.89 / −2.41 | 50.319 | −24.89 / −6.06 | 650.12 | +108.80 |
| +22 | 13.75 | 49.651 | −24.28 / −1.53 | 50.319 | −24.28 / −5.41 | 647.04 | +117.51 |
| +22 | 15 | 48.886 | −23.69 / −0.67 | 50.319 | −23.69 / −4.79 | 644.05 | +125.90 |

Gao et al. modeled free convection around a cryobag with h=15 and ambient 22 °C. At those two inputs the C02 phantom passes the three myocardial criteria but fails the whole-specimen rate screen. Gao's modeled myocardial and surrounding iron concentrations differed from the C02 design. The coefficient was prescribed in the primary study. Applying it here supplies a reference scenario without measuring the environment or reproducing the experiment. Joshi & Rabin also considered external boundary conditions in their heart models. The C02 contribution is the set of schedule-specific switch intervals and their mesh and compartment dependence. [Gao et al.](https://doi.org/10.1002/admt.202100873) · [Joshi & Rabin](https://doi.org/10.1371/journal.pone.0290063)

## Calculation and verification

The locally frozen protocol specifies the grid, search rule, controls and falsifiers. The benchmark uses uniform 1.83 mg Fe/mL muscle loading, 10 mg Fe/mL chamber loading, the released radial margin shell and low-temperature clamped SAR. Both schedules integrate to 66.97622646 duty-seconds. The external term contributes h × face area × (ambient − wall cell temperature) watts on each exposed Cartesian face. All such faces belonged to wall cells on both meshes. Interior faces exchange equal and opposite power. Temperature extrema and rate minima include the chamber contents, solution and muscle. Wall heat capacity contributes to total enthalpy. Labels 0, 1, 2 and 3 in the CSV denote wall, surrounding CPA, myocardium and chamber contents respectively. [Frozen protocol](protocol.md) · [Model source and versions](sources.md)

Before scoring new scenarios, the paired 0.50 mm independent adiabatic baselines reproduced the published ramp and flat rates of 44.43576 and 50.31942 °C/min and all four terminal extrema within 0.00001 °C. The numerical source uses 1000 g/m³ per mg/mL and converts K/s to °C/min by multiplying by 60. The analytical myocardial onset calculation for the ramp agrees within 0.000001 °C/min. The full table contains 137 deterministic integrations including controls, bisections and time checks. Ambient labels at h=0 describe identical physical conditions. The rows have no statistical replication interpretation.

Halving the time step at both ends of the first flat transition for the cold and room extremes preserved all eight decisions. The largest changes across these pairs were 0.000562 °C/min in minimum rate and 0.000050 °C in either terminal extreme. Relative enthalpy imbalance across all runs peaked at 1.78 × 10⁻⁶ of net energy input. All h=0 boundary integrals vanished. A one-cell check with six exposed faces, no iron and h=10 reproduced the independent midpoint boundary-energy calculation (2.579893 µJ) and final temperature (−149.985689 °C). [Boundary-unit fixture](unit_check.py) The implementation interpolates first upward crossings of −100 and −40 °C. All completed crossings were ordered. Fourteen cold scenarios had incomplete crossings, including 15,303 specimen voxels in the fine-grid flat run at −150 °C and h=80. Any transit statistic for such a row covers only completing voxels. Every sampled specimen minimum rate was positive. [All scenario and energy results](results/envelope.csv) · [All component and joint switch intervals](results/brackets.csv)

The model uses an analytic phantom, a fixed initial field, a constant far-field temperature and a uniform effective coefficient on a voxelized exterior. Material correlations and the source law are continued as coded, including in already-failing warm scenarios above 0 °C. Rate checks sample RK2 stages and the final RHS. The finite scenario set cannot exclude unsampled transitions, estimate physical boundary conditions or give a continuum threshold. The intermediate ambients have discovery-grid results only. No ice, stress, fracture or recovery outcomes were calculated. This work awaits coordinator verification and has received no human expert review.

Run the offline verifier with `python3 -B verify.py` in this directory. It checks hashes, reproduces the baseline before scoring, reruns the entire envelope, and compares the two numerical tables. [Executable verifier](verify.py) · [Source and download instructions](sources.md) · [Access record](refusals.md)

Execution note: one empty scratch file was briefly created in the system temporary directory and removed. All research inputs, computations and saved deliverables stayed in the assigned directory or used read-only sources.
