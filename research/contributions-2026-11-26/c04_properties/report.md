# C04 · How much constitutive change erases the high-loading thermal pass?

**A 0.8% decrease in SAR at unchanged thermal capacity removed the constant schedule’s all-time pass. A 0.4% decrease retained it.** At 0.5 mm and a 0.01 s step the two minima were 49.917 and 50.118 °C/min, respectively. In both cases the terminal temperature bounds passed. By contrast, the ramp maintained a higher whole-specimen minimum *after the glass transition* at every point of the 3 × 3 discovery grid. Its advantage ranged from 2.554 to 2.6475 °C/min. The calculations concern the pinned high-loading phantom and the specified perturbations. There has been no human expert review of C04.

## What changed

The accepted benchmark supplies the shell, ramp, equal-integrated-duty constant schedule, high uniform myocardial iron (1.83 mg Fe/mL), and the two diagnostic supports. C04's incremental result is the joint source/storage sensitivity and its numerical pass bracket. The test retained the 70 s horizon from −150 °C, the adiabatic exterior and the published clamped SAR shape. Both schedules commanded 66.97622646 duty-seconds. The prior all-time and window rankings, delivery analyses and schedule designs are background.

The dimensionless source multiplier s scales SAR in W/g Fe. Multiplication by 1000 times iron in mg/mL and local duty yields W/m³. The independent multiplier c scales volumetric heat capacity in J/(m³ K) and the associated enthalpy of the three specimen regions. The conducting polyethylene wall retains its original properties. The heating equation divides net deposited plus conducted power density by capacity. Consequently the uniform-temperature initial source term scales as s/c, while diffusivity within the specimen scales as 1/c at fixed conductivity. The discrete solver checks positive conductivity and capacity at every rate evaluation. Primary studies supply the nominal temperature dependence. The common multipliers are hypothetical. [Sources, versions and download instructions](sources.md) · [Frozen protocol](protocol.md)

## Decision envelope

Each discovery cell below gives **constant all-time minimum / ramp advantage after glass**, in °C/min, for the whole specimen. Bold denotes a constant-schedule pass on rate and both terminal bounds. The ramp failed its all-time rate criterion at all nine settings. A post-glass rate is the minimum sampled local derivative while the voxel is between −123 and 0 °C during the 70 s simulation.

| Capacity c ↓ / SAR s → | 0.95 | 1.00 | 1.05 |
|---|---:|---:|---:|
| 0.95 | **50.326** / +2.604 | **52.973** / +2.580 | **55.621** / +2.554 |
| 1.00 | 47.809 / +2.628 | **50.325** / +2.608 | **52.840** / +2.585 |
| 1.05 | 45.533 / +2.647 | 47.929 / +2.631 | **50.324** / +2.612 |

Equal source and capacity multipliers nearly preserved the constant schedule's minimum rate, yet its coldest endpoint moved from −34.274 °C at (0.95, 0.95) to −34.518 °C at (1.05, 1.05). The hottest endpoint moved from −13.679 to −12.497 °C. Dividing the source by capacity alone misses this change in the temperature field.

The fixed ratio probes supplement the 3 × 3 grid. At c = 0.95 they extend s slightly below the grid edge, to 0.9424 and 0.9462. At each prescribed c, neighboring flat-schedule trials at s/c = 0.992 and 0.996 gave fail and pass, respectively. The high-capacity low-ratio trial also missed the cold limit on the discovery mesh (−35.065 °C). The c = 1 bracket held with both time steps on the finer grid. An independent lower bound follows at the uniform start, when face flux vanishes. The myocardial rate is 60,000 × 1.83 × 691 × 0.956803… × (s/c) / 1,442,178 °C/min. Thus a necessary condition is s/c ≥ 0.9933122. The full-field observations bracket the sampled transition without locating an exact continuous threshold. At stronger heating the limiting voxel can change. At the fine-grid high-source/low-storage corner the flat schedule's minimum occurred in surrounding solution at 70 s, whereas at the nominal point it occurred in myocardium at 0.095 s. [All brackets](results/flat_brackets.csv) · [Full fields summarized by case](results/cases.csv)

## Whole specimen and myocardium

Whole-specimen extrema include chamber, myocardium and surrounding CPA. The wall participates in heat exchange and storage. At 0.5 mm and 0.01 s the three paired checks gave the following results. Parentheses contain myocardial-only values. Cold and hot endpoints are for the constant schedule.

| SAR s, capacity c | Flat all-time min, °C/min | Flat coldest, °C | Flat hottest, °C | Ramp − flat post-glass min, °C/min |
|---|---:|---:|---:|---:|
| 0.95, 1.05 | 45.528 (45.528) | -39.686 (-39.686) | -19.488 (-20.947) | +2.638 (+2.669) |
| 1.0, 1.0 | 50.319 (50.319) | -32.670 (-32.670) | -13.000 (-14.353) | +2.598 (+2.629) |
| 1.05, 0.95 | 55.604 (55.615) | -25.652 (-25.652) | -6.639 (-7.879) | +2.544 (+2.575) |

At 0.75 mm the myocardial post-glass advantage ranged from 2.638 to 2.726 °C/min. The finite grid supports positive paired differences at the sampled values only. The nine combinations are deterministic cases, with no distribution or population frequency attached. [Paired results](results/pairs.csv)

## Numerical checks and scope

- Both nominal 0.5 mm runs first recovered the archived RK2 all-time and −100 to −40 °C rates, cold/hot endpoints and absorbed Joules to within 1e-7 in the respective units. The post-glass minima at the nominal point were 59.014 and 56.416 °C/min for ramp and constant. These baseline recoveries are controls.
- Halving the 0.5 mm step from 0.02 to 0.01 s changed the checked whole-specimen all-time minima by at most 0.000458 °C/min, post-glass minima by 0.000051 °C/min and terminal extrema by 0.000034 °C. Paired ranking signs and pass decisions agreed across checked meshes and time steps. Spatial errors mattered more for absolute endpoints. At the weak-heating corner, the constant cold spot changed from −41.318 °C at 0.75 mm to −39.686 °C at 0.5 mm.
- At that same corner 40 ramp and 30 flat voxels lacked a completed −123 to −40 °C transit at 0.75 mm. At 0.5 mm the corresponding counts were 2 and 0 at both steps. All other cases completed both crossings. Every specimen minimum rate was positive. The transit column gives the slowest average among completers only. The post-glass statistic includes all attained in-window states, including those of incomplete voxels. None of the cases warmed the entire specimen to 0 °C.
- A single adiabatic cell with scaled, temperature-dependent capacity reproduced the constant-source scalar integral after 15 s to 1.3 × 10⁻⁹ °C. With 1 mm³ volume it absorbed 0.01801955 J. An independent quadrature through the capacity joins and SAR knee gave a −123 to −40 °C transit rate of 28.609027 °C/min, compared with 28.609074 by interpolated RK2 crossings. Across all 40 full-field runs the largest absolute relative discrepancy between integrated source and the independently integrated wall-plus-specimen enthalpy rise was 4.16 × 10⁻⁷. [Numerical checks](results/checks.json)

The calculation addresses thermal decisions for one synthetic geometry and model of the material properties. It leaves the uncertainty of actual myocardial SAR and heat capacity unresolved. No experimental evidence here connects the schedule rankings to ice formation, stress, fracture or recovery. Two mesh spacings and selected time refinements also leave the continuum pass boundary unresolved.

## Verification

From the assigned C04 directory run the command below. It checks snapshot hashes, reproduces the nominal comparisons before scoring the new ones, rebuilds the four result files and checks the reported sensitivity bracket and ranking. [Verifier](verify.py) · [Source/download details](sources.md) · [Access record](refusals.md)

```sh
python3 -B verify.py
```
