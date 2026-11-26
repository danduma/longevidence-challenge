# C03 — Localized delivery-deficit geometry: a bounded null

**Executed result:** in the predefined compact-patch family, placement did **not** meet the protocol's materiality criterion. At 0.50 mm, moving an equal-volume, equal-deficit patch among chamber-adjacent, deep and exterior-adjacent seeds changed the terminal myocardial minimum by at most **0.03997 °C** and the longest complete −123→−35 °C transit by **0.02107 s**. The predefined thresholds were 0.10 °C or 0.10 s with refinement support. This is an informative null for these maps, not evidence that localized delivery is generally unimportant.

The patches still have a local thermal effect: at 0.50 mm their average terminal temperature is 0.476–0.528 °C below the temperature of the same voxels in the uniform control, with maximum voxel depressions of 0.699–0.746 °C. Those are descriptive secondary outputs; they do not rescue the primary placement claim. All results are agent-produced and await coordinator verification; no human expert has reviewed this work.

## What was fixed and computed

The question, map rules, comparisons, endpoints and falsification criteria were recorded in [protocol.md](protocol.md) before new scoring (SHA-256 `7198a4a22a21e9ad1edc42aba3841cbd495150160ae1521946f10a93f8bc3611`). This is a local protocol, without an external registration timestamp. The model uses the released Vitrify margin shell/ramp, 70 s exact horizon, initial −150 °C, nominal myocardial mean 1.47 mg Fe/mL, chamber 10 mg Fe/mL, clamped SAR below −80 °C and adiabatic exterior. The portable independent finite-volume midpoint RK2 solver is adapted only to expose voxel diagnostics. Thermal properties follow [Joshi and Rabin 2023](https://doi.org/10.1371/journal.pone.0290063); analytic geometry/design coefficients are pinned to [Vitrify commit 1d498ee](https://github.com/Marc-Dvci/Vitrify/tree/1d498ee5c583e0112a762e60892d576ec2e27f57). These physical and geometric assumptions are not newly validated here.

Three seeds are chosen using geometry alone in a prespecified positive-y myocardial slab: minimize distance to chamber, maximize minimum distance to chamber/surrounding CPA, or minimize distance to surrounding CPA. A compact patch consists of the N nearest myocardial voxel centers, with index tie breaking. All maps are connected under face adjacency. A common local deficit of 0.36 mg Fe/mL gives patch concentration 1.11. All remaining muscle receives a uniform compensating increment of `0.36*N/(M−N)`. Thus placements have identical deficit mass and mean; the two sizes have different deficit masses. This does not estimate real coronary delivery.

Discovery: uniform nominal plus three placements and one lower-mean control at each of two sizes, at 0.75 mm (9 cases). Refinement: uniform nominal plus the large-patch trio at 0.50 mm and the trio again with dt halved from 0.01 to 0.005 s (7 cases). Before new scoring, the archived 0.50-mm uniform-1.83 clamped-ramp baseline was reproduced. No coefficients, seed rules or amplitude were tuned after scoring.

## Endpoint results

Terminal coldest temperatures are identical between myocardium and whole specimen in these cases; the cold spot lies in myocardium. Longest transits and incomplete crossing counts are also identical between the two domains here. The code nevertheless computes these domains separately. Whole specimen excludes the wall.

| Grid mm | Map | Target patch mL | Myocardial/whole cold °C | Longest complete −123→−35 s | Incomplete crossings | Longest −100→−40 s |
|---|---|---:|---:|---:|---:|---:|
| 0.75 | Uniform | 0 | −36.517393 | 55.127279* | 33 | 38.024554 |
| 0.75 | Chamber | .010 | −36.499120 | 55.117675* | 33 | 38.023973 |
| 0.75 | Deep | .010 | −36.494710 | 55.116059* | 33 | 38.022314 |
| 0.75 | Exterior | .010 | −36.494221 | 55.115853* | 33 | 38.023104 |
| 0.75 | Chamber | .040 | −36.437602 | 55.122756* | 31 | 38.019631 |
| 0.75 | Deep | .040 | −36.424050 | 55.116951* | 31 | 38.015812 |
| 0.75 | Exterior | .040 | −36.436205 | 55.123960* | 31 | 38.023060 |
| 0.50 | Uniform | 0 | −34.730521 | 54.711906 | 0 | 37.373669 |
| 0.50 | Chamber | .040 | −34.694556 | 54.708070 | 0 | 37.371568 |
| 0.50 | Deep | .040 | −34.660844 | 54.690692 | 0 | 37.361677 |
| 0.50 | Exterior | .040 | −34.654585 | 54.686999 | 0 | 37.359601 |

`*` Longest among completed crossings only: these are censored summaries, not maximum transit bounds. At 0.75 mm, some voxels never reach −35 °C by 70 s. Consequently the post-glass transit comparison cannot satisfy the protocol's cross-grid confirmation requirement. All −100→−40 crossings are complete in every new case; that interval is motivated by [Gao et al.](https://doi.org/10.1002/admt.202100873), but no ice calculation is performed.

At 0.75 mm, terminal placement ranges are 0.00490 °C for the small patch and 0.01355 °C for the large patch. At 0.50 mm, the large-patch range is 0.03997 °C. Deep minus chamber is +0.01355 °C at 0.75 mm and +0.03371 °C at 0.50 mm. Deep minus exterior changes sign, +0.01216 to −0.00626 °C. Neither primary pair qualifies as material. The complete-window transit placement range is 0.00725 s at 0.75 mm and 0.01197 s at 0.50 mm for the large patch.

Size increases the local depression: at 0.75 mm the small patches' average temperature changes are −0.263 to −0.279 °C, versus −0.509 to −0.534 °C for large patches, relative to the same voxels under uniform loading. Larger patches also remove more iron locally; this is not a size effect at fixed deficit mass. Global cold spots warm slightly with compensated patches because compensation adds iron elsewhere. This illustrates why patch temperatures and a global minimum answer different questions.

The lower-mean controls uniformly distribute the *uncompensated* deficit: mean 1.46210 for the small patch and 1.43871 for the large at 0.75 mm. Their terminal minima are −36.55827 and −36.67930 °C, respectively. All compensated maps instead have exactly the nominal mean, shared with the uniform nominal control. The lower-mean controls are additional total-dose comparisons, not equal-total-iron placement controls.

Myocardial versus whole-specimen terminal maxima differ. At 0.50 mm, the large trio's myocardial maxima range −14.66652 to −14.66465 °C; whole-specimen maxima range −13.11432 to −13.11334 °C. At 0.75 mm, their post-glass minimum rates are 69.87885–69.88982 °C/min in myocardium versus 69.44011–69.45223 for whole specimen. The global endpoint screen is therefore not interchangeable with a tissue-specific summary. See [summary.csv](results/summary.csv) for every domain-specific endpoint.

## Matching and numerical evidence

| Grid mm | Target mL | Patch voxels | Actual mL | Deficit mg Fe | Total muscle voxels |
|---|---:|---:|---:|---:|---:|
| .75 | .010 | 24 | .010125 | .003645 | 1093 |
| .75 | .040 | 95 | .040078125 | .014428125 | 1093 |
| .50 | .040 | 320 | .040000 | .014400 | 3700 |

At each size/grid, all placements match exactly in voxel count and deficit mass. Myocardial mean is 1.47 within 1e-12 mg/mL. Total specimen iron is 0.016993881910123 g at .75 mm and approximately 0.016983686596031 g at .50 mm, unchanged among placements to <1e-12 g. Rounding changes large-patch mass across grids by 0.195%, so the spatial check is a matched construction rule, not identical continuous anatomy.

The central baseline matched six archived metrics within 1e-8. Maximum full-field relative energy imbalance was 1.23e-7. Zero source preserved −150 °C exactly. The unit check independently evaluates `60*1000*C_Fe*SAR*duty/capacity`. A one-cell check against independent scalar quadrature reproduced all four crossing times within 0.000061 s and final temperature within 0.000044 °C. Temporal refinement changes large-patch terminal minima by at most 0.00000744 °C and post-glass transit by 0.0000144 s. These errors are much smaller than the already small placement contrasts. Restricted instantaneous minima can show larger sampling changes (~0.00266 °C/min); they are secondary diagnostics.

Spatial error dominates: the uniform terminal minimum moves by 1.787 °C between grids, and the −35 °C terminal decision changes. No continuum endpoint certificate follows. At 0.50 mm all three large-patch placements clear the terminal temperature bounds; at 0.75 mm all three miss the lower bound. This is a grid effect shared by the controls, not evidence for a placement-induced terminal transition. The below-glass all-time onset floor is not used as a discovery endpoint.

## Limits and contribution boundary

The positive-y slab deliberately isolates nearby depth positions, but constrains how distinct the placements can be. Large patches overlap 71.6–76.8% at .75 mm and 64.7–91.9% at .50 mm. Small patches overlap 29.2–50.0% at .75 mm. Seeds shift with voxelization; deep and exterior labels describe seed-selection criteria, not entirely separated tissue compartments. This overlap and the limited spatial family weaken any general null claim. No whole-organ placement survey or anatomical model is implied.

C03 contributes new executed compact-patch maps, matched-mass checks, post-glass/crossing diagnostics and refined bounded-null results. Baseline reproduction, prior half-muscle perturbations, equal-iron ideas, heating schedules and rate windows are not counted as new discoveries. No measured-heart triplicate failure inference, population estimate, ice, stress, fracture or biological recovery claim is made.

Run the executable [verify.py](verify.py) with the shared Python; its default executes fresh integrations from pinned local snapshots. [sources.md](sources.md) gives licenses, versions and public download instructions; [refusals.md](refusals.md) distinguishes access failures from refusals. Key outputs: [contrasts.csv](results/contrasts.csv), [maps.csv](results/maps.csv), [local_effects.csv](results/local_effects.csv), [time_refinement.csv](results/time_refinement.csv), and [verification.json](results/verification.json). Status: **agent-produced; awaiting coordinator verification**.
