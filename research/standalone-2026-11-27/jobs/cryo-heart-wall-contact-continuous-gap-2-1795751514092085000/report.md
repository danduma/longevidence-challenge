# A heated CPA gap removes the coarse startup dip in the Vitrify phantom

**Conditional numerical inference.** With the same 1.1 mm lattice, phantom and three schedule vectors, adding a thermal state for the CPA missed between adjacent myocardium and wall samples removed the startup-rate depression. Merely correcting the distance-weighted resistance left the depression in place. For the strict uniform witness, the original time step showed two failing phases. At a common 2 ms step only half-y remained below the 51.5 °C/min audit threshold, and inserting the heated gap cleared it.

## Matched comparison

| Phase | Original voxel | Distance-weighted series | Heated gap |
| --- | ---: | ---: | ---: |
| 1.1 mm native | 51.5429 | 51.5305 | 52.5978 |
| 1.1 mm half-x | 51.4010 | 51.3853 | 52.5977 |
| 1.1 mm half-y | 50.9745 | 50.9556 | 52.5975 |
| 1.1 mm half-z | 52.5941 | 52.5941 | 52.5941 |
| 0.75 mm half-y | 52.5933 | 52.5933 | 52.5933 |

Entries are the strict witness's minimum specimen warming rates through the first time level at or beyond 70 s, in °C/min. The heated-gap screen includes the added CPA nodes. The other two schedules gain 1.046–1.611 °C/min in the three affected phases. Across all three schedules, the affected heated-gap minima are within 0.005 °C/min of the 0.75 mm half-y baseline, which has its own native pillow fit and loading origin. The controls have zero myocardium–wall links and produce identical histories across treatments. Every reported minimum is below −149.6 °C. The comparison uses a model audit floor, outside the glass-transition temperature range used for interpreting devitrification risk in the primary rat-heart study. [Primary thermal study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/).

The independent startup solve uses the same step across all treatments. It separates the spatial result from the time-step sensitivity of the original calculation.

| Phase | Original voxel | Distance-weighted series | Heated gap |
| --- | ---: | ---: | ---: |
| 1.1 mm native | 51.7866 | 51.7747 | 52.6048 |
| 1.1 mm half-x | 51.6481 | 51.6345 | 52.6048 |
| 1.1 mm half-y | 51.2139 | 51.1957 | 52.6053 |
| 1.1 mm half-z | 52.6048 | 52.6048 | 52.6048 |
| 0.75 mm half-y | 52.5992 | 52.5992 | 52.5992 |

These are first-second minima for the strict witness using a fully explicit 2 ms step. Halving it changes each value by less than 0.01 °C/min. The half-y result is the threshold reversal that survives time refinement. This independent refinement covers startup, while the full-horizon checks use the principal solver.

## Does the gap have to generate heat?

At 2 ms, disabling heat deposition on the added CPA state while keeping its heat capacity and both interface conductances lowered the half-y minimum on the **original specimen nodes** to **49.372 °C/min**, versus **52.605 °C/min** with gap heating. The half-y values at 1 ms were **49.378** and **52.605 °C/min**. Measuring only the original nodes in both arms avoids automatically assigning zero startup rate to the unheated state. Native and half-x also showed larger dips with the gap unheated, while the two phases without these links were identical in both arms. This intervention isolates the effect of the gap's modeled power term conditional on the new geometry and storage. It leaves total deposited power free to differ by that term. The comparison does not infer how much heat actual CPA deposits in this temperature range.

## What changed in the stencil

The pinned geometry gives an ellipsoid for the heart and a superellipsoid for the inner wall. For each coarse myocardium–wall neighbour pair, intersect the centreline with those surfaces. The three affected lattices each have ten such links, with positive CPA intervals spanning 0.810–1.000 mm over all affected links. The native fitted shape, sampling coordinates, original loading coordinates, extrapolated SAR and duty vectors stay fixed within each comparison. The original coarse solver agrees with the previous phase study. [Geometry](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/src/vitrify/geometry.py) · [Thermal stencil](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/src/vitrify/thermal.py).

The resistance-only control places the change in conductivity at the inner wall crossing and gives the intervening CPA the specimen conductivity law. For the heated-gap case, insert a CPA state at the interval midpoint. Use segment lengths and reciprocal conductivities on either side, with one equal-and-opposite flux per link. Place CPA capacity and scheduled volumetric heating on the added state at its physical location. Use one third of the axial gap fraction times the voxel volume as the storage quadrature, debiting the heart and wall volumes by the corresponding overlaps with their half segments. This directional allocation assigns one sixth of a voxel to each incident half link and keeps total modeled volume fixed. Conductance uses the original face area. The original endpoint and wall states remain in the graph. Both solvers use the precursor’s adiabatic exterior and uniform −150 °C initial condition. The full-horizon solver treats the added state’s incident fluxes implicitly within each otherwise explicit step and applies the same exchanged energy to the endpoint updates. A zero-length CPA interval retains a direct specimen–wall link, although no such links occurred in the tested phases. The wall and specimen power balance is conserved across every internal face, to roundoff in the discrete update.

The resistance-only result leaves the half-y failure in place at both time resolutions. The heated-gap intervention changes sampled heat deposition, storage and the timing of exchange together. The heat-source ablation holds storage and conductance fixed and shows that gap heat deposition is needed for the reversal in this discretization. The 70 s temperature window and the interpolated −35 °C heart event checks retain their original pass decisions. The point of the change is to avoid treating the 1.1 mm half-y startup failure as a schedule defect when a wall-adjacent centreline skips heated CPA.

## Prior work and scope

The coordinator's preceding phase study already documented the startup minima, measured selected centreline gaps and found the fine-grid pass. The first attempt computed this comparison but failed the coordinator check at import time because NumPy was absent. This attempt repeats the calculations with a job-local dependency setup and adds the unheated-gap control. The experiment holds the coarse geometry and schedules fixed and changes the affected links, with resistance-only and no-link controls. The public Vitrify spacing sweep and the primary heart study established the importance of spatial modelling. The cached repository issue and pull request listings contain no entries. Biostasis discusses uniform warming near the glass transition, and a CryoNet message discusses spatially selective warming. Neither excerpt contains this stencil comparison. [Released mesh experiment](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/results/grid_convergence_margin.json) · [Biostasis](https://www.biostasis.com/new-warming-breakthrough-for-cryopreserved-organs/) · [CryoNet](https://www.cryonet.org/cgi-bin/dsp.cgi?msg=18925).

Vitrify deliberately reduces early RF power for its mechanical objective. Later heart modeling also locates mechanical risk around the glass transition and makes clear that endpoint temperature spread cannot substitute for a stress history. [Vitrify README](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/README.md) · [Joshi and Rabin (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10431620/).

The added nodes resolve centreline gaps only. The face areas, directional storage allocation and untouched surrounding links are coarse approximations. This is a local conservative comparison with a finer voxel baseline. A full cut-cell or fitted-mesh solve would require interface areas, material volumes and volumetric heat quadrature from the same continuous shape. The results do not establish continuum convergence, mechanical performance, an improved schedule or biological preservation. The model uses extrapolated low-temperature SAR, as in the preceding runs.

## Reproduction

Run `sh verify.sh` with Python 3 and access to PyPI for the first run. The script installs the pinned NumPy and SciPy versions into a local `.venv` inside this job. A populated environment also supports `.venv/bin/python -B verify.py` directly. The verifier recomputes all 45 70 s runs and checks against the cached results and original phase anchors. For the decisive strict startup comparison it also reconstructs the intersections by root finding, builds its own sparse incidence operator and uses shared fully explicit 2 ms and 1 ms steps. The independent check recomputes and asserts the unheated-gap comparison at both steps using the original measurement nodes. Run `.venv/bin/python -B analysis.py > results.json` to regenerate the principal JSON. The input versions are listed in submission.json. The Vitrify code and released inputs retain Marc Davici's attribution and Apache-2.0 license under sources/upstream. The prior coordinator data and scholarly/community excerpts retain their separate attribution.
