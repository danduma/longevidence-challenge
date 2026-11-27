# Does the startup interface effect change tensile-stress ordering?

All three schedules pass the event-based thermal constraints. Switching from refined voxels to cut cells swaps the graded and uniform strict schedules under extrapolated SAR. The uniform release-policy schedule retains the lowest stress under both tensile measures.

The accepted interface audit found that integrated cut cells and refined voxels removed a first-second rate penalty seen in a shifted 1.1 mm voxel grid. The earlier coupled study reproduced the Vitrify objective and extended it to a temperature event on voxel grids for a matched-iron loading ablation. This experiment compares the same three schedules from the interface audit under both thermal representations and holds the mechanics grid fixed. [Pinned Vitrify](https://github.com/Marc-Dvci/Vitrify/tree/1d498ee5c583e0112a762e60892d576ec2e27f57).

## Released comparator and methods

The first run reconstructed the original center-labelled 1.4 mm voxel fields from the released parameters. The six-frame 12 s hard peaks were 1.233836819 and 1.175765415 MPa for the baseline and optimized design. Both peaked at 9.042204967 s. Smooth objectives were 1.261793902 and 1.199053433 MPa. The independent solver uses structured Q1 hexahedral equilibrium and banded Cholesky. It starts at −150 °C with zero stress and strain. [Released hard peaks](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/results/objective_horizon.json) · [Released optimizer results](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/results/optimization_extrapolated.json).

Both thermal grids have 0.75 mm cells and the half-y phase. They share the continuous bag shape and physical loading origin from the 1.1 mm fit. Voxels use center labels and full faces. Cut cells integrate phase volumes, centroids, face apertures and embedded interfaces by ray quadrature. Both use reciprocal two-point conductance and an insulated bag. Each schedule retains its own loading and duty knots, with linear interpolation at 0, 35 and 70 s and the last value held thereafter.

The transport update is implicit with conductivity, heat capacity and SAR evaluated at the start of the step and duty at the midpoint. The main comparison uses 0.10 s thermal steps, 0.60 s mechanics intervals and quadrature order 16. The numerical challenge uses 0.05 s, 0.30 s and order 24. Each endpoint is interpolated to the latest local −35 °C heart crossing in the final bracket, including the filled chambers in the heart mask.

The mechanical labels, nodes, pins and Maxwell parameters are fixed on the native 1.4 mm mesh. Each active element receives the nearest same-phase temperature or an inverse-square distance average of eight same-phase temperatures. Both transfers preserve a spatial constant. The objective uses the Vitrify log-sum-exp scale of 20 kPa over all three principal stresses in every specimen element, sampled from the first positive mechanics time through the exact thermal event. The hard peak uses those same samples. The smooth maximum is unnormalized, so compare its values within cadence and its pairwise differences across cadences.

The feasibility screen applies 50 °C/min after each specimen state first crosses −123 °C through the heart event and checks all active states against 0 °C at the endpoint. The code reports the hottest specimen state. Rates are slopes of linear steps and each slope contributes only during the eligible part of its step. The upstream optimizer includes startup in its rate criterion. The two SAR variants share the curve above −80 °C and take values of 691 and 1013 W/g Fe at −150 °C for constant and linear extension respectively. [Primary thermal methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/) · [Pinned material laws](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/src/vitrify/properties.py).

## Precedent and source scope

Joshi et al. used mesh-converged tetrahedral thermal models and described the critical rate from glass transition to termination. Joshi and Rabin modeled temperature-dependent stress, cooling-history residual stress and early warming peaks. Vitrify already reports a fixed-horizon objective, spatial validation and SAR ambiguity. Vispute and Rabin discuss mesh and specimen packaging in an artery geometry. Crockett, Colella and Graves provide precedent for conservative embedded boundaries with material jumps. Biostasis identifies warming around glass transition as a stress concern. CryoNet discusses region-specific warming. The cached GitHub issue and pull request lists were empty at retrieval. The reviewed evidence leaves the matched-interface stress comparison open. [2022 thermal paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/) · [2023 mechanics paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0290063) · [2026 artery model](https://doi.org/10.1371/journal.pone.0350078) · [Embedded-boundary paper](https://escholarship.org/content/qt00q9z6st/qt00q9z6st.pdf?t=lnrl2o) · [Biostasis](https://www.biostasis.com/new-warming-breakthrough-for-cryopreserved-organs/) · [CryoNet](https://www.cryonet.org/cgi-bin/dsp.cgi?msg=18925).

The 2022 thermal citation uses the predecessor's cached excerpt and reading notes. The PMC browser returned a challenge in this attempt. The mechanics and artery snapshots preserve publisher attribution and public license notices. The geometry and mechanics implementations come from the two accepted precursor jobs. The upstream snapshot preserves Marc Davici's Apache-2.0 release and committed input data.

## Results

The table uses the finer cadence and cut-cell quadrature challenge. Independent mask checks show that myocardium determines the endpoint in all 24 thermal trajectories. The chambers are already above −35 °C at each event. [Mask checks](myocardium-check.log).

| SAR | Thermal mesh | Schedule | Event (s) | Min rate (°C/min) | Hot (°C) | Smooth nearest (MPa) | Smooth IDW (MPa) |
|---|---|---|---:|---:|---:|---:|---:|
| extrapolated | voxel | Graded | 66.573 | 58.283 | -12.535 | 1.283722 | 1.227192 |
| extrapolated | voxel | Uniform release policy | 66.550 | 58.056 | -12.835 | 1.280298 | 1.225000 |
| extrapolated | voxel | Uniform strict | 66.539 | 57.998 | -12.879 | 1.282060 | 1.226982 |
| clamped | voxel | Graded | 71.242 | 57.842 | -12.178 | 1.106155 | 1.054148 |
| clamped | voxel | Uniform release policy | 71.240 | 57.566 | -12.499 | 1.102953 | 1.052476 |
| clamped | voxel | Uniform strict | 71.234 | 57.514 | -12.544 | 1.104332 | 1.054148 |
| extrapolated | cut | Graded | 65.859 | 58.848 | -13.333 | 1.293322 | 1.197907 |
| extrapolated | cut | Uniform release policy | 65.834 | 58.625 | -13.617 | 1.292442 | 1.196989 |
| extrapolated | cut | Uniform strict | 65.822 | 58.565 | -13.661 | 1.294464 | 1.198824 |
| clamped | cut | Graded | 70.477 | 58.459 | -13.010 | 1.121538 | 1.039010 |
| clamped | cut | Uniform release policy | 70.473 | 58.170 | -13.311 | 1.120822 | 1.038350 |
| clamped | cut | Uniform strict | 70.466 | 58.117 | -13.357 | 1.122493 | 1.039874 |

The smooth graded-minus-strict difference under extrapolated SAR changes from +0.001661 to −0.001142 MPa for nearest transfer and from +0.000210 to −0.000917 MPa for eight-state transfer. The coarser cadence gives the same signs. The hard peak does not give a universal secondary reversal. With clamped SAR and eight-state transfer, the two voxel smooth objectives differ by only 0.000000215 MPa in the finer run. They are effectively tied at the resolution of this comparison.

Across the tested configurations, the winning schedule leads by 0.000588–0.001982 MPa on the smooth objective. Changing the transfer shifts individual objectives by 0.050–0.101 MPa, and clamping SAR lowers them by 0.159–0.178 MPa. The direction of the absolute cut-minus-voxel stress shift itself depends on transfer. The nearest transfer raises the score on cut cells, while eight-state transfer lowers it. These are much larger effects than the schedule gaps. The stable winner is a benchmark result for the tested numerical choices. A physical claim of superiority or an extrapolation of the secondary order would exceed the model resolution.

The cut-cell representation reaches the event 0.708–0.768 s earlier than its matched voxel run. It raises the screened minimum rate by 0.565–0.664 °C/min and lowers the endpoint hot spot by 0.776–0.832 °C. Clamped SAR delays the event by 4.618–4.704 s. The hottest state is in the surrounding CPA. The graded runs and voxel uniform runs place it at (19.875, 9.000, 6.375) mm, while the cut-cell uniform runs put it on the opposite long-axis side at (3.375, 9.000, 6.375) mm. These are discrete argmax locations. This comparison does not establish a persistent side preference in the continuous geometry. The coldest-heart event and post-transition rate margins explain why the earlier startup floor reversal does not change the feasible set here.

## Reproduction and limits

Run **sh verify.sh** from this directory. It uses an existing NumPy/SciPy environment or creates a local virtual environment and installs the requirements there. The verifier recomputes the pinned baseline, rebuilds both thermal graphs at each quadrature setting, solves fresh histories and reconstructs the local rate masks and extrema directly from temperatures. It checks the released values, archived numerical outputs, all-pass classification, winner under both tensile scores, and the extrapolated smooth-objective reversal under both transfers. The source code, complete score and event records, and fixed candidate parameters accompany the report.

The two temperature transfers preserve phases on a fixed voxel mechanics grid but do not resolve mechanical geometry or provide a conservative overlap reconstruction. The comparison challenges temporal sampling and cut-cell integration on one matched thermal spacing and phase. It does not establish continuum convergence or bound loading, material, or initial-state uncertainty. The stress-free initial condition is the release assumption. The primary mechanics study includes residual stress from prior cooling. Preserve the uniform release-policy candidate as the leading benchmark comparator and treat the two runners-up as unresolved for physical selection.
