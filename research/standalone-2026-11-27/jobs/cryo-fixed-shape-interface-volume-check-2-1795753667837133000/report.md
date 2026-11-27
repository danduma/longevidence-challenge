# The half-y startup reversal with integrated interface geometry

The fixed strict schedule clears 51.5 °C/min over the first second when CPA, heart and wall are represented with geometric volumes and interfacial areas. The 1.1 mm voxel calculation falls below the same screen. This is a numerical robustness result for the synthetic Vitrify phantom.

The accepted predecessor used heated CPA states on ten centerline links with nominal face areas and one-third directional heat capacity. Its original fine-grid check refitted the bag exponent and shifted the loading origin. The accepted audit already corrected that comparison, obtaining strict minima of 52.599, 52.593 and 52.584 °C/min at 0.75, 0.55 and 0.40 mm with the shape and physical loading origin fixed. The increment here is a separate test of the interface representation. The cut-cell calculation recovers the pass using integrated phase volumes and areas throughout the bag, without any special allocation for those ten links. [Accepted predecessor audit](sources/predecessor-review.json).

## Geometry and heat balance

The upstream 1.1 mm fit gives a pillow exponent of 2.487834527647254 and an array-derived loading origin of (11.55, 8.80, 6.05) mm. The analytic bag and heart retain the upstream geometric centres, axes and chamber overlap. For the half-y phase each lattice moves by half its own cell width along y. Each comparison starts at −150 °C with an insulated outer bag. All three release-derived schedules keep their original iron loading and duty knots. Startup uses the predecessor's low-temperature SAR extrapolation. [Pinned Vitrify geometry, properties and schedule data](https://github.com/Marc-Dvci/Vitrify/tree/1d498ee5c583e0112a762e60892d576ec2e27f57).

For each cut cell, the code intersects rays with the continuous surfaces, splits the resulting intervals by material and integrates their lengths and first moments. It uses transverse midpoint quadrature in all three axes and averages the three volume estimates. Grid-face aperture is integrated by phase. At embedded interfaces a projected area element receives weight `|n_axis|` for the direction of its ray. The three directions sum to the true area integrand through `sum(n_axis²) = 1`. A reciprocal two-point resistance couples material centroids across face apertures and embedded interfaces. Every flux enters adjacent balances with opposite signs. Source and heat capacity scale with each phase volume. Graded loading is evaluated at the phase centroid, while uniform loading is constant within each phase.

An independent geometric check uses exact superellipsoid volumes, the overlap of the two chamber ellipsoids and a separate surface integral for the myocardial ellipsoid. The myocardium fits inside a copy of the convex inner bag scaled by 0.84598406. The minimum clearance is attained on the short z axis and equals 0.808584 mm. The integrated model consequently has zero direct wall–myocardium interface area.

## Startup rates

Each entry is the lowest specimen-state step increment divided by the step duration, expressed in °C/min over the first second. The main table uses 1 ms steps, with explicit diffusion for voxel and gap-state runs and backward diffusion for cut cells. The integrated calculation includes the external CPA, myocardial and chamber portions of cut cells.

| Method | Graded release | Uniform released policy | Uniform strict |
|---|---:|---:|---:|
| Voxel, 1.1 mm | 50.95362 | 50.98923 | 51.21807 |
| Centerline gap state, 1.1 mm | 52.32685 | 52.36513 | 52.60538 |
| Integrated cut cells, 1.1 mm, q=8 | 52.32689 | 52.36519 | 52.60544 |
| Integrated cut cells, 1.1 mm, q=16 | 52.32687 | 52.36517 | 52.60541 |
| Integrated cut cells, 1.1 mm, q=32 | 52.32688 | 52.36518 | 52.60542 |
| Integrated cut cells, 0.75 mm, q=16 | 52.32047 | 52.35937 | 52.59913 |
| Voxel, 0.75 mm | 52.32080 | 52.35965 | 52.59945 |
| Voxel, 0.55 mm | 52.31445 | 52.35380 | 52.59319 |
| Voxel, 0.375 mm | 52.30457 | 52.34464 | 52.58345 |

Here q is the number of transverse midpoint samples per cell dimension. At q=32 the integrated run gives 876.198 mm² of wall–CPA interface and 388.453 mm² of CPA–myocardium interface. The strict coarse voxel minimum occurs 0.301 s into the run on a myocardial cell near the short-axis pole, at −149.741 °C. The strict 1.1 mm q=32 cut-cell minimum occurs at 0.008 s on a different myocardial state, at −149.993 °C. Halving the step from 2 to 1 ms changes the sampled minimum by at most 0.00417 °C/min over all listed method and schedule pairs.

With a common 0.5 ms trapezoidal diffusion update, the strict minima were 51.222295 °C/min on coarse voxels, 52.605542 for the centerline gap, 52.605168 for the q=24 integrated cut cells and 52.599698 for the 0.75 mm fixed-shape voxels. The classification therefore survives using the same update across all four representations. [Matched-update results](matched-results.json).

## Interpretation and earlier work

The initial temperature is uniform, so the initial physical diffusion flux vanishes. The myocardial source divided by volumetric heat capacity gives 52.606 °C/min. The cut-cell minimum is within 0.001 °C/min of that analytic baseline. Its agreement with the gap-state model shows that the coarse startup penalty disappears under geometry-integrated transport and storage. The startup result supports a benchmark interpretation of the outlying coarse voxel rate. It adds little to the estimate of the rate itself beyond the initial source/capacity ratio. The cut-cell flux uses a two-point centroid approximation on irregular regions, so the observed stability across these meshes applies to the discrete specimen states.

Joshi et al. already used mesh-converged tetrahedra for thermal analysis of rat-heart containers. Vispute and Rabin also describe mesh sensitivity at a vessel/CPA/container interface and modify that geometry to avoid high-aspect-ratio elements. Crockett et al. provide a conservative embedded-boundary precedent for discontinuous heat conduction. Their reported second-order accuracy uses a different flux reconstruction. CryoNet's magnetic-fluid discussion anticipated region-specific warming and Biostasis discussed uniform warming through glass transition. The new comparison concerns one Vitrify lattice defect on a common synthetic geometry. The bibliography and cached evidence document those scopes. [Joshi et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/) · [Vispute and Rabin](https://doi.org/10.1371/journal.pone.0350078) · [Crockett et al.](https://escholarship.org/content/qt00q9z6st/qt00q9z6st.pdf?t=lnrl2o) · [CryoNet](https://www.cryonet.org/cgi-bin/dsp.cgi?msg=18925) · [Biostasis](https://www.biostasis.com/new-warming-breakthrough-for-cryopreserved-organs/).

The 51.5 °C/min floor is a numerical startup diagnostic. The minimum occurs near −150 °C under extrapolated SAR. The primary rat-heart study evaluated local rates from the glass transition through warming termination and stopped once the coldest heart point reached −35 °C. Vitrify throttles initial heating to reduce stress. Any decision about the original schedule objective would require the thermomechanical stress history as well as the warmer-interval rates and temperature-based endpoint. This audit evaluates only first-second heat transport on the synthetic phantom. [Primary thermal study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/) · [Primary stress study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0290063).

## Reproduction

Run `sh verify.sh` from this directory. The script uses an available NumPy/SciPy environment or installs the pinned versions into a job-local virtual environment. The verifier reconstructs the graphs and recomputes the transients with independent sparse incidence fluxes and independently written material polynomials. It checks analytic geometric integrals, heat balance, comparable step sizes and the sign of the threshold crossing. A separate quadrature order and a smaller step with a trapezoidal diffusion update challenge the decisive cut-cell result. This yielded 52.605168 °C/min on the 1.1 mm q=24 graph and 52.598777 °C/min on the 0.75 mm graph. Run `python3 -B analysis.py` in the same dependency environment to regenerate the original rate table. The verifier also recomputes the strict case using one common 0.5 ms trapezoidal diffusion update for the coarse voxel, gap state, off-sequence q=24 cut cells and 0.75 mm voxels. Source and material coefficients are evaluated at the beginning of each step. The check compares these runs with the archived matched results and tests the source-only benchmark.
