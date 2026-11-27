# Does phase-overlap temperature transfer resolve the runner-up stress order?

## Benchmark

The previous contribution compared cut cells with thermal voxels at fixed mechanical spacing and observed a reversal of the graded versus uniform-strict smooth stress order. It also found a stable lead for the uniform release-policy schedule in those configurations. That comparison used nearest and eight-neighbor inverse-distance transfers. This audit tests an overlap average with an affine reconstruction, and measures the transfer sensitivity of the **paired** contrast while keeping the thermal histories fixed. The three schedule vectors are reused, including the tuned uniform release-policy comparator. [Predecessor report](sources/predecessor-report.md) · [Pinned Vitrify](https://github.com/Marc-Dvci/Vitrify/tree/1d498ee5c583e0112a762e60892d576ec2e27f57).

The donor is the predecessor's fixed analytic geometry sampled with 0.75 mm cut cells, three-axis order-24 ray integration, extrapolated SAR and 0.05 s steps. Each schedule is integrated from −150 °C to the interpolated coldest-heart −35 °C event. The thermal screen is 50 °C/min after each specimen state first crosses −123 °C and a maximum endpoint temperature below 0 °C. All three pass. The event times are 65.859, 65.834 and 65.822 s for graded, uniform release policy and uniform strict, respectively. This event and post-transition screen match the predecessor's benchmark. The release's optimizer scores an early 12 s horizon. The primary mechanical study includes stress inherited from cooling, which is absent from this stress-free initial state. [Pinned coupler](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/coupler/workflow.py) · [Primary thermal study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/) · [Primary mechanical study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0290063).

## Crossed comparison

The native 1.4 mm mesh checks continuity with the accepted comparator. For refinement, a single continuous shape determines element-center labels at spacings 1.4, 1.1 and 0.9 mm, with a common origin. The enclosing box is rounded up at each spacing. All transfers at a given spacing use the same Q1 hexahedra, label masks, rigid-mode constraints, Maxwell update and 0.30 s sample cadence, plus the exact thermal event. The smooth score is the release's unnormalized log-sum-exp over all three principal stresses with 20 kPa scale. We compare paired gaps within each grid because the smooth score also depends on the number of spatial samples. The maximum tensile stress is a separate diagnostic.

The IDW8 baseline interpolates the donor centroid states at element centers. The constant overlap method instead averages donors weighted by the volume of the intersection of their thermal cell, the target hex and the target label's analytic phase. The affine method first fits slopes from 16 same-phase donor centroids using distance-weighted linear least squares. Each donor's reconstruction is anchored to its thermal state at its own phase centroid. Integration uses analytic phase-boundary intersections along the ray and midpoint quadrature on each transverse interval of the common donor and target grid. A final run changes the ray direction from x to z and the transverse quadrature order from 6 to 9 at 0.9 mm. This checks quadrature sensitivity without changing the mechanical discretization.

Phase-constant and independent x, y and z affine fields with phase jumps serve as patch tests. The affine target is the field evaluated at the target phase centroid, which equals its phase-volume average. The tests check algebraic reproduction for the computed overlap moments. The change of ray direction and order probes the separate effect of approximate phase geometry integration. Slivers absent from the thermal graph use the nearest donor of the same phase with its reconstruction evaluated on that sliver. The diagnostics record the missing weight. The mechanical model still treats each hex as a whole center-labelled element. The phase-overlap rule estimates an element temperature. Conservation of enthalpy and cut-cell mechanical equilibrium would require a different discretization.

## Prior work and significance

Farrell and Maddison describe integration on donor–target intersections for conservative Galerkin projection. This audit uses that established idea in a narrower setting with same-phase volume averages and the release-style cellwise Maxwell constitutive update. Vitrify's existing spacing experiment rebuilds the native phantom and evaluates the associated models together. Joshi and Rabin model cooling and rewarming stresses; their calculation and the later artery study establish the relevance of mechanics and geometry during rewarming. Specialist discussions in Biostasis and CryoNet identify glass-transition stress and uneven warming as longstanding questions. The checked repository issue and pull lists were empty at the predecessor's retrieval. The new result concerns the numerical reliability of one benchmark schedule distinction. [Farrell and Maddison](https://doi.org/10.1016/j.cma.2010.07.015) · [Pinned grid study](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/experiments/grid_convergence.py) · [Artery study](https://doi.org/10.1371/journal.pone.0350078) · [Biostasis](https://www.biostasis.com/new-warming-breakthrough-for-cryopreserved-organs/) · [CryoNet](https://www.cryonet.org/cgi-bin/dsp.cgi?msg=18925).

## Results

**The runner-up order reverses under the tested refinement.**
With affine overlap, the smooth graded-minus-strict contrast is -0.238, +1.070, -0.483
kPa at 1.4, 1.1 and 0.9 mm, respectively. A positive difference favors
uniform strict. The three-grid span is 1.553 kPa, compared with a
0.483 kPa absolute contrast at 0.9 mm. The 1.1-to-0.9 mm change
alone is 1.553 kPa. At 0.9 mm the three transfers span 1.196
kPa on that same contrast. Switching the affine integrator to the alternate
ray direction and denser quadrature shifts the contrast by 0.003 kPa,
to -0.480 kPa. At that spacing, the hard peak favors uniform strict
under affine overlap while the smooth score favors graded.

All entries below are signed kPa. The two rightmost columns give each
runner-up’s smooth or hard score above the uniform release-policy reference.
Positive values in those columns favor the reference.

| Mechanics and transfer | Graded − strict, smooth | Graded − strict, hard | Graded − reference, smooth / hard | Strict − reference, smooth / hard |
|---|---:|---:|---:|---:|
| 1.4 mm, IDW8 | -0.930 | -0.840 | +0.907 / +0.745 | +1.837 / +1.584 |
| 1.4 mm, constant overlap | +0.351 | -0.215 | +2.126 / +1.993 | +1.776 / +2.208 |
| 1.4 mm, affine overlap | -0.238 | -0.407 | +1.769 / +1.510 | +2.008 / +1.917 |
| 1.1 mm, IDW8 | -0.202 | -0.144 | +1.368 / +1.396 | +1.570 / +1.539 |
| 1.1 mm, constant overlap | -0.248 | -1.390 | +1.500 / +0.586 | +1.747 / +1.976 |
| 1.1 mm, affine overlap | +1.070 | +5.791 | +2.596 / +7.165 | +1.526 / +1.373 |
| 0.9 mm, IDW8 | -1.187 | -1.266 | +0.614 / +0.745 | +1.801 / +2.011 |
| 0.9 mm, constant overlap | +0.009 | +1.075 | +1.635 / +2.600 | +1.626 / +1.526 |
| 0.9 mm, affine overlap | -0.483 | +0.722 | +1.269 / +2.469 | +1.752 / +1.747 |
| 0.9 mm alternate quadrature, affine overlap | -0.480 | +0.722 | +1.272 / +2.473 | +1.752 / +1.751 |

| Transfer | Max phase-constant error (°C) | Max phase-affine error (°C) |
|---|---:|---:|
| IDW8 | 2.84e-14 | 0.0913 |
| constant overlap | 2.44e-12 | 0.0382 |
| affine overlap | 2.4e-12 | 2.49e-12 |

The total fraction covered by the missing-donor fallback is at most
8.46e-10 across the target phase volumes, and at most 1.09e-06
in any individual element. The native IDW cross-check reproduces the
predecessor's cut-cell stress scores within 0.0000023 MPa. The uniform
release-policy schedule has the lowest smooth and hard scores in every tested
configuration. Its smallest lead across the runs and two scores is
0.586 kPa.

## Interpretation and reproduction

Retire the secondary rank for this benchmark and retain the uniform
release-policy schedule as the reference. The signed contrast changes on
mechanical refinement even when each element receives an affine-exact
phase average. The fixed-grid comparison measures transfer sensitivity at
that grid. The fixed-map refinement comparison measures sensitivity to
the mechanical discretization, center-labelled geometry and reaggregation
onto smaller elements together. The observed spans are finite numerical
sensitivity measures. The continuum limit remains unresolved. One thermal mesh,
an extrapolated SAR law, a stress-free start and an unlimited local affine
reconstruction restrict the result to this computational experiment.

With the NumPy and SciPy versions in `requirements.txt`, run
`sh verify.sh` from this directory. The command can install those
dependencies into a local environment when they are absent. It rebuilds
the analytic thermal grid, regenerates all three temperature histories,
checks the event and rate screen directly from those histories, and
solves all of the mechanical paths again. It checks the independent
rectangular overlap integral, phase reproduction, predecessor baseline,
archived stresses, signed gaps and observed ranges. It exits with an
assertion failure on disagreement. Use `sh generate.sh` to regenerate
`results.json` in an existing compatible environment.
