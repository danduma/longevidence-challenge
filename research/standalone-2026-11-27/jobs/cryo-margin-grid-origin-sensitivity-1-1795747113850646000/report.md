# Grid-origin sensitivity of the Vitrify margin designs

Vitrify’s released mesh study changes the spacing while retaining one sampling origin at each spacing. The coordinator’s previous comparison added two jointly tuned uniform surrounds and found a 0.043 °C/min strict-screen margin at 1.1 mm. This audit fixes those three protocols and changes the location of the sample points within each grid cell. It asks whether the same rate decision and mechanical ordering survive a discretization choice that leaves the underlying analytic geometry unchanged. The two optimized uniform witnesses are the relevant comparators because the release’s constant-duty 7.5 mg Fe/mL baseline misses the temperature screen. [Released mesh record](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/results/grid_convergence_margin.json).

## Comparison

The recorded phases are the native origin, a half-cell shift along each of x, y and z, and the joint half-cell shift. Test spacings are 1.1 mm, which set the lowest rate for the previous strict witness, and 0.75 mm, the finest released spacing. At each spacing the native pillow exponent is calibrated once and reused. Chamber shapes, heart dimensions, nominal bag dimensions and wall thickness remain fixed. The graded field is evaluated in physical coordinates around the same native loading centre for each phase. The source uses the native array centre, which differs slightly from the analytic phantom centre. The two spacings use their respective native exponents and loading centres, so the within-spacing contrasts isolate phase. Voxel counts change with phase as expected for centre-sampled labels. [Source geometry](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/src/vitrify/geometry.py), [source loading map](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/coupler/workflow.py), [prespecified phases](protocol.json).

The sweep reuses the coordinator's finite-volume and Q1 Maxwell solver with extrapolated SAR, the released material parameters, the original loading and duty bounds and the prior comparison’s exact vectors. The 70 s whole-specimen screen applies both rate floors and the −35 to 0 °C terminal temperature window. Paired mechanics terminate when the coldest heart cell reaches −35 °C, with interpolation at the crossing. Thermal steps are set by the source stability criterion. Mechanics samples use a nominal 0.24 s cadence through 16 s, nominal 1.8 s thereafter and an exact interpolated endpoint sample. The two stress metrics are peak tension and the original unnormalized smooth maximum of principal stresses with 0.02 MPa smoothing. Compare scores within a phase because voxel count and event time can change across phases. [Solver snapshot and attribution](attribution.md), [rat-heart endpoint](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/).

## Rate decision and limiting location

| Grid | Phase | Graded | Release-policy uniform | Strict uniform |
| --- | --- | ---: | ---: | ---: |
| 1.1 mm | Native | 51.2733 | 51.3103 | 51.5429 |
| 1.1 mm | x | 51.1344 | 51.1702 | 51.4010 |
| 1.1 mm | y | 50.7099 | 50.7463 | 50.9745 |
| 1.1 mm | z | 52.3155 | 52.3545 | 52.5941 |
| 1.1 mm | xyz | 52.3155 | 52.3545 | 52.5941 |
| 0.75 mm | Range over five phases | 52.3134–52.3157 | 52.3527–52.3550 | 52.5921–52.5944 |

All rates are in °C/min. At both spacings, all candidates satisfy the 50 °C/min validation screen, 70 s temperature window and event maximum-temperature check in every tested phase. The 51.5 °C/min screen changes its verdict for the strict uniform witness at 1.1 mm. Its worst sampled margin is −0.5255 °C/min. At 0.75 mm each design’s spread is below 0.003 °C/min, and all clear 51.5 °C/min.

Every selected limiting cell is myocardium, below −149.6 °C and within 0.31 s of startup. On the native, x and y 1.1 mm grids the strict witness’s limiting cell has one face against container wall. On the other two 1.1 mm grids and all five 0.75 mm grids the selected limiter has myocardium on all six faces. At the three wall-adjacent sample lines, direct intersection with the analytic surfaces gives a 0.810–0.833 mm CPA gap in the +z direction. Both surfaces lie between the two sample centres. The coarse stencil therefore skips this gap at those locations. The selected locations vary, and symmetries can produce ties or near-ties. The adjacency result describes the selected cells and cannot isolate the wall contact as the sole cause. This is the release's all-time screen. The early minima precede the glass-transition interval used in the 2022 rat-heart thermal interpretation, so the changed numerical verdict gives no direct evidence of devitrification. [Thermal outputs](thermal.json), [2022 rat-heart analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/).

## Paired mechanical ordering

Entries are uniform minus released graded in MPa. Negative values favour the uniform witness for the named metric.

| Grid | Phase | Policy uniform objective | Policy uniform peak | Strict uniform objective | Strict uniform peak |
| --- | --- | ---: | ---: | ---: | ---: |
| 1.1 mm | Native | -0.001419 | -0.001558 | +0.000315 | +0.000234 |
| 1.1 mm | x | -0.001358 | -0.001221 | +0.000300 | +0.000399 |
| 1.1 mm | y | -0.001057 | -0.001403 | +0.000716 | +0.000518 |
| 1.1 mm | z | -0.000762 | +0.000680 | +0.001125 | +0.002837 |
| 1.1 mm | xyz | -0.000836 | -0.000458 | +0.001059 | +0.001500 |
| 0.75 mm | Native | -0.001516 | -0.002418 | +0.000230 | -0.000659 |
| 0.75 mm | x | -0.001540 | -0.001153 | +0.000170 | +0.000676 |
| 0.75 mm | y | -0.000935 | -0.000688 | +0.000873 | +0.001220 |
| 0.75 mm | z | -0.001590 | -0.001041 | +0.000030 | +0.000674 |
| 0.75 mm | xyz | -0.001004 | -0.000822 | +0.000742 | +0.001191 |

The release-policy uniform design retains its sampled smooth-objective advantage in all ten pairs. Its peak-tension advantage survives all five 0.75 mm phases and four of five 1.1 mm phases. A 1.1 mm z shift reverses the native peak ordering even though both schedules remain feasible under the release validation policy.

At 0.75 mm the strict uniform schedule has lower peak tension at the native origin and higher peak tension in all four shifted phases. Its smooth score stays above the released graded score, by as little as 0.000030 MPa. Both designs pass the strict thermal screen at that spacing. Some of the corresponding 1.1 mm scores compare schedules that fail the stricter screen. The sizes are small and conditional on the sampled solver and cadence. The policy uniform schedule remains the model benchmark for the sampled smooth objective, and its advantage in peak tension depends on phase and spacing. [Paired outputs](results.json).

## Closest precedent and scope

The pinned Vitrify geometry fits the pillow exponent to discretized volume, so its spacing sweep also changes the fitted continuous shape. The released experiment samples a single origin and compares to the fixed Case II schedule. The previous coordinator audit retuned the uniform baseline but also used a single phase at each spacing. Joshi and Rabin used mesh convergence to choose finite-element counts and observed case-dependent stress locations. The later popliteal artery model documented mesh sensitivity at material interfaces and used fitted tetrahedral meshes. CryoNet discussed spatially variable heating and Biostasis discussed uniform warming near glass transition. Those discussions address other aspects of the physical problem. The targeted search and cached Vitrify issue, PR and margin-commit records did not establish an earlier test of these fixed candidates under phase shifts. [Vitrify grid experiment](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/experiments/grid_convergence.py), [2023 mechanics study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0290063), [2026 artery model](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0350078), [CryoNet](https://www.cryonet.org/cgi-bin/dsp.cgi?msg=18925), [Biostasis](https://www.biostasis.com/new-warming-breakthrough-for-cryopreserved-organs/).

The conclusions apply to five deterministic phases per spacing for the released synthetic phantom, extrapolated SAR and fixed schedule vectors. They do not establish the worst phase, continuum convergence, an equal-budget comparison of optimal design classes, or organ-level biological evidence. Weak material is used for nonspecimen cells in the mechanical solve, as in the source, and padding is expanded if an active sample would otherwise fall outside the native array. The verifier reruns voxelization, finite-volume integration, the event interpolation and Maxwell equilibrium solves on all ten lattices from the cached inputs. It checks native thermal anchors against the upstream record and the native finest-grid mechanics against the prior coordinator result. [Verification instructions](attribution.md).
