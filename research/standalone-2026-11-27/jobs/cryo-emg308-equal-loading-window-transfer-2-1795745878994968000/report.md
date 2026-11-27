# Smaller-bag temperature profiles improve window timing on the released 2 L trace

At nominal 4.6 mg Fe/mL, the released S22A 1 L trace predicts the 2 L bag’s seven 20 °C transit times with a mean absolute error (MAE) of 0.87 s. The strongest tested smaller-bag constant gives 1.97 s, and the fixed warm-vial reference gives 12.74 s. The 1 L trace was selected as a comparator by proximity in volume. It requires no fit to the 2 L data. Across the paper’s −100 to 0 °C reporting span, a straight line fitted to the 1 L window durations by least absolute deviations gives 0.70 s, compared with 0.81 s for direct 1 L transfer. The results support using temperature information on this release and leave the choice of temperature-dependent predictor for a new run unresolved.

## Sources and estimand

Gangwar et al. report similar aggregate nanowarming at equal loading. S22A releases one temperature curve per volume for 0.5, 1 and 2 L M22 bags with nominal 4.6 mg Fe/mL EMG308 and common nominal RF conditions. This comparison concerns intermediate transit times within those recorded curves. The preprint reported liter-scale warming and the published supplement documents the equal-loading conditions. The coordinator’s measurement bridge already derived the fixed warm benchmark and measured its Figure 7 late-window residual. The linked COMSOL v1 inventory lists cooling models for 0.5, 1 and 3 L. The model binaries were not inspected. [Version of record](https://pmc.ncbi.nlm.nih.gov/articles/PMC12474906/) · [Preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC11601234/) · [Cooling-model inventory](https://zenodo.org/records/16716098)

Biostasis discusses warming across the glass transition, and Fight Aging! discusses volumetric warming in earlier kidney and heart work. Later porcine work concerns organ recovery under a distinct protocol. These sources provide context for temperature dependence and scale, while the window comparison here uses the released M22 bag data. [Biostasis](https://www.biostasis.com/new-warming-breakthrough-for-cryopreserved-organs/) · [Fight Aging!](https://www.fightaging.org/archives/2022/01/nanowarming-of-vitrified-kidneys-and-hearts/) · [Porcine preprint](https://www.biorxiv.org/content/10.64898/2026.10.07.756990v1)

## Window comparison

Each transit time is the difference between first upward crossings, linearly interpolated between one-second samples. Predictions use smaller-bag durations throughout. Means and medians fit a single duration per 20 °C window. Least absolute deviations (LAD) and ordinary least squares (OLS) fit a line in window midpoint. All empirical fits use the scoring span’s training windows. The fixed warm reference uses mean S19E SAR at 35 kA/m for EMG308 in M22 without LM5 (663.43 W/g Fe), nominal 4.6 mg Fe/mL and warm M22 density and specific heat (1100 kg/m³ and 3500 J/kg/K). This yields 47.56 °C/min, or about 25.23 s per window.

| Window (°C) | 0.5 L (s) | 1 L (s) | Equal blend (s) | 2 L observed (s) |
|---|---:|---:|---:|---:|
| −140 to −120 | 6.77 | 10.23 | 8.50 | 8.37 |
| −120 to −100 | 13.90 | 11.76 | 12.83 | 11.58 |
| −100 to −80 | 12.95 | 12.64 | 12.79 | 11.99 |
| −80 to −60 | 12.39 | 11.84 | 12.11 | 11.20 |
| −60 to −40 | 15.47 | 12.86 | 14.16 | 13.13 |
| −40 to −20 | 13.51 | 14.83 | 14.17 | 14.33 |
| −20 to 0 | 19.02 | 18.82 | 18.92 | 16.84 |

| Predictor | −140 to 0 °C MAE (s) | −100 to 0 °C MAE (s) |
|---|---:|---:|
| 0.5 L windows | 1.630 | 1.497 |
| 1 L windows | 0.867 | 0.805 |
| Equal blend | 0.911 | 1.000 |
| 1 L linear LAD | 0.934 | 0.701 |
| 1 L linear OLS | 0.995 | 1.044 |
| Pooled linear LAD | 1.391 | 1.115 |
| Pooled linear OLS | 1.078 | 1.102 |
| 1 L median constant | 1.973 | 1.652 |
| 1 L mean constant | 2.110 | 1.813 |
| Pooled median constant | 2.011 | 1.620 |
| Pooled mean constant | 2.141 | 1.899 |
| Fixed warm-vial reference | 12.740 | 11.734 |

The 20 °C grid and trend comparisons are exploratory. For a common convex blend weight across windows, MAE is convex in that weight and attains its maximum at an endpoint. The worst blend MAE is 1.630 s over the full span and 1.497 s over the shorter span. Those bounds are 0.343 and 0.122 s below the respective best tested constant, so the advantage over constants on these curves needs no weight selection. A blend and its own mean duration give identical total transit predictions across their calibration span. The tested improvement concerns allocation of time between windows.

## Figure lineage and scope

S22B’s three values at nominal 4.6 mg Fe/mL (89.60000, 87.62792, 86.64333 °C/min) match the three Figure 7c 2 L replicate means across left, center and right to five decimal places. Treat these as reused aggregate evidence. The 95-sample S22A 2 L trace shares at most three consecutive identical temperatures with any of the three 97-sample Figure 7b probes across integer shifts. This rules out a direct copy of a long unchanged probe segment. Run identifiers are absent, so same-run or processed-trace overlap remains unresolved. On the Figure 7b probes, the best full-span predictor is the equal blend for left and center and direct 1 L transfer for right. The probes give a within-release check only.

Each S22A volume has one released curve, and windows and probes are dependent. The fitted predictor ranking is descriptive. The small ranking margins, one-second sampling and unresolved run lineage limit claims of generalization. Temperature rise combines heat capacity, sources and heat exchange. The calculation estimates net transit timing in these bags. Tissue warming, particle SAR and independent repeatability remain untested.

## Reproduction and attribution

Run `python3 -B analysis.py` and `python3 -B verify.py` in this directory. Both use the Python standard library. The verifier shares the workbook parser but uses rational interpolation, a different LAD reconstruction, exact enumeration of blend residual breakpoints and an independent trace-matching algorithm. It checks source hashes, figure labels, rankings, the saved results and aggregate reuse.

Source snapshots were copied from the prior attempt’s archived inputs. The unmodified supplement and workbook are attributed to Lakshya Gangwar et al., DOI 10.1038/s41467-025-63483-2. The version of record reports a change to CC BY 4.0 on 10 November 2025. Other primary and community snapshots retain their own attribution and terms. [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
