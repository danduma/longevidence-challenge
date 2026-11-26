# C10 | How accurately does reporter AUC predict a relative initial rate?

A fixed linear calibration reduced average held-out enzyme-panel RMSE from **0.4302 to 0.3081**. It beat the training mean on each of the six selected TEVp panels. Holding out every substitution at one substrate position across all enzymes was harder. Linear regression beat the training mean in five of seven positions, reducing average fold RMSE from **0.3596 to 0.3144**, but systematically underestimated rates at P5. These are errors in dimensionless initial-rate ratios, where 1 denotes the corresponding canonical control.

## Data and comparison

The matched table contains 786 unique enzyme–motif pairs. Removing the six canonical keys left 780 measured pairs spanning 133 single-substitution motifs. All 170 negative rate estimates remained in the analysis. The assay measured recorder activity as a dimensionless time-averaged area under the flipped-fraction curve. The analysis used its published correction and the prior reconstruction's DMS-parent substitution for panel 0. Both raw and corrected reporter depth had to reach 100. No missing biochemical pair was filled in.

Each fit saw corrected AUC as its sole predictor. In the enzyme comparison it saw zero biochemical labels for the held-out enzyme, though many of the same substrates appeared for other enzymes. In the separate position comparison it saw labels for the same enzymes at the six other positions, with zero instances of a held-out motif. Each pair received one out-of-fold prediction per design and model. The primary comparison was linear regression against the pooled training mean. An increasing isotonic fit and, where available, a within-enzyme training mean were fixed in advance. [Frozen protocol](protocol.md)

## Out-of-sample performance

Mean fold RMSE weights the held groups equally. Pooled errors weight observed enzyme–motif pairs equally. Bias is prediction minus observation. The final two columns describe a regression of observed rates on out-of-fold predictions, with ideal intercept 0 and slope 1. These diagnostics did not alter the fits.

| Held group | Model | Mean fold RMSE | Pooled RMSE | Pooled MAE | Bias | Intercept | Slope |
|---|---|---:|---:|---:|---:|---:|---:|
| Enzyme | Training mean | .4302 | .4302 | .2955 | .0000 | 1.3225 | −5.0363 |
| Enzyme | Linear | .3081 | .3089 | .2024 | +.0023 | .0082 | .9527 |
| Enzyme | Isotonic | .3008 | .3018 | .1585 | +.0045 | .0146 | .9146 |
| Position | Training mean | .3596 | .4849 | .3362 | −.0004 | 1.5131 | −5.9189 |
| Position | Within-enzyme training mean | .3597 | .4849 | .3357 | −.0004 | 1.3281 | −5.0728 |
| Position | Linear | .3144 | .4123 | .2557 | −.0517 | .1102 | .6503 |
| Position | Isotonic | .2934 | .3957 | .2052 | −.0520 | .0914 | .7639 |

The negative pooled calibration slopes for the mean controls reflect complementary training folds. Holding out a high-mean group lowers the mean of the remaining training labels. These pooled slopes have no within-fold predictive interpretation for the constant baseline. The small pooled linear bias in the enzyme design also hides offsetting panel biases, from −0.130 on NX1 to +0.119 on VI.

| Held-out enzyme or position | Pairs | Mean RMSE | Linear RMSE | Isotonic RMSE |
|---|---:|---:|---:|---:|
| 0 | 133 | .413 | .262 | .249 |
| NX1 | 125 | .439 | .357 | .344 |
| NX2 | 130 | .439 | .310 | .315 |
| VI | 131 | .419 | .319 | .308 |
| VI* | 132 | .430 | .295 | .277 |
| VII | 129 | .441 | .306 | .313 |
| P6 | 108 | .207 | .266 | .239 |
| P5 | 113 | 1.145 | .950 | .915 |
| P4 | 113 | .259 | .186 | .174 |
| P3 | 108 | .208 | .275 | .265 |
| P2 | 110 | .228 | .196 | .174 |
| P1 | 114 | .250 | .116 | .043 |
| P1′ | 114 | .222 | .211 | .244 |

The prespecified average-improvement criterion passed in both designs. Linear regression lost to the mean at P3 and P6. At P5 it improved on the baseline yet underpredicted by 0.836 rate units on average. That fold accounts for much of the gap between mean fold RMSE and pooled RMSE, which squares errors across groups before taking the root. Isotonic regression also lost to the mean at P1′. No universally reliable quantitative conversion follows. [Fold errors](results/fold_errors.csv) · [All out-of-fold predictions](results/predictions.csv)

## What normalization allows us to test

The deposited initial estimates identify the reference more precisely than “one canonical rate per enzyme.” For each enzyme and position, dividing the 20 initial estimates by that series' own canonical entry reproduces all 840 deposited relative values to within 6 × 10⁻¹⁵. Each enzyme has seven distinct canonical initial estimates. The largest/smallest within-enzyme reference ratio ranges from 1.19 to 3.16. The resulting 42 canonical relative entries are exactly 1 by construction and collapse to six keys in the matched table. They cannot quantify replicate uncertainty. [Control audit](results/control_audit.csv)

Including those six keys in enzyme training and scoring moves linear mean fold RMSE from 0.3081 to 0.3085 (mean baseline 0.4302 to 0.4339). Including canonical keys only in position training moves it from 0.3144 to 0.3110 on identical tests (baseline 0.3596 to 0.3630). Direct inclusion of unit labels therefore does little to explain the primary advantage. Removing them cannot undo the division already present in every response. Substitutions in a given enzyme–position series share a denominator and any error in it. Both primary designs keep each entire denominator group on one side of each fit. A random pair split could share that denominator across sides. Shared enzyme, substrate or assay context can still induce dependence across groups. The held-panel predictions require only an AUC as input, while evaluating those predictions against the published endpoint requires held-panel canonical measurements. [Sensitivity results](results/summary.csv)

## Scope and verification

Earlier studies reported within-panel rank agreement and tested a target-free time-course score. C10 measures held-out **rate error** for a calibration fitted on other groups. Six related, selected enzymes and one canonical substrate neighborhood support a finite comparison. There is no test here of simultaneous novelty in both axes, catalytic turnover in absolute units, a whole proteome, or cleavage at unseen folded sites. The paired summary estimates supply no replicate variances. No row-independence intervals or human expert review are claimed.

Run `./verify.py` from this directory. It checks frozen hashes, reconciles the source labels and denominators, verifies split coverage, regenerates the out-of-fold predictions and checks all four result tables. [Source versions and download instructions](sources.md) · [Machine-readable contribution](contribution.json)
