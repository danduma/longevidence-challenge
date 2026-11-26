# A conditional concentration-only fidelity limit for SceThg1

Status: independently reproduced and critically reviewed conditional mathematical result. The reviewer supports the proof and finds its biological predictive validity unestablished and its practical value modest. This is not a new enzyme, measured fidelity, or a solution to the 5′ polymerase challenge. [Full skeptical review](../agents/review-thg1-sol/report.md).

## Question

Can one fixed nucleotide mixture make all four Watson–Crick additions highly accurate for the SceThg1 single-addition substrates characterized by [Patel et al., RNA (2021), Table 1](https://doi.org/10.1261/rna.078686.121)? The study's raw table image was separately downloaded and visually checked against the scout's sixteen transcribed efficiencies. Its substrate is preactivated tRNA, not a generic processive polymerase template.

## Model and derivation

Assume mutually exclusive, dilute nucleotide capture, with competing hazards proportional to each published `kpol/KD` value times that nucleotide's concentration. This is an explicit approximation, not an experimentally established competition law for these substrates. One positive concentration vector must serve all four substrate variants.

Index rows by the nucleotide that would make the correct Watson–Crick pair. Let `a[i,j]` be the measured efficiency for incoming nucleotide `j` on the substrate requiring nucleotide `i`, and `x[j]` its concentration. The predicted correct-addition probability is

`p[i] = a[i,i] x[i] / sum_j a[i,j] x[j]`.

Define `B[i,j] = a[i,j]/a[i,i]` off diagonal and zero on diagonal. Then

`p[i] = 1 / (1 + (B x)[i]/x[i])`.

The Collatz–Wielandt characterization of the Perron eigenvalue gives the optimal worst-case probability over all positive mixtures:

`sup_x min_i p[i] = 1 / (1 + rho(B))`.

This uses a standard theorem. Any potential novelty would be the dataset-specific application and checked consequence, not the mathematics itself.

## Result

One correct-incorporation efficiency is lower-censored in the paper. Substituting its published lower bound gives an optimum worst-substrate probability of **92.8349261%**. Increasing that censored entry can improve the result, so treating this number as an unconditional upper bound would be wrong.

A stronger bound removes the row and column corresponding to that censored entry. The spectral radius of this nonnegative principal submatrix is a lower bound on the full matrix's spectral radius. The resulting **94.2961809% upper bound** therefore holds for any value of the censored correct-incorporation efficiency, while retaining the other published rounded point estimates and the model assumptions.

There is also an elementary two-substrate certificate: achieving 99% correctness on both the U73 and G73 substrates requires reciprocal nucleotide ratios whose product is at least **5.46867**. Their product must actually equal one. Thus 99% uniform correctness is infeasible even without using the censored entry or the other two substrate variants.

The scout used a power iteration with Collatz–Wielandt bounds. An independent implementation used an eigensystem and linear-program feasibility bisection; the numerical optima agree within `2.3e-9` in three censoring scenarios. The LP independently rejects the 99% target.

As a deterministic robustness check, allow every uncensored `kpol` and `KD` to vary independently by up to twice its reported ± magnitude, choosing errors in whichever direction most helps fidelity. Lower-bound every mismatch/correct efficiency ratio using the least mismatching and most favorable correct rates in that rectangle. The principal-submatrix argument still yields an upper bound of **98.8857727%**, below 99%. This is an adversarial sensitivity envelope, not a confidence interval; the statistical meaning and correlations of the reported ± quantities have not been established here. Recomputing the unperturbed efficiencies from `kpol/KD` rather than the paper's rounded ratio column gives a slightly different 94.3030154% ceiling.

The independent reviewer obtained exact rational certificates, including a cubic-root bracket, without importing the coordinator's numerical code. The ±2 exclusion has a narrow margin: expanding that particular favorable parameter box to approximately 2.0897 times the reported errors permits a 99% ceiling. The approximate phrase "robust to uncertainty" would conceal this limitation.

## Interpretation and limits

Within this model, adjusting nucleotide abundance alone cannot turn the measured SceThg1 substrate behavior into uniformly high-fidelity copying. That supports investigating different scaffolds or changes to recognition rather than treating mixture tuning as sufficient. It does not rule out changes to the enzyme, substrate, conditions, kinetic regime, proofreading, staged nucleotide delivery, or Thg1-like proteins.

The inputs are rounded point estimates with experimental uncertainty. The upper bound is mathematical conditional on those inputs, not a statistical confidence limit. Processivity, sequence-context dependence, actual multi-substrate competition, yield and extension fidelity remain unmeasured here. The challenge also requests DNA synthesis and other polymerase activities that this RNA/tRNA analysis does not address.

An external diagnostic reinforces the scope restriction. [Rao et al. 2011](https://doi.org/10.1093/nar/gkq976) measured a roughly 5:1 G/U product ratio on one yeast substrate, whereas inserting the 2021 coefficients into this model gives about 21.8:1. Activation state and concentration regime differ, so this is not a matched-condition falsification; it shows why the coefficients cannot be transferred without testing. [2024 extension work](https://doi.org/10.1261/rna.080156.124) and [2026 specificity studies](https://doi.org/10.1261/rna.080815.125) already motivate alternative TLP scaffolds. Poor native SceThg1 Watson–Crick discrimination and the preference for exploring TLPs are established biology.

The bounded potential contribution is the particular common-mixture optimization, censoring-independent certificate and sensitivity analysis. The skeptical reviewer found no matching bound in its inspected paper/preprint/code/patent searches, but this does not establish universal priority. The proof is not being counted as substantial challenge-level progress or as satisfaction of the goal's novelty gate.

## Reproduce

From the workspace root:

```sh
research/millennium/.venv/bin/python research/millennium/polymerase/verify_thg1.py
```

Independent outputs and all caveats are recorded in `results/thg1_independent_verification.json`. The coordinator's [verified table](thg1_table1_verified.csv) and [metadata](table1_metadata.json) label the published ± quantities as reported error magnitudes. The raw scout CSV remains preserved under its original path; its `_sd` headers are not verified statistical metadata. The scout implementation is `../agents/scout-mechanisms-sol/analysis/thg1_competition.py`; the exact-arithmetic independent checker is `../agents/review-thg1-sol/check_bound.py`.
