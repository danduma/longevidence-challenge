**Independent skeptical review — 2026-11-26**

**Verdict: SUPPORTED as conditional mathematics and numerical reproduction; UNDERQUALIFIED if presented as a biological ceiling or a substantial advance toward Millennium problem 11.** I found no transcription, indexing, optimization, censoring-direction, or parameter-rectangle error that invalidates the stated conditional results. An independent exact-arithmetic calculation reproduces both 94.2961809% and 98.8857727%. The major unresolved issue is whether the assumed competing-incorporation law predicts the relevant experiment. The broad biological observations are established prior art; the particular fixed-mixture optimization may be a useful small reanalysis, but this search does not establish priority.

| Claim being assessed | Assessment |
|---|---|
| Numerical reproduction | Supported: independently transcribed inputs; exact rational certificates agree with the supplied output. |
| Conditional global optimum | Supported for each finite positive correct-G efficiency and the specified homogeneous competing-hazard model. |
| Independence from the censored correct-G value | Supported for the 94.2961809% ceiling; the 92.8349261% result is conditional on substituting the lower bound. |
| ±2 rectangle bound | Supported as a deterministic sensitivity result over the chosen physical parameter box. It has no established coverage probability. |
| Prediction of measured mixed-NTP fidelity | Unvalidated. Existing competition data warn against extrapolation across activation state and kinetic regime. |
| General limit for Thg1/TLP reverse polymerases | Unsupported. Enzyme, substrate, context, activation, and conditions matter. |
| New mathematics / discovery of low Thg1 fidelity | Unsupported novelty claims. Both the mathematical machinery and biological background predate this analysis. |
| Progress toward the full challenge | Modest diagnostic value; no new catalyst, processivity demonstration, or challenge-level error measurement. |

I read the specified shared Markdown, Python, JSON, CSV, table image, and 2021 primary-source XML without running the shared writer script. I also read the supplied 2024 primary-source XML and retrieved an identical 2024 XML independently. All review artifacts are in this job directory. No subagents, other model sessions, contacts, posts, or experiments were used. `sources/shared_inputs_sha256.json` records the shared inputs; `refusals.md` and `access_failures.md` distinguish refusals from retrieval limitations.

**The input and indexing audit passed.** Visually checking the original Table 1 image against a fresh manual transcription reproduced all 16 efficiencies and all reported kinetic parameters/error magnitudes in the CSV. The correct incoming nucleotides for template A, C, G, U are U, G, C, A respectively. For incoming-nucleotide order A,C,G,U, the reordered template rows must therefore be U,G,C,A. The proposed verifier makes this permutation correctly; the censored entry is template C / incoming G, not template G / incoming C. [Patel et al., 2021, Table 1](https://doi.org/10.1261/rna.078686.121).

For clarity, the printed efficiency matrix, after this row permutation, is:

| Correct incoming / template | Incoming A | Incoming C | Incoming G | Incoming U |
|---|---:|---:|---:|---:|
| A / U73 | 0.023 | 0.00011 | 0.41 | 0.00069 |
| C / G73 | 0.0035 | 0.030 | 0.16 | 0.00024 |
| G / C73 | 0.00041 | 0.000033 | g ≥ 10 | 0.00027 |
| U / A73 | 0.00046 | 0.000045 | 0.12 | 0.0055 |

Efficiencies have units µM⁻¹ min⁻¹. Dividing each row by its correct diagonal entry and setting the diagonal to zero gives the proposed B. Deleting incoming G means deleting **both** its row and column. The retained matrix, in A,C,U order, is exactly

```text
C = [    0       11/2300      3/100 ]
    [  7/60         0        1/125 ]
    [ 23/275      9/1100        0  ].
```

Keeping a row indexed by the physical template instead of the correct incoming nucleotide would invalidate the eigenvalue argument. That error is absent here. Transposing B alone happens to preserve its spectral radius, so numerical agreement of a radius would not by itself verify indexing or a proposed optimal mixture; the direct probability expression also has to agree.

**The optimum is global and attained for finite g.** Let x be a strictly positive free-NTP concentration vector and let a be the reordered efficiency matrix. Under the assumed hazards,

```text
p_i(x) = a_ii x_i / Σ_j a_ij x_j
       = 1 / [1 + (Bx)_i/x_i].
```

Consequently, maximizing the smallest p_i is equivalent to minimizing `max_i (Bx)_i/x_i`. Here B is nonnegative and irreducible, with every off-diagonal entry positive for finite g. It has positive right and left Perron eigenvectors v,w. For any x>0, putting `r=max_i (Bx)_i/x_i` gives Bx≤rx. Multiplication by wᵀ yields `rho(B) wᵀx≤r wᵀx`, hence r≥rho(B). Taking x=v attains equality and makes all four p_i equal. This proves

```text
max_{x>0} min_i p_i(x) = 1/[1+rho(B)].
```

Scale invariance permits normalization of x; no absolute concentration, speed, yield, or reagent-use constraint is included. Irreducibility explains why “maximum” is justified for this dataset. For a general reducible matrix, only a supremum over strictly positive vectors may be attained in the closure. The formulation is a standard Perron/Collatz–Wielandt optimization, closely equivalent to classical interference-limited max-min SIR balancing after replacing fidelity odds by signal/interference ratios. [Stańczak et al., 2009](https://arxiv.org/abs/0901.0824).

**Censoring is handled in the correct direction.** The paper reports kpol≥10, apparent KD≤1, and efficiency≥10 for correct G on C73. It explains that GTP binding was too tight to reach subsaturating measurements. This is a lower limit, not a measured efficiency of exactly 10. Raising g decreases the off-diagonal entries of the incoming-G row of B, lowering its spectral radius and improving the optimum. Therefore 92.8349261%, calculated at g=10, is a pessimistic endpoint of the possible model optima, not a universal upper ceiling. [Patel et al., 2021](https://doi.org/10.1261/rna.078686.121).

The bound from C survives for any g. For every x>0, removing the G contribution only reduces interference in the retained rows. Thus `max_i (Bx)_i/x_i ≥ max_{i∈ACU}(Cx_ACU)_i/x_i ≥ rho(C)`. No estimate of g enters. For finite g and the other positive entries, the full optimum is strictly below the principal ceiling. As g tends to infinity, B tends to a block triangular matrix with blocks C and 0, so its optimum approaches that ceiling. It is a sharp supremum over the unbounded censored parameter, rather than the attainable optimum of a specified finite-g enzyme.

**Independent numerical checks use exact arithmetic, not the supplied eigensolver/LP route.** `check_bound.py` embeds the fresh image transcription; it does not import the shared verifier, scout calculation, NumPy, SciPy, or an eigenvalue routine. It compares its transcription with the shared CSV only after constructing its own data. For 4×4 matrices, rational iteration proposes a positive vector; exact Collatz–Wielandt residual inequalities certify a tight enclosure regardless of iteration convergence. For each 3×3 principal matrix it independently constructs the characteristic cubic and bisects its unique positive root using rational signs. The saved outputs include exact endpoints and positive-vector certificates.

For C, the cubic is

```text
f(r) = r³ - (118879/37950000) r - 1751/55000000.
```

It has exactly one positive root, by its coefficient signs. Exact rational bisection gives `rho(C)=0.06048833626998886…` and the probability ceiling `0.9429618090070259995…`. The independently recomputed cases are:

| Input choice | Optimal worst probability or ceiling |
|---|---:|
| Rounded efficiencies, g=10 | 92.83492610279933% |
| Rounded efficiencies, g=100 | 94.11478067757066% |
| Rounded efficiencies, g=10000 | 94.29430523275815% |
| Rounded efficiencies, principal ceiling for every g | 94.29618090070260% |
| Recompute efficiencies from central kpol/KD, principal ceiling | 94.30301542645167% |
| Favorable ±1 reported-error rectangle, principal ceiling | 97.19705369886649% |
| Favorable ±2 reported-error rectangle, principal ceiling | 98.88577273224246% |

These reproduce the supplied JSON at the quoted precision. The tiny central-column difference is ordinary rounding: one must specify whether the inputs are the printed efficiency column or ratios recomputed from the printed fitted parameters. Many displayed digits document arithmetic reproducibility; biologically, “about 94.30% under this model” is more appropriate.

**The two-template 99% rejection has an elementary exact certificate.** On U73 the correct incoming nucleotide is A. Retaining only the wrong C contribution, 99% requires

```text
x_A/x_C ≥ 99 × 0.00011/0.023 = 1089/2300.
```

On G73 the correct incoming nucleotide is C. Retaining only the wrong A contribution requires

```text
x_C/x_A ≥ 99 × 0.0035/0.030 = 231/20.
```

The left-hand products equal one, but the right-hand product is `251559/46000 = 5.468673913… > 1`. Contradiction. The value 5.46867 is a lower requirement obtained by omitting additional mismatches, not an equality for the full four-substrate problem. This certificate uses no censored data and no numerical optimizer. It alone excludes 99%; it does not establish the stronger 94.296% ceiling.

**The rectangle calculation is valid, with restricted interpretation.** Let m multiply each reported ± magnitude. For each retained correct reaction, its most favorable efficiency is `(k+m δk)/(K-m δK)`. For each mismatch its smallest allowed efficiency is `max(0,k-m δk)/(K+m δK)`. Taking mismatch-minimum / correct-maximum produces an entrywise lower matrix L_m for every allowed C. Therefore `rho(C)≥rho(L_m)` and fidelity is at most `1/[1+rho(L_m)]`. The choices for the retained entries can be made simultaneously in a Cartesian parameter box; there is no hidden incompatible-corner error.

At m=2, all three required correct-efficiency denominators remain positive. Some mismatch KD intervals cross zero; the parameter domain must be intersected with positive KD and nonnegative rates. That does not harm this upper-fidelity certificate, which uses the positive *upper* mismatch KD endpoints. This is a bound favorable to fidelity, not a worst-fidelity guarantee.

There is also an exact certificate for the rectangle's 99% rejection. Its characteristic cubic evaluated at the required odds `r=1/99` is

```text
f_2(1/99) = -2312052443/8570475248821200 < 0.
```

Its unique positive root therefore exceeds 1/99. The result does not depend on floating-point tolerances. The margin to 99% is only **0.1142273 percentage points**. As an additional diagnostic, the same favorable principal ceiling reaches 99% at `m≈2.08972479182477346`. This does not falsify the ±2 assertion; it shows how little additional expansion of this particular box removes the exclusion.

**The meaning of ± is not established well enough to call these “SDs” or “two sigma.”** The 2021 Figure 2 caption describes pooling three independent experiments into a hyperbolic fit; its methods say parameters derive from at least three independent assays per combination. In the inspected Table 1 image, caption, and methods, I found no explicit identification of its ± numbers as between-assay SD, SEM, or a specified regression standard error. Inferring a statistical definition from their numerical size would be unwarranted. [Patel et al., 2021](https://doi.org/10.1261/rna.078686.121).

The shared CSV headers `kpol_sd` and `KD_sd` consequently overstate the verified metadata. Rename them to reported error magnitudes in a future shared-data revision, unless a primary-source definition is obtained. This review leaves shared files untouched. Even established marginal SDs would not make a joint ±2 rectangle a 95% confidence region. Fit covariance, shared systematic effects, parameter positivity, censored data, and finite rounding are not calibrated by this exercise. Allowing parameters to vary independently means allowing Cartesian choices; it is not evidence for statistical independence. The proposed note's explicit “not a confidence interval” qualification is essential and correct.

**Biological predictive validity remains the weak link.** The measured quantities are apparent parameters from separate-NTP single-turnover assays on full-length, preactivated 5′-triphosphorylated tRNAHis variants. Their low-concentration slopes are plausible effective incorporation efficiencies. They are not direct measurements of four simultaneous independent capture hazards. For a simple common-state, competing-substrate mechanism with compatible kinetic commitments, an efficiency-weighted branching law can be justified. The word “dilute” alone does not establish that mechanism, its common starting state, or the transferability of apparent fitted KD values. Treating KD as an independently measured microscopic equilibrium dissociation constant would be an additional unsupported step. [Patel et al., 2021](https://doi.org/10.1261/rna.078686.121).

The model further assumes constant *free* substrate concentrations, negligible depletion, a single relevant incorporation decision, and no mixture-dependent effects on recognition or commitment. Multi-NTP binding associated with distinct activation and transfer roles, slow conformational states, product-dependent termination, and finite-time selection of products can break the mapping. A shared rapid-equilibrium competitive denominator would cancel from product flux ratios, so merely citing saturation does not by itself prove that the ratio law must change. Conversely, separately saturating four isolated hyperbolas and normalizing them is a different unvalidated model, not a correction established by these data.

**The 2011 competition data provide a concrete external check.** Rao, Maris and Jackman supplied GTP and a competing WC nucleotide simultaneously. On A73, yeast Thg1 produced approximately fivefold more G-addition than U-addition; BtTLP favored U over G by approximately fortyfold. The yeast extrapolation from the 2021 dilute coefficients instead gives G/U=`0.12/0.0055=21.818…`, independent of other competitors in this model. This is a roughly 4.4-fold discrepancy in that product ratio. However, the older assay used monophosphorylated RNA with ATP-dependent activation and high NTP concentrations, rather than the preactivated dilute case. It is evidence against unrestricted extrapolation, not a matched-condition falsification of the narrow theorem. [Rao et al., 2011, Figure 3 and accompanying results](https://doi.org/10.1093/nar/gkq976).

The native G:A reaction is biologically correct for tRNAHis maturation but counts as an error under the proposed WC-copying objective. Calling this a universal “Thg1 fidelity” ceiling obscures that change of objective. The four tRNA variants also do not establish an efficiency matrix for every position of a growing generic RNA. Nor does a single-addition correctness probability establish processivity, the sequence accuracy of recovered full-length products, or the absence of selection against poorly extended mismatch products. An unconstrained mixture can be scaled toward vanishing throughput while preserving its modeled fidelity.

**The relevant recent synthesis literature strengthens the scope restriction.** Jayasinghe, Patel and Jackman demonstrated targeted RNA 5′ extension using short duplexes, with AcaTLP2 and DdiTLP4 markedly more effective than ScThg1. Complete eight-nucleotide extension was demonstrated most clearly with alternating GC templates. Mixed/all-four and AU templates showed poorer extension; neither of those TLPs completed the full target under the tested mixed/AU conditions. These are distinct scaffolds and substrates, and product-length bands are not a four-base, position-resolved fidelity measurement. The paper already recommends these TLPs as candidates for further development. [Jayasinghe et al., 2024](https://doi.org/10.1261/rna.080156.124).

The 2025 metal-ion paper reports enzyme-specific metal dependence and evidence that the nucleotide destined for incorporation can influence the preceding activation step, alongside ATP. This supports coupled recognition beyond a one-decision capture description. It does **not** directly contradict an isolated transfer-step model on preactivated SceThg1. I inspected the accessible primary preprint and indexed journal discussion; the journal version is Biochemistry 64:3559–3569, DOI below, published online August 2025. [Iwaniec et al., 2025](https://doi.org/10.1021/acs.biochem.5c00162), [accessible preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC11957112/).

A 2026 journal article, not merely its 2025 preprint, must also be cited in a current assessment. Johnecheck, Kihm and Jackman describe altered WC/non-WC specificity and an activation-step recognition checkpoint in Dictyostelium enzymes. I verified the final abstract and journal metadata, and inspected the accessible preprint text; the final full text was inaccessible. Its results concern different enzymes and support biochemical specificity as a changeable property, without measuring this four-template mixture optimum. [Johnecheck et al., 2026, RNA 32:843–857](https://doi.org/10.1261/rna.080815.125), [2025 preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC12632743/). Publisher metadata identifies advance publication in February 2026; bibliographic indexing also lists May/June publication dates. A search engine's PDF download date is not a new scientific publication date.

**Prior art substantially limits novelty.** The source paper itself reports the GTP bias, secondary WC preferences, and comparatively poor WC discrimination. Thus, “wild-type SceThg1 is not a high-fidelity general WC copier,” and “consider TLPs or altered recognition instead,” are not discoveries of this calculation. Earlier primary work also demonstrated extended reverse repair and explicitly aimed at improving fidelity. [Patel et al., 2021](https://doi.org/10.1261/rna.078686.121), [Desai et al., 2018](https://doi.org/10.1080/15476286.2017.1372076).

The genuine extra content is quantifier strength: one common positive mixture, optimized over all ratios, cannot meet a specified uniform target in this particular model, even after removing dependence on the censored correct-G efficiency. That is more precise than observing poor equal-mixture or isolated-substrate performance. It is a direct application of established mathematics to a small published table, with a short certificate. Independent reproduction of a supplied result verifies it; it does not itself establish a novel scientific result.

Searches on 2026-11-26 covered exact numerical strings, Thg1 with Perron/Collatz/spectral/minimax terms, concentration/fidelity/optimization wording, bioRxiv-indexed material, GitHub-indexed code/file names and the source DOI, and patent text searches. I found no matching fixed-mixture bound or exact numerical ceiling in the inspected primary literature or returned relevant indexed results. Exact code/file-name searches returned no results. This is **not proof that nobody has posted, patented, or privately derived it**; public web indexing and accessible full text provide incomplete coverage. The unsuccessful exact-title metadata query was also not evidence of absence: a DOI search recovered the 2026 paper.

Patent searches returned a Thg1-mentioning tRNA-fragment application, US20250215504A1, whose inspected disclosure concerns HisGTG tRNAs rather than the minimax optimization, plus irrelevant conventional RNA-capping patents. No matching bound was identified. This was a limited scientific prior-art screen, not a patentability or freedom-to-operate assessment. [US20250215504A1, primary patent text](https://patents.google.com/patent/US20250215504A1/en). Query wording, findings, scope, and retrieval limits are in `sources/search_log.md`.

A defensible novelty statement would be: “We provide a reproducible, dataset-specific max-min analysis of Patel et al. Table 1 under an explicitly assumed competing-incorporation model, including a censoring-independent ceiling, an exact target-infeasibility certificate, and deterministic parameter-box sensitivity.” Add: “We did not identify a prior report of this particular analysis in the sources searched.” Avoid “first proof,” “fundamental reverse-polymerase fidelity limit,” or “new solution toward producing a reverse polymerase.”

**Useful falsification requires testing the kinetic mapping, not rerunning an optimizer.** The following comparisons would discriminate the model without changing nucleotide or enzyme sequences:

1. Reanalyze the published 2011 competition product identities/ratios against the model's predicted pairwise odds, explicitly separating activation and transfer assumptions. The discrepancy quantified above is already a useful diagnostic. A matched, preactivated low-concentration competition dataset would be much stronger evidence; I did not locate a complete public dataset covering all four template variants under those conditions.
2. For those four existing substrate variants, compare measured product identities in common mixtures with the prediction `P_i(j)/P_i(k)=(a_ij/a_ik)(x_j/x_k)`. Systematic failure of linear pairwise odds, mixture-dependent coefficients, or failure to converge on scale-invariant fractions in the dilute regime would falsify the proposed physical law. Agreement should be assessed against measurement uncertainty and the full identity distribution, not only the best-fitting worst probability. Absolute-rate measurements should accompany composition measurements to detect depletion and heterogeneous inactive fractions.
3. A measured uniform WC correctness above the claimed ceiling on the same preactivated substrates, with parameter estimates compatible with the stated box and confirmed model regime, would falsify the model-plus-parameter package. It would not refute the algebraic theorem. A value above 94.296% but below 98.8858% could instead reflect uncertainty in the point estimates. Incompatibility with the box alone would falsify its coverage as a description of the assay, not the rectangle theorem.
4. Use the 2024 public extension profiles to challenge any claim of context-independent transfer to generic duplexes. Their pauses and substrate-dependent outcomes warrant a model with extension/termination states. Those profiles cannot alone validate or refute this single-step probability bound, and their TLP results cannot be inserted into the SceThg1 table. A practical successor analysis would need position-resolved product identities and yield/processivity measurements on the intended scaffold, not a fitted 4×4 matrix borrowed from tRNAHis.

No experiment was performed, and these are comparison criteria rather than experimental or engineering recipes.

**Practical value is real but small.** The calculation rules out a concentration-only rescue of this particular measured behavior under a transparent approximation, exposes a conflict that persists without GTP-related entries, and makes the 99% exclusion auditable without an optimizer. It could support a methods note, a reproducibility supplement, or a baseline for a later mixed-substrate model. It does not materially resolve the dominant challenges of generic substrate recognition, extended synthesis, processivity, errors in full products, or switching template/product chemistries. The official challenge requires four reverse-direction polymerase classes with comparator-level processivity and errors; this analysis concerns one RNA/tRNA transfer system and supplies none of those catalysts. [Millennium problem 11, current challenge text](https://millenniumproblems.bio/).

**Critical corrections before wider use:** retain the model qualifier whenever quoting the ceiling; identify the substrate and enzyme in the result name; state that g=10 is a lower-limit substitution; distinguish printed efficiencies from recomputed ratios; replace unsupported SD metadata; present the ±2 result as deterministic with its narrow numerical margin; cite the 2011 competition evidence and 2024–2026 work; and restrict novelty to the particular analysis with an explicit search limitation. The supplied note already contains many of these cautions. Its mathematical core is stronger than its present biological validation or practical contribution.

**Reproduction and deliverables.** From this job directory, run:

```sh
/path/to/longevidence-checkout/research/millennium/.venv/bin/python check_bound.py
```

The check uses only the Python standard library, performs no experiments or network access, reads the two shared tabular/output files for comparisons, and writes only `outputs/independent_check.json`. Its image-transcribed numbers are embedded in the script. `outputs/check_stdout.txt` records the run. `fetch_sources.py` separately reproduces the literature metadata/retrieval attempts and records successes, hashes, and access failures in `sources/retrieval_log.json`; network availability can differ on rerun. Primary-source snapshots, relevant browser-returned evidence, source links, and search limitations accompany this report.
