# Mechanism-first computational opportunities in the twelve Millennium Problems for Biology

**Research date: 2026-11-26. Independent bounded scout; work confined to this directory.**

Recommend **02 Cryopreservation, 05 The quadruplet cell, and 11 5′ polymerases** for this mechanism-focused computational portfolio. I directly pursued all three. The strongest contribution is a reproducible, conditional kinetic bound derived from published Thg1 measurements: nucleotide-pool adjustment alone cannot raise the modeled worst-template Watson–Crick incorporation probability above **94.30%**, even allowing the censored favorable G:C efficiency to become arbitrarily large. This is a constraint on an explicitly stated model of an existing enzyme, not an experimental polymerase discovery.

The other completed contributions are an audit of actual nanowarming source measurements and an exact enumeration showing why equal mutation sensitivity does not by itself establish four-base amino-acid information or frame discrimination. All three analyses run locally with Python's standard library. **No Millennium challenge was solved, no biological sequence was validated, and no experiments were performed.**

## Scope and exact challenge definitions

I read the [canonical challenge website](https://millenniumproblems.bio/) live. Its exact twelve titles and full definitions are preserved in challenges_exact.json (archival input; see source-downloads.json and the source URLs), extracted from the downloaded primary HTML (archival input; see source-downloads.json and the source URLs). The following brief target sketches are not substitutes for those definitions. In particular, experimental controls, exclusions, preregistration, and ethics provisions remain binding.

| ID | Exact challenge title | Target sketch |
|---|---|---|
| 01 | Origins of life | Unassisted primordial cells; heredity; ≥10⁶ expansion and sustainable division. |
| 02 | Cryopreservation | Intact adult wild-type mice; ≥24 h cryostorage; >99% viability without lasting harm. |
| 03 | The reverse translatase | Untagged peptides→DNA/RNA; 100×≥50-aa pool; ≥90% accuracy; reads ≥25 aa, Q10. |
| 04 | Improve Rubisco | Paired controls: Galdieria partita specificity and maize turnover. |
| 05 | The quadruplet cell | Replicating quadruplet-only translation throughout; four positions matter; triplet decoding undetectable. |
| 06 | Somatic limb regeneration | Adult wild-type mouse limbs regain indistinguishable appearance, motor and sensory function. |
| 07 | Bacterial production of gene therapies | Bacterial infectious, replication-incompetent AAV and lentivirus; packaged genomes; benchmark-quality ratios. |
| 08 | Programmable proteases | Prospective intracellular cleavage of blinded endogenous sites; specificity; >80% of 20. |
| 09 | Cell-penetrating protein binders | Extracellularly administered zero-shot intracellular binders; 80% of 20 targets. |
| 10 | Protein amplification chain reaction | Monomers→arbitrary peptide copies; 100×≥50-aa pool; ≥1000× amplification; ≥90%/residue. |
| 11 | 5′ polymerases | Reverse DNA/RNA polymerase, RT, RDRP; Taq/T7/M-MLV/Φ6 comparators. |
| 12 | New nitrogenases | Nonhomologous sequence and fold; N₂→NH₃ at roughly native rates. |

For 08/09, the design deadline is 24 hours and target-specific laboratory optimization is excluded. The exact source distinguishes partial successes and expressly excludes nucleic-acid-mediated copying from 10.

## Ranking all twelve

I started with **mechanistic bottlenecks**, then asked which admit a useful falsifiable calculation on a laptop. This differs from ranking by the size of available datasets. D = public data relevant to the bottleneck; T = tractability of a meaningful additional computational result after the prior-art search; I = potential impact; V = clarity and accessibility of a path to validation. Scores are ordinal judgments from 1 to 5, with 5 most favorable. V concerns validation of a computational subproblem; whole-challenge validation is much harder.

Composite = **0.25D + 0.40T + 0.15I + 0.20V**. Arithmetic is provided for transparency, not objective precision or probabilities of success. The equal scores for 05/11 are a tie; their ordering is not substantive. A data-first ordering would place Rubisco and regeneration first. Greater weight on near-term experimental validation would favor Rubisco/proteases over this portfolio. [CSV with separate dimension ranks](outputs/challenge_ranking.csv).

| Rank | ID | Challenge | D | T | I | V | Composite | Mechanistic computational opportunity |
|---:|---:|---|---:|---:|---:|---:|---:|---|
| 1 | 02 | Cryopreservation | 4 | 5 | 5 | 3 | 4.35 | Spatial loading uncertainty, thermal margins, and failure certificates. |
| 2= | 05 | The quadruplet cell | 4 | 5 | 4 | 3 | 4.20 | Joint information, substitution, and reading-frame audits. |
| 2= | 11 | 5′ polymerases | 4 | 5 | 4 | 3 | 4.20 | Mixed-substrate selectivity bounds; separate activation, extension, and processivity. |
| 4 | 04 | Improve Rubisco | 5 | 3 | 5 | 4 | 4.00 | Assay-aware frontier uncertainty and phylogenetic leakage controls. |
| 5 | 08 | Programmable proteases | 4 | 3 | 5 | 4 | 3.75 | Separate peptide preference, folded-site accessibility, and cellular off-targets. |
| 6 | 06 | Somatic limb regeneration | 5 | 3 | 5 | 2 | 3.60 | Cross-species positional-memory programs with spatial/lineage controls. |
| 7 | 09 | Cell-penetrating protein binders | 3 | 3 | 5 | 3 | 3.30 | Separate uptake, cytosolic access, target engagement, and cell damage. |
| 8 | 01 | Origins of life | 3 | 3 | 5 | 2 | 3.10 | Coupled replication, strand separation, compartment growth, and dilution. |
| 9 | 07 | Bacterial production of gene therapies | 3 | 2 | 5 | 3 | 2.90 | Assembly/packaging bookkeeping and evidence-gap analysis. |
| 10 | 12 | New nitrogenases | 3 | 2 | 5 | 2 | 2.70 | Mechanistic/cofactor constraints and rigorous exclusion of homologues. |
| 11 | 03 | The reverse translatase | 2 | 2 | 5 | 2 | 2.45 | Information-transfer requirements and recognition/encoding coupling. |
| 12 | 10 | Protein amplification chain reaction | 1 | 1 | 5 | 1 | 1.60 | Error-threshold and product-inhibition bounds, presently weakly parameterized. |

These are rankings of **this worker's computational opportunity**, not rankings of scientific importance. High T for a rigorous constraint does not imply a short route to a functioning biological system.

## Pursuit 1: cryopreservation—local feasibility rather than bulk heating speed

**Mechanism.** Rewarming can fail through ice growth, thermomechanical fracture, or prior cryoprotectant-associated injury. Fast volumetric heating addresses only part of this coupled problem. Nanoparticles follow vasculature; cryoprotectants must also reach extravascular regions. Source uniformity, CPA equilibration, and tolerable exposure therefore cannot be collapsed into a single average temperature slope.

The relevant published progression is substantial. [Han et al. 2023](https://doi.org/10.1038/s41467-023-38824-8) recovered rat kidneys stored for up to 100 days, with life-supporting transplantation. The reported adjacent-solution warming rate was 72.0 ± 8.0 °C/min; that probe is not a measurement of every tissue voxel. [Gangwar et al. 2025](https://doi.org/10.1038/s41467-025-63483-2) demonstrated bulk physical vitrification up to 3 L and nanowarming up to 2 L; their vitrified porcine liver was not nanowarmed. Neither result establishes whole-body recovery.

The **2026-10-14 pig-kidney preprint** is essential recent prior art: [Functional recovery of vitrified and nanowarmed human-scale pig kidneys](https://doi.org/10.64898/2026.10.07.756990). It reports approximately 100-g organs, one stored for 366 days, and acute transplantation with n=3 per group over two hours. Creatinine clearance after transplantation was 1.0 ± 1.0 mL/min/100 g versus 47.8 ± 36.3 for fresh grafts. Thus urine production and immediate reperfusion should not be equated with normal sustained filtration or long-term survival. This remains a preprint. Its posted date was independently verified with the official bioRxiv API, rather than inferred from the DOI's date fragment.

**Real data reanalysis.** I downloaded the original 2025 Source Data XLSX and extracted Fig.7 measurements. Three probes within an experiment are clustered into the same replicate; they are not treated as independent n=9 experiments. [Code](analysis/cryo_source_audit.py), extracted observations (archival input; see source-downloads.json and the source URLs), [output](outputs/cryo_source_audit.json).

| Bulk M22 volume | Fe concentration, mg/mL | Mean rate ± SD of three replicate means, °C/min | Slowest listed probe rate | Mean of six listed ΔT values | Largest listed ΔT |
|---|---:|---:|---:|---:|---:|
| 1 L | 10.7 | 173.52 ± 3.14 | 168.85 | 4.55 °C | 6.3 °C |
| 2 L | 4.6 | 87.96 ± 1.51 | 84.29 | 3.63 °C | 8.3 °C |

The prose's approximate sub-5 °C uniformity description must not become a hard upper bound in an optimizer: individual source observations reach 6.3/8.3 °C. This is an **extreme-versus-average audit**, not evidence that the experiments failed. Rate per iron concentration is 17.9% higher in the 2-L measurements, so volume alone does not explain the rate difference. Concentration, temperatures, and measurement conditions matter. The slowest listed 2-L rate is 210.7 times the nominal 0.4 °C/min bulk-M22 critical warming rate, but this ratio cannot certify underloaded tissue, diluted CPA, or an unmeasured cold spot.

**A concrete feasibility screen.** Start from the heat equation

\[
\rho(c,T)c_p(c,T)\partial_tT
=\nabla\cdot[k(c,T)\nabla T]
+u(t)\,C_{\rm Fe}(x)\,\operatorname{SAR}(T,H,f).
\]

Here c is local CPA composition and u is a common duty control. A useful robust extension is to constrain the minimum local warming margin over uncertain loading fields, while separately checking tensile stress and CPA exposure. A pointwise critical-rate constraint is a screening approximation to a measured thermal history criterion, not a complete ice-nucleation model.

For an exactly solvable limiting screen, neglect conduction, assume constant local heat capacity and source ratios, and write \(T_i=T_0+a_i U\), where \(U=\int u(t)dt\). Requiring a cold region to reach \(T_{\rm exit}\) while a hot region stays below \(T_{\max}\) is feasible only if

\[
\frac{a_h}{a_c}\leq
\frac{T_{\max}-T_0}{T_{\rm exit}-T_0}.
\]

With **illustrative chosen endpoints**, −150, −40, and 0 °C, the ratio is 1.364. A ratio of 1.5 would put the hot region at +15 °C when the cold region first reaches −40 °C. Changing a common scalar duty history cannot repair this limiting case; redistribution or effective heat transport is necessary. The endpoints are not asserted to be M22 material thresholds. Temperature-dependent SAR/heat capacity, conduction, boundaries, and phase changes can change the conclusion. This calculation is a screening certificate under assumptions, not a whole-organ impossibility theorem.

**Prior art and what remains useful.** [Joshi and Rabin 2023](https://doi.org/10.1371/journal.pone.0290063) already coupled loading-dependent thermal fields to Maxwell thermomechanics. The public [Vitrify repository](https://github.com/Marc-Dvci/Vitrify), created August 2026, already implements inverse loading/duty design and resolution margins. I inspected its actual `properties.py` and `thermal.py`. Its own documentation reports an unresolved low-temperature SAR modeling ambiguity. I did not reproduce that repository's optimization or treat its reported gains as experiments. “Optimize nanowarming with differentiable simulation” is therefore **not a new proposal**.

The next defensible contribution is a **formulation-specific uncertainty benchmark**: perturb correlated CPA delivery and iron loading, quantify when nominally feasible protocols become infeasible, and identify which measurement would most reduce uncertainty. Validate first against held-out source thermometry and physical CPA phantoms; validate injury predictions separately against appropriate published biological outcomes. [55-MHz dielectric warming](https://doi.org/10.1016/j.cryobiol.2025.105257) is a substantive alternative because it heats polar molecules rather than relying on particle delivery. It still requires spatial electrical/material-property characterization. Comparing both mechanisms under the same measured loading uncertainty would be more informative than ranking nominal mean heating rates.

## Pursuit 2: quadruplet codes—information, detection, and synchronization are different

**Mechanism.** A four-base word has 256 possible states, but additional states do not automatically provide four-position amino-acid information, tolerance of substitutions, or resistance to frame errors. A functioning code also needs a decoding apparatus. Published evidence includes [four distinct quadruplet codons in one reporter in 2021](https://doi.org/10.1038/s41467-021-25948-y), [improved orthogonal decoding in 2024](https://doi.org/10.1038/s41587-024-02385-y), and [mitochondrial UAGA suppression in 2024](https://doi.org/10.1021/acssynbio.4c00630). These are important partial mechanisms; they do not establish a quadruplet-only cell. The 2024 mitochondrial paper also reports that a similarly modified cysteinyl-tRNA failed, illustrating that anticodon expansion is not uniformly tolerated.

**Audit contribution.** For mapping f on valid sense words S, I computed, under uniform source-codon use and uniform alternative bases,

\[
p_i=\frac{1}{3|S|}\sum_{c\in S}\sum_{b\ne c_i}
1\{f(c)\ne f(c_{i\leftarrow b})\},\qquad
h_i=H(A\mid C_{-i},C\in S).
\]

An invalid/null result is counted in p as a change but is reported separately as an erasure. This convention is deliberately explicit; it is not assumed to be the challenge adjudicator's definition of nonsynonymous. The conditional entropy h measures whether the omitted base supplies additional amino-acid information **among valid codewords**. It is a diagnostic, not a new official acceptance threshold. I also enumerated every ordered sense-codon pair and all three internal four-base windows crossing its boundary.

| Mathematical control | Sense words / amino acids | p₁,p₂,p₃,p₄ | h₁,h₂,h₃,h₄, bits | Shifted-window sense acceptance, offsets 1/2/3 |
|---|---|---|---|---|
| Artificial triplet mapping + padding | 256 / 20 | 1,1,1,0 | 2,2,2,0 | 1 / 1 / 1 |
| Sparse parity subset | 20 / 20 | 1,1,1,1 | 0,0,0,0 | 0.240 / 0.080 / 0.1125 |
| Dense rotation-symmetric initial map | 240 / 20 | 0.9667 at all positions | 1.8170 at all positions | 0.950 / 0.9375 / 0.950 |
| Same dense structure, heuristic optimized assignment | 240 / 20 | 0.8167 at all positions | 1.4604 at all positions | 0.950 / 0.9375 / 0.950 |

[Code](analysis/quadruplet_audit.py), [exact results](outputs/quadruplet_audit.json), [256-entry final map](outputs/quadruplet_dense_map.csv).

The sparse code chooses 20 words from \(\sum_i c_i=0\pmod4\), mapping every other word to null. Every single substitution is detected as an erasure, yet any three bases already identify the amino acid. Thus **100% balanced change sensitivity can coexist with zero additional fourth-base amino-acid information**. A checksum-like fourth base is useful for detection, but would not establish compliance with the explicit rejection of encodings whose first three positions determine the amino acid. This distinction is the principal actionable audit finding.

The dense control uses the 60 aperiodic rotation orbits of length four: assign three whole orbits to each of 20 amino acids and the remaining 16 periodic words to stop. Rotation symmetry guarantees equal p and h across positions under the stated weighting. Swapping orbit labels retains that guarantee and equal 12-word degeneracy per amino acid. With seed 20261126, 120,000 annealing steps reduced p from 0.9667 to 0.8167 while h stayed positive. The result is a heuristic example; it minimizes symbolic label changes, not a validated physicochemical or fitness cost. Crucially, **93.75–95% of shifted windows still look like sense words**. Informational balance has not created frame discrimination. Rotation invariance also makes cyclicly rotated individual words synonymous, another reason this control is not a recommended biological design.

A separate textbook bound rules out one tempting objective: correcting **every** single-base substitution in an isolated four-base word for all 20 amino-acid identities, without outside information. Radius-one Hamming balls contain \(1+4(4-1)=13\) words. Twenty mutually distinguishable identities would need at least 260 received words, exceeding 256. This applies even if one first selects one representative per identity; synonymous alternatives do not rescue that universal guarantee. It does not rule out detection, partial correction, unequal-error protection, or context-dependent decoding. This is established coding theory, not a novelty claim.

**Prior art and validation.** [Calles et al. 2019](https://doi.org/10.1093/nar/gkz745) already designed fail-safe quadruplet codes. I inspected the [actual generator](https://github.com/EndyLab/codon-tables/blob/manuscript/src/ffgen.py); parity-style sparse coding is not new. The [October 2026 triplet-proof v2 archive](https://doi.org/10.5281/zenodo.23125462) also discusses quadruplet robustness and explicitly retracts an earlier null-model interpretation. Its peer-review status was not established and its GitHub URL returned 404; I did not reproduce its numerical claims.

A useful next result would be a certified tradeoff frontier jointly constraining h, position-dependent mutation rates, observed codon–anticodon cross-decoding, and shifted-window rejection. Freeze the weighting before comparisons and include synonymous-block-preserving nulls. Biological validation would require measuring decoder confusion and triplet leakage separately, followed by translation-system-wide tests. Our synthetic map proves only enumerated properties of a mapping; neither a decoder nor a cell was built.

## Pursuit 3: reverse polymerases—a joint selectivity bound from real kinetic measurements

**Mechanism and prior art.** Reverse 3′→5′ synthesis is chemically possible. Thg1/TLP enzymes extend the growing 5′ end: the incoming nucleotide's 3′ hydroxyl attacks the activated end, and the incorporated nucleotide supplies the next terminal triphosphate. Initial activation of a monophosphorylated substrate is distinct from subsequent extensions; ATP activation is **not required anew at every ordinary extension**. [Desai et al.](https://doi.org/10.1080/15476286.2017.1372076) reported up to 13 nucleotides of repair. An endpoint extension length does not by itself measure processivity.

[Jayasinghe et al. 2024](https://doi.org/10.1261/rna.080156.124) already tested simple RNA duplexes, finding AcaTLP2 and DdiTLP4 comparatively capable. Reported single-turnover product-formation rates include 0.63 min⁻¹ for AcaTLP2 on an 8-bp duplex and 3.5 min⁻¹ for DdiTLP4 on a 7-bp duplex. These rates summarize formation of any extension product, not nucleotide-per-second velocities. [Iwaniec et al. 2025](https://doi.org/10.1021/acs.biochem.5c00162) distinguished activation and incorporation metal dependence. A [2026 comparison](https://doi.org/10.1261/rna.080815.125) further narrows novelty around catalytic residue comparisons; its full journal text was subscription-gated here, so I verified metadata/discovery but do not claim a full-text assessment. General RNA-duplex extension or metal/residue screening is already established work.

**Primary dataset.** [Patel et al. 2021](https://doi.org/10.1261/rna.078686.121), Table 1, supplies 16 single-addition efficiencies \(E_{tn}=k_{\rm pol}/K_{D,\rm app}\) for four discriminator templates and four NTPs. I saved the complete primary XML, manually transcribed all table entries, retained reported rounding, and kept the favorable template-C/GTP efficiency as **≥10**, not an exact estimate. The measurements use excess enzyme and preactivated tRNA substrates, rather than generic long templates.

**Conditional model.** Assume mutually exclusive dilute-substrate capture events with effective hazards \(\lambda_{tn}=E_{tn}x_n\), where x is one shared positive NTP pool. This is a testable empirical approximation; single-NTP kinetic fits do not themselves prove that this mixed-NTP competition model holds. With w(t) the Watson–Crick complement,

\[
P_t(x)=\frac{E_{t,w(t)}x_{w(t)}}{\sum_n E_{tn}x_n}.
\]

Reindex rows by correct incoming nucleotide i and define

\[
B_{ij}=\begin{cases}E_{\operatorname{comp}(i),j}/E_{\operatorname{comp}(i),i},&i\ne j\\0,&i=j.\end{cases}
\]

Then \(P_i=1/[1+(Bx)_i/x_i]\). The Collatz–Wielandt characterization gives the global optimum

\[
\boxed{\max_{x>0}\min_t P_t(x)=\frac{1}{1+\rho(B)}}.
\]

The Perron vector supplies the pool proportions; this is a global solution under the model, not a local concentration search. The implementation verifies componentwise spectral bounds with width below 10⁻⁹. [Code](analysis/thg1_competition.py), 16-row dataset (archival input; see source-downloads.json and the source URLs), [full output](outputs/thg1_competition.json).

| Diagnostic | Result |
|---|---|
| Equal-pool correct probability, template A | 4.365% |
| Equal-pool correct probability, template C | ≥99.9929%, because favorable efficiency is censored |
| Equal-pool correct probability, template G | 15.485% |
| Equal-pool correct probability, template U | 5.302% |
| Minimax, setting censored efficiency to its lower bound 10 | 92.8349% worst-template probability |
| Corresponding A/C/G/U pool fractions | 0.229423 / 0.403321 / 0.0002675 / 0.366988 |
| Upper bound independent of favorable G:C efficiency | **94.2962%** |

The last bound uses only the principal submatrix indexed by incoming A,C,U:

\[
B_{ACU}=\begin{pmatrix}
0&0.00478261&0.0300000\\
0.1166667&0&0.0080000\\
0.08363636&0.00818182&0
\end{pmatrix},\quad \rho=0.0604883363.
\]

Its spectral radius lower-bounds the full matrix's radius. Therefore even an arbitrarily favorable censored G:C reaction cannot remove competition among the other three correct incorporations. Holding the other reported point estimates fixed, the possible optimum lies between **92.8349% and 94.2962%** as that censored efficiency varies. This interval is not an experimental confidence interval.

An elementary cross-check uses just two templates. For 99% correctness, template U requires A/C ≥0.47348, while template G requires C/A ≥11.55; their product is 5.4687, impossible for reciprocal ratios. The 99% target is an analyst-chosen diagnostic, not a reinterpretation of the challenge's comparator requirements. The stronger three-component bound accounts for the additional conflicts.

**Interpretation and limitations.** Pool adjustment can substantially improve a badly mismatched substrate environment, but it cannot eliminate the enzyme's residual competing preferences in this model. The result concerns Watson–Crick synthesis on arbitrary templates; SceThg1 normally performs biologically useful non-Watson–Crick G addition to tRNA and is not simply “defective.” Pool ratios are mathematical outputs, not validated operating conditions. Absolute concentrations, saturation, apparent-binding kinetics, activation, depletion, inhibition, and extension-context changes remain unresolved. The tiny modeled GTP fraction may sacrifice useful rate, and the relevant binding constant is censored. No long-copy fidelity or processivity is inferred.

As a sensitivity check, 2,000 independent lognormal perturbations of reported kpol/KD parameters, with the censored row fixed at 10, yielded a 2.5–97.5 percentile range of **0.894–0.948** for the model optimum. This is explicitly **not a confidence interval**: raw replicate time courses and fit covariances were not obtained, and independence/distribution assumptions are imposed.

**Testable contribution and prior-art caveat.** Measure four-template single-addition outcomes under common mixed pools, initially with preactivated substrates so activation is not conflated with selection; compare locked model predictions with held-out mixed-pool counts and dwell times. Then obtain the analogous 4×4 matrix for the more promising TLPs before ranking enzyme candidates. Processivity needs independent dissociation/trap measurements; error rates need actual product sequencing; DNA polymerase/RT capability needs separate evidence. No experiment is proposed as having occurred here.

Competition itself is old prior art: [Rao et al., online 2010 / issue 2011](https://doi.org/10.1093/nar/gkq976) already performed two-NTP competition assays, and the 2021 paper related relative efficiencies to in-vivo incorporation. The contribution here is the joint four-template mathematical reanalysis, explicit censoring-independent bound, and reproducible audit. Searches for Thg1 mixed-pool optimization/Perron analyses did not identify an exact predecessor in this bounded pass; that is **not proof of novelty**.

## Mechanistic assessment of the remaining nine

**04 Rubisco.** This has the strongest readily usable kinetics resource: [Flamholz et al. 2019](https://doi.org/10.1021/acs.biochem.9b00237), approximately 300 enzymes, with [public data/code](https://github.com/flamholz/rubisco). [Sequence-based prediction](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) is already prior art, with a public [Mining-Rubisco-kinetics repository](https://github.com/Iqbalwasim01/Mining-Rubisco-kinetics). A useful alternative is to estimate the probability of frontier dominance after assay/batch corrections and phylogenetic holdout, not simply regress specificity on sequence. The bottleneck is simultaneous catalytic selectivity, turnover, assembly, and activation; a strong correlation is not an immutable physical tradeoff. Paired biochemical controls give a clear validation path. Separate software reuse licenses were not verified, so accessibility is not treated as permission for unrestricted reuse. This is the first reserve choice, displaced mainly because simple data/model analyses are crowded prior art.

**08 Proteases.** Recent evidence substantially improves plausibility: [Choi et al. 2026](https://doi.org/10.1038/s41586-026-11143-w) report 13/69 active de novo cysteine-protease designs, up to 3×10⁷ background rate enhancement, and cellular cleavage. The paper links [design](https://github.com/hojaec/RFD2_MI_cysteine_protease_pipeline) and [analysis](https://github.com/hojaec/Cysteine-protease-data-analysis) repositories. That hit rate does not establish the blinded endogenous-site challenge. [CleaveNet, published 2026-01-06](https://doi.org/10.1038/s41467-025-67226-1), uses approximately 18,500 peptide substrates across 18 MMPs and offers [data/code](https://github.com/microsoft/cleavenet); substrate design is distinct from designing a new enzyme. A laptop project could quantify whether accessibility and competing proteome sites explain failures better than isolated peptide scores. Validate against held-out endogenous folded substrates and off-target cleavage, with prospective design timing. Dataset existence is verified; repository-specific licenses and large models were not audited or downloaded.

**06 Regeneration.** [Zhong et al. 2023](https://doi.org/10.1038/s41467-023-41944-w) combine five vertebrates and spatial evidence, including [GSE243225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243225). Developmental AER programs are not simply replayed during axolotl regeneration. [Molecular basis of positional memory in limb regeneration, 2025](https://doi.org/10.1038/s41586-025-09036-5) further makes positional memory a concrete mechanistic target. Public transcriptomes enable cross-species module tests, but developmental stage, positional identity, cell composition, orthology, and batch must be controlled. Validation should progress from independent perturbation/spatial evidence to functional regeneration; a mouse expression signature is not a regenerated limb. Public sequencing access is established; individual dataset reuse terms must be checked before redistribution. No large atlas was downloaded.

**09 Intracellular binders.** Binding and entry are serial requirements, with endosomal escape a separate failure point. [SLEEQ](https://doi.org/10.1038/s41467-021-23997-x) provides a quantitative escape assay; [Comparing cell penetration of biotherapeutics across human cell lines](https://doi.org/10.1021/acschembio.4c00211), 2024, compares CAPA measurements across seven cell lines plus primary cells. The opportunity is a measurement-aware transfer model that predicts cytosolic engagement from extracellular dose while detecting membrane injury. Uptake fluorescence alone is inadequate. Some apparent delivery methods use membrane destabilization and would need careful challenge-eligibility review; they cannot be assumed to satisfy the allowed-entry mechanism. Validate with independent cytosolic localization and target-engagement measurements at intact-membrane conditions. Public papers/supplements exist; those assay datasets were not downloaded here.

**01 Origins.** Strand separation and product inhibition are mechanistic opportunities. [Trinucleotide substrates under pH–freeze–thaw cycles enable open-ended exponential RNA replication by a polymerase ribozyme](https://doi.org/10.1038/s41557-025-01830-y), 2025, and [QT45](https://doi.org/10.1126/science.adt2760), online 2026-02-12, are relevant recent primary work. The latter is a 45-nt ribozyme using triplet substrates; this is not unassisted emergence of RNA/protein cells. A useful model must jointly track replication, duplex trapping, substrate/energy replenishment, heritable sequence, compartment growth, and dilution. Validation requires predicting independent cycle trajectories and persistence under transfers, not a single favorable rate equation. Author manuscripts are publicly accessible; no claim that all experimental precursors or conditions are prebiotically established is made.

**07 Bacterial gene-therapy production.** [Le et al. 2022](https://doi.org/10.1021/acssynbio.2c00451) demonstrate AAV5 virus-like particle production in E. coli with assembly and cellular internalization evidence. This is not evidence of correctly packaged, infectious gene-delivery vectors. Earlier [defined AAV capsid assembly](https://doi.org/10.1038/s41598-019-54928-y), 2019, is also known. The missing functions differ from making an empty shell; lentiviral envelope formation adds another host-compatibility barrier. A safe, useful computational scope is evidence/stoichiometry accounting that explicitly distinguishes capsids, protected genomes, and transducing units. Validation must measure those separately against mammalian-produced references. No viral constructs, genome designs, or production experiments were generated here. Primary literature is accessible; a standardized public bacterial full-vector performance dataset was not located in this pass.

**12 Nitrogenases.** Breaking N₂ while managing electrons/protons, ATP expenditure, cofactor assembly and competing H₂ production is the core mechanism. [Maslać et al. 2026](https://doi.org/10.1038/s41467-026-77173-0) provide a structurally simplified hyperthermophilic nitrogenase, still clearly homologous to known nitrogenases. Reported ammonia-production activity at 50 °C is 280 ± 18.2 nmol/min/mg; higher-temperature activity in that paper is extrapolated, not measured. Structural databases support homology exclusion and electron-transfer/cofactor geometry analyses, but those do not identify an active unrelated scaffold. [Ancestral resurrection](https://doi.org/10.7554/eLife.85003) likewise remains within the known family. Future validation needs direct N₂-derived ammonia with stringent contamination/isotope controls and matched turnover. Sequence/structure exclusions alone cannot validate nitrogen fixation. No general nonhomologous-enzyme training dataset was established.

**03 Reverse translatase.** A protein reader must couple twenty-way residue discrimination, processive movement, and covalent codon output. [Protein-templated synthesis of dinucleotide repeat DNA by an antiphage reverse transcriptase](https://doi.org/10.1126/science.aed1656), online 2026-04-16, is genuine recent primary evidence, but its protein side chains enforce a repetitive product rather than read an arbitrary external peptide. It should not be represented as solving reverse translation. A computational subproblem is an information/kinetic error budget for a coupled reader/writer, benchmarked against measured recognition confusion where available. For a uniform 20-letter residue source, each residue carries log₂20 ≈4.322 bits before coding redundancy; recognition and writing errors consume the same budget. Validation needs sequence-disjoint arbitrary-peptide pools and actual encoded reads; database retrieval, peptide barcoding, or generating an enzyme sequence is not that validation. Appropriate end-to-end training pairs were not located.

**10 Protein amplification.** Template-directed ligation, foldamer replication, and fibril seeding are important neighboring mechanisms, but sequence-faithful arbitrary copying from amino-acid monomers requires different information transfer. [Light-Driven Competitive Selection in a Protein-Catalyzed Dissipative Peptide Replication](https://doi.org/10.1002/anie.202518911), 2026, demonstrates a foldamer/disulfide-based replication network with calmodulin; it is not the target monomer-to-arbitrary-polypeptide copying reaction. Public primary data can parameterize specific autocatalytic networks, not a general copying enzyme. A growth/error-threshold analysis is useful only once substrate-specific copying and product release are measured. Validate covalent sequence and new-monomer incorporation, distinguishing exponential synthesis from aggregation or templated assembly. No dataset of arbitrary faithful peptide copying was found; negative search evidence is limited.

## Deliverables, evidence limits, and next decision

The recommendation is to prioritize the **Thg1 mixed-pool bound for further computational/mixed-substrate testing**, the **quadruplet audit as a mandatory design-screening tool**, and the **cryopreservation source audit as the starting point for robust uncertainty benchmarks**. Their value is in ruling out misleading interpretations and identifying measurements that would distinguish mechanisms.

The closest next computational work is to fit a mixed-substrate competition model with genuine multi-NTP measurements, or to extend the quadruplet audit to empirically observed cross-decoding matrices. Both have sharper failure criteria than a new unconstrained generative sequence. For cryopreservation, do not launch a new nominal inverse-design project before testing whether loading/material uncertainty overwhelms its claimed gain.

Available artifacts:

- [Reproduction instructions](analysis/README.md), three analysis scripts, deterministic JSON outputs, and [rankings CSV](outputs/challenge_ranking.csv).
- Data/access/license notes (archival input; see source-downloads.json and the source URLs), downloaded small primary workbook, transcribed kinetics, and extracted cryopreservation observations.
- [Durable research log](notes.md), exact canonical definitions, inspected public source-code snapshots, and SHA-256 manifest (archival input; see source-downloads.json and the source URLs).
- [Refusal log](refusals.md): **none observed**. HTTP errors, access checks, DNS failure, GitHub 404, and subscription access limits are recorded separately as access errors.

This was a bounded search of primary papers, preprints, code, and datasets through 2026-11-26, not an exhaustive patent or literature review. Key recent 2026 papers were verified against primary records; DOI years and search-engine relative dates were not used alone to date publications. Full texts were read where available; gated or metadata-only items are identified. Numerical calculations are original exploratory analyses under explicit assumptions; known mathematical tools and prior biological mechanisms are credited. **No blanket novelty claim, experimental validation claim, or challenge-completion claim follows from this report.**

## Exact-title reference register for the directly pursued work

All links below were checked on 2026-11-26 through primary records or their public author manuscripts. Access qualifications and licenses are detailed in data/README.md.

1. **Vitrification and nanowarming enable long-term organ cryopreservation and life-sustaining kidney transplantation in a rat model**. Han et al.; 2023-06-09; Nature Communications. DOI [10.1038/s41467-023-38824-8](https://doi.org/10.1038/s41467-023-38824-8). Primary XML downloaded; CC BY 4.0.
2. **Physical vitrification and nanowarming at liter-scale CPA volumes: toward organ cryopreservation**. Gangwar et al.; 2025-09-26; Nature Communications. DOI [10.1038/s41467-025-63483-2](https://doi.org/10.1038/s41467-025-63483-2). Primary XML and original XLSX downloaded; CC BY 4.0. This is the journal successor to preprint 10.1101/2024.11.08.622572.
3. **Thermomechanical stress analyses of nanowarming-assisted recovery from cryopreservation by vitrification in human heart and rat heart models**. Joshi and Rabin; 2023-08-16; PLOS ONE. DOI [10.1371/journal.pone.0290063](https://doi.org/10.1371/journal.pone.0290063). Primary full text read; CC BY.
4. **Functional recovery of vitrified and nanowarmed human-scale pig kidneys**. Han et al.; posted 2026-10-14, v1; bioRxiv preprint. DOI [10.64898/2026.10.07.756990](https://doi.org/10.64898/2026.10.07.756990). Indexed primary text and official metadata read; CC BY-NC-ND 4.0; not peer reviewed.
5. **55 MHz constant field dielectric warming of kidneys and ovaries cryopreserved by vitrification**. Cryobiology; online 2025-05-08. DOI [10.1016/j.cryobiol.2025.105257](https://doi.org/10.1016/j.cryobiol.2025.105257). Public [author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC12161498/) inspected; separate reuse license not asserted.
6. **Fidelity of base-pair recognition by a 3′–5′ polymerase: mechanism of the Saccharomyces cerevisiae tRNAHis guanylyltransferase**. Patel et al.; RNA, June 2021. DOI [10.1261/rna.078686.121](https://doi.org/10.1261/rna.078686.121). Primary XML read; Table 1 transcribed; CC BY-NC 4.0 after embargo.
7. **Minimal requirements for reverse polymerization and tRNA repair by tRNAHis guanylyltransferase**. Desai et al.; online 2017-09-29, issue 2018; RNA Biology. DOI [10.1080/15476286.2017.1372076](https://doi.org/10.1080/15476286.2017.1372076). Public primary full text inspected; no raw data download or blanket reuse-license claim.
8. **Thg1 family 3′–5′ RNA polymerases as tools for targeted RNA synthesis**. Jayasinghe, Patel and Jackman; RNA, October 2024. DOI [10.1261/rna.080156.124](https://doi.org/10.1261/rna.080156.124). Primary XML read; CC BY-NC 4.0 after embargo.
9. **Metal Ion Requirement for Catalysis by 3′-5′ RNA Polymerases**. Iwaniec, Allegretti and Jackman; online 2025-08-07, issue 2025-08-19; Biochemistry. DOI [10.1021/acs.biochem.5c00162](https://doi.org/10.1021/acs.biochem.5c00162). Published record verified against preprint [10.1101/2025.03.21.644660](https://doi.org/10.1101/2025.03.21.644660); primary preprint full text read. Do not confuse that preprint with the journal version.
10. **A comparison of Dictyostelium discoideum 3′-5′ RNA polymerases reveals a conserved tRNAHis guanylyltransferase residue that plays a dual role in catalysis**. Johnecheck, Kihm and Jackman; journal advance publication 2026-02-04, June 2026 issue. DOI [10.1261/rna.080815.125](https://doi.org/10.1261/rna.080815.125). Journal metadata checked; full text gated. PubMed displays an additional May date; the publisher's advance/issue dates are distinguished here.
11. **tRNA 5′-end repair activities of tRNAHis guanylyltransferase (Thg1)-like proteins from Bacteria and Archaea**. Rao, Maris and Jackman; online 2010-11-03, 2011 issue; Nucleic Acids Research. DOI [10.1093/nar/gkq976](https://doi.org/10.1093/nar/gkq976). Primary competition methods inspected.
12. **Fail-safe genetic codes designed to intrinsically contain engineered organisms**. Calles et al.; online 2019-09-12, issue 2019-11-04; Nucleic Acids Research. DOI [10.1093/nar/gkz745](https://doi.org/10.1093/nar/gkz745). Primary article CC BY 4.0; accompanying code's separate license unverified.
13. **Multiplex suppression of four quadruplet codons via tRNA directed evolution**. DeBenedictis et al.; 2021-09-29; Nature Communications. DOI [10.1038/s41467-021-25948-y](https://doi.org/10.1038/s41467-021-25948-y). Primary paper and current-document record checked.
14. **Efficient genetic code expansion without host genome modifications**. Costello et al.; online 2024-09-11, issue July 2025; Nature Biotechnology. DOI [10.1038/s41587-024-02385-y](https://doi.org/10.1038/s41587-024-02385-y). Primary results and source-data availability inspected; correction dated 2024-10-07 noted. Publisher retains exclusive article rights; an accessible author manuscript is not an unrestricted license.
15. **Towards a Quadruplet Codon Mitochondrial Genetic Code**. Pigula, Ban and Schultz; online 2024-12-04, issue 2024-12-20; ACS Synthetic Biology. DOI [10.1021/acssynbio.4c00630](https://doi.org/10.1021/acssynbio.4c00630). Public author manuscript includes supporting sequence/control data; separate unrestricted reuse license not verified.
16. **triplet-proof v2.0.0: code and data for “Error robustness of doublet and triplet genetic codes”**. Ribbeck, Paudyal and Kiran; 2026-10-03; Zenodo repository/manuscript record. DOI [10.5281/zenodo.23125462](https://doi.org/10.5281/zenodo.23125462). Corrected v2 statement read; numerical archive not reproduced.
17. **Error Detecting and Error Correcting Codes**. R. W. Hamming; April 1950; Bell System Technical Journal. DOI [10.1002/j.1538-7305.1950.tb00463.x](https://doi.org/10.1002/j.1538-7305.1950.tb00463.x). [Primary paper copy](https://bioinfo.uib.es/~joemiro/comptrans/HammingIntro.pdf); used for established coding-theory context. The quaternary sphere-packing count is derived explicitly above.

Additional bibliographic records used in the nine-challenge assessment follow. Dates are publication dates where verified; year-only entries do not assert an unverified day. No secondary news claim is used as biological evidence.

| Exact title | Date | DOI / primary access | Data/access qualification |
|---|---|---|---|
| Revisiting Trade-offs between Rubisco Kinetic Parameters | Online 2019-07-22; issue 2019-08-06 | [10.1021/acs.biochem.9b00237](https://doi.org/10.1021/acs.biochem.9b00237) | Public author paper, supplement and repository; separate software license unverified. |
| Predicting plant Rubisco kinetics from RbcL sequence data using machine learning | Online 2022-09-12; issue 2023-01-11 | [10.1093/jxb/erac368](https://doi.org/10.1093/jxb/erac368) | Public paper and named GitHub notebooks/data; primary record inspected; no code license assumption. |
| De novo design of cysteine proteases | 2026-10-14 | [10.1038/s41586-026-11143-w](https://doi.org/10.1038/s41586-026-11143-w) | Primary author-paper PDF read; public design/analysis repositories; large models not downloaded. |
| Deep learning guided design of protease substrates | 2026-01-06 | [10.1038/s41467-025-67226-1](https://doi.org/10.1038/s41467-025-67226-1) | Public article, source data and code; 2025 in DOI is not publication year. |
| Multi-species atlas resolves an axolotl limb development and regeneration paradox | 2023-10-10 | [10.1038/s41467-023-41944-w](https://doi.org/10.1038/s41467-023-41944-w) | Public source datasets, GEO accessions and paper; atlas not downloaded. |
| Molecular basis of positional memory in limb regeneration | 2025-05-21 | [10.1038/s41586-025-09036-5](https://doi.org/10.1038/s41586-025-09036-5) | Primary article methods/results inspected; no new transcriptomic reanalysis. |
| Unravelling cytosolic delivery of cell penetrating peptides with a quantitative endosomal escape assay | 2021-06-17 | [10.1038/s41467-021-23997-x](https://doi.org/10.1038/s41467-021-23997-x) | Public primary SLEEQ paper; not a general binder-engagement dataset. |
| Comparing cell penetration of biotherapeutics across human cell lines | 2024-06-05 | [10.1021/acschembio.4c00211](https://doi.org/10.1021/acschembio.4c00211) | Public author manuscript; distinct assays/cargos cannot be pooled without harmonization. |
| Trinucleotide substrates under pH–freeze–thaw cycles enable open-ended exponential RNA replication by a polymerase ribozyme | 2025-05-28 | [10.1038/s41557-025-01830-y](https://doi.org/10.1038/s41557-025-01830-y) | Public primary author manuscript; no claim of primordial cells. |
| A small polymerase ribozyme that can synthesize itself and its complementary strand | Online 2026-02-12; issue 2026-03-05 | [10.1126/science.adt2760](https://doi.org/10.1126/science.adt2760) | Public [author manuscript](https://discovery.ucl.ac.uk/id/eprint/10222270/1/EMS212593.pdf); QT45, not an unassisted cell. |
| Synthesis and Concomitant Assembly of Adeno-Associated Virus-like Particles in Escherichia coli | Online 2022-10-24; issue 2022-11-18 | [10.1021/acssynbio.2c00451](https://doi.org/10.1021/acssynbio.2c00451) | Primary abstract/metadata verified; full experimental dataset not obtained. |
| Adeno-associated virus capsid protein expression in Escherichia coli and chemically defined capsid assembly | 2019-12-09 | [10.1038/s41598-019-54928-y](https://doi.org/10.1038/s41598-019-54928-y) | Public primary full text; assembly evidence is distinct from gene delivery. |
| Molecular basis of N₂ fixation in a hyperthermophilic archaeon | 2026-09-08 | [10.1038/s41467-026-77173-0](https://doi.org/10.1038/s41467-026-77173-0) | Open primary paper; measured versus extrapolated activities distinguished. |
| Nitrogenase resurrection and the evolution of a singular enzymatic mechanism | 2023-02-17 | [10.7554/eLife.85003](https://doi.org/10.7554/eLife.85003) | CC BY primary article; ancestral homologues remain homologues. |
| Protein-templated synthesis of dinucleotide repeat DNA by an antiphage reverse transcriptase | Online 2026-04-16; issue 2026-06-18 | [10.1126/science.aed1656](https://doi.org/10.1126/science.aed1656) | Primary record/abstract and paper text inspected; no arbitrary-peptide reverse-translation dataset. |
| Light-Driven Competitive Selection in a Protein-Catalyzed Dissipative Peptide Replication | 2026; issue 2026-04-27 | [10.1002/anie.202518911](https://doi.org/10.1002/anie.202518911) | Public author paper, CC BY 4.0; foldamer/disulfide reaction, not arbitrary polypeptide copying. |
