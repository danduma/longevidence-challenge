# Millennium biology: public-data prioritization and three executed analyses

**Research date: 2026-11-26. Independent computational worker; this job directory only.**

The strongest public-data opportunities are **04 Improve Rubisco**, **08 Programmable proteases**, and **06 Somatic limb regeneration**. I directly pursued all three, downloaded primary data, inspected associated code, and wrote runnable analyses. The useful outputs are an audited Rubisco evidence table and frontier sensitivity analysis, a protease extrapolation/assay-portability benchmark, and an animal-aware regeneration contrast. **No experimental challenge was solved, no new enzyme was verified, and no sequence was generated or represented as a discovery.**

The clearest next project is a **provenance-corrected, assay-aware Rubisco frontier update incorporating the new 2026 paired kinetic measurements**, followed by phylogenetic blocking and explicit uncertainty. The old speed–specificity narrative is insufficient: my results depend strongly on taxonomic composition; the five-strain 2026 primary cohort has a positive, inconclusive within-cohort association; and newer fast-turnover libraries lack specificity. Importantly, I verified misdirected citation links at influential historical frontier points. These findings make data curation a necessary scientific step, not clerical work.

## Deliverables and evidence standard

- [Source inventory with 40 entries, exact titles, DOIs, dates, licenses/access and comparability](source_inventory.md); [machine-readable source inventory](source_inventory.csv).
- Downloaded-data sizes, checksums and provenance (archival input; see source-downloads.json and the source URLs); download request log (archival input; see source-downloads.json and the source URLs); access errors (archival input; see source-downloads.json and the source URLs).
- [Analysis protocol and scientific nulls](notes/analysis_protocol.md); [durable work log](notes/worklog.md).
- [Rankings and sensitivity](results/challenge_ranking.csv); [ranking method](results/ranking_method.json).
- [Rubisco results](results/rubisco_summary.json), [protease results](results/protease_summary.json), [regeneration results](results/regeneration_summary.json).

The scan used live searches, primary full text or author manuscripts, repository trees/code, source spreadsheets, and archive metadata. Some lower-priority sources were accessible only through abstracts, indexed primary text or supplementary-data statements; those limits are identified in the inventory. “Public,” “downloaded,” “licensed for reuse,” “measured,” and “validated against this challenge” are separate attributes. Not every underlying numerical transcription in a historical compilation has been rechecked against its original paper. Every paper actually cited here has an identifiable primary record; resolving a DOI alone was insufficient when its title disagreed with the citation.

## The twelve definitions used

The names below exactly follow the [challenge site](https://millenniumproblems.bio/), retrieved today and saved in the original HTML snapshot (archival input; see source-downloads.json and the source URLs). The endpoint descriptions are compressed summaries, not substitutes for the full definitions. In particular, I did not reinterpret an adjacent assay as the challenge itself.

| ID | Exact challenge name | Required endpoint and key exclusion |
|---|---|---|
| 01 | Origins of life | Unassisted emergence of heritable, bounded RNA/protein cells from plausible precursors and energy; ≥10⁶ population increase, with plausibly indefinite division rather than progressive shrinkage. |
| 02 | Cryopreservation | Intact adult wild-type mice recovered after ≥24 h whole-body freezing/vitrification, >99% viability, no permanent damage; ethics approval required. |
| 03 | The reverse translatase | A purified protein catalyst/fixed complex reads untagged peptide into amplifiable nucleic acid; no template, barcode, residue-specific operator cycle or lookup. Preregistered pool: ≥100 random peptides ≥50 residues; ≥90% accuracy, mean read ≥25 residues, mean Q≥10. |
| 04 | Improve Rubisco | Specificity at least that of G. partita and turnover at least that of maize, demonstrated with the specified paired enzyme controls. Natural discovery, engineering and de novo design are eligible. |
| 05 | The quadruplet cell | A replicating cell with every coding sequence, including translation machinery, encoded in genuine nonoverlapping quadruplets; no detectable triplet decoding and comparable nonsynonymous-mutation probability across all four positions. |
| 06 | Somatic limb regeneration | Adult wild-type mouse limbs regrown reproducibly, matching controls in motor and sensory function and blinded noninvasive appearance; ethics approval required. |
| 07 | Bacterial production of gene therapies | Bacteria produce infectious replication-incompetent AAV and lentivirus with specified, nuclease-protected genomes and mammalian-comparable particle/genome and infectious-unit/genome ratios; AAV alone is partial success. |
| 08 | Programmable proteases | Prospectively designed enzymes cleave blinded accessible sites in endogenous folded proteins in cells; reference-level efficiency/off-target performance; >80% success over 20 preregistered sites, designs within 24 h, no target-specific experimental optimization. Peptide-only cleavage is partial success. |
| 09 | Cell-penetrating protein binders | Zero-shot extracellularly dosed binders engage intracellular targets at supported concentrations: 80% success over 20 preregistered targets, design within 24 h, no subsequent evolution. Transfection, intrabodies, electroporation and membrane disruption do not qualify; antibody-like scaffolds are preferred. |
| 10 | Protein amplification chain reaction | Single-pot input-dependent covalent copying from amino-acid monomers: ≥100 random ≥50-residue peptides pooled, ≥1000× amplification and ≥90% residue accuracy. No nucleic-acid template/intermediate, cognate scaffold or explicit sequencing. |
| 11 | 5′ polymerases | Complete 3′→5′ DNA, RNA, reverse-transcriptase and RdRP toolkit, with processivity/error performance comparable to Taq, T7 RNA polymerase, M-MLV RT and Phi6 RdRP respectively. |
| 12 | New nitrogenases | N₂→ammonia at roughly natural-enzyme rates, while below sequence/structure similarity thresholds against all known nitrogenase and nitrogenase-like families. |

## Transparent ranking

Scores are subjective, evidence-informed priorities for **meaningful bounded computational progress**, not probabilities of solving the entire challenge. Each dimension runs from 1 (poor) to 5 (strong), with half-points. Data scores reflect relevant labels, provenance, usable downloads and licenses, not simply sequence volume. Tractability rewards falsifiable computations that can alter a decision. Validation scores reward a short path from computational output to the actual endpoint.

Primary weights: **data 35%, tractability 30%, impact 20%, validation path 15%**. Ties are broken by validation path, then tractability. Impact is prospective scientific/social value, not a quantitative forecast.

| Rank | Challenge | Data | Tractability | Impact | Validation | Weighted score |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 04 Improve Rubisco | 5 | 4.5 | 5 | 4.5 | 4.775 |
| 2 | 08 Programmable proteases | 5 | 4 | 5 | 4 | 4.550 |
| 3 | 06 Somatic limb regeneration | 5 | 3.5 | 5 | 2 | 4.100 |
| 4 | 05 The quadruplet cell | 4 | 3.5 | 4 | 3.5 | 3.775 |
| 5 | 09 Cell-penetrating protein binders | 3.5 | 3 | 5 | 3 | 3.575 |
| 6 | 02 Cryopreservation | 3.5 | 3.5 | 5 | 2 | 3.575 |
| 7 | 12 New nitrogenases | 4 | 2.5 | 5 | 2 | 3.450 |
| 8 | 01 Origins of life | 3.5 | 2.5 | 5 | 1.5 | 3.200 |
| 9 | 07 Bacterial production of gene therapies | 3 | 2.5 | 4.5 | 2.5 | 3.075 |
| 10 | 11 5′ polymerases | 2.5 | 3 | 4 | 3 | 3.025 |
| 11 | 03 The reverse translatase | 3 | 1.5 | 4.5 | 2 | 2.700 |
| 12 | 10 Protein amplification chain reaction | 2 | 1 | 4.5 | 1.5 | 2.125 |

Rubisco and proteases remain first and second with equal weights, stronger data weighting, and stronger validation weighting. **The third choice is less robust:** regeneration remains third under equal/data-first weights, but falls to fifth when validation receives 45%; quadruplet decoding then becomes third. Thus regeneration is recommended for computational discovery and reproducibility work, while quadruplet decoding is a reasonable replacement if near-term experimental validation dominates the project brief.

## Independent investigation of all twelve

**01 — Origins of life.** There are real RNA activity measurements and generative-model code rather than merely proposed mechanisms. The 2025 paper *Exploring the space of self-reproducing ribozymes using generative models* supplies a [CC BY archive](https://zenodo.org/records/16531362) of about 108 MB; its correction history must be retained. The 2025 [trinucleotide replication study](https://doi.org/10.1038/s41557-025-01830-y) and 2026 [45-nt QT45 primary paper](https://doi.org/10.1126/science.adt2760) are substantial advances in RNA copying. A worthwhile benchmark is whether sequence–activity models generalize across ribozyme families and environmental cycles, rather than recognizing close relatives. Supplied activated substrates, selected RNA, environmental forcing and compartment conditions are major domain differences. Neither an isolated ribozyme activity predictor nor RNA amplification establishes emergent RNA/protein cellular life.

**02 — Cryopreservation.** Downloaded the 2023 rat-kidney source workbook (archival input; see source-downloads.json and the source URLs): 17 sheets include functional outcomes and eGFR. The [primary paper](https://doi.org/10.1038/s41467-023-38824-8) reports life-sustaining transplantation after storage up to 100 days. A [2026 pig-kidney preprint](https://doi.org/10.64898/2026.10.07.756990), posted October 14, extends scale to 84–132 g and reports acute recovery, with impaired filtration relative to fresh controls. It is not a long-term intact-animal result. Thermal/transport uncertainty analysis is executable, but this workbook does not provide complete spatial CPA, nanoparticle and temperature fields. The validation bottleneck is simultaneous toxicity, cooling, warming and organ-system recovery; fitting a heat equation does not establish the >99% whole-mouse criterion.

**03 — The reverse translatase.** Two 2026 primary papers are especially relevant. [*Single-molecule peptide sequencing through reverse translation of peptides into DNA*](https://doi.org/10.1038/s41587-026-03061-z), March 18, has MIT [analysis code](https://github.com/whulwzheng-source/smPeptideSeq) and SRA BioProjects PRJNA1420480/PRJNA1423337. It uses iterative Edman chemistry and barcodes, so the title must not be mistaken for fulfillment of the enzyme challenge. [*From peptides to DNA: All required steps can be catalyzed*](https://doi.org/10.1073/pnas.2534387123), February 26, tests separate steps using aminopeptidase, aptazymes and ligation. It does not demonstrate integrated sequence-order-preserving processivity across an arbitrary peptide pool. Error-channel benchmarking is possible, but there is no verified training set of successful challenge-conforming enzymes.

**04 — Improve Rubisco.** Highest priority because source-linked kinetic tables, errors, multiple assay studies, exact biochemical outputs and specified controls make a falsifiable frontier audit possible. Historical data and normalization code are [public](https://github.com/flamholz/rubisco); 2025 [corrected rates](https://doi.org/10.1038/s44318-025-00419-y), the [tree-of-life library](https://doi.org/10.1073/pnas.2501433122), and new [2026 paired cyanobacterial measurements](https://doi.org/10.1093/plphys/kiag259) materially change the evidence base. Computational progress is most credible as an uncertainty/provenance/coverage result, not as an inferred new enzyme. Details below.

**05 — The quadruplet cell.** The [2021 tRNA-evolution paper](https://doi.org/10.1038/s41467-021-25948-y) has public source sequences and reporter measurements, including up to 80-fold improvements. The [2025 QCODE paper](https://doi.org/10.1093/nar/gkae1292) supplies supplementary data for engineered quadruplet insertions and growth phenotypes. The [mitochondrial cox3 study](https://doi.org/10.1021/acssynbio.4c00630) is a 2024 publication, despite its manuscript becoming public in 2025. All are partial decoding systems, not whole quadruplet cells. A useful computational null is that reported “quadruplet” codebooks are effectively triplet codes with a fourth-position frameshift signal: enumerate nonsynonymous mutation probabilities at each position and compare with balanced four-position codebooks. Reporter contexts, tRNA charging, off-frame decoding and genome-wide triplet machinery prevent extrapolation from suppression efficiency to the challenge endpoint.

**06 — Somatic limb regeneration.** Public single-cell, spatial and bulk data are abundant, including the [2023 multi-species atlas](https://doi.org/10.1038/s41467-023-41944-w), [2025 retinoic-acid positional-identity study](https://doi.org/10.1038/s41467-025-59497-5), and [October 2026 tissue-plane study](https://doi.org/10.1038/s41467-026-78607-5). Downloaded normalized feature counts and GEO metadata for the last study. [GSE267553](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267553) provides adult mouse WT/KO digit-tip material, but only two sample records require careful biological-replication auditing; its [Mest paper](https://doi.org/10.1038/s41536-024-00376-w) is from 2024. Development versus regeneration, distal digit versus limb, species, orthology, age, injury plane and cell composition are strong confounders. Expression conservation is a prioritization signal, not evidence of causal mouse limb restoration.

**07 — Bacterial production of gene therapies.** [*Synthesis and Concomitant Assembly of Adeno-Associated Virus-like Particles in Escherichia coli*](https://doi.org/10.1021/acssynbio.2c00451), 2022, establishes VLP assembly/internalization, not infectious genome-containing therapeutic particles. A downloaded author-related URL turned out to be a 2021 dissertation and is explicitly classified that way. The [2025 HEK293 multi-omics/GSMM paper](https://doi.org/10.1016/j.ymben.2025.03.011) links code/data and reports 2.5× higher capsid production with fewer vector genomes: a useful warning against optimizing capsids alone. Exact repository access was not resolved in this run. A host-comparison mass-balance/QC analysis is executable, but mammalian expression labels cannot establish bacterial assembly, genome protection, packaging and infectivity. No challenge-conforming bacterial lentiviral dataset was verified in this bounded scan.

**08 — Programmable proteases.** Particularly strong data: [ProtRec/MLDEEP](https://doi.org/10.1038/s41467-025-60622-7) provides measured enzyme–substrate activities, controls, in-vitro validation, split code and licensed archives. I downloaded a 14.7 MB compressed DMS table. Recent progress is real: the October 14, 2026 papers [*De novo design of cysteine proteases*](https://doi.org/10.1038/s41586-026-11143-w) and [*Computational design of metalloproteases*](https://doi.org/10.1038/s41586-026-11167-2) report 13/69 active initial cysteine designs and 36% of 135 active initial metalloprotease designs, respectively. Optimized variants and target examples are not a blinded 20-site success rate. The [January 2026 CleaveNet paper](https://doi.org/10.1038/s41467-025-67226-1) designs substrates for enzymes, the inverse task. Access, folding, endogenous-cell context and proteome-wide off-target cleavage remain the decisive benchmark gaps.

**09 — Cell-penetrating protein binders.** [CAPA across cell types](https://doi.org/10.1021/acschembio.4c00211) is a 2024 paper with public manuscript release in 2025. It distinguishes productive cytosolic/nuclear access from total uptake. The April 17, 2026 [cell-permeable CFTR nanobody paper](https://doi.org/10.1038/s41589-026-02199-w) provides an endogenous functional target example, but begins with a selected binder and carrier chemistry; raw imaging is upon request. An assay-stratified delivery benchmark using cargo, concentration, time, cell type, viability and compartment readout is tractable. There is no verified broad dataset equating fluorescence uptake with zero-shot, therapeutically plausible intracellular engagement. Membrane-disruptive delivery would also fail the specified challenge, even if uptake were impressive.

**10 — Protein amplification chain reaction.** The February 2026 [light-driven dissipative peptide replication study](https://doi.org/10.1002/anie.202518911) has real kinetic evidence and supporting data for 78 replicator types assembled from 12 foldamer building blocks. Those precursors and disulfide exchange are far from arbitrary polypeptide copying from amino-acid monomers. A mechanistic-model comparison can test whether autocatalytic kinetics are distinguishable from ordinary assembly, but would not establish sequence-faithful full-length covalent copying. Seeded aggregation assays and nucleic-acid-dependent cell-free translation cannot be counted as challenge training labels. This was lowest priority because positive data matching the requested transformation are essentially absent from the sources verified here.

**11 — 5′ polymerases.** The actual [2026 Thg1/TLP comparison](https://doi.org/10.1261/rna.080815.125), May 18, has a verified [2025 preprint](https://doi.org/10.1101/2025.10.12.681950), whose full text was retrieved through its correct PMID/PMCID. Short 3′→5′ additions, tRNA recognition, activation and fidelity can be analyzed using existing structures and transient kinetics. The relevant baseline is template-dependent tRNA repair, not ordinary 5′→3′ reverse transcription. A falsifiable subproblem is whether substrate-independent nucleotide addition can be distinguished from tRNA-recognition effects. There is no verified complete long-template reverse-direction toolkit benchmark comparable to all four named reference polymerases. The journal PDF returned 403; the underlying preprint and indexed journal abstract were read.

**12 — New nitrogenases.** [*Structural evolution of nitrogenase over 3 billion years*](https://doi.org/10.7554/eLife.105613.4), September 2025, supplies 5,378 predicted complexes, experimental structures, [code](https://github.com/kacarlab/nif-structures) and an archived structural database. The [2025 NFEmbed paper](https://doi.org/10.1093/bioadv/vbaf204) has downloadable real tables: 402 processed strain samples after removing metadata headers, plus source activities and original units. These mainly represent whole-strain acetylene reduction, not purified N₂-to-ammonia turnover. I downloaded the small CSVs and inspected their schemas. A homology/structure out-of-distribution audit and species-held-out predictor are executable. They cannot establish a genuinely unrelated nitrogenase: all verified positive families are homologous, and a distant fold or predicted metal site does not demonstrate catalysis. A direct isotope-supported activity measurement would ultimately be necessary, with contamination controls.

## Pursuit 1: Rubisco source-linked kinetics and the frontier

### Inputs, baseline and null

The 2019 original workbook (archival input; see source-downloads.json and the source URLs) retains 634 source rows and the original `use`, primary-study, mutant, temperature and error fields. The source CSV drops the `use` flag, and the supplied merged table can include imputed errors and cross-pH joins. I therefore used the original workbook, retaining `use=True`, primary observations, nonmutants, Form I and 25°C. This gave **459 rows**, **130 with both endpoints on the same original row**, and **226 species–study pairs from 27 studies / 198 species** after joining within study. Within-study joins are explicitly distinguished from same-row measurements; neither proves the same preparation or the required paired controls.

Baseline: species-median log specificity versus log turnover, plus a point-estimate nondominated set. Scientific null: the relationship and frontier are stable to source blocking, pH restrictions, study deletion and plausible measurement perturbation. I also examined joint study/taxonomy adjustment; this is an exploratory confounding diagnostic, not a phylogenetic causal model.

### Quantitative results

| Quantity | Computed result |
|---|---:|
| Species-median log–log Pearson r | −0.534, n=198 |
| Source-paired log–log Pearson r | −0.556, n=226 |
| Same-original-row log–log Pearson r | −0.658, n=130 |
| Within-study log slope | −0.307 |
| Study-cluster bootstrap 95% slope interval | [−0.590, −0.122] |
| Within-taxonomy log slope | −0.050 |
| Joint study + coarse-taxonomy log slope | −0.036 |
| Leave-one-study-out within-study slope range | [−0.442, −0.225] |
| Point-estimate frontier | 14 source pairs, 13 distinct species, 8 studies |
| Pairs not certainly dominated under ±10% bounds on both endpoints | 145 / 226 |

The ±10% analysis is **a sensitivity assumption**, not measured confidence limits: a row is certainly dominated only if another row's lower bounds exceed its upper bounds. It demonstrates that sharp frontier membership can be fragile, not that 145 enzymes are truly optimal. The SD perturbation output similarly assumes independent marginal errors and is not a posterior probability. Reported SDs were never silently replaced with zero; the historical compilation’s `_SD` labels still require per-paper SD/SE verification. Broad taxonomy combines ecology and evolutionary history, so the adjusted slope cannot separate them, establish no chemical constraint, or substitute for a tree-based analysis.

Using only diagnostic historical values from the workbook, **S≥238.1 and kcat≥4.05 s⁻¹ produces zero qualifying pairs**. Among pairs above the maize turnover reference, the highest specificity is 111.5; above the Galdieria specificity reference, the highest turnover is 1.6 s⁻¹. These are a historical-data diagnostic, not the official challenge's universal numerical thresholds: actual success requires contemporaneous paired control assays. The original Galdieria paper's abstract rounds specificity to 238; its numerical table has not been fully rechecked here.

The new 2026 primary cyanobacterial table (archival input; see source-downloads.json and the source URLs) provides five strains and a wheat control at 25°C. Across the five cyanobacteria, log–log **r=+0.684, exploratory p=0.203; Spearman ρ=0.80, p=0.104**. The highest-specificity strain also has the highest turnover in that cohort (74.6 and 11.8 s⁻¹). This is a small, selected cohort with reported **SE**, not SD. It contradicts an assumption that every within-cohort relationship must be inverse; it does not demonstrate escape from the global frontier or the Galdieria requirement. The 2026 aquatic 103-species article is a review explicitly reporting no new data, so I did not count it as 103 independent new assays.

### What the newer rate library can and cannot show

The 2025 tree-of-life screened table (archival input; see source-downloads.json and the source URLs) has **469 rows, 261 numeric rate estimates, 26 explicit below-0.5 flags and no specificity column**. Unreported/failed screens must not become zero-rate training labels. Its rates are corrected coupled-assay values at 30°C, scaled by a shared factor approximately 2.1 (published 95% CI 1.7–2.8). The largest numeric row, 33.00 s⁻¹, is a reconstructed ancestor; the largest naturally sampled row is 28.53 s⁻¹. Neither should be mixed uncritically with historical 25°C direct assays. Sequence-derived rate predictions are not specificity measurements.

The [April 2025 correction](https://doi.org/10.1038/s44318-025-00419-y) documents both a missed dilution and inaccurate active-site quantification. Correcting only a factor of two is insufficient. The paper's repeated direct assays and the later library's shared calibration are different evidence layers; common calibration uncertainty is correlated across rows.

### A concrete source-audit result

I resolved the identifiers for all 27 studies used in the paired audit and compared titles with source citations. The downloaded compilation contains influential bibliographic mistakes:

| Source label | Observed mismatch | Verified correct reference |
|---|---|---|
| Uemura 1997 | Galdieria/red-algal rows link to the 1996 assay-method DOI. | [1997 red-algal kinetics paper, DOI 10.1006/bbrc.1997.6497](https://pubmed.ncbi.nlm.nih.gov/9144578/) |
| Uemura 1996 | Specificity rows carry the 1997 PMID/title. | [1996 anion-exchange method, DOI 10.1093/oxfordjournals.pcp.a028949](https://academic.oup.com/pcp/article-abstract/37/3/325/1930998) |
| Whitney 2001 | PMID 26662726 resolves to Occhialini 2015. | [*Form I Rubiscos from non-green algae are expressed abundantly but not assembled in tobacco chloroplasts*, DOI 10.1046/j.1365-313x.2001.01056.x](https://onlinelibrary.wiley.com/doi/10.1046/j.1365-313x.2001.01056.x) |
| Bainbridge 1995 | DOI suffix increments across rows; the wild-type row points to nonresolving page 1272. | [*Engineering Rubisco to change its catalytic properties*, DOI 10.1093/jxb/46.special_issue.1269](https://academic.oup.com/jxb/article-abstract/46/special_issue/1269/613044) |

Original files were preserved. [The reviewed audit](results/rubisco_citation_audit_reviewed.csv), [24-row bibliographic overlay](results/rubisco_bibliographic_overlay.csv), and [226 pairs with corrected source links](results/rubisco_source_pairs_linked.csv) distinguish identifier correction from numerical verification. Not all 24 overlay rows had distinct errors; the overlay covers all rows in the affected source groups. No kinetic measurement was altered on the strength of a citation correction.

**Prior art and novelty limit.** [Flamholz 2019](https://doi.org/10.1021/acs.biochem.9b00237) already weakened the simple tradeoff interpretation. [Bouvier 2021](https://doi.org/10.1093/molbev/msab079) explicitly analyzed phylogenetic constraint; the interpretation was contested and [a 2023 response](https://doi.org/10.1016/j.jplph.2023.154021) defended it. Thus “the tradeoff is weaker after correction” is not a novel discovery. A 2026 [Rubisco-design preprint](https://doi.org/10.64898/2026.09.04.749267) explicitly distinguishes endpoint apparent specificity from conventional kinetic specificity; those values cannot be placed directly on this frontier. A potentially new contribution is an updated provenance/assay/phylogeny/uncertainty analysis with the post-2024 data, provided it is checked against those analyses and subsequent supplements. This run establishes feasibility and concrete curation targets, not novelty of that future study.

## Pursuit 2: Protease benchmark generalization and assay portability

**Inputs.** The public ProtRec DMS file contains **413,779 measured pairs, 6,199 enzyme labels and 163 substrates**, mostly local variants of a seven-residue TEV motif. I also downloaded 3,482 control records, 840 in-vitro relative-rate rows and the 29,040-row published MLDEEP result table, and read its preprocessing/splitting code. The publication's approximately 600,000 pairs refer to a broader collection of screens, not the denominator of this downloaded DMS file.

**Baseline/null.** Can substrate identity or a simple additive model explain random-pair performance, and does that prediction transfer to unseen combinations of substrate substitutions? I compared a global mean, substrate-only mean, enzyme-only mean, and a regularized additive model with enzyme identity plus substrate-position/residue effects. The audit retained ≥100 minimum reads, standard seven-residue substrates and single-substitution enzyme labels, leaving **294,557 pairs / 4,411 enzyme labels / 161 substrates**. No joint enzyme–substrate interaction features were fitted. This is a deliberately simple baseline audit, not a rerun or claimed defeat of MLDEEP.

| Evaluation | Mean baseline RMSE | Additive RMSE | Additive Spearman | Additive average precision / positive prevalence |
|---|---:|---:|---:|---:|
| Random 20% pairs, n=59,054 | 0.02241 | 0.01671 | 0.644 | 0.205 / 0.0120 |
| Unseen multi-mutant substrates, n=10,018 across 27 motifs | 0.02150 | 0.02906 | 0.356 | 0.0188 / 0.0221 |

The blocked test has **no previously unseen nonzero sequence feature**: the model can combine substitution features observed individually. Yet its RMSE is about 35% worse than the mean predictor, and average precision falls below prevalence. A bootstrap over the 27 held-out motifs gives a 95% RMSE difference, additive minus mean, of **[+0.00171,+0.01252]**. This rejects adequacy of this additive extrapolation baseline under the chosen split. It does not prove every deep model fails, establish biological epistasis independently of assay background, or estimate prospective endogenous-site success. Positives here use an explicit alternative threshold—the 99th percentile of 467 noncatalytic controls, AUC 0.136788—not the authors' p-value hysteresis labels.

**Independent biochemical check.** I matched **802 distinct selected enzyme–substrate pairs** between reporter controls and in-vitro initial-rate measurements. Per-enzyme Spearman correlations range **0.620–0.816** before background subtraction and **0.680–0.808** after it. For the parental label, subtraction improves ρ from 0.698 to 0.739. Correlation is substantial but imperfect; background subtraction is not uniformly beneficial. These are selected validation constructs, not an unbiased screen of arbitrary proteases. Negative initial-rate estimates were retained as noise rather than presented as negative catalysis.

**Limitations/validation path.** Reporter AUC depends on expression, binding, recording and read depth and is not kcat/Km. Raw and matched-background metrics should accompany a future sequence-cluster/substrate-held-out MLDEEP rerun. The 2026 cysteine-protease repository contains analysis notebooks and placeholders but no original plate-reader XLSX inputs in its inspected tree; the metalloprotease archive is CC BY but about 2.31 GB and was not downloaded. Those availability differences matter when planning compute. Ultimate validation needs preregistered endogenous folded sites, cellular cleavage, reference efficiency and proteome-wide off-target measurements. Neither the local baseline nor the published optimized examples supplies that evidence.

## Pursuit 3: Regeneration injury-plane interaction

**Inputs/design.** Extracted the primary 2026 Supplementary Data 5 workbook into **139,698 feature rows × 18 biological sample columns**. Three intact animals contribute paired epiphysis/diaphysis regions; regenerated samples are independent animals, three per plane/time at 3 and 5 days after amputation. Two forelimbs pooled from an animal count as one sample. The paper uses young axolotls and specialized experimental lines; it is not adult wild-type mouse material.

**Baseline/null.** The common mistake would be to attribute every difference between regenerating planes to a new regenerative program. I estimated, on log₂(CPM+1),

`interaction(day) = [diaphysis(day) − epiphysis(day)] − [diaphysis(intact) − epiphysis(intact)]`.

The null is zero injury-by-plane interaction after accounting for pre-existing differences. Uncertainty uses the paired intact difference and independent injured-group variances, with a Welch/Satterthwaite approximation. I retained 20,052 rows with CPM≥1 in at least three samples and applied BH separately at each time. Normalized CPM was not forced to sum to one million: effective-library normalization can change column sums.

**Results.** At 3 days, 1,570 expressed rows have an unadjusted plane difference of at least one log₂ unit; **412 (26.2%) fall below that threshold after intact subtraction**. At 5 days, **78/235 (33.2%)** do so. **No interaction passes BH q<0.05 at either time** under this conservative analysis. Effect rankings are relatively insensitive to changing the pseudocount from 1 to 0.1 (Spearman 0.987 and 0.981). The two times share the same intact baseline; their effect correlation of 0.423 is therefore not independent replication.

The already-published Loc138491483 chemokine has estimated interactions +2.15 at 3 days and −0.19 at 5 days (BH q≈0.38 and 0.97). These whole-tissue contrasts do not contradict the paper's spatial localization, HCR or perturbation data: spatial enrichment, absolute injury induction and interaction after baseline subtraction are different estimands. Merged-feature annotation is incomplete; 3,162 expressed rows have blank symbols. I did not infer cross-species causal targets from them.

**Falsification and next validation.** The stronger version of “all large injury-plane differences are new regeneration-specific programs” is unsupported by the threshold attrition. The interaction null is not rejected after multiplicity correction; that is not evidence that all true interactions vanish. This small-n log-CPM test is not an edgeR raw-count reanalysis. The next meaningful computation would obtain unnormalized count/annotation files and independent mouse regenerating versus fibrosing digit cohorts, then test ortholog-matched, animal-level effects with cell-composition controls. Whole-limb function in adult mice remains a separate experimental challenge.

## Reproduction and verification

From this directory, with the included environment or the versions in requirements.txt (archival input; see source-downloads.json and the source URLs):

```sh
.venv/bin/python scripts/extract_cyano2026.py
.venv/bin/python scripts/analyze_rubisco.py
.venv/bin/python scripts/apply_bibliographic_overlay.py
.venv/bin/python scripts/build_assay_ledger.py
.venv/bin/python scripts/analyze_protease.py
.venv/bin/python scripts/extract_regeneration.py
.venv/bin/python scripts/analyze_regeneration.py
.venv/bin/python scripts/rank_challenges.py
.venv/bin/python scripts/build_inventory.py
```

The extraction scripts assert sample order, feature counts and table headers; analysis outputs retain source IDs. The separate [citation resolver](scripts/audit_rubisco_citations.py) uses live PubMed/Crossref and can recreate the underlying identifier-target audit. It is not required for a network-free numerical rerun because the checked audit is saved. Bootstrap/model randomness is fixed. Download checksums, content validation, direct reruns and output-schema/numerical consistency checks are recorded in [verification notes](notes/verification.md). This is verification of computations and provenance, not independent experimental replication.

**Refusals: none observed.** Network/access errors, including HTTP403/500/503, browser checks, an invalid link and Python TLS installation errors, are recorded separately. No model change, other model session, subagent, external posting, contact, deployment, or human/animal experiment was performed.

## Most promising specific next action

A [22-record initial assay ledger](results/rubisco_assay_ledger.csv) already combines the historical frontier/control rows with the six 2026 records and marks unverified fields explicitly. **Start the next Rubisco session by rechecking the original numerical tables and assay conditions for the Galdieria/Uemura1997, Griffithsia/Whitney2001 and maize/Kubien2008 anchor measurements, using the corrected links already delivered.** Then append the six extracted 2026 Table 1 records as a separate assay cohort, reconcile exact strain identities, and add gas convention, pH, active-site quantification, SD-versus-SE and shared-calibration fields. Fit a preregistered source/phylogeny-aware uncertainty frontier only after that ledger is complete. Select any missing-specificity measurements from the natural fast-rate library using coverage and uncertainty reduction, not a predicted discovery claim. A future qualifying candidate would still need the specified paired controls in actual enzyme assays. This is a concrete data-and-analysis task with immediate inputs, identifiable errors and a falsifiable endpoint; it is the best next use of effort from this scan.
