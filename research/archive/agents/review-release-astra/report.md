# Cross-species Rubisco kinetics release audit · 26 November 2026

**Decision.** The release supports a retrospective specificity ranking signal. It does not establish performance on unseen RbcL sequences or the ability to select an enzyme meeting both Millennium thresholds. The strongest evidence is a combination of directly traceable target duplication, a confirmed wrong-organism sequence association, and poor absolute turnover predictions. The separate mutant active-learning campaigns were outside this audit. [Young, paper and Table 3](https://www.mdpi.com/3042-6723/1/2/7)

## What the release actually contains

I compared the 12 supplied kinetic package files at GitHub ade87f3 with the corresponding members of the paper-linked Version3 archive (archive prefix 92ce749). All 12 match byte-for-byte. Version3 also supplies four nested-CV prediction files, their summary and chosen parameters. Its alternate UniProt FASTA is empty and its alternate metadata file is one byte. Neither supplied source contains the Flamholz t12 mean embedding array named in the README. The t33 embedding discussed elsewhere in that README is a different input. File digests are in manifest.csv. [Pinned repository](https://github.com/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/tree/ade87f378aa25c24a24c19915722aa04be92936e) · [Version3 archive](https://doi.org/10.5281/zenodo.19546367)

The source has 634 rows, including 504 primary wild-type measurements under 315 literal species strings. Fourteen identifiers repeat beyond their first occurrence, but no full WT record repeats after removing the exported row index. Repeated identifiers alone are insufficient grounds for deleting measurements. The FASTA and metadata agree on record order, labels and lengths. Their 214 rows comprise 207 labels, 207 literal names, 200 accessions, 200 taxon IDs and **169 exact amino-acid sequences**. Seven metadata rows are exact repeats beyond the first occurrence. Six sequence records contain nonstandard amino-acid letters. Lengths range from 16 to 490, including a 16-residue Capsicum annuum record. Identity here means equality of the entire deposited string, including ambiguous letters. Closely related, nonidentical sequences can still span folds.

**One organism mismatch can be checked against the original experiment.** The source row dilatatum_carmo_silva has S 88, KC 19.9 and KO 415. The experimental paper identifies the grass *Paspalum dilatatum* and gives those three values. The release assigns the abbreviated name “P. dilatatum” to *Orthomnion dilatatum*, accession A0A2Z6FKH0. This establishes a wrong metadata pairing without having to infer identity from predictions. The missing embedding matrix prevents checking which vector was actually used. Other genus differences can represent synonyms and were not counted as errors. [Original experimental Table 3](https://academic.oup.com/view-large/83477399)

## From measurements to deposited targets

Every deposited y_true equals log10 of the **median of the primary WT values grouped by the literal species spelling printed in that prediction file**, for both analyses. The arithmetic survives the audit. The upstream join does not. Both shipped scripts aggregate before joining and map each embedding label to a single metadata species using a Python dictionary. The supplied metadata has no trailing-space species, while 17 WT source rows have trailing-space names. With the supplied FASTA order and metadata, a literal implementation cannot emit those names. Below, “literal n” uses the supplied files, positive median filter and one embedding ID per FASTA record.

|Trait|Deposited n|Literal n|Distinct deposited names|Distinct exact sequences inferred|Test rows with same literal name in train|Test rows with identical sequence in train|
|---|---:|---:|---:|---:|---:|---:|
|S|183|177|176|142|12|57|
|vC|181|173|174|141|14|53|
|KC|197|191|190|155|10|58|
|KO|169|161|164|131|8|53|
|vO|38|42|38|32|0|8|
|KRuBP|34|24|28|19|10|20|

A specific reconstruction explains all six ordered output name lists. For each stripped metadata species, substitute the lexicographically last corresponding **raw** WT spelling, then run the original literal-name median join. This is consistent with a whitespace-preserving precursor and a last-value-wins lookup. The historical precursor is unavailable. The primary evidence is the mismatch between the supplied inputs and outputs, not a claim about the authors' runtime environment.

For example, the clean-spelling maize rows have median S 78.6, and the single trailing-space row has S 74.9 and vC 4.05. The deposited specificity file uses 74.9 twice under the trailing-space name. The Carella cultivar also supplies observations. Pooling all these raw maize records gives S 78.6 and vC 4.075. In Flaveria angustifolia the two spellings have S medians 83.2 and 90.4, and the deposit repeats 90.4. Nicotiana tabacum and two W38 spellings all map to P00876. Across the specificity file, 19 of 183 row targets change under the baseline's species pooling, and 17 of 181 change for vC. This is a change of estimand when distinct cultivars are pooled. It is a data-handling defect when trailing whitespace alone determines which measurements survive.

## Scores and fold reconstruction

All ten deposited prediction files reproduce their summary R², Spearman rho and log10 RMSE to floating-point precision. The paper's Table 3 rank scores for the four traits correspond to Version3's nested analysis. The fixed model uses five shuffled row folds, fold-local scaling and PCA (16 components), then ridge (alpha 10). Nested tuning selects PCA dimension and ridge penalty by inner pooled rank correlation. The archived summary records ordinary KFold and WT median, with no native-only or temperature/pH restrictions. The optional taxonomy split in the code was not selected. Train-only preprocessing does not ensure biological independence.

|Trait|Fixed rho|Fixed R²|Fixed RMSE|Nested rho|Nested R²|Nested RMSE|
|---|---:|---:|---:|---:|---:|---:|
|S|.642|.296|.144|.753|.381|.135|
|vC|.339|−1.157|.273|.375|−2.211|.333|
|KC|.598|.414|.242|.605|.424|.240|
|KO|.197|−.197|.340|.133|−.274|.350|
|vO|.267|−1.016|.388| | | |
|KRuBP|−.014|−.863|.253| | | |

Overlap counts above assume the saved prediction row order was the modeling order, five shuffled sklearn folds with default script seed 0, and mapping each trimmed deposited name to its metadata sequence. All rows map, and duplicates for the same name agree in sequence. Each test row is counted once if any training row matches. Under these assumptions 31.1% of S rows and 29.3% of vC rows have an exact sequence in training. Even different species labels can share RbcL. Without saved fold indices or embeddings these are reconstructed counts. Dependence can affect weighting and generalization, but the data do not identify the direction or magnitude of change in original-model performance under grouped splits. Negative turnover R² means larger squared error than the overall evaluation-set mean, even though its rank correlation is positive.

## Independent sequence baseline

The design was saved in baseline-plan.md before scoring. Whitespace, explicit cultivar, variety, authority and selected strain suffixes were normalized using the inspectable function in audit.py. One record represents each normalized species, with a median computed on original-scale WT measurements before logging. All 202 normalized metadata names have one distinct sequence each. This operational grouping does not constitute a comprehensive taxonomic revision. Remaining names sharing a taxon ID also share the same sequence here and are held together. Five folds of shuffled unique sequence hashes, repeated over 20 seeds, protect both exact-sequence and normalized-species boundaries. The model uses length and 20 amino-acid fractions, train-only scaling and fixed ridge alpha 10. Fold-specific training-target means provide a matched control.

|Trait|Species / sequences|Median ridge rho|Median ridge R² [min, max]|Median ridge / mean RMSE|
|---|---:|---:|---:|---:|
|S|174 / 142|.456|.071 [−.286, .353]|.169 / .177|
|vC|173 / 142|.274|−.189 [−.719, .056]|.204 / .188|
|KC|188 / 156|.402|.165 [−.208, .311]|.278 / .306|
|KO|163 / 133|.260|−.294 [−1.487, .094]|.251 / .222|
|vO|41 / 35|.482|−.201 [−1.271, .082]|.289 / .268|
|KRuBP|28 / 19|−.100|−1.512 [−2.285, −.274]|.321 / .210|

Ranges measure split sensitivity on one dataset. Each repeat pools all out-of-fold predictions before scoring. The fold mean is constant within a fold, so pooled mean-control ranks carry no within-fold ranking power. In a **post hoc sensitivity analysis**, excluding the verified organism mismatch and retaining lengths 400–550 improved median vC R² to .078 (range −.081 to .277, ridge/mean RMSE .173/.182). S R² was .102 (−.117 to .363). This illustrates sensitivity to sequence curation. It neither validates each remaining accession nor measures ESM performance. Near homologues, phylogenetic and assay differences, small groups, cultivar variation and unmodeled small subunits remain relevant.

## Decision against the Millennium criterion

One enzyme must jointly equal or exceed *Galdieria partita* specificity and maize turnover in paired assays with the reference enzymes. The source has a G. partita Tokara row (S 238 and vC 1.60), but that organism is absent from the metadata and prediction benchmark. The maize turnover label depends on spelling and cultivar handling. Those historical values are not substitutes for the assay controls. Marginal rank correlation across diverse WT enzymes does not measure joint threshold precision, calibration in the extreme tail, or extrapolation to engineered variants. The modest specificity signal may inform exploratory ranking, but neither the deposited turnover predictions nor this baseline establish a reliable candidate-selection benchmark for the challenge. [Millennium acceptance criterion](https://millenniumproblems.bio/)

## Search for prior reporting and reproducibility

Searched the paper and supplied code, the public repository issue/PR listing in all states (four records, no issue comments returned), and web-indexed literature using the title/DOI with Flaveria, duplication, whitespace, leakage and Orthomnion. The paper acknowledges weaker kinetic traits and assay context. Earlier Rubisco prediction research explicitly tested held-out genera, and later work examined dependence on sequence identity. Neither establishes that the particular metadata mismatch, whitespace selection or reconstructed overlap counts here had previously been documented. I found no such report in the inspected material. This search cannot establish universal novelty. [Issue/PR listing](https://api.github.com/repos/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/issues?state=all) · [Prior genus-holdout study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9833099/) · [Identity-dependent rate prediction](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/) · [Original kinetic compilation](https://doi.org/10.1021/acs.biochem.9b00237)

From this directory run `../../.venv/bin/python audit.py` and `../../.venv/bin/python audit.py --quality`. The former produces the predeclared baseline and the latter adds the exploratory quality sensitivity. identities.csv, source_rows.csv, literal_join_*.csv and provenance.csv expose the record-level joins. deposited.csv and overlap.csv contain the original-output checks. baseline*.csv contain replicate metrics and held-out predictions with fold membership. issues.json and issue_comments.json preserve the API results. refusals.md and network.md record refusals and retrieval failures separately.
