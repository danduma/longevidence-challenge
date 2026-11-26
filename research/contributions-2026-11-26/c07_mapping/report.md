# C07 — Source-verified Rubisco linkage repairs

**Executed result:** a version1.0.0 mapping patch replaces one additional, demonstrably wrong-strain native sequence link and explicitly quarantines two additional unresolved nonfragment links. It also operationalizes quarantine of seven fragments already recognized in the prior audit. Those seven fragment observations and the established Paspalum-to-Orthomnion error are not counted as new findings. No original embeddings or fitted model were reconstructed or refitted. Status: **agent-produced; awaiting coordinator verification**. No human expert review.

## Additional repair and explicit quarantines

| Released assay label | Released accession | C07 operation | Evidence and boundary |
|---|---|---|---|
| Rhodopseudomonas palustris | Q07N61, BisA53, 460 aa | Replace native core with Q6N0W9, CGA009, 461 aa | Original assay explicitly uses CGA009 cbbM; its WT kinetic tuple matches the compendium record. Authors' 4LF1 deposition identifies CGA009 and its canonical sequence equals the full Q6N0W9 native core plus the documented 20-aa tag. The deposited BisA53 sequence differs. |
| Avena sativa | A0A510B6V9, Avena nuda, 477 aa | Quarantine pending cultivar/sequence bridge | Original crop assay names Avena sativa ‘Forridena’. Last UniSave entry v17 names Avena nuda, taxid4497. No evidence bridges this accession to Forridena. This is an unresolved preparation link, not a categorical claim that all oat taxonomic treatments distinguish the species or that their RbcL proteins must differ. |
| Glycine max cv. Wayne | P27066, Glycine max, 475 aa | Quarantine strict cultivar link | The source labels four records Wayne. The accession's available sequence references do not identify Wayne; one identifies PI437654. Primary Laing text identity could not be recovered beyond the scanned-paper landing page/abstract. This is missing provenance, not a proven wrong protein or cultivar effect. |

The R. palustris replacement is supported jointly by [the original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4118107/?pdf=1), [4LF1](https://www.rcsb.org/structure/4LF1), and [the accession](https://www.uniprot.org/uniprotkb/Q6N0W9/entry). It corrects the native encoded core for the single matched WT source record. The paper includes expression/construct contexts; the replacement does not assert identical tags across kinetic and crystallographic preparations. A future representation must be regenerated from the replacement sequence; simply relabeling the original vector would be invalid.

The oat identity follows the original crop study's Materials and Methods and Table I, cross-checked against the archived accession. The [crop paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4972260/?pdf=1) identifies Forridena and its published kinetic values match both eligible oat records in the compendium. A replacement oat accession is not guessed.

The patch is `results/mapping_patch_v1.json` and `.csv`, with exact old/new accessions, actions, rationales, local evidence paths, source identifiers and record counts. `replacement_native_core_v1.fasta` contains only the supported replacement. `mapping_with_status_v1.csv` preserves the full mapping with audit status; quarantined rows remain visible. Consumers must filter status rather than interpret an unresolved retained link as certified.

## Selection, taxonomy and unresolved cases

`protocol.md` was written before accession retrieval/outcome scoring. Suspects were the union of sequence length below400 aa, literal species/taxon-name difference after fixed whitespace/case/minus normalization, and explicit cultivar/strain/numeric qualifiers. No kinetic value selected a suspect or threshold. There were60 selected unique species links covering58 accessions; one is the excluded prior Paspalum finding. All58 selected current accession endpoints were retrieved. Thirty-four returned Inactive/deleted status. Last-entry histories were recovered for all34, and every selected released sequence matched its active or historical accession sequence exactly. Deletion alone did not cause quarantine: it reflects current reference-proteome policy, not biological invalidity.

The complete60-row ledger is `results/link_audit.csv`:46 remain unresolved, three are retained with explicit cultivar-reference support, nine are quarantined, one is repaired, and one is the prior excluded error. These statuses concern mapping evidence, not biological independent replicates. Not every original assay full text was available; source identities for unresolved rows are retained from the compendium with citations, explicitly distinguished from independently retrieved primary evidence.

The negative controls matter. O03042 contains Columbia nucleotide-sequence provenance; P00876 contains a Wisconsin-38 structural reference. Thus the Columbia label and both W38 spellings are retained with reference support instead of being called wrong solely because the organism field lacks a cultivar. This is accession-reference support, not complete verification of every historical assay preparation.

Flaveria australasica→trinervia remains unresolved. The accession names trinervia, whereas the primary Kubien abstract treats australasica and trinervia as distinct assay labels. Taxonomic synonym reports cannot establish identical assayed specimens or sequences. No new taxonomic error or replacement is claimed for this link. Other punctuation, historical genus, subspecies and strain-name differences are retained unresolved rather than transformed automatically into errors.

The seven fragment quarantines are Atriplex glabriuscula, Capsicum annuum, Medicago sativa, Mercurialis annua, Pistacia lentiscus, Polypodium aureum and Urtica membranacea. Their lengths were already known. This patch adds auditable status and accession evidence; it does not rebrand their discovery. The <400-aa rule concerns full-sequence eligibility, not enzyme inactivity or taxonomy.

## Eligibility impact and apparent train/test exposure

The comparison excludes the prior Paspalum error from **both** arms. The baseline has206 unique linked labels and366 eligible source records under `primary==1`, `mutant==0` and fixed whitespace-normalized matching. The additional C07 operations retain204 labels and360 source records: six records are quarantined, and the repaired R. palustris link supports one retained source record. Four Wayne observations represent multiple source records, not four independent enzymes. The full patch including prior fragment quarantine retains197 labels and348 source records; the extra12 removals belong to the prior fragment stratum.

Exposure uses the fixed released prediction order and KFold(5, shuffle=True, random_state=0), keeping original folds after removal. Three conventions were fixed: same normalized species, exact full sequence, and connected components linking exact sequences plus the first two species-name tokens. The component and exact-sequence exposure counts coincide in these data. These are conditional apparent exposures, not original recovered fold manifests. Published predictions/targets are not recalculated.

| Target | Baseline prediction records | Additional C07 only | Full patch including prior fragments | Exact-sequence exposed baseline → additional → full |
|---|---:|---:|---:|---|
| S |182|181|175|57 →57 →57|
| KC |196|194|189|58 →56 →56|
| vC |181|180|176|53 →53 →53|
| KO |168|167|165|53 →51 →51|
| vO |38|38|37|8 →8 →8|
| KRuBP |34|34|33|20 →20 →20|

The additional operations lower absolute exact-sequence exposure by two records for KC and KO, with no reduction for the other targets. Removal can raise the exposed fraction: specificity changes from57/182 (31.32%) to57/181 (31.49%) under the additional operations, and57/175 (32.57%) under the full patch. No general exposure improvement or model-performance benefit follows. Same-normalized-species exposed counts remain12,10,14,8,0,10 respectively. `exposure_impact.csv` gives every numerator/denominator, `exposure_records.csv` the row-level decisions, and `fixed_fold_*.csv` the fixed folds. Counts across targets overlap and must not be summed as independent observations.

## Verification, provenance and limitations

Run `./verify.sh` from this directory. It uses the requested shared Python, no network, pinned local snapshots and no model calls. It checks hashes, reruns the audit, checks fresh result hashes, and verifies shared-input hashes when those files are present. Assertions verify released FASTA order/lengths, all selected active/historical sequences, primary-assay identity text and the WT tuple, native-core equality to the primary structure deposition, patch counts and eligibility arithmetic. `verification.log` records fresh execution; `provenance_hashes.json` pins evidence, protocol, code and outputs. Exact versions, retrieval logs, failure classifications and public download instructions are in `sources.md` and `refusals.md`.

This is a bounded source-link curation contribution. The strongest new positive finding is the source-supported BisA53→CGA009 native-core repair. The two nonfragment quarantines expose unresolved provenance without asserting different enzyme performance. Forty-six selected links and all unselected links remain incompletely validated. Small-subunit identity, assay temperature/context, synonyms, cultivar equivalence and historical preparations require further curation. The operational fragment threshold is conservative and not experimentally validated. No statistical confidence interval, improved enzyme, original-model score change, laboratory finding, or universal novelty claim is made. Local third-party snapshots are for verification; redistribution permissions are identified or explicitly withheld in sources.md.
