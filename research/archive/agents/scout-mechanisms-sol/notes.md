# Research log

Date: 2026-11-26. Scope: current job directory only; no subagents, external writes, or experiments.

## Mechanism-first starting hypotheses
- Cryopreservation: delivery heterogeneity and locally competing warming/stress constraints, beyond bulk mean warming rate.
- Quadruplet cell: distinguish balanced informational dependence on four bases from triplet code plus padding; quantify error detection, frame synchronization, and triplet shadow limitations.
- Reverse polymerases: activated growing 5-prime end; distinguish initial activation from subsequent extension, mismatch selectivity and substrate recognition. Thg1 is a starting mechanism, not a general polymerase.

Exact challenge definitions are preserved in sources/millenniumproblems_2026-11-26.html (raw downloaded primary site); report will summarize.

## Access issues (not refusals)
- Initial Python urllib request for canonical website: HTTP 403. Standard curl subsequently succeeded, without authentication or circumvention.
- Web open returned browser/recaptcha checks for PMC11792677, PMC6103667, PMC8127993. Indexed primary text and journal/Europe PMC copies used where available. PMC11957112 text opened normally.

## Early prior art
- Han et al., Nature Communications 2023, 10.1038/s41467-023-38824-8: kidney recovery in rats; not whole adult mice.
- Pigula et al., ACS Synthetic Biology 2024, 10.1021/acssynbio.4c00630: mitochondrial partial quadruplet decoding; not quadruplet-only cell.
- Patel et al., RNA 2021, 10.1261/rna.078686.121: 16 incorporation efficiencies available in Table 1; SceThg1 is biologically selected for G addition, not arbitrary templated synthesis.
- Iwaniec et al., 2025 preprint, 10.1101/2025.03.21.644660: activation and elongation differ in metal requirements; initial activation is not required every extension.
- Vitrify repository created 2026-08-09, latest push 2026-08-10, Apache-2.0: inverse nanowarming design is already public prior art; inspect implementation before proposing novelty.

## Additional access and source-identity checks
- Europe PMC REST fullTextXML succeeded for the 2021 Thg1 paper, 2024 short-duplex study, 2023 rat kidney study, and 2025 liter-scale cryopreservation study. Initial requests for the 2018 and 2025-metal papers timed out; indexed primary texts remain readable.
- Corrected a guessed PMCID before use: PMC10256762 is an unrelated blood-pressure paper, excluded. Correct rat kidney record verified by DOI/title is PMC10256770.
- Nature web pages sometimes redirect to inaccessible idp.nature.com. Official static-content.springer.com Source Data workbook downloaded successfully (845 KiB). Initial springer-cdn hostname failed DNS; this is not a refusal.
- RNA journal PDF/image requests returned HTTP 403; PMC table image is represented as a graphic in downloaded XML. Table 1 numbers verified against indexed author-paper PDF text, with all 16 values and the censored G:C entry retained.
- GitHub khatvangi/triplet-proof returned 404 through connector; its public Zenodo v2 archive record was readable. Its corrected v2 statement is treated as repository/manuscript prior art, not independently reproduced or peer-reviewed evidence.
- RNA 2026 Johnecheck full text is subscription-gated. Primary journal metadata and abstract/preprint discovery inspected; not represented as a fully read paper.

## Computation completed
- analysis/thg1_competition.py: Perron-Frobenius minimax analysis of published 4-by-4 single-addition efficiencies. A numerical stopping issue for extremely small GTP fractions was fixed by using componentwise Collatz-Wielandt residuals instead of an absolute vector-change cutoff. Fresh run succeeds; no biological inference beyond explicit model.
- analysis/quadruplet_audit.py: exhaustively enumerates mutations and shifted windows of synthetic maps. Sparse parity design shows equal mutation sensitivity but zero conditional amino-acid entropy; dense rotation control has balanced dependence and high shifted-window acceptance.
- analysis/cryo_source_audit.py: extracts exact Fig. 7 source cells, clusters probes by experimental replicate, and reports individual temperature differences rather than asserting all are below 5 C. Two-region duty-schedule bound is an idealized feasibility screen.

## Final source/date corrections and verification scope
- Europe PMC initial 2018/2025-metal requests encountered timeout, followed by HTTP 500; mitochondrial full-text XML also returned HTTP 500. Indexed primary manuscripts/journal records were used; no raw data claimed for those items.
- The unrelated guessed-PMCID download was removed after checking its title; it never entered any computation.
- Costello paper is online 2024-09-11 / July 2025 issue, with publisher rights reserved; not assumed CC BY.
- Iqbal exact title is “Predicting plant Rubisco kinetics from RbcL sequence data using machine learning”, DOI 10.1093/jxb/erac368; online 2022-09-12 / issue 2023-01-11.
- Verification-before-completion skill read and applied to fresh analysis runs, explicit mathematical checks, full source IDs/licenses, all twelve canonical titles, and local artifact links.
