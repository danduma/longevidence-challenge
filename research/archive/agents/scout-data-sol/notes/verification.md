# Verification —2026-11-26

Fresh numerical reruns: extract_cyano2026.py, analyze_rubisco.py, apply_bibliographic_overlay.py, analyze_protease.py, extract_regeneration.py, analyze_regeneration.py, rank_challenges.py, build_assay_ledger.py, build_inventory.py.

- All local analysis/acquisition scripts parse successfully.
- Ranking contains exactly all twelve challenges and reproduces primary/sensitivity selections.
- Independent brute-force Pareto calculation yields14 pairs/13 species; all226 paired records have source links, known citation overlays are correct; 2026 SE/strain values and natural/ancestor distinction checked.
- Protease join is unique802 pairs; filtered denominator, blocked failure, observed-feature coverage and motif-bootstrap direction checked.
- Regeneration count/filter,18-sample design, multiplicity outputs and threshold-attrition percentages checked.
- All downloaded XLSX/ZIP data files are valid archives; known HTML masquerading as ZIP excluded.
- 86 current downloaded artifacts match recorded SHA256 checksums; two renames and one rejected unrelated source are explicit.
- All local links in report resolve to existing deliverables.

Limitations: calculations and source integrity are verified, not original experiments; historical numeric transcriptions, phylogenetic models, all external data licenses and inaccessible raw files remain subject to the stated limits.
