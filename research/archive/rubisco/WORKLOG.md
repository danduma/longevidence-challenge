# Rubisco investigation — 2026-11-26

## Scope and preregistration of the next checks

We are evaluating public data and prediction methods as a route to selecting Rubisco candidates against the challenge's two-objective criterion. No new enzyme or experimental improvement is claimed.

The initial input audit found seven duplicated sequence headers in the 2026 Young et al. public release at commit `ade87f378aa25c24a24c19915722aa04be92936e`. The deposited prediction files independently confirm duplicated species/targets. This is an observed data-quality issue, not yet a verified novel contribution. The source paper already discusses small-subunit effects and measurement heterogeneity; those general observations are not novel.

Next checks, set before fitting any new predictive model:

1. Verify duplicated records across source data, metadata, FASTA and deposited predictions without relying on positional alignment. Record ambiguous and missing matches explicitly.
2. Reconstruct exposure to identical species/sequences under the published K-fold procedure and its default seed. Describe it as a reconstruction where the original fold assignments are not deposited. Do not infer the direction or size of score bias without refitting.
3. Check source target values against the deposited targets and identify any unreproducible preprocessing. Preserve cultivar distinctions and do not combine incompatible assay temperatures.
4. Build a corrected, accession-linked and species-deduplicated evaluation dataset, documenting the inability to identify assay-sample genotypes where necessary.
5. Compare simple sequence baselines under row, exact-sequence-group and genus holdout, with all fitted transforms restricted to training data. Include a mean-prediction control. Group identity and data overlap must be checked independently of scores.
6. Search the 2026 paper, repository issues/history and wider literature for this precise release-specific defect and analysis. Ask independent Codex reviewers to reproduce or falsify the finding before drafting a bounded novelty statement.

## Evidence obtained

- `python3 -m unittest discover -s research/millennium/rubisco -p 'test_*.py'`: four tests pass, covering whitespace/cultivar identity, conflicting FASTA IDs, incompatible source/temperature pairings and missing assay-condition handling.
- `python3 research/millennium/rubisco/input_audit.py`: 634 source rows; 504 primary wild-type rows; 315 raw species keys vs 308 whitespace-normalized keys; 214 sequence records / 207 headers / 169 sequence strings.
- The six deposited prediction tables contain `(rows, unique species)` of KC `(197,190)`, KO `(169,164)`, KRuBP `(34,28)`, S `(183,176)`, vC `(181,174)` and vO `(38,38)`.
- With the shipped default shuffled five-fold seed 0 applied to deposited row order, identical-species exposure is KC 10/197, KO 8/169, KRuBP 10/34, S 12/183 and vC 14/181. These are reconstructed exposures, not recovered original fold assignments. Cross-species sequence duplication may add further exposure.
- The original Flamholz source and the redistributed source each contain 634 rows and 19 untrimmed species labels. Source-to-output reconciliation remains required; the deposited predictions cannot yet be called reproduced.

## Limitations

The small-subunit sequence is not included. The downloaded large-subunit sequence may not be the genotype assayed in the historical publication. Identical downloaded RbcL does not establish identical full enzymes. The particular wild-type ESM embeddings referenced in the reproduction command are absent from the GitHub tree; the original prediction metrics have not been independently regenerated.
