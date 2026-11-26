# C08 | Which high-rate comparisons survive missing Rubisco measurements?

**At 30 °C, 183 of 261 numeric rates (70.1%) meet the paper’s 10 s⁻¹ cutoff.** Including the explicitly low records gives 183/287 (63.8%). Across the 469 eligible records, the screen confirmed 183 hits (39.0%). If every unresolved record has a comparable rate, the sharp bounds on the full-library fraction are **183/469 to 365/469 (39.0–77.8%)**. A majority at or above 10 requires at least 52 of the 182 unresolved records to meet the cutoff. Both majority outcomes remain possible. These calculations concern the selected library at the deposited calibration. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/)

## Where the form comparison holds

| Form | Eligible | Numeric | Flagged low | Unresolved | ≥10 among numeric | ≥10 among eligible, lower–upper |
|---|---:|---:|---:|---:|---:|---:|
| I | 144 | 98 | 8 | 38 | 57/98 (58.2%) | 39.6–66.0% |
| I′ | 18 | 10 | 1 | 7 | 0/10 (0%) | 0–38.9% |
| II | 178 | 142 | 0 | 36 | 124/142 (87.3%) | 69.7–89.9% |
| II/III | 24 | 6 | 2 | 16 | 1/6 (16.7%) | 4.2–70.8% |
| III | 105 | 5 | 15 | 85 | 1/5 (20.0%) | 1.0–81.9% |

At ≥10 s⁻¹, three of the ten pairwise strict orderings are forced. Form II exceeds I and I′, and I exceeds I′. Even the minimum for II (124/178) exceeds the maximum for I (95/144) by 3.69 percentage points. For each of the other seven comparisons, either form can lead. For example, if all 85 unresolved form III records met the cutoff and none of the 36 unresolved form II records did, III would lead 86/105 to 124/178. Reversing those assignments gives II the lead. These are endpoint assignments proving identifiability bounds. The data supply no probability for either assignment. [All signed contrasts and attainable ties](form_comparisons.csv)

At the secondary cutoff of ≥15 s⁻¹, the numeric, ascertained and full-library fractions are respectively 115/261 (44.1%), 115/287 (40.1%) and 115–297/469 (24.5–63.3%). Only II over I and II over I′ remain forced. For II versus I the minimum gap shrinks to 0.53 percentage points (85/178 versus 68/144). This narrow margin conditions on the reported means and their common calibration. It does not establish robustness to measurement error.

The overall numeric median is 14.05 s⁻¹. Adding the low flags identifies the ascertained median as 13.20 s⁻¹ without assigning any flagged value. Including the unresolved records yields a full median between 2.98 and 19.59 s⁻¹. The form II median remains above the form I median even at their closest admissible endpoints (14.69 and 14.16 s⁻¹). For III, 15 low flags among 20 ascertained records put the ascertained median below the censoring ceiling, but the 85 unresolved records allow the full median to be arbitrarily high. [Medians and denominator results](denominators.csv)

## Selection and ascertainment

The 180 entries labeled new to the study have a higher ≥10 fraction among numeric rates than the 289 previously screened entries (50/66, 75.8%, versus 133/195, 68.2%). Including explicit low flags **reverses the observed ordering** (50/82, 61.0%, versus 133/205, 64.9%). Confirmed yields are 27.8% versus 46.0%. Full-library bounds overlap (27.8–82.2% versus 46.0–75.1%). The new cohort includes all 105 form III entries. An additional descriptive form-by-origin cross-tab shows that 49 of its 50 numeric hits come from form II. The study expanded two fast clades after its initial screen. This composition precludes interpreting a cohort difference as a change in biological rates over time. Nine source records carry the eukaryotic label. All three numeric rates there meet 10, leaving the eligible fraction between 3/9 and 9/9. [Study methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/) · [Stratified results](denominators.csv)

## Source reconciliation and assumptions

The calibration notebook selects the master CSV’s 469 lab entries from 2,260 rows, corrects the numeric lab means and SDs by a common multiplier, and carries the low-rate flag forward. The verifier joins by unique internal ID, checks identifiers, sequences, form, taxonomic label, flag, numeric mean, SD and assay N against the master, and checks every Supplementary Data 3 field and reconstructed state against the earlier ledger. The paper’s 221 numeric measurements for its representative survey refer to a narrower selection. The present denominator includes the expanded collection, controls and reconstructed entries. A row receives one vote. The authors’ sequence-identity coverage has a different denominator from the endpoint coverage here. [Pinned notebook](https://gitlab.com/milo-lab-public/rubisco_is_slow/-/blob/e1c9ed945864a0744e81b9cfdbfeb609ccf67907/1_regression_measured_lit.ipynb)

Figure 2 groups unresolved experimental outcomes as insoluble or inactive. The released lab table supplies a numeric value, a below-0.5 flag, or neither. It has no row-level field separating insolubility from inactivity. Neither an expression count in the paper nor a blank lab value establishes an intrinsic turnover rate. The primary analysis uses the labeled 0.5 s⁻¹ ceiling. Because the notebook rescales the numeric values while passing the flag through, the verifier also tests a ceiling of 0.5 times the deposited factor (1.065 s⁻¹) if that flag originally used the raw scale. Both ceilings are below both hit thresholds. The overall medians and their bounds are also unchanged. The form III ascertained median’s upper limit follows whichever ceiling applies.

For a group with N entries, H observed hits and U unresolved rates, every integer hit count from H through H+U is compatible with the observed states. Groups in a form contrast are disjoint, so their endpoints can be attained jointly. The verifier tests ties on the exact integer grids. For medians it orders all lower endpoints, then all upper endpoints, averaging the two central values for even N. Monotonicity gives sharp limits by moving unresolved values simultaneously toward the corresponding endpoints. Zero denotes the nonnegative floor or infimum. A flag ceiling denotes a supremum, and infinity denotes the absence of a finite upper bound. This calculation imputes no individual kinetic value. If an unresolved enzyme lacks a well-defined comparable endpoint, confirmed yield remains interpretable and the latent whole-library frequency requires that extra assumption.

**Falsifiable scope.** Holding eligibility and existing rates fixed, completing the missing records cannot put form I level with or ahead of II at either threshold. Such a result would falsify the source reconciliation or endpoint assumptions. Both signs are attainable for II versus III. Neither the finite-library bounds nor the median ordering estimates frequencies across natural Rubisco diversity. The previous synthesis established the three state totals and calibration limitation. C08 adds the denominator-specific frequencies, order bounds and selection reversal. No human expert has reviewed C08.

## Reproduction

From the C08 directory, run the command below. The default verifier reads pinned snapshots, checks SHA-256 hashes, and compares regenerated result tables byte-for-byte. The calculation uses the Python standard library and ran with Python 3.14.10.

```sh
../../../millennium/.venv/bin/python -B verify.py
```

For public reproduction, download the three files listed under `download_url` in [sources.json](sources.json) into one directory, keeping their original basenames. Run `python3 -B verify.py --public-inputs /path/to/downloads` alongside the distributed result tables and manifest. This mode checks primary-input hashes and regenerates the same tables without the local prior-ledger cross-check. The manifest pins the 2025 paper, the GitLab commit and blob IDs, snapshot hashes and access date. Primary inputs stay in their existing locations. [Protocol](protocol.md) · [Row-state audit keyed by supplement CSV line](states.csv) · [Verifier](verify.py)
