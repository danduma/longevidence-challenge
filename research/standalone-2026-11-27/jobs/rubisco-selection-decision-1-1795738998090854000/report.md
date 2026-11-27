# Which Rubisco would get one confirmation slot?

For the released 30 °C coupled-carboxylation screen, **ANC_3 is the defensible single confirmation target if reconstructed candidates are eligible**. The paper already identified it as the highest measured point estimate. The new check shows it also leads on the individual assay readings before and after the authors' reference normalization, remains first after omitting each of 131 assay labels in turn, and minimizes worst-case regret under both tested dispersion scenarios. This is evidence for allocating a confirmation slot on the screened endpoint. Confirmation would still be needed to establish a common-format kinetic winner.

**For extant candidates, the released evidence does not identify an equally stable single choice.** The expanded extant point leader is Form II RBC2_36 (*Mariprofundus ferrooxydans*). Form I RBC8_39 (*Synechococcus* sp. BS56D) is the representative-survey leader. The published gap is about 0.21 s⁻¹ at the deposited multiplier, based on raw-scale means rounded to one decimal. Re-averaging the published individual reference-normalized readings gives 13.3714 versus 13.3000 s⁻¹ before global calibration, a 0.0714 s⁻¹ gap. The highest two estimates have no shared assay label. In the two top candidates, dropping CA105 or CA99 displaces RBC2_36. Across the entire extant numeric pool, eight label omissions change its leader, including those two. Ranking by the readings before reference normalization puts RBCSeed_6 first in the extant pool. That replay is a normalization sensitivity check. The authors' normalized ranking is the stronger point baseline because it uses the within-assay reference.

## Reconstructing the eligible pools

| Cohort | Released rows | Numeric | Numeric with SD | Low flag | Unresolved |
|---|---:|---:|---:|---:|---:|
| Representative screen as recovered from the release | 423 | 221 | 144 | 26 | 176 |
| Targeted clade expansion | 44 | 38 | 36 | 0 | 6 |
| Reconstructed ancestors | 2 | 2 | 2 | 0 | 0 |

We used the histogram notebook's exclusions: rows with `expanded=1` in the expansion table and IDs starting with `ANC_`. The remaining 221 numeric rows match the paper's representative count. The supplement's `new_to_this_study` flag marks a different distinction. We ranked only the authors' numeric candidates. The S02 workbook also contains readings for five low-flag rows. We preserved that published classification. S02 counts match the master table's N for all 261 numeric entries. Means and sample SDs recalculated from the one-decimal S02 normalized readings agree with deposited raw-scale summaries within 0.1 s⁻¹. The 79 singleton numeric entries lack an SD. They participate in the assay-level replay and are excluded from the interval exercise.

## Stress test of the best observed candidate

For each eligible candidate, we used a nonnegative interval centered on the deposited mean with radius equal to its SD, or to SD divided by √N. Each radius can be scaled by a common width parameter. This is an exploratory sensitivity model. The S02 assay labels do not establish independent biological preparations or the covariance of the estimates. SD/√N is therefore a narrower scenario, with no confidence-coverage interpretation. Under the rectangular model, each row can move independently within its interval. Worst-case regret is the largest admissible shortfall against any rival under the same realization. At the first tie width, the point leader's lower endpoint meets the first rival's upper endpoint. Reported regrets use the deposited units.

| Pool | Point leader | First tie in SD units | Minimax at one SD | Minimax at one SD/√N | S02 normalized leadership changes |
|---|---|---|---|---|---:|
| Representative | RBC8_39 | 0.03125, RBCSeed_3 | RBC8_39 | RBC8_39 | 6/131 |
| Extant with expansion | RBC2_36 | 0.05263, RBC8_39 | RBC8_39 | RBC2_36 | 8/131 |
| Extant Form II only | RBC2_36 | 0.05405, RBCSeed_3 | RBC2_36 | RBC2_36 | 5/131 |
| Expanded with ancestors | ANC_3 | 0.46269, RBC5_23 | ANC_3 | ANC_3 | 0/131 |

In the extant pool, switching from RBC2_36 to RBC8_39 reduces worst-case regret from 9.79 to 8.94 s⁻¹ at one SD. At one SD/√N the choice switches back to RBC2_36. Even the more comparable Form II-only pool permits a tie after just 0.054 SD of opposite perturbation to RBC2_36 and RBCSeed_3. ANC_3 is the unique possible leader at one SD/√N on the dispersion-eligible subset, and one of 12 possible leaders at one SD. Its minimax choice survives in both cases.

In the assay replay, all observations with the omitted label are removed together. Means are recalculated from remaining readings, with candidates having no remaining readings unavailable for that iteration. Counts of changed winners are descriptive, with no interpretation as independent trials or as probabilities of being best.

## What can transfer to a follow-up decision?

A common positive multiplier scales all means, widths and regrets together. Across 1.7–2.8 it cannot change the choices above. The range describes the fitted common factor, leaving candidate-specific bias unconstrained by that interval. For the deposited extant top two, the raw means differ by 13.4/13.3. A relative factor for RBC8_39 just over 1.00752 times the factor for RBC2_36 reverses the point ranking. Form I used soluble extracts and Form II purified enzyme. This makes a common-format ranking especially uncertain for that comparison, while the within-Form II result shows that the problem persists after restricting preparation class.

An unresolved entry may lack a comparable turnover endpoint in this pipeline. Even if all unresolved entries have such endpoints, the release gives no ranking upper bound for them. No uncertainty radius was assigned to a singleton. This limits the choice to confirmation within the measured collection and rules out a statement of superiority throughout the release. The assay measures saturating-CO₂ turnover at 30 °C. Selection for air-level catalysis, oxygen discrimination, or performance in a host requires other endpoints.

## Prior work and source provenance

C08 and the reassessment bound form-level hit fractions on the expanded screen. The 2025 paper and Fig. 4 already name ANC_3 as the point leader. The archive's fastest-candidate scout tabulates specificity coverage. The increment here is a release-linked replay of assay-label influence and a comparison of individual-risk decisions with the observed-rank baseline, including a matched preparation-class restriction. Neither cohort counts nor the identity of the published point leader is claimed as a discovery. The 2026 Kehl et al. study revisits sparsely sampled clades and turnover but does not provide a direct comparison of the leaders here. Plantae's specialist commentary on Davidi et al. describes the earlier Gallionella screening result. Sciety listed zero evaluations for the 2025 paper at retrieval, and the GitLab issues query returned an empty list. These checks have limited coverage.

Data and notebooks are the authors' releases, kept as local snapshots with attribution to de Pins et al. The PNAS article is CC BY 4.0. Its supplement is preserved as supplied. GitLab inputs retain their upstream status and are not relicensed. See the source entries in submission.json for URLs and retrieval scopes.

## Reproduce

Run `python3 -B analysis.py` to regenerate results.json and `python3 -B verify.py` to verify it. Python's standard library is sufficient. The verifier reconciles S02 observations to the pinned master data, recomputes label-omission leaders via grouped sums, and reconstructs interval decisions on the uncorrected scale with adverse and favorable corners. It finds first ties by bisection independently of the analysis formula.
