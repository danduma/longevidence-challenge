# Codex CLI research configurations — observed results

Date: 26 November 2026. **Use Sol 6.1 medium for substantial analysis/code and Astra 6.2 low as a fast independent critic.** Sol high and Astra medium are useful for deeper discovery and targeted extensions. This is a provisional workflow recommendation from the actual runs below, not a general model leaderboard. No global model configuration was changed.

## Identical-task comparison

Two independent workers received byte-identical prompts and the same pinned source data. They were prohibited from reading the coordinator's Rubisco analysis or each other's reports. Both completed their actual audit and fit a predeclared lightweight baseline, rather than returning a plan. Their generated analyses were subsequently rerun by the coordinator; selected numerical output hashes were unchanged. [Verification evidence](archive/verification/release-and-prior-art-reproduction.json).

| Configuration | Wall time | Input minus cached tokens | Reported output tokens | Useful distinctions |
|---|---:|---:|---:|---|
| Sol 6.1 medium | 531.88 s | 103,659 | 20,941 | Richer composition/dipeptide baseline; exact-sequence/species groups plus a ≥99% sequence-component stress test; training-mean controls and component bootstrap. |
| Astra 6.2 low | 324.09 s | 46,842 | 7,982 | Independently recovered the core audit and found a source-confirmed grass-to-moss metadata mismatch; added a clearly labeled post hoc fragment/taxon sensitivity. |

Both recovered 214 sequence records / 207 headers / 169 distinct sequence strings, the repeated deposited targets, the input/output whitespace discrepancy, fixed versus nested output branches, and conditional reconstructed overlap. Both identified missing wild-type embeddings and avoided pretending to refit the original model. They also correctly distinguished CO₂ affinity from CO₂/O₂ specificity and marginal retrospective scores from joint experimental challenge success.

Their baseline scores cannot be directly compared as model-quality outcomes: feature sets, ridge penalties, cultivar grouping, split methods and data counts differ. Sol's richer baseline found stronger signal; that is not a fair demonstration that its scientific algorithm beats Astra's. Neither worker established prospective candidate improvement or universal novelty.

## Evidence against the predeclared audit rubric

The ten criteria are in EVALUATION.md (archival input; see source-downloads.json and the source URLs). Both delivered the core raw-identity/prediction counts, explicit split assumptions, metric/reweight/refit distinctions, target provenance, executable grouped baselines, mean controls, source checks and challenge-aware interpretation. Both scripts passed the coordinator's numerical reproduction. Sol's ≥99% stress test covers a stronger sequence-shift condition than Astra's exact-sequence-only primary split. Astra's wrong-organism check provides a concrete primary-source curation finding beyond the common defect list.

One minor Sol report error is preserved for audit: it says the all-state GitHub listing contained three PRs, while its own `web_evidence/issues_api.json` contains four records, as do the coordinator's and Astra's snapshots. All four were inspected by the coordinator and none documents these specific defects. This does not change the source-data result, but it prevents treating either narrative as an infallible authority.

No numeric overall quality score is assigned. One run per model/effort, different analytical choices, shared local skill availability, and coordinator-only adjudication cannot establish general superiority or reliability. Astra's shorter runtime on this task is measured; a universal speed or cost claim is not.

## Different-task discovery round

| Configuration / starting point | Wall time | Input minus cached tokens | Reported output tokens | Actual artifacts |
|---|---:|---:|---:|---|
| Sol 6.1 high, data-first | 1,730.87 s | 288,862 | 71,668 | Twelve-problem inventory, source/citation audits, executed Rubisco/protease/regeneration investigations. |
| Sol 6.1 high, mechanism-first | 1,221.03 s | 260,473 | 50,735 | Twelve-problem assessment; executable cryopreservation, quadruplet-code and Thg1 calculations. |
| Astra 6.2 medium, skeptical prior art | 1,031.15 s | 103,516 | 23,829 | Twelve-problem ranking; independent source calculations; counterexamples and recent prior-art checks. |

These prompts differ deliberately; this table does not compare equal work. Full/cache/raw usage, commands, prompts, failures and terminal status are preserved per agent, with a generated [run summary](runs.json). Tokens are reported as the CLI records them; reasoning-output counts remain separate raw fields and are not added to reported output totals. No dollar costs were inferred.

## Practical division of work

Use separate data-first and mechanism-first Sol workers to avoid one favored hypothesis controlling the search. Give Astra a specific falsification/source-provenance task, and retain concrete counterexamples even when they contradict a Sol result. Use identical packets when comparing configurations. Evaluate usable verified artifacts before verbosity, speed or refusal count. Increasing effort is justified by an unresolved scientific or implementation question, rather than a policy stop being treated as model difficulty.

Five completed research runs in this comparison reported **no model refusals**. The subsequent Thg1 and protease reviews also completed without model refusals, bringing the overall completed count to seven; they are not part of the matched comparison. HTTP errors, robot checks, unavailable files and TLS issues occurred and were logged separately. The prior Opus session's safeguard stops are historical events outside these counts. The corrected-population Sol investigation and current Rubisco Sol investigation subsequently completed, as did a deeper Astra medium thermal investigation. The independent Sol medium thermal review subsequently completed as well: all eleven completed jobs report no model refusals.

The targeted reviews added value by narrowing claims: Sol's exact arithmetic verified the Thg1 proof while rejecting its interpretation as substantial challenge-level progress, and Astra reproduced the protease RMSE result while showing that the below-prevalence AP conclusion changes under eligible controls. Both reviews were freshly rerun by the coordinator with identical selected outputs (`verification/critical-review-reproduction.json`). Their critical conclusions are retained rather than counted as a failure to answer biology questions.

## Completed targeted extensions

| Configuration | Wall time | Result and verification state |
|---|---:|---|
| Sol 6.1 medium, corrected protease population | 937.51 s | Exact final data/source-figure reconstruction; weakens the earlier raw-error conclusion. Full coordinator rerun passes, including three arithmetic tests and source reconciliation. |
| Sol 6.1 high, current Rubisco controls | 1,401.49 s | Source/control ledger, calibration identifiability and uncertainty-provenance audit. All 56 checks pass on fresh coordinator rerun; selected hashes unchanged. |
| Astra 6.2 medium, thermal robustness | 833.09 s | New thermal benchmark/counterexample; full coordinator rerun passes, selected hashes match, and independent Sol implementation supports the precise result. |
| Sol 6.1 medium, independent thermal review | 829.59 s | Independently implements the equations, finds prior art narrowing the concept claim, and rejects an unsupported measured-triplicate inference while accepting the specific numerical contribution. |

These different tasks further demonstrate usable contributions from both model families, not comparative superiority. The protease extension is particularly valuable because it overturns an attractive interpretation rather than preserving it. The coordinator's additional fixed hazard-score experiment failed and is retained as a negative result.

The thermal pairing demonstrates complementary strengths: Astra medium produced a useful new computed comparison; Sol medium verified the equations independently and materially narrowed the empirical and priority claims. The coordinator accepted those corrections and freshly reran both the original full pipeline and four central cases in the separate solver. See [accepted findings](FINDINGS.md); the original uncorrected worker report remains an explicitly labeled research record.
