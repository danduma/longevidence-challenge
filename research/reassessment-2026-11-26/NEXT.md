# What next

This is a prioritized research brief, not a record of completed work or permission to contact authors. No upstream submission, new model training or experiment is claimed here.

## First priority: finish one useful data correction

**Question:** Does repairing the CGA009 sequence–measurement pairing produce a trustworthy corrected benchmark, and does it change any reported prediction or interpretation?

1. Recheck the source assay and exact construct context against the released row. Keep strain identity and crystallography tag context separate.
2. Prepare a minimal correction to the released metadata/sequence mapping, with immutable source links. Keep the previously known grass/moss issue and fragments separate so this one repair can be reviewed on its own.
3. Locate the exact released feature-generation configuration and cached feature row. Regenerate the replacement embedding with that checkpoint/pooling convention; never substitute a label on the old embedding.
4. Hold the existing folds, preprocessing and target definitions fixed. Compare original and corrected feature paths and report paired changes in affected predictions, aggregate metrics and conclusions. A negligible change is a legitimate result.

**Success:** A source-supported, independently runnable correction that preserves provenance. Claim a performance benefit only if the controlled rerun demonstrates it. **Stop or narrow:** source identity cannot be bridged, the public version is already repaired, or the reproduction needs unavailable inputs. Document that outcome rather than expanding the claim. Contacting maintainers or opening an upstream PR needs separate authorization; publishing our own evidence here does not imply acceptance there.

## Second priority: make the library inference decision-relevant

**Question:** Does C08's surviving≥ 10 s⁻¹ bound help allocate work within the actual collection?

Separate the representative survey, deliberately oversampled clades and unsuccessful-expression outcomes. Introduce enzyme-specific uncertainty only where the source supplies defensible data. Compare whether a choice among already released candidates changes under the bounds and calibration assumptions. Avoid interpreting form labels as a causal property or inferring natural prevalence.

**Success:** A clearly named selection or interpretation remains justified under the relevant uncertainties. **Stop or narrow:** individual errors or endpoint noncomparability make that choice unidentifiable. Do not fill missing measurements with invented rates.

## Deprioritize broad thermal sweeps and weak-baseline prediction

- **Thermal work:** pursue a focused benchmark/objective audit only if it implements the original stress objective and temperature-based termination, with physically justified temperature support and adequate numerical resolution. A higher all-time warming rate alone is insufficient. Reassess utility before another parameter grid.
- **Protease work:** retain the stronger motif baseline. Test a clearly defined genuinely unseen-motif or unseen-enzyme question against feasible sequence/context baselines, preserving training-only normalization and grouped dependence. Do not call a small gain meaningful without defining the decision or calibration accuracy it improves. The original AUC-only result can remain as a restricted comparator.

## Shared first gate

Before scaling any of these paths: inspect primary/SI/code evidence and later work, search relevant specialist/community archives, state the remaining novelty precisely, and run one adversarial comparator or analytic check. Use the [skill](../../skills/auditing-research-contributions/SKILL.md). Keep work bounded to what the user has authorized; no automatic laboratory work, outreach or open-ended agent campaigns follows from this brief.
