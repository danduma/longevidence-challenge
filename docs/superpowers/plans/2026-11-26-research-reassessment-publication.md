# Research reassessment publication

**Goal:** Publish the surviving findings and disconfirming evidence, correct the public interpretation, and install a reusable research-audit skill.
**Architecture:** Preserve original numerical artifacts; add a controlling reassessment with executable checks and update entry points and status metadata.
**Tech Stack:** Python, Markdown, existing static website and GitHub publication.
**Current Milestone:** A published reassessment, actionable next-work brief, and installed/validated skill. No upstream contacts or new biological experiments.

### Task 1: Publish inspectable evidence and lessons

**Files:** `research/reassessment-2026-11-26/{README.md,verify.py,results.json,LESSONS.md,NEXT.md,requirements.txt}`.

- [x] Reproduce the C07 accession/deposition check, C08 common-factor sensitivity, and stronger C10 held-enzyme baselines; document exact inputs and scope.
- [x] Write the claim-by-claim reassessment, sources/forum coverage, lessons, and prioritized next-work brief with success and stop criteria.

**Verification:** Run the calculations against hash-pinned sources and compare generated results; read every conclusion against its actual check.

### Task 2: Correct public interpretation without erasing history

**Files:** Batch README/AUDIT/ledger and bundle assessment notices; research/root entry points; research and process website pages; bundle/global manifests.

- [x] Distinguish ten reproduced artifacts from established scientific advances; mark C07/C08 scope and C10's omitted comparator.
- [x] Put the reassessment before historical claims in the website and repository; retain original reports and attribution.

**Verification:** Check links, status consistency, unchanged original numerical files, and all public hashes; inspect final prose.

### Task 3: Create and install the reusable skill

**Files:** `skills/auditing-research-contributions/SKILL.md`, focused case-study reference and UI metadata; matching installed skill in the current Codex skill directory.

- [x] Encode source triangulation, field-specific forum searches, strong baselines, actual objectives, uncertainty, and significance-based reporting.
- [x] Validate discovery/frontmatter and apply the skill to realistic audit cases; verify installed/public copies agree.

**Verification:** Skill validator, trigger/non-trigger examples, and worked C10/C07/C08 assessments. No agent is credited as a human reviewer.

### Task 4: Publish and verify

- [x] Commit/push the public source and reassessment; deploy the existing static site.
- [x] Read back GitHub and live pages, check applicable CI, and record the final publication evidence.

**Verification:** Exact published content/hash checks and successful relevant CI. Checkboxes stay open until their requirements pass.

Publication evidence: implementation commit `c4fc53471daf29918d656af8f70deaeb4a093601` pushed; exact GitHub results/skill readback and exact live research/process content match. Deployment helper also copied pre-existing unrelated homepage/style edits; those two live files were restored from its backup, preserving the workspace edits. Remote diff is restricted to the two authorized pages. Applicable implementation CI is still being checked.

Final implementation verification: GitHub run `45847064901` succeeded (API/database tests and existing thermal benchmark). Reassessment checks passed both fresh public-download and subsequent offline modes. Skill validator and public/installed byte parity passed. All409 original numerical/code/protocol/report files remain unchanged. Public input downloads are ignored, not staged. Final live source matches include research/process and the restored homepage/style.
