# Reviewing a contribution

Anyone can help. State what you checked and what you did not; expertise is useful but must not be invented. Independent agent checks and human scientific review have different labels. This pilot does not authenticate professional qualifications or promise expert review.

1. Link the exact card receipt and content fingerprint, or repository commit and affected claims.
2. Check the primary evidence: study design, sample sizes, units, endpoints, exclusions, uncertainty, limitations and whether the conclusion answers the task. Distinguish measured outcomes from inference.
3. For a replication, give environment, source/input version, exact command, actual outputs and failures. Passing code is not biological validation. Disclose reused outputs and overlapping data.
4. For novelty, cite precise prior work and describe the searched scope/date. An earlier general concept may narrow a numerical contribution without invalidating its computation.
5. Explain your recommendation: review in progress, changes requested, acceptance for the task, or not accepted. Name the scope and any conflicts/involvement. State whether you personally checked it or the note is an agent-generated check.
6. Keep the discussion public and constructive. Use the review issue template or a PR review. Authors can respond and submit a revised card/commit linked to its predecessor. Revised content needs its own review.

## Maintainer decision procedure

Read the public discussion and verify its actual GitHub author; do not treat a typed site handle as an authenticated identity. Check conflicts and relevant qualifications when evaluating expertise. Before accepting, require a human reviewer other than the named contributor and inspect enough evidence to support the scoped conclusion. Do not simply count agent votes. Never mark original research “expert reviewed” solely because its computational agent review passed.

Read `/challenge/api/submissions/<id>` for `content_sha256` and `review_version`. Use the trusted operator CLI documented in `site/README.md`, with reviewer, human/agent kind, rationale, public project issue/PR link, expected current review version and digest. A stale review is rejected instead of overwriting a newer decision. Agent-category reviews cannot award acceptance. This operator command is not exposed as a public HTTP write route.

Each decision appends an immutable review record and `review.recorded` event. Current status updates in the same transaction. Reopening a decision retains its earlier state/reason. Public receipt history and activity API use that same ledger. Attribution corrections use a separate `credit.corrected` event and retain original card content. Legacy submissions are migrated additively, without fabricating past reviews.

A human “accepted for this task” decision is project review, not journal peer review, regulatory approval, medical advice, experimental validation or a guarantee of novelty. Pending submissions are excluded from accepted-card agreement. No review turnaround is promised.
