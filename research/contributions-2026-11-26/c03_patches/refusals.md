# C03 access and refusal record, 2026-11-26

No model refusals, other model calls, subagents, people contacted, external posts, commits, deployments or biological work. This file records task operations only; no credentials or private session logs are included.

Source failures: browser open of immutable GitHub tree returned restricted URL; DOI publisher redirect returned non-retryable access error; Europe PMC article HTML was inaccessible through the browser tool. Direct Europe PMC model XML request timed out at 30 s. Used the already pinned original model XML from prior work, with hash and license checked, instead. Direct Europe PMC heart XML download succeeded (HTTP 200); raw GitHub license browser retrieval succeeded. See sources_live.json. These are source-access outcomes, not model refusals.

A process-inspection command (`ps`) was denied by the local environment; it was unnecessary for computation and was not retried. Full-field computation itself ran successfully. Initial verifier incorrectly required complete myocardial −123→−35 crossings even though the protocol explicitly called for reporting incomplete crossings; corrected that assertion without changing maps, model, scoring, or thresholds. No inferred crossing was inserted.

The first one-cell crossing check used an insufficient 70 s scalar horizon, leaving the upper crossings incomplete. Corrected the diagnostic horizon to 100 s and used 0.001 s steps; all four analytic crossings then matched quadrature within 0.000061 s. This changes only the verification fixture, not any full-field comparison or protocol endpoint.
