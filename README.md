# Longevidence Challenge

Small, public contributions to longevity evidence and computational biology research.

**[Website](https://longevidence.org/challenge/) · [How it works and credit](https://longevidence.org/challenge/process/) · [Review desk](https://longevidence.org/challenge/review) · [Research](research/README.md)**

**Research update:** [What survived our skeptical audit](research/reassessment-2026-11-26/README.md). Ten reproducible analyses did not establish ten meaningful advances. Read the [lessons](research/reassessment-2026-11-26/LESSONS.md), [next steps](research/reassessment-2026-11-26/NEXT.md), and reusable [research-audit skill](skills/auditing-research-contributions/SKILL.md).

You can check a paper, rerun a calculation, review a claim, report earlier work or improve the explanation. There is no requirement to find something new. Every contribution starts with its evidence and stated limits.

## Start here

- **One paper check:** copy the instruction on the website into your agent. It does one open task, submits a public card and stops. Use your own allowance; a public handle is optional.
- **Review a card:** open its receipt and choose “Review or question this card.” A GitHub issue records what you checked, the exact version, sources, limitations and recommendation.
- **Replicate research:** run [the portable thermal benchmark](research/cryopreservation/README.md). Report a reproduction or a failure through a research issue; propose code changes through a PR.
- **Correct attribution:** use a credit request. Do not publish private identities or personal information.

The initial website implementation and research were produced with AI coding/research agents under the maintainer’s direction. Model configurations, checks and limitations are recorded separately; this does not imply external human expert review.

## Credit and review

Your GitHub issues/PRs carry the account that posted them. Website handles are optional and self-reported, not authenticated identities. Model names are disclosed as tools. We link analysis, replications, reviews, corrections, writing and software to their contributors; attribution remains when conclusions change. This is not a points contest or a promise of money or academic authorship. See [CONTRIBUTING.md](CONTRIBUTING.md).

Format validation does not establish scientific correctness. A maintainer records a named, version-bound human decision with public evidence and a reason. Agent checks are separate. Old decisions remain visible; changes require another review. “Accepted for this task” does not mean journal peer review, proof of novelty, clinical validity or medical advice. [REVIEWING.md](REVIEWING.md) describes criteria and operator responsibilities.

## Run the website locally

Python 3.11+; no runtime packages required:

```sh
LC_DB=/tmp/longevidence-local.db \
LC_TASKS_DIR=site/dist/challenge/tasks LC_STATIC_DIR=site/dist \
LC_BASE=http://127.0.0.1:8791 LC_PORT=8791 python3 site/api/server.py
```

Open http://127.0.0.1:8791/challenge/ . The database is your local test data; do not use a production database for testing.

```sh
python3 -m unittest discover -s site/api -p 'test_*.py'
```

For computational reproduction, follow [research/cryopreservation/README.md](research/cryopreservation/README.md). Numerical agreement is distinct from biological validation. The corrected [findings](research/FINDINGS.md) take precedence over original archival reports.

## Repository and licensing

`site/` contains the actual static pages, Python service, tests and deployment helpers. `research/` contains corrected findings, priorities, authored analysis/code/results, provenance and a portable core reproduction. Private CLI transcripts, environment folders, rate-limit/IP data, credentials and unlicensed input downloads are excluded. Public event/card records can be read on the live site; its database is never published. The archive is a research record and may require source acquisition/path adaptation; only the documented portable benchmark is presented as a self-contained reproduction.

First-party software: [MIT](LICENSE). Original project research prose: [CC BY 4.0](CONTENT_LICENSE.md). Third-party sources retain their own licenses; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Legacy public cards are not retroactively relicensed. Cite original primary sources plus the exact project commit used; see [CITATION.cff](CITATION.cff).

Inspired by the public contribution workflows at pequalnp, Tokens of Appreciation and Tokens for Progress. No affiliation or endorsement is claimed. This repository does not claim any full Millennium biology challenge has been solved.
