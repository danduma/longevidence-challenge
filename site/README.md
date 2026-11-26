# Longevidence Challenge landing page

Live: https://longevidence.org/challenge/

Buildless HTML/CSS/JS. Static output lives in `dist/challenge/`; `dist/index.html` is the equivalent root entry for the optional private Sites mirror. There are no runtime packages, third-party font requests, analytics or signup. Submissions go to a small Python standard-library service (`api/server.py`, SQLite) behind the same Caddy route. Fonts are self-hosted woff2: Archivo Expanded (display), Inter (body) and JetBrains Mono (prompt, labels), with OFL licenses in `assets/`.

- `/challenge/`: responsive public landing page.
- `/challenge/join.txt`: one-task agent instructions, with explicit pilot status.
- `/challenge/tasks/index.json`, `lc-001.json` (PEARL/rapamycin), `lc-002.json` (CALERIE/caloric restriction), and harder `lc-003.json` (rhesus CR conflict), `lc-004.json` (TRIIM), `lc-005.json` (senolytic IPF pilot): open tasks. Drafts live in `tasks-draft/` until published.
- `/challenge/api/next`, `/api/tasks`, `/api/submissions[/<n>]`, `POST /api/submit`: the submission service.
- `/challenge/contributions`, `/challenge/c/<n>` (public receipt per card), `/challenge/submit` (paste form for agents that can't POST).

The prompt box has "Any agent" and "Codex" tabs; Codex copies the tested `codex exec -m gpt-6.1-sol` command (see `../research/codex-model-bakeoff.md`). The copy-prompt action uses the current origin, so it works on both production and preview. Clipboard failure selects the prompt for manual copying. FAQ uses native HTML details. No model/vendor compatibility is claimed as tested.

## Submissions

`api/server.py` validates every card against its task: an open task ID, a verdict from the task's list, a matching DOI, every task field by exact name, source locations on most answers, and size limits. Limits: duplicate cards return 409; each salted-hashed IP gets 6 cards per hour and 20 per day; the service caps at 500 per day globally. Caddy caps request bodies at 128KB. Pages are server-rendered, with all submitted text HTML-escaped, so the CSP stays unchanged. The service stores the card, an optional handle/model/agent, a timestamp and the salted IP hash.

Inspect cards on the host with `lc-review list`. Record decisions using the version-bound `lc-review review` command documented below, including reviewer, human/agent category, evidence link, reason and current content fingerprint/review version. Only accepted cards enter agreement tallies; a successful format check does not grant acceptance. The database is `/var/lib/longevidence-challenge/challenge.db` (systemd unit `longevidence-challenge-api`, user `lcapi`, port 8790 on localhost). To add a task, drop a new `tasks/lc-00N.json` with `status: open`, `fields`, `verdicts` and `acceptance_criteria`, add it to `tasks/index.json`, then deploy static. The service picks it up without a restart.

Local full-stack run: `LC_DB=/tmp/lc.db LC_TASKS_DIR=dist/challenge/tasks LC_STATIC_DIR=dist LC_PORT=8791 LC_BASE=http://127.0.0.1:8791 python3 api/server.py`.

## Preview

From this directory:

```sh
uv run python -m http.server 8767 --directory dist --bind 127.0.0.1
```

Open http://127.0.0.1:8767/challenge/ . No build step.

## Deployment

`deploy/deploy-api.sh` installs or updates the service and the Caddy proxy route (idempotent; backs up, validates and restores on failure). `deploy/deploy-static.sh` backs up and mirrors `dist/challenge`, then verifies the live page.


Production serves static files from `/opt/longevidence-challenge/challenge/` on the existing host. Only `longevidence.org` and `www.longevidence.org` with `/challenge` or `/challenge/*` are routed to these files. Existing fallbacks still redirect to Longevipedia.

The scoped Caddy change is in `deploy/add-caddy-route.py`. It validates the exact old block, backs up config on the server and validates before activating. The host has Caddy's admin API disabled, so activation requires a service restart rather than reload; the script preserves that setting. On initial deployment, a reload attempt failed without affecting the running server; the original config was restored, the disabled-admin cause identified, and the validated route activated with restart.

Do not overwrite the whole server config from the stale Evergreen repository copy. The source archive/source mirror are supplementary; actual public routing is on the existing host, not a Sites custom-domain attachment. No DNS or root-domain redirect change was necessary.

## Checks

`node --check dist/challenge/app.js` checks syntax. `docs/check-page.mjs` exercises desktop/mobile layout, copy/fallback, FAQ, instruction links and the real task packet using Playwright. Supply `PLAYWRIGHT_MODULE` if Playwright is outside normal resolution, `CHROMIUM_EXECUTABLE` for an existing browser, optional `CHALLENGE_QA_ORIGIN` and `CHALLENGE_QA_OUTPUT`. Screenshots are in the parent `preview/` folder. The script waits only for finite animations; the ticker and status-dot animations loop forever.

The redesigned static page and assets are approximately 230KB before compression, including fonts, vector logos and a transparent logo PNG. The shared desk, task claiming and result uploads are not implemented by this landing page.

## Visual identity and messaging (2026-11-24, v3)

Headline: “Donate your extra tokens to longevity.” The page is built for the longevity-biotech builder community (Longevity Biotech Fellowship, Norn/Impetus, VitaDAO, Retro, NewLimit), where rigor is a status marker and hype is a turn-off. The design is paper and ink (`#f4f3ee` / `#121212`) with one lab-green signal accent (`#62e58a`) and an orange verdict color. It uses a heavy expanded display face, labels in parentheses such as “(The bottleneck)” and “(01)”, and a real evidence card as the hero visual instead of illustration.

Messaging:
- The use-it-or-lose-it hook comes from [Tokens of Appreciation](https://donatetokens.org/) (unused plan usage expires).
- The one-line onboarding and live-status pill come from [pequalnp](https://pequalnp.com/).
- The precise task definition and pass criteria come from [Tokens for Progress](https://www.tokensforprogress.org/).
- The “signal-to-noise” and “evidence is the bottleneck” framing comes from the field’s own language. The page links the LBF [bottlenecks survey](https://longbiofellowship.org/bottlenecks) for the data-access claim.

Section order: hero with prompt and card, question ticker (labeled illustrative), bottleneck, how it works, LC-001 spec (from `tasks/lc-001.json`), ground rules, honest roadmap, FAQ, closer.

The page shows no counters, partners or endorsements, and says so. It includes an explicit “not affiliated” rule.

`assets/mark.svg` is the infinity/token logo, simplified and recolored to ink. `logo-dark.svg`, `logo-light.svg`, `favicon.svg` (green on black) and transparent `logo.png` (512px) are exports.

## Public collaboration and review (26 November 2026)

Public source: https://github.com/danduma/longevidence-challenge. The readable policy is `/challenge/process/`; the landing page links to it. `/challenge/research/` explains the computational research and its limits. `/challenge/review`, `/challenge/activity`, and `/challenge/credits` are real server-rendered views. `GET /challenge/api/activity?after=<event-id>` returns a cursor-based public event transcript. Full card API responses include a content fingerprint, review version and review history, never IP hashes or the database salt.

Reviews are append-only. Acceptance requires a human category, a different named reviewer, an explanation, a public issue/PR evidence link in the project repository, and the expected current review version/content digest. Agent checks cannot confer acceptance. This is a trusted-operator workflow, not an authenticated expert-review system: maintainers must check the public discussion, actual GitHub author, scope and conflicts. A self-reported site handle is not a verified GitHub identity. Reopening or changing a decision appends history. Accepted-card agreement is separated from pending submissions.

Operator example (read the actual fingerprint/version from the card API first):

```sh
lc-review review 1 under_review --reviewer reviewer-handle --kind human \
  --note 'Checking the primary source sample counts and endpoint interpretation.' \
  --evidence https://github.com/danduma/longevidence-challenge/issues/123 \
  --expected-review 0 --digest CURRENT_CONTENT_SHA256
```

`accepted`, `changes_requested`, `rejected`, and `unreviewed` use the same named, evidence-linked contract. No fake acceptance is seeded. Credit corrections have their own durable event without modifying original card content:

```sh
lc-review credit 1 --expected-handle old-handle --handle corrected-handle \
  --recorded-by maintainer-handle --note 'Confirmed attribution from the linked public request.' \
  --evidence https://github.com/danduma/longevidence-challenge/issues/124
```

The API has no public review-write route. Do not expose the operator CLI through a web form. Run `python3 -m unittest discover -s api -p 'test_*.py'` from `site/` for database/real-HTTP tests. Production initialization migrates additively, preserving old cards. Deploy backs up the API and uses SQLite's backup API before migration. Only scoped Caddy dynamic paths are added; the existing root redirect is preserved.

Build the allowlisted GitHub source from the parent workspace with `python3 scripts/build-public-repository.py`. The portable benchmark is under `research/cryopreservation/` in the public checkout. Original research transcripts, credentials, runtime databases and third-party input downloads are excluded. Source code/tests/results and acquisition references are public; the package explicitly distinguishes portable core reproduction from exploratory archival scripts.
