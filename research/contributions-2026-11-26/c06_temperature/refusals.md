# C06 — refusals and source failures

No model refusal occurred. No subagents, other model calls, people contacted, external posts, commits, deployments or laboratory/animal/human work were used. No credentials or private session logs were used as evidence.

Source/tool failures, distinct from refusals:

- PMC supplement endpoints returned HTTP200 **HTML browser-check pages**, rather than the requested PDF/DOCX. The saved bytes and content types are logged; these responses are not counted as data.
- Europe PMC supplementary-file endpoints returned HTTP200 XML `errorBean`, not a ZIP archive. Full-text XML requests timed out. These are transport/content failures, not evidence that the articles lack data.
- Orr supplement alternatives returned HTTP403, HTTP404 or timed out; original 75-species numerical Table S1 remains unavailable in this bundle. Its known existence is preserved.
- An initial Hermida author-repository URL guessed document `/2/` and returned HTTP404; the repository page linked `/3/`, which successfully supplied the original PDF. The failed request is retained in the log.
- Jordan/Gubernator publisher pages expose subscription previews; original full numerical tables were not obtained. Their potential microbial-temperature evidence remains unresolved.
- The shared Python lacks bs4/lxml. The executable extraction uses the standard-library HTML parser and vendored pypdf instead; this had no effect on the retained numerical scope.

See `inputs/retrievals.jsonl` for exact URL/status/timestamp/hash evidence. None of these failures warrants a claim of biological absence or that Q10=2.2 is wrong.
