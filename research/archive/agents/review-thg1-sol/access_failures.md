Observed access limitations during the 2026-11-26 review (not model refusals):

- Europe PMC fullTextXML: PMC3061083 and PMC8127993 timed out on initial and retry requests; PMC12369627 returned HTTP 500 and later timed out. Exact URL/error records are in sources/retrieval_log.json. Metadata searches worked; the 2024 fullTextXML retrieval succeeded. The user-supplied 2021 primary XML and table image were readable and used.
- Publisher/PMC browser access was inconsistent: multiple rnajournal.cshlp.org full/abstract/PDF URLs returned 403 or inaccessible errors; PMC3061083, PMC11404450, and PMC12369627 returned browser-check/CAPTCHA pages. The wrong guessed PMC3064775 was also a browser-check page and was not used as evidence.
- The 2011 OUP PDF initially returned searchable text, including its figure/results evidence through indexed search; subsequent requests redirected to inaccessible Silverchair URLs (500). One initially guessed OUP article path was inaccessible; the correct article is /nar/article/39/5/1833/2409527. Primary indexed results and the Figure 3 caption were usable.
- The 2026 journal full-text page required a subscription. Final abstract/metadata and the accessible 2025 preprint were inspected; the final full paper was not read. A bare PubMed open returned only an HHS footer, while indexed PubMed search returned the final abstract and bibliographic record.
- The 2025 ACS article direct open returned 403. Its indexed discussion, bibliographic metadata, and accessible primary preprint were usable.
- Alternative europepmc.org article-page browser opens were inaccessible. No paywall bypass, authentication changes, or contacts were attempted.
- A local `ps` command used to check retrieval progress was denied by the environment (`operation not permitted: ps`). Retrieval completion was checked via the existing command session instead.

A zero-result exact-title metadata search was a query limitation, not an access failure or evidence of absence: the 2026 paper was found by DOI. A missing retrieval log while the retrieval script was still running was not loss of source access. No source-access failure affected the independent arithmetic certificate.
