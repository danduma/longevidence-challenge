2026-11-26 access limitations, separate from model refusals:

- Python urllib: MDPI article and version-notes URLs returned HTTP 403 Forbidden. Both pages were accessible through the web tool; substantive article evidence was available there.
- Web tool: initial MDPI /xml open showed an automated robot check rather than article XML. Main HTML was accessible; some XML content appeared in indexed search results.
- Web tool: Zenodo record HTML returned an inaccessible/internal-error response. Zenodo record API succeeded, allowing checksum verification of the supplied archive.
- Web tool: Journal of Experimental Botany article 72/22/7846 redirected to an inaccessible silverchair CDN (reported 500). Relevant indexed excerpts were accessible; full-text claim is not made.
- Web tool: PMC8233502 returned a browser/reCAPTCHA check. PubMed and an Oxford-hosted primary manuscript were available through search.

Machine-readable direct HTTP attempts: web_evidence/retrieval_log.json; command log: retrieval_run.log. No failure was interpreted as a model refusal.
