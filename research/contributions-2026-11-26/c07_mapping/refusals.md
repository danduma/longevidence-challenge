# C07 refusals and source-access failures

No model refusal occurred. No subagents, other LLM calls, contacts, posts, commits, deployments or laboratory work were performed. This report is not human expert reviewed.

Source-access and tooling failures, 2026-11-26:

- The shared scientific Python lacks `requests`; the initial fetch failed with ModuleNotFoundError. Fetching was changed to standard-library urllib. This is a dependency failure, not a refusal.
- Web-reader PMC pages returned browser/recaptcha interstitials; the OUP reader route failed with a redirected-page error. Shell requests to the PMC `?pdf=1` HTML routes succeeded for the crop and palustris papers and were parsed for identity evidence. The suffix requests HTML, not a PDF conversion.
- NCBI BioC requests returned an HTTP429 or HTTP200 HTML error saying no result. Those HTTP200 files are semantic failures, not primary-text evidence. Europe PMC fullTextXML returned HTTP500. OUP Kubien HTML returned HTTP403. The PubMed-indexed primary abstract supports only the described separation of assay labels; full specimen identity remains unresolved.
- PubMed requests for Laing and Walker returned HTTP203 and were not saved by the downloader. Laing PMC accessible HTML provides an abstract/scanned-paper entry; cultivar identity was not independently recovered from full text. Wayne is explicitly unresolved, not a proven sequence mismatch.
- The UniProt license/help and API documentation shell responses were JavaScript fallback shells. The live indexed primary UniProt API documentation explicitly identifies CC BY 4.0; citation is in sources.md. A guessed documentation GitHub path returned404.
- `pymupdf` and `pdftotext` were absent; the downloaded crop PDF was not used for scored evidence. Accessible original-paper HTML supplied the relevant identity evidence.
- All 58 selected-accession UniProt requests succeeded. 34 were Inactive/deleted records, which is database status, not access failure or evidence of wrong biology. All 34 last-version histories were subsequently retrieved (68 successful requests). Historic sequences were checked byte-for-byte against the release.

Machine retrieval logs preserve actual HTTP outcomes and URLs; successful HTTP status alone was never treated as evidence of valid assay text. No credentials or private session logs are included.
