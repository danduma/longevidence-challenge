**Review date: 26 November 2026.** These are observed access limitations, separate from scientific findings.

- PMC browser access to PMC9164386 and subsequent PMC8823202/PMC9248894 requests returned a reCAPTCHA page. Fresh Europe PMC heart XML and the original heart figure/supplement succeeded; the earlier 2022 cardiac browser text was accessible and archived.
- Direct Europe PMC XML requests for PMC8498880 and PMC8823202 initially timed out. Kidney retry succeeded; cardiac retry returned HTTP 500. Correct liter-scale PMC12474906, DSC PMC11611618, and ice-kinetics PMC9248894 XML attempts timed out. The relevant primary indexed/browser text supplied the stated evidence; no claim of successful full-text XML retrieval is made for these.
- Browser GitHub REST opens returned an internal inaccessible-URL error. Connected GitHub GET tools and direct Python HTTP requests succeeded for the history, issues and branch endpoints. Public repository audit was completed.
- Direct browser bioRxiv v1 open returned HTTP 403; indexed primary text, direct Python HTML retrieval and the bioRxiv metadata API succeeded.
- Direct browser DOI access to the tissue calorimetry paper failed after an Elsevier redirect. Its primary PMC Table 1 was accessible. The ferumoxytol paper was inspected only through its primary publisher abstract.
- Nature's liter-scale article browser open failed on its identity-service redirect. Correct PMC indexed primary text and figure definition were retrieved instead.
- A guessed Skow PDF URL failed; following the actual article download link succeeded. Direct Python downloaded the PDF and its sections were inspected.
- The heart supplement browser open was inaccessible, but direct S3 download succeeded. Local PDFKit text extraction succeeded and was inspected.
- An initial guessed liter-scale PMC identifier, PMC12326946, returned HTTP 200 for an unrelated diabetes paper. Title verification caught it; it is excluded from evidence. Correct identifier PMC12474906 is documented in sources.md.
- Shell process inspection (`ps`) was denied with “operation not permitted.” It was unnecessary for numerical work.
- Environment inspection found no JAX or pytest. The original JAX/Tesseract/FEniCSx container stack was not run. Original-code reproduction here means compiling the reviewed author's C++ transcription and comparing its outputs to pinned baseline JSONs; the independent calculation uses separate reviewer code. No installation was attempted.

Exact HTTP URLs, errors, successful fallbacks, retrieval dates and source hashes are retained in sources/retrieval*.json, sources/supplement-retrieval.json, sources/primary*.json and fetch*.log. A reviewer harness JSON serialization error for NumPy int32 was corrected with `.item()` and the reproduction rerun successfully; that was a local test-harness error, not a failure in the contribution.
