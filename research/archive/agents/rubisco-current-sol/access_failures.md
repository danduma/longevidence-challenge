# Access/retrieval failures

- Initial web DOI open: 10.1006/bbrc.1997.6497 returned redirect/safety retrieval error (500 in web tool).
- Wiley full-text 10.1046/j.1365-313x.2001.01056.x returned HTTP 403 in web tool.
- Initially guessed Oxford URL /jxb/article/59/7/1767/642023 unavailable; correct URL to be resolved from primary DOI.
- Shared scientific Python lacks beautifulsoup4, pypdf, requests and networkx. This is an environment limitation, not a model refusal; use stdlib HTML/XML/urllib and available packages or local dependencies as needed.

- {"utc": "2026-11-26T12:20:21.938215+00:00", "url": "https://academic.oup.com/jxb/article/59/7/1767/639023", "file": "raw/kubien2008.html", "http": "403", "stderr": "", "exit": 0, "bytes": 5437, "sha256": "aed23a9cb5917252e14ef8a505d6a38e1669cc9547fdd4759799af1ed26bac3c"}

- {"utc": "2026-11-26T12:20:22.237669+00:00", "url": "https://www.sciencedirect.com/science/article/pii/S0006291X97964974", "file": "raw/uemura1997.html", "http": "403", "stderr": "", "exit": 0, "bytes": 832981, "sha256": "fc7c694f766582bb10f1c76afb5cc6cd969bf1557bc228327132983ec4818e19"}

- {"utc": "2026-11-26T12:23:23.591275+00:00", "url": "https://academic.oup.com/jxb/article/59/7/1767/639023", "file": "raw/kubien2008.html", "http": "403", "stderr": "", "exit": 0, "bytes": 5437, "sha256": "a7202771d324ec3b4ac2367cc867799363c4569560ea81b9765d870526e05e30"}

- {"utc": "2026-11-26T12:23:23.853783+00:00", "url": "https://www.sciencedirect.com/science/article/pii/S0006291X97964974", "file": "raw/uemura1997.html", "http": "403", "stderr": "", "exit": 0, "bytes": 832981, "sha256": "6d066bd48c749a70eae740aa768d887dadd671b58fbd74b40980703ec7e67899"}

- {"utc": "2026-11-26T12:23:25.285170+00:00", "url": "https://www.sciencedirect.com/science/article/pii/S0006291X97964974/pdf", "file": "raw/uemura1997_publisher.pdf", "http": "403", "stderr": "", "exit": 0, "bytes": 832981, "sha256": "8b6d4660c4e5e9705449732088292b06b20f84133945d6799aa36eab387f6c8d"}

- {"utc": "2026-11-26T12:23:25.388225+00:00", "url": "https://jxb.oxfordjournals.org/content/59/7/1767.full.pdf", "file": "raw/kubien2008_alternate.pdf", "http": "403", "stderr": "", "exit": 0, "bytes": 93, "sha256": "0d3e98ca727fc1201b436170af5a63f23348aaf146a3ac6234f6c4da283e8b34"}

- {"utc": "2026-11-26T12:23:25.007409+00:00", "url": "https://academic.oup.com/pcp/article/37/3/325/1930998", "file": "raw/uemura1996.html", "http": "403", "stderr": "", "exit": 0, "bytes": 5437, "sha256": "d37b233cfd09ec35de75ab62bb670ba2e92e5ca74901280ba30950676e33ff1a"}

- {"utc": "2026-11-26T12:24:15.048681+00:00", "url": "https://www.publish.csiro.au/FP/PP9940449", "file": "raw/kane1994.html", "http": "000", "stderr": "curl: (28) Connection timed out after 50030 milliseconds\n", "exit": 28, "bytes": 0, "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"}

- {"utc": "2026-11-26T12:27:49.680725+00:00", "url": "https://openresearch-repository.anu.edu.au/server/api/core/bitstreams/f5135dd3-7a7b-4a69-94c9-e74f0fe7d324/content", "file": "raw/kubien2008_primary.pdf", "http": "401", "stderr": "", "exit": 0, "bytes": 202, "sha256": "de899d73d4eea2dc80d767cbee614e4d68460c2a170fc35814be464b2253337c"}

- {"utc": "2026-11-26T12:27:48.435791+00:00", "url": "https://static-content.springer-cdn.com/esm/art%3A10.1007%2Fs11120-023-01008-y/MediaObjects/11120_2023_1008_MOESM1_ESM.pdf", "file": "raw/cyano2023_supp.pdf", "http": "000", "stderr": "curl: (6) Could not resolve host: static-content.springer-cdn.com\n", "exit": 6, "bytes": 0, "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"}

- Shared Python also lacks pip; pypdf7.1.1 was extracted from the already downloaded wheel into this job's vendor directory. No shared environment was modified.
- Elsevier API for Uemura1997 returned HTTP200 with bibliographic metadata only, not numerical tables/full text. This was a content-access failure despite successful HTTP.
- Oxford's cached/indexed Kubien2008 PDF text was readable through web search; direct token download returned HTTP400/session expired. ANU's public metadata resolved the original PDF but content returned HTTP401. No restricted access was bypassed.
- PMC find on Orr2016 later returned a browser/reCAPTCHA page; the earlier successfully retrieved primary HTML is retained.
- The initial Matplotlib run could not write its default home cache and automatically used the permitted temporary area. Plotting now explicitly sets MPLCONFIGDIR inside this job directory; no analysis inputs were written outside the job.
