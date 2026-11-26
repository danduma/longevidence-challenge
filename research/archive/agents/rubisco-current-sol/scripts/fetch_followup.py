from pathlib import Path
import json,concurrent.futures
from fetch_sources import fetch,ROOT
urls={
 'orr2016.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5047088/',
 'hermida2016.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4972260/',
 'gallionella2026.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC13544237/',
 'editing2026.html':'https://www.nature.com/articles/s41467-026-73783-w',
 'kubien2008_anu.html':'https://openresearch-repository.anu.edu.au/handle/1885/58448',
 'cyano2023.pdf':'https://link.springer.com/content/pdf/10.1007/s11120-023-01008-y.pdf',
}
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,urls.items()))
with (ROOT/'raw/retrievals.jsonl').open('a') as f:
 for r in rows:f.write(json.dumps(r)+'\n');print(r['file'],r['http'],r['bytes'])
with (ROOT/'access_failures.md').open('a') as f:
 for r in rows:
  if r['http']!='200' or r['exit']:f.write('\n- '+json.dumps(r)+'\n')
