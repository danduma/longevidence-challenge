"""Read-only public source retrieval; archive bytes, dates, failures and versions locally."""
from pathlib import Path
import urllib.request, json, hashlib, datetime, xml.etree.ElementTree as ET
OUT=Path('sources');OUT.mkdir(exist_ok=True)
urls={
 'heart.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9164386/fullTextXML',
 'jr23.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10431620/fullTextXML',
 'kidney2021.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8498880/fullTextXML',
 'heart2022.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8823202/fullTextXML',
 'liter2025-correct.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12474906/fullTextXML',
 'github-repo.json':'https://api.github.com/repos/Marc-Dvci/Vitrify',
 'github-commits.json':'https://api.github.com/repos/Marc-Dvci/Vitrify/commits?until=2026-11-26T23:59:59Z&per_page=100',
 'github-issues.json':'https://api.github.com/repos/Marc-Dvci/Vitrify/issues?state=all&per_page=100',
 'github-branches.json':'https://api.github.com/repos/Marc-Dvci/Vitrify/branches?per_page=100',
 'github-tags.json':'https://api.github.com/repos/Marc-Dvci/Vitrify/tags?per_page=100',
 'pig-v1.html':'https://www.biorxiv.org/content/10.64898/2026.10.07.756990v1.full',
 'pig-api.json':'https://api.biorxiv.org/details/biorxiv/10.64898/2026.10.07.756990',
 'dsc-crossref.json':'https://api.crossref.org/works/10.1016/j.cryobiol.2024.104977',
}
manifest=[]
for name,url in urls.items():
 record=dict(name=name,url=url,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'IndependentComputationalReview/1.0','Accept':'application/json' if name.endswith('.json') else '*/*'})
  with urllib.request.urlopen(req,timeout=45) as response:
   content=response.read();record.update(status=response.status,headers=dict(response.headers),final_url=response.url)
  (OUT/name).write_bytes(content);record.update(bytes=len(content),sha256=hashlib.sha256(content).hexdigest())
  if name.endswith('.xml'):
   root=ET.fromstring(content)
   paragraphs=[' '.join(e.itertext()) for e in root.iter() if e.tag in ('p','table','article-title','caption')]
   (OUT/(name+'.txt')).write_text('\n\n'.join(paragraphs))
  print(name,record['status'],len(content),flush=True)
 except Exception as e:
  record['error']=repr(e);print(name,record['error'],flush=True)
 manifest.append(record)
 (OUT/'retrieval.json').write_text(json.dumps(manifest,indent=2)+'\n')
