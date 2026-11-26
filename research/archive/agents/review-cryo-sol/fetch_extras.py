from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request,json,hashlib,datetime,xml.etree.ElementTree as E
urls={
 'kidney2021-retry.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8498880/fullTextXML',
 'heart2022-retry.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8823202/fullTextXML',
 'liter2025-correct.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12474906/fullTextXML',
 'dsc.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11611618/fullTextXML',
 'icekinetics.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9248894/fullTextXML',
 'skow2026.pdf':'https://truthinresearch.pw/index.php/TIR/article/download/83/66/269',
 'heart_fig1.jpg':'https://s3.amazonaws.com/pmc-oa-opendata/PMC9164386.1/nihms-1745591-f0002.jpg',
}
def fetch(item):
 name,url=item;record=dict(name=name,url=url,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 try:
  with urllib.request.urlopen(url,timeout=30) as res:data=res.read();record.update(status=res.status,headers=dict(res.headers))
  Path('sources',name).write_bytes(data);record.update(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
  if name.endswith('.xml'):
   root=E.fromstring(data);record['title']=''.join(root.find('.//article-title').itertext())
   Path('sources',name+'.txt').write_text('\n\n'.join(''.join(p.itertext()) for p in root.iter('p')))
 except Exception as e:record['error']=repr(e)
 print(json.dumps(record),flush=True);return record
with ThreadPoolExecutor(max_workers=7) as pool:records=list(pool.map(fetch,urls.items()))
Path('sources/retrieval-extra.json').write_text(json.dumps(records,indent=2)+'\n')
