#!/usr/bin/env python3
"""Read-only public literature retrieval; all artifacts stay in this directory."""
import json, urllib.request, urllib.parse, hashlib, datetime
from pathlib import Path
import xml.etree.ElementTree as ET
BASE=Path(__file__).resolve().parent
queries={
 'rao2011':'DOI:10.1093/nar/gkq976',
 'patel2021':'DOI:10.1261/rna.078686.121',
 'jayasinghe2024':'DOI:10.1261/rna.080156.124',
 'metal2025':'DOI:10.1021/acs.biochem.5c00162',
 'comparison2026':'DOI:10.1261/rna.080815.125',
 'recent_thg1':'(Thg1 OR "Thg1-like" OR "reverse RNA polymerase") AND FIRST_PDATE:[2024-01-01 TO 2026-11-26]',
}
log=json.loads((BASE/'sources/retrieval_log.json').read_text()) if (BASE/'sources/retrieval_log.json').exists() else []
def get(url,name):
 if (BASE/'sources'/name).exists():
  return (BASE/'sources'/name).read_bytes()
 try:
  data=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'IndependentLiteratureAudit/1.0'}),timeout=15).read()
  (BASE/'sources'/name).write_bytes(data)
  log.append({'url':url,'file':'sources/'+name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'status':'ok'})
  return data
 except Exception as ex:
  log.append({'url':url,'file':name,'status':'access failure','error':repr(ex)})
  (BASE/'sources/retrieval_log.json').write_text(json.dumps(log,indent=2)+'\n')
  return None
for name,q in queries.items():
 url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=100&query='+urllib.parse.quote(q)
 data=get(url,name+'_metadata.json')
 if not data: continue
 records=json.loads(data)['resultList']['result']
 print(name,len(records))
 for rec in records:
  if name=='recent_thg1' and not any(w in rec.get('title','').lower() for w in ['polymeras','guanylyltransferase']): continue
  print(rec.get('id'),rec.get('pmcid'),rec.get('firstPublicationDate'),rec.get('doi'),rec.get('title'))
  pmcid=rec.get('pmcid')
  if name!='recent_thg1' and pmcid:
   xml=get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmcid+'/fullTextXML',name+'_'+pmcid+'.xml')
   if xml:
    root=ET.fromstring(xml)
    paragraphs=[' '.join(''.join(n.itertext()).split()) for n in root.iter() if n.tag in ['article-title','p','table-wrap','title']]
    (BASE/'sources'/(name+'_'+pmcid+'.txt')).write_text('\n\n'.join(paragraphs))
(BASE/'sources/retrieval_log.json').write_text(json.dumps(log,indent=2)+'\n')
