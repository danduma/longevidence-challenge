#!/usr/bin/env python3
"""Optional primary-source refresh. Not required for offline numerical verification."""
import datetime,hashlib,json,urllib.request
from pathlib import Path
R=Path(__file__).resolve().parent
sources=[('model.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10431620/fullTextXML','CC BY 4.0'),('heart.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9164386/fullTextXML','CC BY 4.0'),('thermal2022.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8823202/fullTextXML','inspect article license; download only, do not redistribute without checking'),('vitrify_readme.txt','https://raw.githubusercontent.com/Marc-Dvci/Vitrify/1d498ee5c583e0112a762e60892d576ec2e27f57/README.md','Apache-2.0')]
records=[]
for name,url,permission in sources:
 rec=dict(name=name,url=url,permission=permission,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 try:
  with urllib.request.urlopen(url,timeout=40) as f:data=f.read();rec['http_status']=f.status
  (R/'inputs'/name).write_bytes(data);rec.update(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data),status='source retrieval succeeded')
 except Exception as e:rec.update(status='source access failure',error=str(e))
 records.append(rec)
(R/'inputs/retrievals.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(records,indent=2))
