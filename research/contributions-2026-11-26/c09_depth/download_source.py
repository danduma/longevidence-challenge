#!/usr/bin/env python3
"""Optional public primary-source check. Offline C09 verification uses inputs/."""
from pathlib import Path
import urllib.request,urllib.parse,json,hashlib,datetime
ROOT=Path(__file__).resolve().parent
commit='8b9a1865108fa44c3e0a24ac20a59754222847ca'
p='2 Data analysis (Lukas Huber)/In vitro data by Luzius Pestalozzi/data_vinirel_Luzi.txt'
u='https://raw.githubusercontent.com/JeschekLab/ProtRec/'+commit+'/'+urllib.parse.quote(p)
record={'id':'C09','url':u,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'commit':commit}
try:
 with urllib.request.urlopen(u,timeout=30) as r:b=r.read();record.update(status=r.status)
 record['sha256']=hashlib.sha256(b).hexdigest();record['git_blob']=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
 assert record['git_blob']=='17a2b0225edafcda1fed7463161b4d8bb70389df'
 assert b==(ROOT/'inputs/protease_vitro_relative.tsv').read_bytes()
 record['local_snapshot_exact']=True
except Exception as e:record['source_failure']=str(e)
(ROOT/'sources/live_download_check.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
