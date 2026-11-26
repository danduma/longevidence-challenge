"""Read-only public retrieval; all outputs are confined to this job directory."""
from pathlib import Path
import subprocess,json,hashlib,datetime,concurrent.futures
ROOT=Path(__file__).resolve().parents[1]
URLS={
 'whitney2001.pdf':'https://ripe.illinois.edu/sites/ripe.illinois.edu/files/2018-06/Whitney_et_al-2001-The_Plant_Journal.pdf',
 'kubien2008.html':'https://academic.oup.com/jxb/article/59/7/1767/639023',
 'uemura1997.html':'https://www.sciencedirect.com/science/article/pii/S0006291X97964974',
 'uemura1997_crossref.json':'https://api.crossref.org/works/10.1006/bbrc.1997.6497',
 'uemura1997_openalex.json':'https://api.openalex.org/works/https://doi.org/10.1006/bbrc.1997.6497',
}
def fetch(item):
 name,url=item; dest=ROOT/'raw'/name
 r=subprocess.run(['curl','-L','--max-time','50','--silent','--show-error','-w','%{http_code}',url,'-o',str(dest)],capture_output=True,text=True)
 b=dest.read_bytes() if dest.exists() else b''
 out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'url':url,'file':str(dest.relative_to(ROOT)),'http':r.stdout,'stderr':r.stderr,'exit':r.returncode,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
 return out
if __name__ == "__main__":
 with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
  rows=list(ex.map(fetch,URLS.items()))
 with (ROOT/'raw'/'retrievals.jsonl').open('a') as f:
  for r in rows:f.write(json.dumps(r)+'\n');print(r['file'],r['http'],r['bytes'])
 with (ROOT/'access_failures.md').open('a') as f:
  for r in rows:
   if r['http']!='200' or r['exit']:f.write('\n- '+json.dumps(r)+'\n')
