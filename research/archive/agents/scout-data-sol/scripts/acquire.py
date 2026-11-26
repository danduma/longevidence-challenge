"""Verified-TLS public downloads, bounded size, durable source provenance."""
import subprocess,json,hashlib,datetime
from pathlib import Path

def fetch(url,path,limit=30_000_000):
 p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
 r=subprocess.run(['curl','--fail','--location','--silent','--show-error','--max-time','90','--max-filesize',str(limit),'-w','%{http_code}\t%{content_type}\t%{url_effective}',url,'-o',str(p)],capture_output=True,text=True)
 row={'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'url':url,'path':str(p),'curl_exit':r.returncode,'http_metadata':r.stdout,'error':r.stderr.strip()}
 if r.returncode==0:
  row.update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 else:
  if p.exists(): p.unlink()
 with open('sources/download_manifest.jsonl','a') as f: f.write(json.dumps(row)+'\n')
 print(json.dumps(row),flush=True)
 return p if r.returncode==0 else None

if __name__=='__main__':
 import sys
 fetch(sys.argv[1],sys.argv[2])
