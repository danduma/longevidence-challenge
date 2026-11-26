"""Download verified HTTPS wheels with system curl (Python TLS unavailable here)."""
import json,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path('.venv/lib').resolve()/('python'+'.'.join(map(str,sys.version_info[:2]))) / 'site-packages'))
from pip._vendor.packaging.tags import sys_tags
from pip._vendor.packaging.utils import parse_wheel_filename
priority={str(t):i for i,t in enumerate(sys_tags())}
Path('sources/wheels').mkdir(exist_ok=True)
for pkg in ['numpy','pandas','scipy','openpyxl','et-xmlfile','python-dateutil','six','tzdata','pypdf']:
 j=json.loads(subprocess.check_output(['curl','-fsSL','https://pypi.org/pypi/'+pkg+'/json']))
 candidates=[]
 for f in j['urls']:
  if f['filename'].endswith('.whl'):
   _,_,_,tags=parse_wheel_filename(f['filename']); ranks=[priority[str(t)] for t in tags if str(t) in priority]
   if ranks: candidates.append((min(ranks),f))
 if not candidates: raise RuntimeError(pkg)
 f=min(candidates,key=lambda x:x[0])[1]; p=Path('sources/wheels')/f['filename']
 if not p.exists(): subprocess.run(['curl','-fsSL',f['url'],'-o',str(p)],check=True)
 import hashlib
 assert hashlib.sha256(p.read_bytes()).hexdigest()==f['digests']['sha256']
 print(pkg,f['filename'],flush=True)
subprocess.run(['.venv/bin/pip','install','--no-index','--find-links=sources/wheels','numpy','pandas','scipy','openpyxl','pypdf'],check=True)
