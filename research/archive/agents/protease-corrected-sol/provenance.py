from pathlib import Path
import hashlib,json,csv,importlib.metadata,platform,sys
from reconstruct import ROOT,DATA,O,S
rows=[]
for p in sorted(S.rglob('*')):
 if p.is_file() and '__pycache__' not in str(p) and p.name!='input_hashes.tsv':rows.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),role='local_primary_or_review_source'))
shared=[DATA/x for x in ['protease_DMS.zip','protease_controls.tsv','protease_vitro_relative.tsv','protease_vitro_absolute.tsv','file_inventory.csv']]
review=ROOT.parent/'review-protease-astra'
shared += [review/x for x in ['report.md','protocol.md','audit.py','sources/screen3.tsv']]
shared += sorted((review/'results').glob('*'))
for p in shared:rows.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),role='read_only_shared_input'))
for p in sorted(ROOT.glob('*.py'))+sorted(ROOT.glob('*.sh'))+[ROOT/'protocol.md']:
 rows.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),role='analysis_code_or_protocol'))
with (S/'input_hashes.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','sha256','role'],delimiter='\t');w.writeheader();w.writerows(rows)
packages={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','openpyxl']}
(ROOT/'requirements.txt').write_text('\n'.join(k+'=='+v for k,v in packages.items())+'\n')
(O/'environment.json').write_text(json.dumps(dict(python=sys.version,executable=sys.executable,platform=platform.platform(),packages=packages),indent=2)+'\n')
