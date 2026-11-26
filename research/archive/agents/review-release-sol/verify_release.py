#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
import pandas as pd
S=Path('/path/to/longevidence-checkout/research/millennium/sources/rubisco-active-learning')
Z=Path('/path/to/longevidence-checkout/research/millennium/sources/rubisco-active-learning-Version3.zip')
tree=json.load(open('web_evidence/pin_tree.json'));assert not tree.get('truncated');lookup={v['path']:v['sha'] for v in tree['tree'] if v['type']=='blob'}
rows=[]
for r in pd.read_csv('outputs/release_manifest.csv').itertuples():
 b=(S/r.path).read_bytes();gitsha=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest();match=gitsha==lookup.get(r.path)
 rows.append(dict(path=r.path,local_git_blob=gitsha,pinned_git_blob=lookup.get(r.path),matches_pin=match));assert match,r.path
pd.DataFrame(rows).to_csv('outputs/pinned_git_verification.csv',index=False)
z=json.load(open('web_evidence/zenodo.json'));md5=hashlib.md5(Z.read_bytes()).hexdigest();assert 'md5:'+md5==z['files'][0]['checksum']
Path('outputs/archive_verification.json').write_text(json.dumps(dict(archive_md5=md5,zenodo_checksum=z['files'][0]['checksum'],matches=True,archive_sha256=hashlib.sha256(Z.read_bytes()).hexdigest(),archive_size=Z.stat().st_size),indent=2))
print(f'All {len(rows)} provided files match pinned Git blobs. Archive MD5 matches Zenodo Version3 record.')
