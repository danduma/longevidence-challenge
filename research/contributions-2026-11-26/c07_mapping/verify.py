#!/usr/bin/env python3
"""C07 offline verifier. Run audit and compare frozen outputs/evidence, with shared-input guard."""
import hashlib,json,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parent
manifest=json.loads((R/'provenance_hashes.json').read_text())
checks=0
for rec in manifest['files']:
 p=R/rec['file'];assert hashlib.sha256(p.read_bytes()).hexdigest()==rec['sha256'],str(p);checks+=1
subprocess.run([sys.executable,'-B',str(R/'audit.py')],cwd=R,check=True)
for rec in manifest['files']:
 p=R/rec['file'];assert hashlib.sha256(p.read_bytes()).hexdigest()==rec['sha256'],str(p);checks+=1
for rec in json.loads((R/'inputs/manifest.json').read_text()):
 p=Path(rec['shared_path'])
 if p.exists():assert hashlib.sha256(p.read_bytes()).hexdigest()==rec['sha256'];checks+=1
s=json.loads((R/'results/summary.json').read_text())
assert s['new_replaced_links']==1 and s['new_nonfragment_quarantines']==2
assert s['baseline_wt_source_records']-s['additional_only_wt_source_records']==6
assert s['additional_only_wt_source_records']-s['after_wt_source_records']==12
assert s['repaired_wt_source_records']==1
print(f'C07 PASS: {checks} hash/shared-input checks; source, sequence, primary-assay, patch and exposure assertions passed. No model refit.')
