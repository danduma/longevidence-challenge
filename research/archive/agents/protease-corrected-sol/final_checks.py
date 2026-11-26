"""Assertions against primary deposited source data and complete provenance ledgers."""
import json,hashlib
from pathlib import Path
import pandas as pd,numpy as np
from reconstruct import ROOT,O,PAIR
v=json.loads((O/'source_validation.json').read_text())
assert v['S8e_final_singles']['common']==358058
assert v['S8e_final_singles']['max_AUC_error']==0
assert v['S8e_final_singles']['only_deposited']==v['S8e_final_singles']['only_reconstructed']==0
assert all(x['max_error']==0 for x in v['S7a_controls'].values())
assert v['S7c_manual_parent_branch']['discrepancies_over_1e10']==0
assert v['S7c_manual_parent_branch']['common_nonmissing']==836
for row in pd.read_csv(ROOT.parent/'scout-data-sol/data/file_inventory.csv').itertuples(index=False,name=None):
 p=ROOT.parent/'scout-data-sol/data'/Path(row[0]).name
 if p.name.startswith('protease_'):assert hashlib.sha256(p.read_bytes()).hexdigest()==row[2]
r=json.loads((O/'reconstruction.json').read_text())
assert r['script05_counts']==dict(all_singles_including_stops=358058,fig2_pairs=355097,fig2_enzymes=4437,fig2_substrates=134)
for name,nrows,npairs in [('DMS',1725713,532935),('controls',24103,3715)]:
 x=pd.read_csv(O/f'{name}_all_upstream_pair_ledger.tsv.gz',sep='\t')
 assert not x.duplicated(PAIR).any();assert len(x)==npairs
 assert x.raw_barcodes.sum()==nrows
 assert not x.reason.eq('UNEXPLAINED').any()
 assert x.reason.eq('retained').sum()==r[name]['final_pairs']
 for p in [O/f'{name}_corrected.tsv.gz',O/f'{name}_barcode_ledger.tsv.gz']:
  assert p.stat().st_size>0
m=pd.read_csv(O/'models.csv');x=m[(m.population=='original_raw')&(m.split=='blocked')&(m.reference=='all')].iloc[0]
assert abs(x.RMSE-.02906313509057658)<1e-12
assert abs(x.AP-.018764184159272827)<1e-12
k=pd.read_csv(O/'model_common_keys.csv.gz');assert len(k)==282646;assert not k.duplicated(PAIR).any()
rows=[]
for p in sorted(O.iterdir()):
 if p.is_file() and p.name not in ['output_hashes.tsv','final_checks.txt','run.log']:
  rows.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(rows).to_csv(O/'output_hashes.tsv',sep='\t',index=False)
(O/'final_checks.txt').write_text('PASS: exact S8e keys/AUCs, S7a controls, S7c manual-parent branch; headline counts; complete unique-key upstream exclusion ledgers; original review score replication; common-key count/uniqueness.\n')
print((O/'final_checks.txt').read_text())
