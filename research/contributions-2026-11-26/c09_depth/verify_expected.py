#!/usr/bin/env python3
"""C09 frozen result checks after independent regeneration."""
from pathlib import Path
import json, hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parent
expected=json.loads((ROOT/'sources/expected_scientific_hashes.json').read_text())
for name,digest in expected.items():
 assert hashlib.sha256((ROOT/'results'/name).read_bytes()).hexdigest()==digest, name
b=pd.read_csv(ROOT/'results/biochemical_trajectories.csv')
assert len(b)==1608 and not b.duplicated(['branch','TEVp','TEVs']).any()
assert b.vinirel.min()<0 # signed noisy rates must not be silently clipped
x=pd.read_csv(ROOT/'results/candidate_summary.csv')
z=x[(x.cutoff==500)&(x.selection=='within_motif_decile')].iloc[0]
assert z.pairs==95963 and z.baseline_selected_still_selected==3946
print(f'C09 PASS: {len(expected)} scientific output hashes reproduced; key and signed-rate checks passed.')
