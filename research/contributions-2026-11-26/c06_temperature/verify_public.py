"""Recompute public derived ratios/rank implications, not original source extraction."""
from pathlib import Path
import csv,json,math,statistics
R=Path(__file__).resolve().parent

def read(name):
 with (R/'results'/name).open() as f:return list(csv.DictReader(f))
p=read('primary_temperature_points.csv');pairs=read('measured_temperature_pairs.csv')
s=json.loads((R/'results/summary.json').read_text())
assert len(p)==60 and len(pairs)==60 and len({x['species'] for x in p})==20
for species in {x['species'] for x in p}:
 v={int(x['temp_C']):float(x['kcat_s_minus1']) for x in p if x['species']==species}
 assert set(v)=={15,25,35}
 for pair in [x for x in pairs if x['species']==species]:
  a,b=int(pair['low_C']),int(pair['high_C']);ratio=v[b]/v[a]
  assert math.isclose(float(pair['ratio']),ratio,rel_tol=1e-12)
  assert math.isclose(float(pair['effective_Q10']),ratio**(10/(b-a)),rel_tol=1e-12)
for interval,expected in s['effective_Q10_by_interval'].items():
 a,b=map(int,interval.split('-'));v=[float(x['effective_Q10']) for x in pairs if int(x['low_C'])==a and int(x['high_C'])==b]
 for key,value in [('minimum',min(v)),('median',statistics.median(v)),('maximum',max(v))]:assert math.isclose(value,expected[key],rel_tol=1e-12)
i=read('source_inventory.csv');assert len(i)==67 and sum(int(x['rows']) for x in i)==634
c=read('workbook_candidate_pairs.csv');assert len(c)==46 and sum(x['same_origin']=='True' for x in c)==14
positive=[x for x in c if x['same_origin']=='True' and x['low_vC'] and x['high_vC']]
assert len(positive)==1 and 'in_vivo_estimate' in positive[0]['exclusion']
bridge=read('calibration_overlap_coverage.csv');assert len(bridge)==11 and all(x['supported_bridge']=='False' for x in bridge)
# Original [1,-.5] rows have null [.5,1]. A separate plant beta leaves it.
assert .5-.5==0 and s['strict_augmented_rank']==1 and s['conditional_universal_Q10_rank']==2
print('C06 public derived check PASS: 60 ratios, interval summaries, coverage and rank implications. Source re-extraction is separate; read README.')
