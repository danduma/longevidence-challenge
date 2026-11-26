"""Run reviewed author's thermal implementation locally; no writes/import cache upstream.
This is reproduction, separate from independent.cpp and run_review.py.
"""
import sys,json
import numpy as np
from pathlib import Path
sys.dont_write_bytecode=True
TARGET=Path(__file__).resolve().parent.parent/'cryo-robust-astra'
code=(TARGET/'investigate.py').read_text().replace("U = Path('upstream')",'U = Path('+repr(str(TARGET/'upstream'))+')')
scope={'__name__':'reviewed_original'}
exec(compile(code,str(TARGET/'investigate.py'),'exec'),scope)
run,read=scope['run'],scope['read']
records=[]
for mode,values in read('thermal_validation')['readings'].items():
 for case,expected in values.items():
  actual=run(.5,'case'+case,sar=mode,duration=30)['glass_rate']
  records.append(dict(case=case,mode=mode,mm=.5,metric='glass_rate',actual=actual,expected=expected['mean_rate_c_per_min']))
for file,design in [('grid_convergence','original'),('grid_convergence_margin','margin')]:
 for level in read(file)['levels']:
  for key,d in [('baseline','caseII'),('optimized',design)]:
   actual=run(level['grid_spacing_mm'],d)
   for metric,field in [('min_rate','min_warming_rate_c_per_min'),('cold','final_min_temperature_c'),('hot','final_max_temperature_c'),('mean_rate','mean_warming_rate_c_per_min')]:
    records.append(dict(case=d,mode=key,mm=level['grid_spacing_mm'],metric=metric,actual=actual[metric],expected=level['designs'][key][field]))
for record in records:record['error']=record['actual']-record['expected']
Path('results/original_baselines.json').write_text(json.dumps(records,indent=2)+'\n')
print('baseline_metrics',len(records),'max_abs_error',max(abs(x['error']) for x in records),flush=True)
rows=[]
for d,s,mode in [('margin','nominal','extrapolated'),('equal','nominal','extrapolated'),('margin','high','clamped'),('equal','high','clamped')]:
 result=run(.5,d,s,mode);rows.append(result);print(json.dumps(result,default=lambda x:x.item()),flush=True)
Path('results/original_comparison.json').write_text(json.dumps(rows,indent=2,default=lambda x:x.item())+'\n')
