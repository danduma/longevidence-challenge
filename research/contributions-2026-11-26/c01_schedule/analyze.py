#!/usr/bin/env python3
"""C01 deterministic tables and sampled conditional nondominance."""
import sys
sys.dont_write_bytecode=True
import csv,json
from pathlib import Path
R=Path(__file__).resolve().parent
DOMAINS=['specimen','myocardium'];WINDOWS=['all','glass','glass35','window']
def writecsv(name,rows):
 if not rows:return
 with (R/'results'/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def dominates(a,b):return all(x>=y-.02 for x,y in zip(a,b)) and any(x>y+.02 for x,y in zip(a,b))
def main():
 rows=json.loads((R/'results/all_runs.json').read_text()); flat=[];minima=[]
 for r in rows:
  for d in DOMAINS:
   v=r[d];base={k:r[k] for k in ['id','phase','mm','mean','mode','schedule']};base['domain']=d
   flat.append(dict(base,**{f'{w}_rate':v[w]['rate'] for w in WINDOWS},cold=v['cold'],hot=v['hot'],terminal_pass=v['terminal_pass'],all_pass=v['terminal_pass'] and v['all']['rate']>=50,glass_pass=v['terminal_pass'] and v['glass']['rate']>=50,source_J=r['source_J'],enthalpy_J=r['enthalpy_J'],duty_s=r['duty_s'],energy_rel_error=r['energy_rel_error'],missing_100_40=v['missing_100_40'],transit_100_40=v['transit_100_40']))
   for w in WINDOWS:minima.append(dict(base,window=w,**{k:v[w][k] for k in ['rate','time','temp','label','cell','crossing_123_time']},ijk=json.dumps(v[w]['ijk']),xyz_mm=json.dumps(v[w]['xyz_mm'])))
 writecsv('summary.csv',flat);writecsv('minima.csv',minima)
 pareto=[];ranking=[]
 for phase in ['discovery','spatial','time']:
  for mean in [1.47,1.83]:
   for mode in ['clamped','extrapolated']:
    rs=[r for r in rows if r['phase']==phase and r['mean']==mean and r['mode']==mode]
    if not rs:continue
    for d in DOMAINS:
     for eligible in [False,True]:
      pool=[r for r in rs if not eligible or r[d]['terminal_pass']]
      for r in pool:
       pair=lambda x:[x[d][w]['rate'] for w in ['all','glass']]
       better=[s['schedule'] for s in pool if dominates(pair(s),pair(r))]
       pareto.append(dict(id='C01',phase=phase,mean=mean,mode=mode,domain=d,terminal_required=eligible,schedule=r['schedule'],all_rate=r[d]['all']['rate'],glass_rate=r[d]['glass']['rate'],nondominated=not better,dominated_by=';'.join(better)))
      if pool:
       bestall=max(r[d]['all']['rate'] for r in pool);bestglass=max(r[d]['glass']['rate'] for r in pool)
       for r in pool:
        if bestall-r[d]['all']['rate']<=.02:
         ranking.append(dict(id='C01',phase=phase,mean=mean,mode=mode,domain=d,terminal_required=eligible,all_best_tie=r['schedule'],all_rate=r[d]['all']['rate'],all_time=r[d]['all']['time'],all_temp=r[d]['all']['temp'],glass_rate=r[d]['glass']['rate'],glass_regret=bestglass-r[d]['glass']['rate'],best_glass=bestglass,best_glass_schedules=';'.join(s['schedule'] for s in pool if bestglass-s[d]['glass']['rate']<=.02)))
 writecsv('nondominance.csv',pareto);writecsv('objective_choices.csv',ranking)
 checks=[]
 for r in rows:
  if r['phase'] not in ['spatial','time']:continue
  prior='discovery' if r['phase']=='spatial' else 'spatial'
  ref=next(s for s in rows if s['phase']==prior and all(s[k]==r[k] for k in ['mean','mode','schedule']))
  for d in DOMAINS:
   checks.append(dict(id='C01',comparison='space' if prior=='discovery' else 'time',mean=r['mean'],mode=r['mode'],schedule=r['schedule'],domain=d,**{f'delta_{w}':r[d][w]['rate']-ref[d][w]['rate'] for w in WINDOWS},delta_cold=r[d]['cold']-ref[d]['cold'],delta_hot=r[d]['hot']-ref[d]['hot'],terminal_same=r[d]['terminal_pass']==ref[d]['terminal_pass']))
 writecsv('refinement.csv',checks)
 secondary=[]
 for phase in ['discovery','spatial','time']:
  for mean in [1.47,1.83]:
   for mode in ['clamped','extrapolated']:
    rs=[r for r in rows if r['phase']==phase and r['mean']==mean and r['mode']==mode]
    for d in DOMAINS:
     for support in ['glass35','window']:
      for required in [False,True]:
       pool=[r for r in rs if not required or r[d]['terminal_pass']]
       for r in pool:
        better=[s['schedule'] for s in pool if dominates([s[d]['all']['rate'],s[d][support]['rate']],[r[d]['all']['rate'],r[d][support]['rate']])]
        secondary.append(dict(id='C01',phase=phase,mean=mean,mode=mode,domain=d,rate_support=support,terminal_required=required,schedule=r['schedule'],nondominated=not better,dominated_by=';'.join(better)))
 writecsv('secondary_nondominance.csv',secondary)
 print('C01 tables written:',len(flat),'domain rows, not independent replicates')
if __name__=='__main__':main()
