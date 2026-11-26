#!/usr/bin/env python3
"""C01 offline result, protocol, analytic and refinement verification."""
import sys
sys.dont_write_bytecode=True
import hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parent

def main():
 assert hashlib.sha256((R/'protocol.md').read_bytes()).hexdigest()==(R/'inputs/protocol.sha256').read_text().split()[0]
 if (R/'provenance.json').exists():
  provenance=json.loads((R/'provenance.json').read_text())
  for item in provenance['files']:
   assert hashlib.sha256((R/item['path']).read_bytes()).hexdigest()==item['sha256'], item['path']
 baselines=json.loads((R/'results/baseline_verified.json').read_text())
 assert abs(baselines[0]['min_rate']-44.4358)<.03 and abs(baselines[1]['min_rate']-50.3155)<.03
 assert all(abs(b['energy_rel_error'])<1e-5 and b['missing']==0 for b in baselines)
 rows=json.loads((R/'results/all_runs.json').read_text())
 assert len(rows)==66 and len({(r['phase'],r['mean'],r['mode'],r['schedule']) for r in rows})==66
 assert sum(r['phase']=='discovery' for r in rows)==36
 target=66.97622645671459
 assert all(abs(r['duty_s']-target)<1e-10 for r in rows)
 assert all(abs(r['energy_rel_error'])<1e-5 and r['boundary_J']==0 for r in rows)
 assert all(r[d]['missing_100_40']==0 for r in rows for d in ['specimen','myocardium'])
 assert all(r[d]['all']['rate']>0 for r in rows for d in ['specimen','myocardium'])
 assert all(r[d]['missing_crossings'][4]>0 for r in rows for d in ['specimen','myocardium'])
 for r in rows:
  for d in ['specimen','myocardium']:
   v=r[d]
   assert v['cold']<=v['hot']<0
   for w,lo,hi in [('glass',-123,0),('glass35',-123,-35),('window',-100,-40)]:assert lo-1e-9<=v[w]['temp']<=hi+1e-9
  # Every initial myocardium has this rate; a sampled trajectory minimum cannot exceed it.
  assert r['myocardium']['all']['rate']<=r['analytic_myocardium_onset_rate']+1e-8
 max_time=0;max_temp=0
 for r in rows:
  if r['phase']!='time':continue
  ref=next(s for s in rows if s['phase']=='spatial' and all(s[k]==r[k] for k in ['mean','mode','schedule']))
  for d in ['specimen','myocardium']:
   max_time=max(max_time,*[abs(r[d][w]['rate']-ref[d][w]['rate']) for w in ['all','glass','glass35','window']])
   max_temp=max(max_temp,abs(r[d]['cold']-ref[d]['cold']),abs(r[d]['hot']-ref[d]['hot']))
 assert max_time<.02 and max_temp<.01
 # Report the finite primary null and the narrower-window reversal without asserting global optimality.
 for mean in [1.47,1.83]:
  for mode in ['clamped','extrapolated']:
   group=[r for r in rows if r['phase']=='discovery' and r['mean']==mean and r['mode']==mode]
   mid=next(r for r in group if r['schedule']=='deficit_mid')
   assert all(mid['specimen']['all']['rate']>=s['specimen']['all']['rate']-.02 and mid['specimen']['glass']['rate']>=s['specimen']['glass']['rate']-.02 for s in group)
 def get(n):return next(r for r in rows if r['phase']=='time' and r['mean']==1.83 and r['mode']=='clamped' and r['schedule']==n)
 mid=get('deficit_mid');orig=get('original')
 assert mid['specimen']['all']['rate']>orig['specimen']['all']['rate']+1
 assert mid['specimen']['glass']['rate']>orig['specimen']['glass']['rate']+.1
 assert mid['specimen']['glass35']['rate']<orig['specimen']['glass35']['rate']-.1
 assert mid['specimen']['window']['rate']<orig['specimen']['window']['rate']-.1
 assert mid['specimen']['all']['temp']< -123 and mid['specimen']['terminal_pass']
 proof=dict(id='C01',runs=len(rows),continuous_command_duty_s=rows[0]['duty_s'],max_relative_energy_imbalance=max(abs(r['energy_rel_error']) for r in rows),max_time_refinement_rate_change_C_min=max_time,max_time_refinement_endpoint_change_C=max_temp,primary_family_result='deficit_mid weakly dominates entire nine-schedule discovery family within 0.02 C/min in all four whole-specimen scenarios; nominal clamped discovery has no terminal-eligible schedule',secondary_result='narrower-window tradeoffs persist in time-refined high-loading clamped comparison',status='agent-produced; awaiting coordinator verification')
 (R/'results/verification.json').write_text(json.dumps(proof,indent=2)+'\n')
 print(json.dumps(proof,indent=2))
if __name__=='__main__':main()
