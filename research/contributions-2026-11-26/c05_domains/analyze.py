#!/usr/bin/env python3
"""C05 recompute domain attribution and decision concordance from local tables."""
import csv,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
O=ROOT/'results'
rows=list(csv.DictReader((O/'domains.csv').open()))
def same(a,b,keys):return all(a[k]==b[k] for k in keys)
def truth(s):return s=='True'
comparison=[]
for m in [r for r in rows if r['domain']=='myocardium']:
 w=next(r for r in rows if r['case']==m['case'] and r['domain']=='whole_specimen')
 comparison.append({k:m[k] for k in ['case','mm','mean','mode','design','dt_requested']}|{s+'_decision_differs':truth(m[s+'_pass'])!=truth(w[s+'_pass']) for s in ['all','window','glass']}|{s+'_myocardium_minus_whole':float(m[s+'_rate'])-float(w[s+'_rate']) for s in ['all','window','glass']})
with (O/'domain_comparison.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(comparison[0]));wr.writeheader();wr.writerows(comparison)
paired=list(csv.DictReader((O/'paired.csv').open())); rankings=[]
for m in [r for r in paired if r['domain']=='myocardium']:
 w=next(r for r in paired if same(r,m,['mm','mean','mode','dt_requested']) and r['domain']=='whole_specimen')
 rankings.append({k:m[k] for k in ['mm','mean','mode','dt_requested']}|{s+'_ranking_differs':np.sign(float(m[s+'_rate_flat_minus_original']))!=np.sign(float(w[s+'_rate_flat_minus_original'])) for s in ['all','window','glass']})
with (O/'domain_rankings.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(rankings[0]));wr.writeheader();wr.writerows(rankings)
summary=dict(contribution_id='C05',simulations=len(comparison),paired_comparisons=len(rankings),decision_differences={s:sum(r[s+'_decision_differs'] for r in comparison) for s in ['all','window','glass']},ranking_differences={s:int(sum(r[s+'_ranking_differs'] for r in rankings)) for s in ['all','window','glass']},discovery_decision_differences={s:sum(r[s+'_decision_differs'] for r in comparison if r['mm']=='0.75') for s in ['all','window','glass']})
attrs=[]
for r in rows:
 if r['domain']!='whole_specimen':continue
 attrs.append({k:r[k] for k in ['case','mm','mean','mode','design','dt_requested']}|{s+'_setting_label':json.loads(r[s+'_location'])['label'] for s in ['all','window','glass','cold','hot']}|{s+'_rate':float(r[s+'_rate']) for s in ['all','window','glass']})
with (O/'attribution.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(attrs[0]));wr.writeheader();wr.writerows(attrs)
summary['whole_setting_labels']={s:{str(l):sum(a[s+'_setting_label']==l for a in attrs) for l in [1,2,3]} for s in ['all','window','glass','cold','hot']}
summary['discovery_setting_labels']={s:{str(l):sum(a[s+'_setting_label']==l for a in attrs if a['mm']=='0.75') for l in [1,2,3]} for s in ['all','window','glass','cold','hot']}
# Isolate adding surrounding bath: compare heart contents (muscle + chamber) to whole.
bath=[]
for h in [r for r in rows if r['domain']=='heart_contents']:
 w=next(r for r in rows if r['case']==h['case'] and r['domain']=='whole_specimen')
 bath.append({'case':h['case']}|{s+'_decision_differs':truth(h[s+'_pass'])!=truth(w[s+'_pass']) for s in ['all','window','glass']})
with (O/'bath_addition.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(bath[0]));wr.writeheader();wr.writerows(bath)
summary['bath_addition_decision_differences']={s:sum(r[s+'_decision_differs'] for r in bath) for s in ['all','window','glass']}
summary['max_myocardium_minus_whole_all_rate']=max(r['all_myocardium_minus_whole'] for r in comparison)
summary['max_myocardium_minus_whole_window_rate']=max(r['window_myocardium_minus_whole'] for r in comparison)
(O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
# Convenient unrounded subset for report review; all values remain available in domains.csv.
for r in rows:
 if r['domain'] in ['myocardium','whole_specimen']:
  print(r['case'],r['domain'],*[f'{float(r[k]):.8f}' for k in ['all_rate','window_rate','glass_rate','cold','hot','window_transit_rate_min','window_transit_rate_max']],r['all_pass'],r['window_pass'],r['glass_pass'],r['glass_missing_exit'])
