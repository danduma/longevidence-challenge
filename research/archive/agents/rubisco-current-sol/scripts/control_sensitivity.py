from pathlib import Path
import pandas as pd,numpy as np,json,itertools
R=Path(__file__).resolve().parents[1]
m=pd.read_csv(R/'derived/historical_measurements.csv');ed=pd.read_csv(R/'results/graph_edges.csv');rr=pd.read_csv(R/'results/control_rectangles.csv');rr=rr[rr['mode'].eq('condition')]
notes=m.set_index('measurement_id').note.fillna('')
flags={}
for _,r in ed[ed['mode'].eq('condition')].iterrows():
 ids=r.members.split(';');g=m[m.measurement_id.isin(ids)]
 flags[(r.trait,r.source,r.enzyme)]={'single':len(g)==1,'no_suspect':all('suspect' not in notes[i].lower() for i in ids),'native':all(str(v) in ['False','native'] for v in g.heterologous_expression),'known_pH':all(v!='unknown' for v in g.pH_compiled)}
rows=[]
for sensitivity in ['all_predeclared','exclude_suspect_notes','single_member_edges','native_only','equal_known_pH']:
 for trait in ['S','vC']:
  z=rr[rr.trait.eq(trait)].copy();keep=[]
  for _,r in z.iterrows():
   fs=[flags[(trait,source,enzyme)] for source in [r.s,r.t] for enzyme in [r.a,r.b]]
   b=True
   if sensitivity=='exclude_suspect_notes':b=all(f['no_suspect'] for f in fs)
   if sensitivity=='single_member_edges':b=all(f['single'] for f in fs)
   if sensitivity=='native_only':b=all(f['native'] for f in fs)
   if sensitivity=='equal_known_pH':b=all(f['known_pH'] for f in fs) and r.s.split(' | pH=')[-1]==r.t.split(' | pH=')[-1]
   keep.append(b)
  z=z[keep]
  if not len(z):continue
  sp=z.groupby(['s','t'])[['raw_squared_log_error','calibrated_squared_log_error']].mean()
  rows.append({'sensitivity':sensitivity,'trait':trait,'rectangles':len(z),'source_pairs':len(sp),'median_factor':z.factor.median(),'max_factor':z.factor.max(),'above1.1_n':int((z.factor>1.1).sum()),'above1.25_n':int((z.factor>1.25).sum()),'source_pair_balanced_raw_rms':np.sqrt(sp.raw_squared_log_error.mean()),'source_pair_balanced_calibrated_rms':np.sqrt(sp.calibrated_squared_log_error.mean())})
pd.DataFrame(rows).to_csv(R/'results/control_sensitivity.csv',index=False)
# Graph consequences of explicit source-origin corrections; broad species remains an assumption.
z=pd.read_csv(R/'derived/measurement_ledger.csv');summ=[];comps=[]
for regime in ['paper_label','endpoint_origin']:
 for trait in ['S','vC']:
  v=z[z.trait.eq(trait)].copy();v['source']=v.study if regime=='paper_label' else v.origin_source
  edges=v[['species_key','source']].drop_duplicates();adj={}
  for e,s in edges.itertuples(index=False):
   e='enzyme:'+e;s='source:'+s;adj.setdefault(e,set()).add(s);adj.setdefault(s,set()).add(e)
  seen=set();cc=[]
  for node in sorted(adj):
   if node in seen:continue
   c=set();q=[node]
   for n in q:
    if n in c:continue
    c.add(n);q+=list(adj[n]-c)
   cc.append(c);seen|=c
  for i,c in enumerate(cc):
   for node in sorted(c):comps.append({'regime':regime,'trait':trait,'component':i,'node':node})
  def connected(a,b):return any('enzyme:'+a in c and 'enzyme:'+b in c for c in cc)
  summ.append({'regime':regime,'trait':trait,'edges':len(edges),'nodes':len(adj),'components':len(cc),'cycle_rank':len(edges)-len(adj)+len(cc),'Galdieria_Griffithsia_connected':connected('Galdieria partita Tokara','Griffithsia monilis'),'Griffithsia_maize_connected':connected('Griffithsia monilis','Zea mays'),'Galdieria_maize_connected':connected('Galdieria partita Tokara','Zea mays')})
pd.DataFrame(summ).to_csv(R/'results/origin_graph_summary.csv',index=False);pd.DataFrame(comps).to_csv(R/'results/origin_graph_components.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(summ).to_string(index=False))
# Exact calibration/temperature confounding in the 2025 deposited overlap design.
a=pd.read_csv(R/'results/library2025_calibration_pairs.csv')
assert len(a)==11 and all(set(str(s).split('/'))=={'25.0'} for s in a.Tdeg_lit)
X=np.column_stack([np.ones(len(a)),np.repeat(.5,len(a))]); _,sv,_=np.linalg.svd(X)
c=float(a.deposited_factor.iloc[0]);q=2.2
s={'n_overlaps':len(a),'n_original_30C':0,'design_columns':['intercept','(30-T_literature)/10'],'rank':int(np.linalg.matrix_rank(X)),'columns':2,'singular_values':sv.tolist(),'author_Q10':q,'deposited_factor':c,'temperature_independent_identifiable_ratio_c_over_sqrtQ':c/np.sqrt(q),'limitation':'Only log(c)-0.5*log(Q10) is identified by this calibration table. This does not invalidate correction conditional on the stipulated Q10 or additional information outside the table. No rates converted here.'}
(R/'results/calibration_identifiability.json').write_text(json.dumps(s,indent=2))
print(json.dumps(s,indent=2))
