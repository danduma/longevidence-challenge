"""Frozen-protocol source/control portability benchmark; no network or mutations to shared inputs."""
from pathlib import Path
import itertools,json,hashlib,datetime,platform,sys,importlib.metadata
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
SHARED=ROOT/'raw/shared_inputs'
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
def write(x,name):pd.DataFrame(x).to_csv(OUT/name,index=False)
def norm(s):return ' '.join(str(s).split())
def ph(s):
 if pd.isna(s):return 'unknown'
 try:return '%g'%float(s)
 except (TypeError,ValueError):return norm(s)
inputs=[SHARED/'data/rubisco2019_original.xlsx',SHARED/'sources/rubisco_cyano2026_pmc.html',SHARED/'data/rubisco2025_library.csv',SHARED/'data/rubisco2025_supp3.csv',SHARED/'results/rubisco_bibliographic_overlay.csv']
manifest=[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in inputs]
write(manifest,'input_manifest.csv')
x=pd.read_excel(inputs[0]);x.insert(0,'excel_row',np.arange(len(x))+2)
mask=x['use'].eq(True)&x.primary.eq(True)&x.mutant.eq(False)&x.temp_C.eq(25)
x['accepted']=mask;x['species_key']=x.species.map(norm);x['pH_key']=x.pH.map(ph)
x.to_csv(ROOT/'derived/workbook_with_selection.csv',index=False)
x=x[mask].copy();x['study']=x.short_ref.map(norm)
# Bookkeeping preserves all measurements and flags; SD fields remain compiled error values.
measurement=[]
for _,r in x.iterrows():
 for trait in ['S','vC']:
  if pd.notna(r[trait]) and r[trait]>0:
   measurement.append({'measurement_id':f"2019:row{r.excel_row}:{trait}",'excel_row':r.excel_row,'identifier':r.identifier,'species_raw':r.species,'species_key':r.species_key,'study':r.study,'trait':trait,'value':r[trait],'unit':'mol/mol dissolved CO2/O2 convention' if trait=='S' else 's^-1 per catalytic site','compiled_error':r[trait+'_SD'],'error_kind':'unverified compilation _SD field','temp_C':r.temp_C,'pH_compiled':r.pH_key,'heterologous_expression':r.heterologous_expression,'note':r.note,'source_identifier':r.pmid_or_doi,'cohort':'historical2019','strict_shared_preparation':'not established'})
m=pd.DataFrame(measurement);m.to_csv(ROOT/'derived/historical_measurements.csv',index=False)
all_summary=[];edges_all=[];rect_all=[];spair_all=[];loso_all=[];components_all=[];paths=[]
for mode in ['condition','study_only']:
 for trait in ['S','vC']:
  z=m[m.trait==trait].copy();z['source']=z.study+(' | pH='+z.pH_compiled if mode=='condition' else '')
  edges=[]
  for (e,s),g in z.groupby(['species_key','source'],sort=True):
   edges.append({'mode':mode,'trait':trait,'enzyme':e,'source':s,'study':g.study.iloc[0],'value':float(np.exp(np.log(g.value).mean())),'members':';'.join(g.measurement_id),'n_rows':len(g),'pH':';'.join(sorted(g.pH_compiled.unique())),'expression':';'.join(sorted(set(g.heterologous_expression.astype(str))))})
  ed=pd.DataFrame(edges);edges_all+=edges
  adj={}
  for r in edges:
   e='enzyme:'+r['enzyme'];s='source:'+r['source'];adj.setdefault(e,set()).add(s);adj.setdefault(s,set()).add(e)
  seen=set();components=[]
  for node in sorted(adj):
   if node in seen:continue
   todo=[node];comp=set()
   while todo:
    v=todo.pop()
    if v in comp:continue
    comp.add(v);todo+=sorted(adj[v]-comp)
   seen|=comp;components.append(comp)
  components.sort(key=lambda c:(-len(c),sorted(c)))
  for ci,c in enumerate(components):
   for node in sorted(c):components_all.append({'mode':mode,'trait':trait,'component':ci,'node':node,'component_size':len(c)})
  anchors=['Galdieria partita Tokara','Griffithsia monilis','Zea mays','Spinacia oleracea','Nicotiana tabacum']
  for a,b in itertools.combinations(anchors,2):
   start='enzyme:'+a;end='enzyme:'+b;q=[start];prev={start:None}
   for node in q:
    if node==end:break
    for v in sorted(adj.get(node,set())):
     if v not in prev:prev[v]=node;q.append(v)
   path=[]
   if end in prev:
    cur=end
    while cur is not None:path.append(cur);cur=prev[cur]
    path.reverse()
   paths.append({'mode':mode,'trait':trait,'a':a,'b':b,'connected':bool(path),'edges':len(path)-1 if path else None,'path':json.dumps(path)})
  values={(r['enzyme'],r['source']):r['value'] for r in edges}
  studies={r['source']:r['study'] for r in edges}
  bysource={s:set(g.enzyme) for s,g in ed.groupby('source')}
  rect=[]
  for s,t in itertools.combinations(sorted(bysource),2):
   if studies[s]==studies[t]:continue
   common=sorted(bysource[s]&bysource[t])
   for a,b in itertools.combinations(common,2):
    a1,b1,a2,b2=[values[k] for k in [(a,s),(b,s),(a,t),(b,t)]]
    da,db=np.log(a2/a1),np.log(b2/b1);d=da-db
    rect.append({'mode':mode,'trait':trait,'s':s,'t':t,'study_s':studies[s],'study_t':studies[t],'a':a,'b':b,'a_s':a1,'b_s':b1,'a_t':a2,'b_t':b2,'log_offset_a':da,'log_offset_b':db,'cycle_log':d,'factor':float(np.exp(abs(d))),'raw_abs_log_error':float((abs(da)+abs(db))/2),'calibrated_abs_log_error':float(abs(d)),'raw_squared_log_error':float((da*da+db*db)/2),'calibrated_squared_log_error':float(d*d),'order_reversal':bool((a1>b1)!=(a2>b2)),'strict_identical_control':'not certified'})
  rd=pd.DataFrame(rect);rect_all+=rect
  sp=[]
  if len(rd):
   for (s,t),g in rd.groupby(['s','t']):
    sp.append({'mode':mode,'trait':trait,'s':s,'t':t,'study_s':g.study_s.iloc[0],'study_t':g.study_t.iloc[0],'rectangles':len(g),'median_factor':g.factor.median(),'max_factor':g.factor.max(),'fraction_above_1.1':(g.factor>1.1).mean(),'fraction_above_1.25':(g.factor>1.25).mean(),'raw_mae_log':g.raw_abs_log_error.mean(),'calibrated_mae_log':g.calibrated_abs_log_error.mean(),'raw_mse_log':g.raw_squared_log_error.mean(),'calibrated_mse_log':g.calibrated_squared_log_error.mean()})
  spd=pd.DataFrame(sp);spair_all+=sp
  if len(spd):
   for omit in sorted(set(spd.study_s)|set(spd.study_t)):
    v=spd[~(spd.study_s.eq(omit)|spd.study_t.eq(omit))]
    if len(v):loso_all.append({'mode':mode,'trait':trait,'omit_study':omit,'source_pairs':len(v),'raw_rms_log':float(np.sqrt(v.raw_mse_log.mean())),'calibrated_rms_log':float(np.sqrt(v.calibrated_mse_log.mean())),'raw_mae_log':v.raw_mae_log.mean(),'calibrated_mae_log':v.calibrated_mae_log.mean()})
  summ={'mode':mode,'trait':trait,'raw_measurements':len(z),'edges':len(ed),'enzymes':ed.enzyme.nunique(),'source_nodes':ed.source.nunique(),'studies':ed.study.nunique(),'components':len(components),'cycle_rank':len(ed)-len(adj)+len(components),'rectangles':len(rd),'source_pairs':len(spd)}
  if len(rd):summ.update({'median_cycle_factor':rd.factor.median(),'max_cycle_factor':rd.factor.max(),'fraction_above_1.1':(rd.factor>1.1).mean(),'fraction_above_1.25':(rd.factor>1.25).mean(),'raw_mae_log':rd.raw_abs_log_error.mean(),'calibrated_mae_log':rd.calibrated_abs_log_error.mean(),'raw_rms_log':float(np.sqrt(rd.raw_squared_log_error.mean())),'calibrated_rms_log':float(np.sqrt(rd.calibrated_squared_log_error.mean())),'source_pair_balanced_raw_rms_log':float(np.sqrt(spd.raw_mse_log.mean())),'source_pair_balanced_calibrated_rms_log':float(np.sqrt(spd.calibrated_mse_log.mean())),'source_pair_balanced_raw_mae_log':spd.raw_mae_log.mean(),'source_pair_balanced_calibrated_mae_log':spd.calibrated_mae_log.mean()})
  all_summary.append(summ)
write(edges_all,'graph_edges.csv');write(rect_all,'control_rectangles.csv');write(spair_all,'source_pair_metrics.csv');write(loso_all,'leave_one_source_out.csv');write(components_all,'graph_components.csv');write(paths,'anchor_paths.csv');write(all_summary,'graph_summary.csv')
(OUT/'summary.json').write_text(json.dumps(all_summary,indent=2))
env={'executed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'versions':{p:importlib.metadata.version(p) for p in ['numpy','pandas','openpyxl','scipy','matplotlib']},'protocol_sha256':hashlib.sha256((ROOT/'protocol.md').read_bytes()).hexdigest()}
(OUT/'environment.json').write_text(json.dumps(env,indent=2))
print(pd.DataFrame(all_summary).to_string(index=False))
