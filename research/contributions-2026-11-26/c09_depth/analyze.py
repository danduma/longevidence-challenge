#!/usr/bin/env python3
"""C09: offline final-pair depth sensitivity; no fitting or network calls."""
import os
os.environ['PYTHONDONTWRITEBYTECODE']='1'
from pathlib import Path
import json, hashlib, platform, sys
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr, rankdata
ROOT=Path(__file__).resolve().parent
I=ROOT/'inputs'; O=ROOT/'results'; O.mkdir(exist_ok=True)
CUTS=[20,50,100,200,500]; KEY=['TEVp','TEVs']; WT='ENLYFQS'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(d,n): pd.DataFrame(d).to_csv(O/n,index=False, float_format='%.12g',compression={'method':'gzip','mtime':0} if n.endswith('.gz') else None)
def intact(d): return d.TEVs.str.fullmatch('[ACDEFGHIKLMNPQRSTVWY]{7}')
def load(n):
 d=pd.read_csv(I/n,sep='\t',dtype={'TEVp':str})
 assert not d.duplicated(KEY).any()
 totals=np.column_stack([d[f'N{i}']+d[f'F{i}'] for i in range(1,6)])
 assert np.array_equal(totals.min(axis=1),d['min'])
 assert np.array_equal(totals.sum(axis=1),d['sum'])
 return d

def rho(x,y):
 if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
 return float(spearmanr(x,y).statistic)
def selection(d):
 return d.AUC.ge(d.groupby('TEVs').AUC.transform('quantile',q=.9))

def main():
 manifest=json.loads((ROOT/'sources/input_manifest.json').read_text())
 for x in manifest: assert sha(ROOT/x['local'])==x['sha256'],x['local']
 d=load('DMS_script05_singles.tsv.gz'); full=load('DMS_corrected.tsv.gz'); controls=load('controls_corrected.tsv.gz')
 assert (d.n_mut_TEVp<=1).all() and (d.n_mut_TEVs<=1).all()
 d=d[intact(d)].copy().reset_index(drop=True)
 assert len(d)==355242 and d.TEVs.nunique()==134 and d.TEVp.nunique()==4440
 # Check final intended population is an exact keyed subset, not an off-panel surrogate.
 sub=d[KEY+['AUC','min']].merge(full[KEY+['AUC','min']],on=KEY,validate='one_to_one',suffixes=('_single','_full'))
 assert len(sub)==len(d) and (sub.AUC_single==sub.AUC_full).all() and (sub.min_single==sub.min_full).all()
 base=selection(d); top=d.AUC.eq(d.groupby('TEVs').AUC.transform('max')); glob=d.AUC.ge(d.AUC.quantile(.9))
 ledger=d[KEY+['min','AUC','source']].copy(); ledger['rank20_within_motif']=d.groupby('TEVs').AUC.rank(ascending=False,method='average')
 rows=[]; erows=[]
 common=d[d['min']>=500].copy(); common_sel=selection(common); common_top=common.AUC.eq(common.groupby('TEVs').AUC.transform('max'))
 fixed=[]
 for c in CUTS:
  mask=d['min']>=c; a=d[mask]; s=selection(a); t=a.AUC.eq(a.groupby('TEVs').AUC.transform('max')); g=a.AUC.ge(a.AUC.quantile(.9))
  sets=[('within_motif_decile',set(d.index[base]),set(a.index[s])),('within_motif_top',set(d.index[top]),set(a.index[t])),('global_decile',set(d.index[glob]),set(a.index[g]))]
  counts=a.groupby('TEVp').TEVs.nunique()
  for e,n in counts.items():erows.append(dict(cutoff=c,TEVp=e,motifs=n,complete_134=n==134))
  for label,b,z in sets:
   survive=b&set(a.index); lost=b-set(a.index)
   rows.append(dict(cutoff=c,selection=label,pairs=len(a),enzymes=a.TEVp.nunique(),motifs=a.TEVs.nunique(),complete_enzyme_profiles=int(counts.eq(134).sum()),selected=len(z),baseline_selected=len(b),baseline_selected_eligible=len(survive),baseline_selected_lost_depth=len(lost),baseline_selected_still_selected=len(b&z),baseline_selected_demoted=len(survive-z),newly_selected=len(z-b),jaccard=len(b&z)/len(b|z)))
  ledger[f'eligible{c}']=mask
  ledger[f'selected{c}']=False;ledger.loc[a.index,f'selected{c}']=s
  ledger[f'top{c}']=False;ledger.loc[a.index,f'top{c}']=t
  ledger[f'rank{c}']=a.groupby('TEVs').AUC.rank(ascending=False,method='average').reindex(d.index)
  k=a.loc[common.index]
  assert np.array_equal(k.AUC,common.AUC)
  assert np.array_equal(selection(k),common_sel)
  fixed.append(dict(cutoff=c,pairs=len(k),motifs=k.TEVs.nunique(),auc_max_change=float(np.max(np.abs(k.AUC-common.AUC))),rho=rho(k.AUC,common.AUC),within_motif_rank_changes=int((k.groupby('TEVs').AUC.rank()!=common.groupby('TEVs').AUC.rank()).sum()),selected=int(common_sel.sum()),top=int(common_top.sum())))
 save(rows,'candidate_summary.csv');save(erows,'enzyme_coverage.csv');save(fixed,'fixed_common.csv');save(ledger,'candidate_trajectories.csv.gz')
 # Off-panel population never enters primary candidate scores.
 off=[]
 for c in CUTS:
  a=full[intact(full)&(full['min']>=c)]
  for lab,m in [('intended',a.n_mut_TEVp.le(1)&a.n_mut_TEVs.le(1)),('off_panel_motif',a.n_mut_TEVp.le(1)&a.n_mut_TEVs.gt(1)),('other_enzyme',a.n_mut_TEVp.gt(1))]:
   z=a[m];off.append(dict(cutoff=c,population=lab,pairs=len(z),enzymes=z.TEVp.nunique(),motifs=z.TEVs.nunique()))
 save(off,'off_panel_coverage.csv')
 # Compact identity ledgers for review; full ledger retains all keys.
 save(ledger[ledger.top20], 'top_candidates20.csv')
 save(ledger[ledger.top20 & ledger.top500], 'top_candidates_stable.csv')
 v=pd.read_csv(I/'protease_vitro_relative.tsv',sep='\t',dtype={'TEVp':str});pos=['P6','P5','P4','P3','P2','P1',"P1'"]
 v['TEVs']=[WT[:pos.index(p)]+aa+WT[pos.index(p)+1:] for p,aa in zip(v.pos,v.TEVs_AA)]
 raw_n=len(v);v=v.groupby(KEY,as_index=False).vinirel.mean()
 panels=sorted(v.TEVp.unique()); bio=[]; bl=[]; boot=[]
 for branch in ['DMS_wt_replacement','manual_parent']:
  ref=controls[controls.TEVp.isin(panels)&intact(controls)].copy();ref['parent_provenance']='corrected_controls'
  if branch=='DMS_wt_replacement':
   p=d[d.TEVp.eq('wt')].copy();p.TEVp='0';p['parent_provenance']='corrected_DMS_wt'
   ref=pd.concat([ref[ref.TEVp.ne('0')],p],ignore_index=True)
  assert not ref.duplicated(KEY).any()
  j=v.merge(ref[KEY+['AUC','min','source','parent_provenance']],on=KEY,how='left',validate='one_to_one')
  j['branch']=branch
  for e,g in j.groupby('TEVp'):
   b=g[g['min']>=20]; bs=set(b.index[b.AUC.ge(b.AUC.quantile(.9))]);bh=set(b.index[b.vinirel.ge(b.vinirel.quantile(.9))])
   j.loc[g.index,'rate_high20']=g.index.isin(bh);j.loc[g.index,'recorder_selected20']=g.index.isin(bs)
   fixedg=g[g['min']>=500]
   for c in CUTS:
    eligible=g['min']>=c;a=g[eligible];s=a.AUC.ge(a.AUC.quantile(.9));h=a.vinirel.ge(a.vinirel.quantile(.9))
    ix=set(a.index[s]);av=set(a.index)
    for population,z in [('available',a),('fixed500',fixedg)]:
     sz=z.AUC.ge(z.AUC.quantile(.9));hz=z.vinirel.ge(z.vinirel.quantile(.9))
     bio.append(dict(branch=branch,enzyme=e,cutoff=c,population=population,n=len(z),rho=rho(z.AUC,z.vinirel),selected=int(sz.sum()),rate_high=int(hz.sum()),overlap=int((sz&hz).sum()),baseline_selected_lost=len(bs-av) if population=='available' else None,baseline_selected_retained=len(bs&ix) if population=='available' else None,baseline_selected_demoted=len((bs&av)-ix) if population=='available' else None))
    j.loc[g.index,f'eligible{c}']=eligible
    j.loc[g.index,f'selected{c}']=g.index.isin(ix)
    j.loc[g.index,f'rank{c}']=a.AUC.rank(ascending=False).reindex(g.index)
  bl.append(j)
  # Paired motif bootstrap across six panels, allowing eligibility populations to differ.
  good=j[j['min']>=20];motifs=sorted(good.TEVs.unique());mapping={m:i for i,m in enumerate(motifs)}
  arr=[]
  for e,g in good.groupby('TEVp'):arr.append((g.TEVs.map(mapping).to_numpy(),g.AUC.to_numpy(),g.vinirel.to_numpy(),g['min'].to_numpy()))
  rng=np.random.default_rng(20261126);values=[]
  for _ in range(2000):
   multiplicity=np.bincount(rng.integers(len(motifs),size=len(motifs)),minlength=len(motifs));ds=[]
   for ix,x,y,dep in arr:
    rep=np.repeat(np.arange(len(ix)),multiplicity[ix]);hi=rep[dep[rep]>=500]
    ds.append(rho(x[hi],y[hi])-rho(x[rep],y[rep]))
   values.append(np.nanmean(ds))
  btab=pd.DataFrame(bio);bb=btab[(btab.branch==branch)&(btab.population=='available')].pivot(index='enzyme',columns='cutoff',values='rho')
  boot.append(dict(branch=branch,comparison='500 minus 20; changing eligible populations',mean_delta=float((bb[500]-bb[20]).mean()),lo=float(np.nanquantile(values,.025)),hi=float(np.nanquantile(values,.975)),draws=2000,motifs=len(motifs),panels=len(arr),seed=20261126))
 save(bio,'biochemical_summary.csv');joined=pd.concat(bl,ignore_index=True)
 save(joined,'biochemical_trajectories.csv');save(boot,'biochemical_uncertainty.csv')
 save(joined[joined.recorder_selected20.eq(True)],'biochemical_selected20.csv')
 save(joined[joined.recorder_selected20.eq(True)&joined.rate_high20.eq(True)&joined.eligible500.eq(False)],'biochemical_high_calls_lost500.csv')
 # Fixed-panel null checked numerically, with the existing depth100 source result
 # serving only as a reconstruction sanity check (not counted as new work).
 bdf=pd.DataFrame(bio)
 fixed_bio=bdf[bdf.population.eq('fixed500')]
 assert fixed_bio.groupby(['branch','enzyme']).rho.nunique().eq(1).all()
 assert fixed_bio.groupby(['branch','enzyme']).overlap.nunique().eq(1).all()
 expected={'0':(134,.811599010326),'NX1':(126,.742630538653),'NX2':(131,.838313268253),'VI':(132,.772200510345),'VI*':(133,.706759131610),'VII':(130,.670064306798)}
 for e,(n,r) in expected.items():
  z=bdf[(bdf.branch=='DMS_wt_replacement')&(bdf.enzyme==e)&(bdf.cutoff==100)&(bdf.population=='available')].iloc[0]
  assert z['n']==n and abs(z.rho-r)<1e-11
 # Assert a direct keyed set reconciliation for every depth table.
 for row in rows:
  assert row['baseline_selected']==row['baseline_selected_lost_depth']+row['baseline_selected_demoted']+row['baseline_selected_still_selected']
  assert row['selected']==row['baseline_selected_still_selected']+row['newly_selected']
 environment=dict(id='C09',python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,platform=platform.platform(),cutoffs=CUTS,biochemical_input_rows=raw_n,biochemical_distinct_keys=len(v),checks='input SHA256, pooled depth/sum, unique keys, exact intended subset, fixed common scores and calls invariant: PASS')
 (O/'verification.json').write_text(json.dumps(environment,indent=2)+'\n')
 outputs={p.name:sha(p) for p in sorted(O.iterdir()) if p.is_file() and p.name!='output_hashes.json'}
 (O/'output_hashes.json').write_text(json.dumps(outputs,indent=2)+'\n')
 print(json.dumps(environment,indent=2))
if __name__=='__main__':main()
