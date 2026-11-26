#!/usr/bin/env python3
from pathlib import Path
import json,zipfile,importlib.util
import numpy as np
import pandas as pd
from scipy import stats,sparse
from scipy.sparse.linalg import spsolve
from reconstruct import ROOT,O,DATA,NEG,PAIR,WT
REVIEW=ROOT.parent/'review-protease-astra'
spec=importlib.util.spec_from_file_location('review_audit',ROOT/'sources/review_audit.py');review=importlib.util.module_from_spec(spec);spec.loader.exec_module(review)
CANON=review.CANON
RNG=np.random.default_rng(20261126)
def save(df,name):df.to_csv(O/name,index=False)
def elig(d,depth=100,maxmut=7):return d[d.TEVs.str.fullmatch(CANON)&(d['min']>=depth)&(d.TEVs.map(review.hamming)<=maxmut)]
def load():
 c=pd.read_csv(DATA/'protease_controls.tsv',sep='\t',dtype={'TEVp':str})
 with zipfile.ZipFile(DATA/'protease_DMS.zip') as z:d=pd.read_csv(z.open('Data_DMS.txt'),sep='\t')
 return d,c,pd.read_csv(O/'DMS_corrected.tsv.gz',sep='\t'),pd.read_csv(O/'controls_corrected.tsv.gz',sep='\t',dtype={'TEVp':str})
def transfer(c,k):
 records=[]; predictions=[]
 for panel,depth,m in [('primary',100,7),('singles_and_wt',100,1),('all_depths',0,7)]:
  a=elig(c,depth,m);b=elig(k,depth,m)
  for pop in ['each','common']:
   raw=a[a.TEVp.isin(NEG)].pivot(index='TEVs',columns='TEVp',values='AUC').reindex(columns=NEG).dropna()
   corr=b[b.TEVp.isin(NEG)].pivot(index='TEVs',columns='TEVp',values='AUC').reindex(columns=NEG).dropna()
   if pop=='common':keys=raw.index.intersection(corr.index);raw=raw.loc[keys];corr=corr.loc[keys]
   for stage,t in [('raw',raw),('corrected',corr)]:
    for held in NEG:
     x=t.drop(columns=held).mean(axis=1);y=t[held];selected=x>=x.quantile(.9);high=y>=y.quantile(.9)
     boot=[]
     for _ in range(3000):
      ix=RNG.integers(len(t),size=len(t));boot.append(stats.spearmanr(x.iloc[ix],y.iloc[ix]).statistic)
     lo,hi=np.nanquantile(boot,[.025,.975])
     records.append(dict(panel=panel,population=pop,stage=stage,held=held,n=len(t),rho=x.corr(y,method='spearman'),rho_lo=lo,rho_hi=hi,selected=int(selected.sum()),high=int(high.sum()),overlap=int((selected&high).sum()),precision=high[selected].mean(),base=high.mean()))
     if panel=='primary' and pop=='common':
      predictions.extend(dict(stage=stage,held=held,TEVs=key,predictor=x[key],observed=y[key],selected=bool(selected[key]),high=bool(high[key])) for key in t.index)
 save(pd.DataFrame(records),'transfer.csv');save(pd.DataFrame(predictions),'transfer_keys.csv')
 # independent fourth construct, complete cases unchanged depth
 rows=[]
 for stage,d in [('raw',c),('corrected',k)]:
  t=elig(d);t=t[t.TEVp.isin(NEG+['MBP'])].pivot(index='TEVs',columns='TEVp',values='AUC').reindex(columns=NEG+['MBP']).dropna()
  x=t[NEG].mean(axis=1);y=t.MBP;sel=x>=x.quantile(.9);high=y>=y.quantile(.9)
  rows.append(dict(stage=stage,n=len(t),rho=x.corr(y,method='spearman'),selected=sel.sum(),high=high.sum(),overlap=(sel&high).sum()))
 save(pd.DataFrame(rows),'MBP.csv')
 # screen3 pre-existing protocol
 old=pd.read_csv(REVIEW/'sources/screen3.tsv',sep='\t');lib=old['library.ID']
 old=old[lib.isin(['libLH015','libLH016','libLH019'])&~(old.TEVs.eq(WT)&lib.isin(['libLH016','libLH019']))&old.TEVp.isin(NEG)&old.TEVs.str.fullmatch(CANON)&(old['min']>=100)]
 rows=[]
 for stage,d in [('raw',c),('corrected',k)]:
  a=elig(d);a=a[a.TEVp.isin(NEG)].pivot(index='TEVs',columns='TEVp',values='AUC').dropna().mean(axis=1)
  for rep,g in old.groupby('rep'):
   b=g.pivot_table(index='TEVs',columns='TEVp',values='AUC4').reindex(columns=NEG).dropna().mean(axis=1)
   p=pd.concat([a.rename('x'),b.rename('y')],axis=1).dropna();sel=p.x>=p.x.quantile(.9);high=p.y>=p.y.quantile(.9)
   rows.append(dict(stage=stage,rep=rep,n=len(p),rho=p.x.corr(p.y,method='spearman'),selected=sel.sum(),high=high.sum(),overlap=(sel&high).sum(),adequate=len(p)>=10))
 save(pd.DataFrame(rows),'cross_screen.csv')

def model(d,f,c,k):
 d=elig(d);d=d[(d.n_mut_TEVp<=1)&d.TEVp.str.fullmatch(r'wt|[A-Z]\d+[A-Z]')].copy()
 # assign once in original eligible order; corrected keys inherit masks.
 d['random_test']=np.random.default_rng(20261126).uniform(size=len(d))<.2
 common=d.merge(elig(f)[PAIR+['AUC','min','source']],on=PAIR,suffixes=('_raw','_corrected'),validate='one_to_one')
 save(common[PAIR+['AUC_raw','AUC_corrected','min_raw','min_corrected','random_test','source']],'model_common_keys.csv.gz')
 thresholds=[]
 for stage,t in [('raw',c),('corrected',k)]:
  t=t[t.TEVp.isin(NEG)]
  for label,z in [('all',t),('eligibility_matched',elig(t))]:
   thresholds.append(dict(stage=stage,reference=label,n=len(z),cutoff=z.AUC.quantile(.99)))
 save(pd.DataFrame(thresholds),'thresholds.csv')
 records=[]; cis=[];preds=[]
 for population in ['original_raw','common']:
  if population=='original_raw':p=d.copy();stages=['raw']
  else:p=common.copy();stages=['raw','corrected']
  X=review.design(p)
  for stage in stages:
   y=p['AUC' if population=='original_raw' else 'AUC_'+stage].to_numpy()
   for split,test in [('random',p.random_test.to_numpy()),('blocked',p.n_mut_TEVs.to_numpy()>1)]:
    T=X[~test];mu=y[~test].mean();coef=spsolve((T.T@T).tocsc()+sparse.eye(T.shape[1],format='csc'),T.T@(y[~test]-mu));pr=X[test]@coef+mu
    mean=np.full(test.sum(),mu)
    for ref in ['all','eligibility_matched']:
     cutoff=next(t['cutoff'] for t in thresholds if t['stage']==stage and t['reference']==ref)
     labels=y[test]>cutoff
     records.append(dict(population=population,stage=stage,split=split,reference=ref,train=int((~test).sum()),test=int(test.sum()),motifs=p.loc[test].TEVs.nunique(),positives=int(labels.sum()),cutoff=cutoff,prevalence=labels.mean(),AP=review.ap(labels,pr),RMSE=np.sqrt(np.mean((pr-y[test])**2)),mean_RMSE=np.sqrt(np.mean((mu-y[test])**2)),rho=stats.spearmanr(pr,y[test]).statistic))
    if split=='blocked':
     g=pd.DataFrame(dict(motif=p.loc[test].TEVs.to_numpy(),n=1,ridge=(pr-y[test])**2,baseline=(mean-y[test])**2)).groupby('motif').sum()
     ix=RNG.integers(len(g),size=(3000,len(g)));n,r,b=g.to_numpy()[ix].sum(axis=1).T
     lo,hi=np.quantile(np.sqrt(r/n)-np.sqrt(b/n),[.025,.975]);cis.append(dict(population=population,stage=stage,motifs=len(g),delta=np.sqrt(np.mean((pr-y[test])**2))-np.sqrt(np.mean((mu-y[test])**2)),lo=lo,hi=hi,worse=int((g.ridge>g.baseline).sum())))
     if population=='common':preds.extend(dict(stage=stage,TEVp=a,TEVs=b,observed=v,predicted=z) for a,b,v,z in zip(p.loc[test].TEVp,p.loc[test].TEVs,y[test],pr))
 save(pd.DataFrame(records),'models.csv');save(pd.DataFrame(cis),'blocked_uncertainty.csv');save(pd.DataFrame(preds),'blocked_predictions.csv.gz')
 # candidate label changes at each stage's eligible threshold, fixed candidates
 tr=next(t['cutoff'] for t in thresholds if t['stage']=='raw' and t['reference']=='eligibility_matched');tf=next(t['cutoff'] for t in thresholds if t['stage']=='corrected' and t['reference']=='eligibility_matched')
 common['positive_raw']=common.AUC_raw>tr;common['positive_corrected']=common.AUC_corrected>tf
 save(common.groupby(['positive_raw','positive_corrected']).size().rename('n').reset_index(),'label_transition.csv')

def biochemical(d,c,f,k):
 pos=['P6','P5','P4','P3','P2','P1',"P1'"]
 v=pd.read_csv(DATA/'protease_vitro_relative.tsv',sep='\t',dtype={'TEVp':str});v['TEVs']=[WT[:i]+aa+WT[i+1:] for i,aa in zip(v.pos.map(pos.index),v.TEVs_AA)]
 v=v.groupby(PAIR).vinirel.mean().reset_index()
 # author parent replacement, and same common eligible keys
 def replace(c,d):
  parent=d[(d.TEVp=='wt')&(d.n_mut_TEVs<=1)].copy();parent.TEVp='0'
  return pd.concat([c[c.TEVp!='0'],parent],ignore_index=True)
 cr=replace(c,d);cf=replace(k,f)
 j=v.merge(elig(cr)[PAIR+['AUC','min']],on=PAIR,validate='one_to_one').merge(elig(cf)[PAIR+['AUC','min','dnAUC']],on=PAIR,suffixes=('_raw','_corrected'),validate='one_to_one')
 bg=elig(k);bg=bg[bg.TEVp.isin(NEG)].groupby('TEVs').AUC.mean();j['background']=j.TEVs.map(bg);j['subtracted']=j.AUC_corrected-j.background
 save(j,'biochemical_common_keys.csv')
 rows=[];selection=[];delta=[]
 for name,g in j.groupby('TEVp'):
  high=g.vinirel>=g.vinirel.quantile(.9)
  for score in ['AUC_raw','AUC_corrected','dnAUC','subtracted']:
   sel=g[score]>=g[score].quantile(.9)
   rows.append(dict(enzyme=name,score=score,n=len(g),rho=g[score].corr(g.vinirel,method='spearman'),selected=sel.sum(),high=high.sum(),overlap=(sel&high).sum(),precision=high[sel].mean(),base=high.mean()))
   selection.extend(dict(TEVp=name,TEVs=a,score=score,selected=bool(s),biochem_high=bool(h)) for a,s,h in zip(g.TEVs,sel,high))
  delta.append(dict(enzyme=name,rho_delta=g.AUC_corrected.corr(g.vinirel,method='spearman')-g.AUC_raw.corr(g.vinirel,method='spearman')))
 save(pd.DataFrame(rows),'biochemical.csv');save(pd.DataFrame(selection),'biochemical_selection.csv');save(pd.DataFrame(delta),'biochemical_panel_deltas.csv')
 # Joint motif bootstrap preserving dependence across all six enzyme panels.
 motifs=j.TEVs.unique();vals=[]
 for _ in range(3000):
  sample=pd.DataFrame({'TEVs':RNG.choice(motifs,size=len(motifs),replace=True)}).merge(j,on='TEVs',how='inner')
  ds=[]
  for _,g in sample.groupby('TEVp'):
   ds.append(g.AUC_corrected.corr(g.vinirel,method='spearman')-g.AUC_raw.corr(g.vinirel,method='spearman'))
  vals.append(np.mean(ds))
 (O/'biochemical_uncertainty.json').write_text(json.dumps(dict(mean_delta=float(np.mean([x['rho_delta'] for x in delta])),joint_motif_bootstrap_ci=np.quantile(vals,[.025,.975]).tolist(),enzymes=len(delta),motifs=len(motifs),interpretation='Technical panel sensitivity at fixed rates and fits; no independent cultivation or biochemical replicate variance. Enzyme panels share motifs.'),indent=2)+'\n')

def main():
 d,c,f,k=load()
 print('transfer',flush=True);transfer(c,k)
 print('models',flush=True);model(d,f,c,k)
 print('biochemical',flush=True);biochemical(d,c,f,k)
 print('done',flush=True)
if __name__=='__main__':main()
