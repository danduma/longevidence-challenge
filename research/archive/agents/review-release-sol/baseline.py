#!/usr/bin/env python3
"""Execute baseline_protocol.md. Inputs and grouping are independent of scores."""
from pathlib import Path
import json,warnings
import numpy as np,pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr
from sklearn.metrics import r2_score,mean_squared_error
O=Path('outputs');T=['S','vC','KC','KO','vO','KRuBP'];AA='ACDEFGHIKLMNPQRSTVWY'
def feat(s):
 return [s.count(a)/len(s) for a in AA]+[sum(s[i:i+2]==a+b for i in range(len(s)-1))/max(1,len(s)-1) for a in AA for b in AA]+[len(s)]
def metric(y,p):
 with warnings.catch_warnings():
  warnings.simplefilter('ignore');rho=float(spearmanr(y,p).statistic)
 return {'r2':r2_score(y,p),'rmse':np.sqrt(mean_squared_error(y,p)),'spearman':rho}
def groups(d,pairs,threshold):
 parent=list(range(len(d)))
 def root(i):
  while i!=parent[i]:parent[i]=parent[parent[i]];i=parent[i]
  return i
 def join(i,j):parent[root(i)]=root(j)
 for col in ['base_species','sequence']:
  for _,inds in d.groupby(col).groups.items():
   inds=list(inds)
   for i in inds[1:]:join(inds[0],i)
 if threshold:
  for r in pairs.itertuples():join(r.i,r.j)
 return np.array([root(i) for i in range(len(d))])
def folds(g,seed):
 rng=np.random.default_rng(seed);unique=np.unique(g);order=rng.permutation(unique);sizes=np.zeros(5,int);f=np.empty(len(g),int)
 for group in order:
  k=int(np.argmin(sizes));inds=np.where(g==group)[0];f[inds]=k;sizes[k]+=len(inds)
 return f

def main():
 d=pd.read_csv(O/'baseline_data.csv');pairs=pd.read_csv(O/'sequence_pairs_99.csv');X=np.array([feat(s) for s in d.sequence]);allmetrics=[];boots=[];groupinfo=[]
 for regime,threshold in [('species_exact',False),('species_99pct',True)]:
  g=groups(d,pairs,threshold);counts=pd.Series(g).value_counts();groupinfo.append(dict(regime=regime,groups=len(counts),largest_group=int(counts.max()),species=len(d)))
  pd.DataFrame({'species':d.species,'group':g}).to_csv(O/f'groups_{regime}.csv',index=False)
  for seed in range(10):
   f=folds(g,seed)
   for t in T:
    valid=d[t].gt(0)&np.isfinite(d[t]);ix=np.where(valid)[0];y=np.log10(d.loc[valid,t].to_numpy());xx=X[valid];ff=f[valid];gg=g[valid];yp=np.full(len(y),np.nan);mu=yp.copy()
    for k in range(5):
     tr=ff!=k;te=ff==k
     if not te.any():continue
     assert not set(gg[tr])&set(gg[te]);assert not set(d.loc[ix[tr],'base_species'])&set(d.loc[ix[te],'base_species']);assert not set(d.loc[ix[tr],'sequence'])&set(d.loc[ix[te],'sequence'])
     model=make_pipeline(StandardScaler(),Ridge(alpha=100));model.fit(xx[tr],y[tr]);yp[te]=model.predict(xx[te]);mu[te]=y[tr].mean()
    assert np.isfinite(yp).all()
    for name,p in [('ridge',yp),('mean',mu)]:allmetrics.append(dict(regime=regime,seed=seed,target=t,model=name,n=len(y),groups=len(set(gg)),**metric(y,p)))
    pred=pd.DataFrame({'species':d.loc[valid,'species'].values,'group':gg,'fold':ff,'y_true':y,'ridge':yp,'mean':mu})
    pred.to_csv(O/f'baseline_pred_{regime}_{t}_seed{seed}.csv',index=False)
    if seed==0:
     rng=np.random.default_rng(20261126);u=np.unique(gg);indices={v:np.where(gg==v)[0] for v in u};diff=[]
     for _ in range(1000):
      b=np.concatenate([indices[v] for v in rng.choice(u,len(u),replace=True)])
      diff.append(np.sqrt(np.mean((y[b]-yp[b])**2))-np.sqrt(np.mean((y[b]-mu[b])**2)))
     boots.append(dict(regime=regime,target=t,n=len(y),rmse_delta=np.sqrt(np.mean((y-yp)**2))-np.sqrt(np.mean((y-mu)**2)),ci025=np.quantile(diff,.025),ci975=np.quantile(diff,.975),draws=1000))
 r=pd.DataFrame(allmetrics);r.to_csv(O/'baseline_metrics.csv',index=False)
 summary=r.groupby(['regime','target','model']).agg(n=('n','first'),groups=('groups','first'),r2_mean=('r2','mean'),r2_min=('r2','min'),r2_max=('r2','max'),rmse_mean=('rmse','mean'),rmse_min=('rmse','min'),rmse_max=('rmse','max'),spearman_mean=('spearman','mean'),spearman_min=('spearman','min'),spearman_max=('spearman','max')).reset_index()
 summary.to_csv(O/'baseline_summary.csv',index=False);pd.DataFrame(boots).to_csv(O/'baseline_bootstrap.csv',index=False);pd.DataFrame(groupinfo).to_csv(O/'baseline_group_summary.csv',index=False)
 print(summary.to_string(index=False));print(pd.DataFrame(boots).to_string(index=False));print(groupinfo)
if __name__=='__main__':main()
