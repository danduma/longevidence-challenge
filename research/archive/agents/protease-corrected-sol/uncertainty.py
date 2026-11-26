import numpy as np,pandas as pd,json
from scipy.stats import spearmanr
from reconstruct import O,PAIR
rng=np.random.default_rng(20261126);t=pd.read_csv(O/'transfer_keys.csv');rows=[]
for held,g in t.groupby('held'):
 r=g[g.stage=='raw'].set_index('TEVs');f=g[g.stage=='corrected'].set_index('TEVs').reindex(r.index)
 a=r[['predictor','observed']].to_numpy();b=f[['predictor','observed']].to_numpy();delta=[];prec=[]
 for _ in range(3000):
  ix=rng.integers(len(a),size=len(a));aa=a[ix];bb=b[ix]
  delta.append(spearmanr(bb[:,0],bb[:,1]).statistic-spearmanr(aa[:,0],aa[:,1]).statistic)
  sel=bb[:,0]>=np.quantile(bb[:,0],.9);high=bb[:,1]>=np.quantile(bb[:,1],.9);prec.append(high[sel].mean())
 rows.append(dict(held=held,n=len(a),rho_delta=spearmanr(b[:,0],b[:,1]).statistic-spearmanr(a[:,0],a[:,1]).statistic,delta_lo=np.quantile(delta,.025),delta_hi=np.quantile(delta,.975),corrected_precision_lo=np.quantile(prec,.025),corrected_precision_hi=np.quantile(prec,.975)))
pd.DataFrame(rows).to_csv(O/'transfer_paired_uncertainty.csv',index=False)
h=pd.read_csv(O/'blocked_fixed_fit_keys.csv.gz');l=pd.read_csv(O/'DMS_all_upstream_pair_ledger.tsv.gz',sep='\t');j=h.merge(l[PAIR+['reason','min']],on=PAIR,validate='one_to_one');j['exclusion']=np.where(j.survives,'retained_depth100',np.where(j.reason=='retained','corrected_min_below100',j.reason))
j.groupby('exclusion').agg(n=('TEVp','size'),excess_sq=('excess_sq','sum')).reset_index().to_csv(O/'blocked_exclusion_reasons.csv',index=False)
