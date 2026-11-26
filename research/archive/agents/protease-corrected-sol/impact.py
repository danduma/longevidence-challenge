import numpy as np,pandas as pd,json
from scipy import sparse
from scipy.sparse.linalg import spsolve
from evaluate import load,elig,review,PAIR,O,save,RNG

d,c,f,k=load();d=elig(d);d=d[(d.n_mut_TEVp<=1)&d.TEVp.str.fullmatch(r'wt|[A-Z]\d+[A-Z]')].copy()
X=review.design(d);test=d.n_mut_TEVs.to_numpy()>1;y=d.AUC.to_numpy();mu=y[~test].mean();T=X[~test]
coef=spsolve((T.T@T).tocsc()+sparse.eye(T.shape[1],format='csc'),T.T@(y[~test]-mu));pred=X[test]@coef+mu
hold=d.loc[test].copy();hold['prediction']=pred;hold['ridge_sq']=(pred-y[test])**2;hold['mean_sq']=(mu-y[test])**2
fin=elig(f);hold=hold.merge(fin[PAIR+['AUC']],on=PAIR,how='left',suffixes=('_raw','_corrected'));hold['survives']=hold.AUC_corrected.notna();hold['excess_sq']=hold.ridge_sq-hold.mean_sq
save(hold,'blocked_fixed_fit_keys.csv.gz')
g=hold.groupby('TEVs').agg(n=('TEVp','size'),retained=('survives','sum'),raw_mean=('AUC_raw','mean'),corrected_mean=('AUC_corrected','mean'),excess_sq=('excess_sq','sum'),ridge_sq=('ridge_sq','sum'),mean_sq=('mean_sq','sum')).reset_index();save(g,'blocked_motif_impact.csv')
g=hold.groupby('survives').agg(n=('TEVp','size'),ridge_sq=('ridge_sq','sum'),mean_sq=('mean_sq','sum'),excess_sq=('excess_sq','sum')).reset_index();g['RMSE']=np.sqrt(g.ridge_sq/g.n);g['mean_RMSE']=np.sqrt(g.mean_sq/g.n);save(g,'blocked_survival_decomposition.csv')
j=d.merge(fin[PAIR+['AUC']],on=PAIR,suffixes=('_raw','_corrected'));sr=j.AUC_raw>=j.AUC_raw.quantile(.9);sf=j.AUC_corrected>=j.AUC_corrected.quantile(.9)
info=dict(common=len(j),raw_selected=int(sr.sum()),corrected_selected=int(sf.sum()),overlap=int((sr&sf).sum()),removed_from_raw_selected=int((sr&~sf).sum()),added=int((~sr&sf).sum()))
(O/'candidate_decile_impact.json').write_text(json.dumps(info,indent=2)+'\n')
