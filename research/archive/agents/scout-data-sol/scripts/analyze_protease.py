"""Audit measured reporter-to-biochemical portability and simple sequence baselines.
No new proteases are designed. This is not a reproduction of MLDEEP training.
"""
import json,re,zipfile
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats,sparse
from scipy.sparse.linalg import lsqr
rng=np.random.default_rng(20261126); out=Path('results')
with zipfile.ZipFile('data/protease_DMS.zip') as z:
 with z.open('Data_DMS.txt') as f:df=pd.read_csv(f,sep='\t')
controls=pd.read_csv('data/protease_controls.tsv',sep='\t',dtype={'TEVp':str})
vitro=pd.read_csv('data/protease_vitro_relative.tsv',sep='\t',dtype={'TEVp':str})
WT='ENLYFQS';AA='ACDEFGHIKLMNPQRSTVWY';pos={'P6':0,'P5':1,'P4':2,'P3':3,'P2':4,'P1':5,"P1'":6}
vitro['TEVs']=vitro.apply(lambda r:WT[:pos[r.pos]]+r.TEVs_AA+WT[pos[r.pos]+1:],axis=1)
v=vitro.groupby(['TEVp','TEVs'],as_index=False).vinirel.mean()
c=controls.groupby(['TEVp','TEVs'],as_index=False).agg(AUC=('AUC','mean'),min_reads=('min','min'))
neg_names=['(1-13)-mCherry','I C151A 30*31*','no TEVp']
neg=controls[controls.TEVp.isin(neg_names)].copy()
neg.to_csv(out/'protease_negative_controls.csv',index=False)
background=neg.groupby('TEVs').AUC.mean()
c['background']=c.TEVs.map(background);c['AUC_background_subtracted']=c.AUC-c.background
paired=v.merge(c,on=['TEVp','TEVs'],validate='one_to_one')
paired.to_csv(out/'protease_reporter_vitro_pairs.csv',index=False)
ports=[]
for enz,d in paired.groupby('TEVp'):
 ports.append({'enzyme_label':enz,'n_distinct_substrates':len(d),'rho_raw_AUC_vs_initial_rate':float(stats.spearmanr(d.AUC,d.vinirel).statistic),'rho_background_subtracted_vs_initial_rate':float(stats.spearmanr(d.AUC_background_subtracted,d.vinirel).statistic),'rho_background_vs_initial_rate':float(stats.spearmanr(d.background,d.vinirel).statistic)})
pd.DataFrame(ports).to_csv(out/'protease_assay_portability.csv',index=False)
# Conservatively use the same population of noncatalytic controls as preprocessing code.
# This alternative empirical threshold is our audit definition, not MLDEEP hysteresis labels.
neg99=float(neg.AUC.quantile(.99))
keep=(df['min']>=100)&(df.n_mut_TEVp<=1)&(df.TEVs.str.len()==7)&~df.TEVs.str.contains('[^'+AA+']',regex=True)
d=df[keep].reset_index(drop=True)
# Drop protein truncations represented by '*' in variant labels; wt and single substitutions only.
valid=d.TEVp.eq('wt')|d.TEVp.str.fullmatch(r'[A-Z]\d+[A-Z]')
d=d[valid].reset_index(drop=True)
# Sequence-additive substrate effects + categorical enzyme effect, L2-regularized.
# Substrate one-hot features differ from WT; an unseen double mutant can combine learned effects.
enzymes={x:i for i,x in enumerate(sorted(d.TEVp.unique()))}
rows=[];cols=[]
for i,r in enumerate(d.itertuples()):
 rows.extend([i,i]);cols.extend([0,1+enzymes[r.TEVp]])
 for p,aa in enumerate(r.TEVs):
  if aa!=WT[p]:rows.append(i);cols.append(1+len(enzymes)+p*20+AA.index(aa))
X=sparse.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(d),1+len(enzymes)+140))
y=d.AUC.to_numpy()

def average_precision(labels,scores):
 # Tie-safe average precision, equivalent to integrating precision at distinct thresholds.
 order=np.argsort(-scores,kind='stable'); yy=labels[order].astype(int);ss=scores[order]
 if yy.sum()==0:return float('nan')
 end=np.r_[np.flatnonzero(ss[1:]!=ss[:-1]),len(ss)-1];tp=np.cumsum(yy)[end];fp=end+1-tp
 recall=tp/yy.sum();precision=tp/(tp+fp)
 return float(np.sum(np.diff(np.r_[0,recall])*precision))
def metrics(actual,pred):
 labels=actual>neg99
 return {'n':len(actual),'RMSE_AUC':float(np.sqrt(np.mean((actual-pred)**2))),'spearman':float(stats.spearmanr(actual,pred).statistic) if np.std(pred)>1e-12 else None,'positive_fraction_audit_threshold':float(labels.mean()),'average_precision':average_precision(labels,pred)}

results=[]
blocked_uncertainty=[]
random_test=rng.random(len(d))<.2
splits={'random_pairs':random_test,'holdout_double_or_more_mutant_substrates':d.n_mut_TEVs.to_numpy()>=2}
for split,test in splits.items():
 train=~test;yt=y[train];trainframe=d[train].copy();trainframe['y']=yt
 baseline=np.repeat(yt.mean(),test.sum())
 sm=trainframe.groupby('TEVs').y.mean();spred=d.loc[test,'TEVs'].map(sm).fillna(yt.mean()).to_numpy()
 em=trainframe.groupby('TEVp').y.mean();epred=d.loc[test,'TEVp'].map(em).fillna(yt.mean()).to_numpy()
 # Center targets to avoid penalizing intercept; ridge damp=1, fixed a priori.
 fit=lsqr(X[train],yt-yt.mean(),damp=1.0,atol=1e-9,btol=1e-9,iter_lim=500)
 pred=np.asarray(X[test]@fit[0]).ravel()+yt.mean()
 if split.startswith('holdout'):
  evaldf=pd.DataFrame({'substrate':d.loc[test,'TEVs'].to_numpy(),'y':y[test],'pred':pred,'global_mean':baseline})
  evaldf['SSE_additive']=(evaldf.y-evaldf.pred)**2
  evaldf['SSE_mean']=(evaldf.y-evaldf.global_mean)**2
  tab=evaldf.groupby('substrate').agg(n=('y','size'),SSE_additive=('SSE_additive','sum'),SSE_mean=('SSE_mean','sum'))
  tab['RMSE_additive']=np.sqrt(tab.SSE_additive/tab.n);tab['RMSE_global_mean']=np.sqrt(tab.SSE_mean/tab.n)
  tab.to_csv(out/'protease_heldout_substrate_errors.csv')
  arr=tab[['n','SSE_additive','SSE_mean']].to_numpy()
  ix=rng.integers(len(arr),size=(3000,len(arr)));bb=arr[ix].sum(axis=1)
  difference=np.sqrt(bb[:,1]/bb[:,0])-np.sqrt(bb[:,2]/bb[:,0])
  trainactive=np.asarray(X[train].getnnz(axis=0)).ravel()>0
  unseen=np.asarray(X[test][:,~trainactive].getnnz(axis=1)).ravel()>0
  blocked_uncertainty.append({'split':split,'substrate_groups':len(arr),'cluster_bootstrap_RMSE_additive_minus_global_CI95':np.quantile(difference,[.025,.975]).tolist(),'test_rows_with_any_unseen_nonzero_feature':int(unseen.sum())})
 for model,p in [('global_mean',baseline),('substrate_identity_mean',spred),('enzyme_identity_mean',epred),('additive_substrate_sequence_plus_enzyme',pred)]:
  results.append({'split':split,'baseline':model,'train_n':int(train.sum()),'unseen_test_substrates':len(set(d.loc[test,'TEVs'])-set(trainframe.TEVs)),'lsqr_stop':int(fit[1]) if model.startswith('additive') else None,**metrics(y[test],p)})
pd.DataFrame(results).to_csv(out/'protease_baselines.csv',index=False)
# Negative-control variation; substrate means are descriptive, not independent biological repeats.
negsummary=neg.groupby('TEVs').AUC.agg(['count','mean','std'])
negsummary.to_csv(out/'protease_substrate_background.csv')
# Published results audit: no new confidence tests because run IDs not included.
bench=pd.read_csv('sources/mldeep_experiments_mldeep_source_data.csv')
res={'original_DMS_rows':len(df),'original_unique_enzyme_labels':df.TEVp.nunique(),'original_unique_substrates':df.TEVs.nunique(),'modeled_rows_after_read_and_single_substitution_filters':len(d),'modeled_enzyme_labels':d.TEVp.nunique(),'modeled_substrates':d.TEVs.nunique(),'negative_control_rows':len(neg),'negative_control_AUC_99th_percentile':neg99,'negative_substrate_mean_AUC_min_max':[float(negsummary['mean'].min()),float(negsummary['mean'].max())],'independent_vitro_matched_distinct_pairs':len(paired),'portability':ports,'baseline_metrics':results,'blocked_uncertainty':blocked_uncertainty,'published_benchmark_output_rows':len(bench),'published_benchmark_columns':bench.columns.to_list(),'caveats':['Recorder AUC is not kcat/Km or proof of endogenous protein cleavage.','In-vitro pairs are selected validation constructs; portability is not an unbiased prospective success rate.','Simple baselines and alternative control-quantile labels differ from the authors MLDEEP evaluation; no SOTA ranking claimed.','The double-substrate split tests local combinatorial extrapolation, not arbitrary substrate or enzyme generalization.','No biological replicate variance estimated from the one-record-per-pair screen.']}
(out/'protease_summary.json').write_text(json.dumps(res,indent=2,allow_nan=False));print(json.dumps(res,indent=2))
