"""Animal-aware exploratory contrast of normalized primary CPM data.
Uses paired intact planes, independent injured groups, n=3 in each.
No edgeR raw-count model or adult mouse regeneration claim is made.
"""
from pathlib import Path
import json
import pandas as pd,numpy as np
from scipy import stats
out=Path('results')
d=pd.read_csv('data/regeneration2026_CPM_extracted.tsv.gz',sep='\t')
ids=d.iloc[:,0].astype(str)
a=d.iloc[:,1:19].apply(pd.to_numeric,errors='raise').to_numpy(float)
assert a.shape==(139698,18) and np.isfinite(a).all() and (a>=0).all()
# Input order verified against workbook: intact E,D pairs, 3d E then D, 5d E then D.
E0=[0,2,4];D0=[1,3,5];E3=[6,7,8];D3=[9,10,11];E5=[12,13,14];D5=[15,16,17]
expressed=(a>=1).sum(axis=1)>=3

def bh(p):
 order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q

def contrast(log,injE,injD):
 b=log[:,D0]-log[:,E0]
 diff0=b.mean(axis=1);ei=log[:,injE];di=log[:,injD]
 raw=di.mean(axis=1)-ei.mean(axis=1);effect=raw-diff0
 var_parts=np.column_stack([b.var(axis=1,ddof=1)/3,di.var(axis=1,ddof=1)/3,ei.var(axis=1,ddof=1)/3])
 se=np.sqrt(var_parts.sum(axis=1));den=(var_parts**2/2).sum(axis=1)
 df=np.divide(var_parts.sum(axis=1)**2,den,out=np.full(len(log),np.inf),where=den>0)
 t=np.divide(effect,se,out=np.zeros(len(log)),where=se>0)
 p=2*stats.t.sf(abs(t),df);p[(se==0)&(effect!=0)]=0
 return raw,diff0,effect,se,df,p

l=np.log2(a+1);x=contrast(l,E3,D3);y=contrast(l,E5,D5)
r=pd.DataFrame({'transcript_id':ids,'annotated_symbol':ids.str.split('|').str[-1].str.removeprefix('gene-'),'expressed':expressed,'intact_log2_D_minus_E':x[1]})
for day,v in [(3,x),(5,y)]:
 raw,_,eff,se,df,p=v
 r[f'raw_plane_log2_difference_{day}d']=raw;r[f'injury_plane_interaction_{day}d']=eff;r[f'SE_interaction_{day}d']=se;r[f'df_{day}d']=df;r[f'p_{day}d']=p
 q=np.full(len(r),np.nan);q[expressed]=bh(p[expressed]);r[f'q_BH_{day}d']=q
 r[f'CI95_low_{day}d']=eff-stats.t.ppf(.975,df)*se;r[f'CI95_high_{day}d']=eff+stats.t.ppf(.975,df)*se
 for label,inds in [('intact_E',E0),('intact_D',D0),(f'{day}d_E',E3 if day==3 else E5),(f'{day}d_D',D3 if day==3 else D5)]:
  r['mean_CPM_'+label]=a[:,inds].mean(axis=1)
r.to_csv(out/'regeneration_interactions.tsv.gz',sep='\t',index=False)
markers=r[r.annotated_symbol.str.fullmatch('TNFSF11|TNFRSF11A|LOC138491483|CTSK|MMP9|SALL4|PRRX1',case=False,na=False)]
markers.to_csv(out/'regeneration_prior_art_markers.csv',index=False)
# Check sensitivity to arbitrary pseudocount (0.1 versus 1); full reanalysis is exploratory.
l01=np.log2(a+.1);x01=contrast(l01,E3,D3)[2];y01=contrast(l01,E5,D5)[2]
large3=expressed&(abs(x[0])>=1);large5=expressed&(abs(y[0])>=1)
q3=r.q_BH_3d.to_numpy();q5=r.q_BH_5d.to_numpy()
res={'input_rows':len(d),'samples':a.shape[1],'expressed_rows_CPM1_in_at_least3':int(expressed.sum()),'samples_per_group':3,'intact_samples_paired_by_animal':True,'regenerating_samples_independent_animals':True,'both_forelimbs_pooled_per_injured_sample':True,'CPM_column_sum_range':[float(a.sum(axis=0).min()),float(a.sum(axis=0).max())],'raw_large_plane_difference_count_3d':int(large3.sum()),'raw_large_plane_difference_count_5d':int(large5.sum()),'large_raw_plane_difference_falls_below1_after_intact_subtraction_3d':int((large3&(abs(x[2])<1)).sum()),'large_raw_plane_difference_falls_below1_after_intact_subtraction_5d':int((large5&(abs(y[2])<1)).sum()),'interaction_BH05_3d':int((q3<.05).sum()),'interaction_BH05_5d':int((q5<.05).sum()),'interaction_BH05_both_times':int(((q3<.05)&(q5<.05)).sum()),'interaction_BH05_both_and_same_direction':int(((q3<.05)&(q5<.05)&(x[2]*y[2]>0)).sum()),'interaction_time_spearman_shared_baseline_not_independent':float(stats.spearmanr(x[2][expressed],y[2][expressed]).statistic),'pseudocount_sensitivity_spearman_3d':float(stats.spearmanr(x[2][expressed],x01[expressed]).statistic),'pseudocount_sensitivity_spearman_5d':float(stats.spearmanr(y[2][expressed],y01[expressed]).statistic),'blank_annotation_rows_expressed':int((expressed&r.annotated_symbol.eq('')).sum()),'marker_records':json.loads(markers[['transcript_id','annotated_symbol','injury_plane_interaction_3d','injury_plane_interaction_5d','q_BH_3d','q_BH_5d']].to_json(orient='records')),'caveats':['Welch/Satterthwaite tests on log CPM with small n are exploratory; do not reproduce edgeR raw-count tests.','Same intact baseline in both times induces correlation; these are not independent replication cohorts.','Gene annotations may include merged loci; each row retains its source identifier.','Differences can reflect cell composition, maturity, dissection, and injury geometry rather than causal regeneration switches.','All work is computational reanalysis of existing axolotl data; no adult mouse limb evidence produced.']}
(out/'regeneration_summary.json').write_text(json.dumps(res,indent=2,allow_nan=False));print(json.dumps(res,indent=2))
