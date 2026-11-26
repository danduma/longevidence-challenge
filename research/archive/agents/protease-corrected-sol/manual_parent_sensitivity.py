import pandas as pd,numpy as np
from evaluate import load,elig,WT,PAIR,save
from reconstruct import DATA
_,c,_,k=load();pos=['P6','P5','P4','P3','P2','P1',"P1'"]
v=pd.read_csv(DATA/'protease_vitro_relative.tsv',sep='\t',dtype={'TEVp':str});v['TEVs']=[WT[:i]+aa+WT[i+1:] for i,aa in zip(v.pos.map(pos.index),v.TEVs_AA)];v=v.groupby(PAIR).vinirel.mean().reset_index()
j=v.merge(elig(c)[PAIR+['AUC']],on=PAIR).merge(elig(k)[PAIR+['AUC','dnAUC']],on=PAIR,suffixes=('_raw','_corrected'))
rows=[]
for e,g in j.groupby('TEVp'):
 high=g.vinirel>=g.vinirel.quantile(.9)
 for score in ['AUC_raw','AUC_corrected','dnAUC']:
  sel=g[score]>=g[score].quantile(.9)
  rows.append(dict(enzyme=e,score=score,n=len(g),rho=g[score].corr(g.vinirel,method='spearman'),selected=sel.sum(),high=high.sum(),overlap=(sel&high).sum()))
save(pd.DataFrame(rows),'biochemical_manual_parent_sensitivity.csv')
