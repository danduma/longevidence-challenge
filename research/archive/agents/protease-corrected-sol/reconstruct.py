#!/usr/bin/env python3
"""Port of pinned ProtRec scripts 01.1, 02, 03, 04, 05 (min20DNA branch).
See protocol.md; reads shared inputs without modification.
"""
from pathlib import Path
import zipfile,json,hashlib,platform,math
import numpy as np
import pandas as pd
from scipy.stats import norm
ROOT=Path(__file__).resolve().parent
S=ROOT/'sources'; O=ROOT/'results';O.mkdir(exist_ok=True)
DATA=ROOT.parent/'scout-data-sol/data'
PAIR=['TEVp','TEVs']; C=[f'{a}{i}' for a in ['N','F'] for i in range(1,6)]
NEG=['(1-13)-mCherry','I C151A 30*31*','no TEVp'];WT='ENLYFQS'

def rround(x,digits):
 # Literal finite, small-nonnegative-digit branch of R 4.2 src/nmath/fround.c.
 x=np.asarray(x,dtype=float);sgn=np.sign(x);z=np.abs(x);p=10.**digits;t=z*p;i=np.floor(t);down=i/p;up=np.ceil(t)/p
 return sgn*np.where((up-z<z-down)|((i%2==1)&(up-z==z-down)),up,down)

def measures(d, auc_digits=4, fraction_round=None):
 d=d.copy(); n=d[[f'N{i}' for i in range(1,6)]].to_numpy();f=d[[f'F{i}' for i in range(1,6)]].to_numpy()
 with np.errstate(invalid='ignore',divide='ignore'):q=f/(n+f)
 d['min']=(n+f).min(axis=1);d['sum']=(n+f).sum(axis=1)
 for i in range(5):d[f't{i+1}']=rround(q[:,i],3 if fraction_round is None else fraction_round)
 if fraction_round is not None:q=rround(q,fraction_round)
 # Preserve original minute scaling for reaggregated >=3 groups, and R's
 # extended-precision rowSums accumulation via correctly-rounded math.fsum.
 if fraction_round is None:
  pieces=q[:,:-1]+(q[:,1:]-q[:,:-1])*.5
  auc=np.array([math.fsum(row) if np.isfinite(row).all() else np.nan for row in pieces])/4
 else:
  pieces=(60*q[:,:-1]+60*(q[:,1:]-q[:,:-1])*.5)/240
  auc=np.array([math.fsum(row) if np.isfinite(row).all() else np.nan for row in pieces])
 d['AUC']=rround(auc,auc_digits)
 return d

def pool(b, fraction_round=None):
 d=measures(b.groupby(PAIR,sort=True)[C].sum().reset_index(),fraction_round=fraction_round)
 cols=PAIR+[x for x in ['n_mut_TEVp','n_mut_TEVs'] if x in b]
 anno=b[cols].drop_duplicates();assert not anno.duplicated(PAIR).any()
 return d.merge(anno,on=PAIR,validate='one_to_one')

def noise(d):
 a=d.loc[d.TEVp.isin(NEG)&(d['min']>=20),'AUC'].to_numpy()
 return dict(n=len(a),mu=float(a.mean()),sigma=float(a.std(ddof=0)),cutoff=float(a.mean()+norm.ppf(.95)*a.std(ddof=0)))

def correct(b, cutoff, posname, name, save=True):
 # b has already had mutated backbones removed. Retain original order for which.max.
 b=b.copy();b['plasmidID']=b.BC_TEVp+'_'+b.BC_TEVs
 b['below20DNA']=b['min']<20
 e=b.loc[~b.below20DNA].copy()
 assert np.isfinite(e.AUC).all()
 orig=pool(e); orig=orig[orig['min']>=20]
 g=e.groupby(PAIR,sort=True).AUC.agg(['size','min','max'])
 div=g[(g['size']==2)&(g['max']>=cutoff)&(g['max']/g['min']>=3)].reset_index()[PAIR]
 active=orig[orig.AUC>=cutoff][PAIR]
 reps=e.merge(active,on=PAIR,how='inner',sort=False)
 reps=reps[reps.groupby(PAIR).AUC.transform('size')>=3]
 top=reps.loc[reps.groupby(PAIR,sort=True).AUC.idxmax(),PAIR+['plasmidID']]
 wo=reps.merge(top.assign(top=True),on=PAIR+['plasmidID'],how='left')
 wo=pool(wo[wo.top.isna()])
 comp=orig[PAIR+['AUC']].merge(wo[PAIR+['AUC','min']],on=PAIR,suffixes=('','_wo_max')).merge(top,on=PAIR)
 comp['fc']=comp.AUC/comp.AUC_wo_max
 bad=comp[(comp.AUC_wo_max!=0)&(comp['min']>=20)&(comp.fc>=2)]
 # R removes by plasmidID globally, as in the original source.
 retained_reps=reps[~reps.plasmidID.isin(bad.plasmidID)]
 adjusted=pool(retained_reps,fraction_round=4)
 adjusted=adjusted[adjusted['min']>=20]
 adjusted=adjusted.merge(bad[PAIR].assign(outlier=True),on=PAIR,how='left')
 adjusted['source']=np.where(adjusted.outlier.fillna(False),'repl_fc_corrected','repl_uncorrected')
 unadj=orig.merge(adjusted[PAIR].assign(reaggregated=True),on=PAIR,how='left')
 unadj=unadj[unadj.reaggregated.isna()].drop(columns='reaggregated');unadj['source']='original'
 final=pd.concat([adjusted.drop(columns='outlier'),unadj],ignore_index=True)
 ledger=orig[PAIR+['AUC','min']].rename(columns={'AUC':'AUC_pre_outlier','min':'min_pre_outlier'})
 ledger=ledger.merge(div.assign(diverging_duplicate=True),on=PAIR,how='left')
 final=final.merge(div.assign(diverging=True),on=PAIR,how='left');final=final[final.diverging.isna()].drop(columns='diverging')
 pos=final[(final.TEVp==posname)&(final.TEVs=='ENLYFQ*')].AUC
 assert len(pos)==1 or posname=='none'
 positive=float(pos.iloc[0]) if len(pos) else np.inf
 above=final[~final.TEVs.str.contains(r'\*')&(final.AUC>=positive)][PAIR]
 ledger=ledger.merge(above.assign(above_positive_control=True),on=PAIR,how='left')
 final=final.merge(above.assign(above=True),on=PAIR,how='left');final=final[final.above.isna()].drop(columns='above')
 ledger=ledger.merge(final[PAIR+['AUC','min','source']],on=PAIR,how='left')
 ledger['retained']=ledger.AUC.notna()
 b['removed_outlier']=b.plasmidID.isin(bad.plasmidID)
 b=b.merge(div.assign(diverging_pair=True),on=PAIR,how='left').merge(above.assign(above_positive_control=True),on=PAIR,how='left')
 meta=dict(input_barcodes=len(b),eligible_barcodes=len(e),pooled_pairs=len(orig),diverging_pairs=len(div),active_threeplus=len(top),outlier_barcodes=len(bad),above_positive_pairs=len(above),positive_control=positive,final_pairs=len(final),noise_cutoff=cutoff)
 if save:
  final.to_csv(O/f'{name}_corrected.tsv.gz',sep='\t',index=False)
  ledger.to_csv(O/f'{name}_pair_ledger.tsv.gz',sep='\t',index=False)
  b[['raw_row']+PAIR+['BC_TEVp','BC_TEVs','min','AUC','below20DNA','removed_outlier','diverging_pair','above_positive_control']].to_csv(O/f'{name}_barcode_ledger.tsv.gz',sep='\t',index=False)
  comp.to_csv(O/f'{name}_outlier_comparisons.tsv.gz',sep='\t',index=False)
 return final,ledger,b,meta

def load_barcodes(path, digits):
 if path.suffix=='.zip':
  with zipfile.ZipFile(path) as z:d=pd.read_csv(z.open('Data_DMS_unmapped.txt'),sep='\t',dtype={'TEVp':str})
 else:d=pd.read_csv(path,sep='\t',dtype={'TEVp':str})
 d['raw_row']=np.arange(len(d));return measures(d,auc_digits=digits)

def main():
 info={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}
 d=load_barcodes(S/'Data_DMS_unmapped.txt.zip',4)
 # Validate total pooling against deposited mapped Data_DMS input (before correction).
 with zipfile.ZipFile(DATA/'protease_DMS.zip') as z: raw=pd.read_csv(z.open('Data_DMS.txt'),sep='\t')
 pooled=pool(d); match=raw.merge(pooled,on=PAIR,suffixes=('_deposited','_pooled'),validate='one_to_one')
 info['raw_pool_check']={'unmapped_rows':len(d),'all_pooled_pairs':len(pooled),'deposited_pairs':len(raw),'matched':len(match),'count_mismatches':int((match[[x+'_deposited' for x in C]].to_numpy()!=match[[x+'_pooled' for x in C]].to_numpy()).any(axis=1).sum()),'max_AUC_difference':float(abs(match.AUC_deposited-match.AUC_pooled).max())}
 # R group_by string ordering is C locale under modern dplyr; count tied boundary explicitly.
 eligible=d[(d['min']>=20)&(d.BC_TEVs!='ACGTCGCTGA')]
 med=eligible.groupby('BC_TEVs',sort=True).AUC.agg(['median','size']).reset_index().sort_values('median',ascending=False,kind='stable')
 mapping=d[['BC_TEVs','TEVs']].drop_duplicates();assert not mapping.BC_TEVs.duplicated().any()
 med=med.merge(mapping,on='BC_TEVs',sort=False);med=med[~med.TEVs.str.contains(r'\*')]
 med.to_csv(O/'backbone_medians.tsv',sep='\t',index=False)
 assert med['median'].iloc[10]>med['median'].iloc[11], 'Manual top11 boundary tied'
 badhigh=med.BC_TEVs.iloc[:11].tolist()
 manual=[x.strip() for x in (S/'TEVs_BCs_suspicious.txt').read_text().splitlines() if x.strip()]
 correctbc=['AAGTCGCTGA','ACGTCGCTTA','ACTTCGCTGA','CCCCCCCCTG','GAAAAGGTTA','GCCCCTCCCC','GGGGGTGGC','CCGCCCCACG','TCGTGGGATC']
 badmanual=[x for x in manual if x not in correctbc]+['GAATAAGAAA']
 bad=pd.DataFrame([(x,'author_high_median_top11') for x in badhigh]+[(x,'author_manual_list_minus_confirmed_correct') for x in badmanual],columns=['BC_TEVs','reason'])
 bad.to_csv(O/'mutated_backbones.tsv',sep='\t',index=False)
 d['mutated_backbone']=d.BC_TEVs.isin(bad.BC_TEVs)
 d.loc[d.mutated_backbone,['raw_row']+PAIR+['BC_TEVp','BC_TEVs','min','AUC']].merge(bad,on='BC_TEVs').to_csv(O/'DMS_mutated_backbone_rows.tsv.gz',sep='\t',index=False)
 info['backbones']={'listed':len(bad),'unique':bad.BC_TEVs.nunique(),'raw_rows_removed':int(d.mutated_backbone.sum()),'top11_boundary_medians':med['median'].iloc[10:12].tolist()}
 ic=pd.read_csv(S/'Data_DMS_internal_controls.txt',sep='\t');inoise=noise(ic);info['internal_noise']=inoise
 c=load_barcodes(S/'Data_DMS_controls_unmapped.txt',3)
 cp=pool(c);rc=pd.read_csv(DATA/'protease_controls.tsv',sep='\t',dtype={'TEVp':str})
 check=rc.merge(cp,on=PAIR,suffixes=('_deposited','_pooled'),validate='one_to_one')
 info['raw_control_pool_check']={'unmapped_rows':len(c),'pooled_pairs':len(cp),'deposited_pairs':len(rc),'matched':len(check),'count_mismatches':int((check[[x+'_deposited' for x in C]].to_numpy()!=check[[x+'_pooled' for x in C]].to_numpy()).any(axis=1).sum()),'max_AUC_difference':float(abs(check.AUC_deposited-check.AUC_pooled).max())}
 c['mutated_backbone']=c.BC_TEVs.isin(bad.BC_TEVs)
 c.loc[c.mutated_backbone,['raw_row']+PAIR+['BC_TEVp','BC_TEVs','min','AUC']].merge(bad,on='BC_TEVs').to_csv(O/'controls_mutated_backbone_rows.tsv.gz',sep='\t',index=False)
 info['control_backbone_removed']=int(c.mutated_backbone.sum())
 finalc,_,_,info['controls']=correct(c[~c.mutated_backbone],inoise['cutoff'],'0 C151A no ssrA','controls')
 fit=noise(finalc[(finalc.n_mut_TEVs<=1)&(finalc.TEVs!='ENLYFQ*')]);info['final_noise']=fit
 finald,_,_,info['DMS']=correct(d[~d.mutated_backbone],fit['cutoff'],'C151A','DMS')
 # Published blanking; separate raw corrected AUC remains primary.
 for name,frame in [('DMS',finald),('controls',finalc)]:
  frame['dnAUC']=np.where(frame.AUC<=fit['cutoff'],0,frame.AUC-fit['mu'])
  frame.to_csv(O/f'{name}_corrected.tsv.gz',sep='\t',index=False)
 singles=finald[(finald.n_mut_TEVp<=1)&(finald.n_mut_TEVs<=1)]
 singles.to_csv(O/'DMS_script05_singles.tsv.gz',sep='\t',index=False)
 # Figure2 count definitions include stops on enzyme axis but not substrate.
 pos=singles.TEVp.str.extract(r'^.([0-9]+).$',expand=False).astype(float)
 fig=singles[(pos.between(0,234)|singles.TEVp.eq('wt'))&~singles.TEVs.str.contains(r'\*')]
 info['script05_counts']={'all_singles_including_stops':len(singles),'fig2_pairs':len(fig),'fig2_enzymes':fig.TEVp.nunique(),'fig2_substrates':fig.TEVs.nunique()}
 (O/'reconstruction.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2))
if __name__=='__main__':main()
