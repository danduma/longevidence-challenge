from pathlib import Path
import json
import openpyxl,numpy as np,pandas as pd
from reconstruct import O,S,PAIR,WT,NEG
w=openpyxl.load_workbook(S/'41467_2025_60622_MOESM4_ESM.xlsx',read_only=True,data_only=True)
def sheet(name):
 it=w[name].values;cols=next(it);return pd.DataFrame(it,columns=cols)
info={}
a=sheet('Fig. S8e');a.columns=['key','AUC_deposited','category'];a.AUC_deposited=pd.to_numeric(a.AUC_deposited,errors='coerce')
b=pd.read_csv(O/'DMS_script05_singles.tsv.gz',sep='\t');b['key']=b.TEVp+'_'+b.TEVs
j=a.merge(b[['key','AUC']],on='key',how='outer',indicator=True,validate='one_to_one')
err=abs(j.AUC-j.AUC_deposited)
info['S8e_final_singles']={'deposited':len(a),'reconstructed':len(b),'common':int((j._merge=='both').sum()),'only_deposited':int((j._merge=='left_only').sum()),'only_reconstructed':int((j._merge=='right_only').sum()),'max_AUC_error':float(err.max()),'AUC_discrepancies_over_1e12':int((err>1e-12).sum())}
j[j._merge.ne('both')|(err>1e-12)].to_csv(O/'source_S8e_discrepancies.csv',index=False)
# S7a is exactly the fitted negative-control population, with abbreviated construct names.
c=pd.read_csv(O/'controls_corrected.tsv.gz',sep='\t',dtype={'TEVp':str})
c=c[c.TEVp.isin(NEG)&(c['min']>=20)&(c.n_mut_TEVs<=1)&(c.TEVs!='ENLYFQ*')]
a=sheet('Fig. S7a');a.columns=['TEVp','AUC'];a.TEVp=a.TEVp.replace({'mCherry':'(1-13)-mCherry','TEVp I (truncated)':'I C151A 30*31*','no TEVp':'no TEVp'})
info['S7a_controls']={}
for name,g in a.groupby('TEVp'):
 x=np.sort(g.AUC.to_numpy());y=np.sort(c[c.TEVp==name].AUC.to_numpy())
 info['S7a_controls'][name]={'deposited_n':len(x),'reconstructed_n':len(y),'max_error':float(np.max(abs(x-y))) if len(x)==len(y) else None}
# Validate normalized reconstructed control data including author parent replacement.
f=pd.read_csv(O/'DMS_corrected.tsv.gz',sep='\t');c=pd.read_csv(O/'controls_corrected.tsv.gz',sep='\t',dtype={'TEVp':str})
p=f[(f.TEVp=='wt')&(f.n_mut_TEVs<=1)].copy();p.TEVp='0';c=pd.concat([c[c.TEVp!='0'],p],ignore_index=True)
den=c[c.TEVs==WT].set_index('TEVp').dnAUC;c['dnAUCrel']=c.dnAUC/c.TEVp.map(den)
a=sheet('Fig. S7c');a.columns=['TEVp','pos','aa','assay','relative'];a=a[a.assay=='DNA recorder (in vivo)'].copy();a.TEVp=a.TEVp.astype(str)
pos=['P6','P5','P4','P3','P2','P1',"P1'"]
a['TEVs']=[WT[:i]+aa+WT[i+1:] for i,aa in zip(a.pos.map(pos.index),a.aa)]
a.relative=pd.to_numeric(a.relative,errors='coerce');j=a.merge(c[PAIR+['dnAUCrel']],on=PAIR,how='left',validate='many_to_one')
both=j.relative.notna()&j.dnAUCrel.notna();err=abs(j.loc[both,'relative']-j.loc[both,'dnAUCrel'])
info['S7c_normalized_controls']={'deposited_rows':len(a),'deposited_nonmissing':int(j.relative.notna().sum()),'reconstructed_nonmissing':int(j.dnAUCrel.notna().sum()),'common_nonmissing':int(both.sum()),'max_error':float(err.max()),'discrepancies_over_1e10':int((err>1e-10).sum())}
j.to_csv(O/'source_S7c_comparison.csv',index=False)
# The deposited figure uses manually cloned parent 0, unlike the current script.
manual=pd.read_csv(O/'controls_corrected.tsv.gz',sep='\t',dtype={'TEVp':str})
den=manual[manual.TEVs==WT].set_index('TEVp').dnAUC
manual['manual_relative']=manual.dnAUC/manual.TEVp.map(den)
t=j.merge(manual[PAIR+['manual_relative']],on=PAIR,how='left',validate='many_to_one')
both=t.relative.notna()&t.manual_relative.notna();err=abs(t.loc[both,'relative']-t.loc[both,'manual_relative'])
info['S7c_manual_parent_branch']={'deposited_nonmissing':int(t.relative.notna().sum()),'reconstructed_nonmissing':int(t.manual_relative.notna().sum()),'common_nonmissing':int(both.sum()),'max_error':float(err.max()),'discrepancies_over_1e10':int((err>1e-10).sum())}
t.to_csv(O/'source_S7c_comparison.csv',index=False)
(O/'source_validation.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2))
