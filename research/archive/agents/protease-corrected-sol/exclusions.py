"""Complete upstream key-level ledger, including pairs lost before pooling."""
import pandas as pd,numpy as np,json,zipfile
from reconstruct import S,O,PAIR,WT
records=[]
for name,path in [('DMS',S/'Data_DMS_unmapped.txt.zip'),('controls',S/'Data_DMS_controls_unmapped.txt')]:
 if name=='DMS':
  with zipfile.ZipFile(path) as z:d=pd.read_csv(z.open('Data_DMS_unmapped.txt'),sep='\t',usecols=PAIR+['BC_TEVs']+[f'{a}{i}' for a in ['N','F'] for i in range(1,6)])
 else:d=pd.read_csv(path,sep='\t',dtype={'TEVp':str})
 bad=pd.read_csv(O/'mutated_backbones.tsv',sep='\t').BC_TEVs
 d['raw_barcodes']=1;d['after_backbone']=~d.BC_TEVs.isin(bad)
 n=d[[f'N{i}' for i in range(1,6)]].to_numpy();f=d[[f'F{i}' for i in range(1,6)]].to_numpy()
 d['eligible_barcodes']=d.after_backbone&((n+f).min(axis=1)>=20)
 g=d.groupby(PAIR)[['raw_barcodes','after_backbone','eligible_barcodes']].sum().reset_index()
 detail=pd.read_csv(O/f'{name}_pair_ledger.tsv.gz',sep='\t')
 g=g.merge(detail,on=PAIR,how='left',validate='one_to_one')
 g['reason']=np.select([g.after_backbone.eq(0),g.eligible_barcodes.eq(0),g.diverging_duplicate.eq(True),g.above_positive_control.eq(True),g.retained.eq(True)],['all_barcodes_author_backbone_excluded','no_barcode_depth20_after_backbone','discordant_duplicate_ge3','at_or_above_positive_control','retained'],default='UNEXPLAINED')
 assert not g.reason.eq('UNEXPLAINED').any()
 g.to_csv(O/f'{name}_all_upstream_pair_ledger.tsv.gz',sep='\t',index=False)
 records.extend(dict(dataset=name,reason=k,n=int(v)) for k,v in g.reason.value_counts().items())
pd.DataFrame(records).to_csv(O/'exclusion_waterfall.csv',index=False)
