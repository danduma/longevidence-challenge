"""Extract the 18 sample CPM columns, rejecting the extra metadata header row."""
import pandas as pd
from pathlib import Path
p='data/41467_2026_78607_MOESM9_ESM.xlsx'
d=pd.read_excel(p,usecols=range(19))
assert d.shape==(139699,19)
assert d.iloc[0,0]=='Full ID'
assert 'Animal 1: Intact Epiphysis'==d.columns[1]
assert 'Animal 1: Intact Diaphysis'==d.columns[2]
assert '3dpa Epiphysis (Replicate 1)'==d.columns[7]
assert '3dpa Diaphysis (Replicate 1)'==d.columns[10]
assert '5dpa Epiphysis (Replicate 1)'==d.columns[13]
assert '5dpa Diaphysis (Replicate 1)'==d.columns[16]
d=d.iloc[1:].copy()
d.to_csv('data/regeneration2026_CPM_extracted.tsv.gz',sep='\t',index=False)
print('Extracted',len(d),'rows x 18 biological sample columns')
