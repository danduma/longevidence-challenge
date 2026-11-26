"""Extract Table 1 of the 2026 cyanobacterial kinetic primary paper.
Preserves source strings and separates reported standard errors from 2019 SDs.
"""
from html.parser import HTMLParser
from pathlib import Path
import csv,re
import pandas as pd
class Tables(HTMLParser):
 def __init__(self):super().__init__();self.tables=[];self.table=None;self.row=None;self.cell=None
 def handle_starttag(self,t,a):
  if t=='table':self.table=[]
  elif t=='tr' and self.table is not None:self.row=[]
  elif t in ['td','th'] and self.row is not None:self.cell=[]
 def handle_endtag(self,t):
  if t in ['td','th'] and self.cell is not None:self.row.append(' '.join(''.join(self.cell).split()));self.cell=None
  elif t=='tr' and self.row is not None:self.table.append(self.row);self.row=None
  elif t=='table' and self.table is not None:self.tables.append(self.table);self.table=None
 def handle_data(self,s):
  if self.cell is not None:self.cell.append(s)
p=Tables();p.feed(Path('sources/rubisco_cyano2026_pmc.html').read_text());table=p.tables[0]
assert table[0][0]=='Species' and 'Sc/o' in table[0][1] and len(table)==7
with open('data/rubisco_cyano2026_table_0.csv','w') as f:csv.writer(f).writerows(table)
x=pd.read_csv('data/rubisco_cyano2026_table_0.csv')
for c,n in [(x.columns[1],'S'),(x.columns[5],'vC')]:
 x[n]=x[c].map(lambda s:float(re.search(r'^([0-9.]+)',s).group(1)))
 x[n+'_SE']=x[c].map(lambda s:float(re.search(r'±\s*([0-9.]+)',s).group(1)))
x['source_doi']='10.1093/plphys/kiag259';x['temp_C']=25
x[['Species','S','S_SE','vC','vC_SE','source_doi','temp_C']].to_csv('data/rubisco_cyano2026_kinetics.csv',index=False)
print('Extracted five cyanobacteria + wheat control; errors are SE, not SD.')
