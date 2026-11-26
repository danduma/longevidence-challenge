"""Exhaustive two-strain reuse audit; manual 2023 SI transcription, parsed 2026 HTML."""
from pathlib import Path
from html.parser import HTMLParser
from decimal import Decimal
import re,json,hashlib
import pandas as pd
R=Path(__file__).resolve().parents[1];S=R/'raw/shared_inputs'
class Tables(HTMLParser):
 def __init__(self):super().__init__();self.tables=[];self.table=None;self.row=None;self.cell=None;self.italic=False
 def handle_starttag(self,t,a):
  if t=='table':self.table=[]
  elif t=='tr' and self.table is not None:self.row=[]
  elif t in ['td','th'] and self.row is not None:self.cell=[];self.italic=False
  elif t in ['i','em'] and self.cell is not None:self.italic=True
 def handle_endtag(self,t):
  if t in ['td','th'] and self.cell is not None:self.row.append({'text':' '.join(''.join(self.cell).split()),'italic':self.italic});self.cell=None
  elif t=='tr' and self.row is not None:self.table.append(self.row);self.row=None
  elif t=='table' and self.table is not None:self.tables.append(self.table);self.table=None
 def handle_data(self,s):
  if self.cell is not None:self.cell.append(s)
p=Tables();p.feed((S/'sources/rubisco_cyano2026_pmc.html').read_text());t=p.tables[0]
(R/'derived/cyano2026_table1_cells.json').write_text(json.dumps(t,indent=2,ensure_ascii=False))
# Exact printed 2023 SI Table 1 pairs, mean and SD; includes unit conversion in efficiency only.
old={'Synechococcus sp.':[('48.2','5.8'),('147.2','4.6'),('215.9','27.7'),('631.6','62.6'),('8.7','1.2'),('0.06','0.01')], 'C. thermalis':[('66.0','2.2'),('87.3','1.3'),('106.9','4.2'),('1162.6','182.4'),('9.1','0.2'),('0.105','0.002')]}
(R/'raw/cyano2023_SI_table1_transcription.json').write_text(json.dumps({'source':'raw/cyano2023_supp_primary.pdf; page 1, Supplementary Table 1','columns':['S','Kc','Kc21','Ko','vC','vC/Kc'],'error_kind':'SD','n_min':3,'n_max':6,'data':old},indent=2))
def interval(s,scale=1):
 d=Decimal(s);half=Decimal(10)**d.as_tuple().exponent/2
 return float(max(Decimal(0),d-half))*scale,float(d+half)*scale
rows=[];new=[]
for row in t[1:]:
 name=row[0]['text'];name=name.replace(' (control)','')
 for j,c in enumerate(row[1:],1):
  m=re.search(r'^([0-9.]+)\s*±s?\s*([0-9.]+)',c['text'])
  if m:new.append({'species_short':name,'column':t[0][j]['text'],'value_string':m[1],'error_string':m[2],'value':float(m[1]),'error':float(m[2]),'error_kind':'SE as 2026 table caption','reused_italic':c['italic'],'source_doi':'10.1093/plphys/kiag259','source_location':f'Table 1 row {name} column {j}','cell_text':c['text']})
 if name not in old:continue
 for k,j in enumerate([1,2,3,4,5,7]):
  a,b=old[name][k];c=row[j];m=re.search(r'^([0-9.]+)\s*±s?\s*([0-9.]+)',c['text']);u,v=m[1],m[2]
  scale=1000 if k==5 else 1
  al,ah=interval(a,scale);ul,uh=interval(u)
  bl,bh=interval(b,scale);vl,vh=interval(v)
  nlo=(bl/vh)**2;nhi=(bh/vl)**2 if vl else float('inf')
  plausible=[n for n in range(3,7) if nlo<=n<=nhi]
  rows.append({'species':name,'trait':['S','Kc','Kc21','Ko','vC','vC/Kc'][k],'old_mean_string':a,'old_error_string':b,'old_error_kind':'SD','old_units_scale_to_new':scale,'old_printed_unit':'s^-1 (denominator unit omitted in SI header)' if k==5 else ('mol/mol' if k==0 else 's^-1' if k==4 else 'uM'),'unit_scale_basis':'1000 inferred from quotient definition and Kc in uM; SI printed unit incomplete' if k==5 else 'same printed units','new_mean_string':u,'new_error_string':v,'new_error_kind':'SE (caption)','new_marked_reused':c['italic'],'means_compatible_under_rounding':bool(al<=uh and ul<=ah),'implied_n_min':nlo,'implied_n_max':nhi,'compatible_n_3_to_6':';'.join(map(str,plausible)) or 'none','compatible_n_3_to_5':';'.join(str(n) for n in plausible if n<=5) or 'none','uncertainty_transport_status':'compatible_under_rounding' if plausible else 'unresolved_incompatible_as_pure_SD_to_SE_conversion'})
a=pd.DataFrame(rows);a.to_csv(R/'results/reused2023_2026_audit.csv',index=False)
pd.DataFrame(new).to_csv(R/'derived/cyano2026_all_parameters.csv',index=False)
print(a.to_string(index=False))
