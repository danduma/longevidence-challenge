"""Verify identifiers for all source studies used in the 2019 paired audit.
Preserves original source IDs and provides a separate reviewed correction table.
"""
import re,sys,json,urllib.parse,concurrent.futures,xml.etree.ElementTree as E
from pathlib import Path
import pandas as pd
sys.path.insert(0,'scripts');from acquire import fetch
x=pd.read_csv('results/rubisco_source_pairs.csv')
refs=x[['short_ref','pmid_or_doi','citation']].drop_duplicates('short_ref')
pmids=[str(v) for v in refs.pmid_or_doi if str(v).isdigit()]
p=fetch('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id='+','.join(pmids)+'&retmode=xml','sources/rubisco_source_pubmed_metadata.xml')
lookup={}
if p:
 root=E.parse(p).getroot()
 for art in root.findall('.//PubmedArticle'):
  pmid=art.findtext('.//MedlineCitation/PMID');title=art.find('.//ArticleTitle')
  doi=next((q.text for q in art.findall('.//ArticleId') if q.get('IdType')=='doi'),None)
  lookup[pmid]={'resolved_title':''.join(title.itertext()) if title is not None else '', 'resolved_doi':doi,'identifier_checked':True}

def check(row):
 key=str(row.pmid_or_doi).strip();rec={'study':row.short_ref,'original_identifier':key,'source_citation':row.citation}
 doi=re.search(r'10\.\d{4,9}/\S+',key)
 if doi:
  path='sources/crossref_'+re.sub('[^A-Za-z0-9]','_',row.short_ref)+'.json'
  p=fetch('https://api.crossref.org/works/'+urllib.parse.quote(doi.group(0),safe=''),path)
  if p:
   j=json.loads(p.read_text()).get('message',{});rec.update(resolved_title='; '.join(j.get('title',[])),resolved_doi=j.get('DOI'),identifier_checked=bool(j.get('title')))
 elif key in lookup:rec.update(lookup[key])
 else:rec['identifier_checked']=False
 # Known title/identifier mismatch separately verified against primary publisher and PubMed.
 if row.short_ref=='Uemura 1997':
  rec.update(reviewed_correct_doi='10.1006/bbrc.1997.6497',reviewed_correct_pmid='9144578',audit_note='Original DOI resolves to the 1996 method paper; cited 1997 red-algal paper is PMID 9144578.')
 return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:result=list(pool.map(check,refs.itertuples()))
pd.DataFrame(result).to_csv('results/rubisco_citation_audit.csv',index=False)
# All source rows using the reversed Uemura identifiers, not only paired rows.
raw=pd.read_excel('data/rubisco2019_original.xlsx',sheet_name='kinetic_data')
u=raw[raw.short_ref.isin(['Uemura 1996','Uemura 1997'])][['species','identifier','short_ref','pmid_or_doi','citation']].copy()
u['reviewed_correct_doi']=u.short_ref.map({'Uemura 1996':'10.1093/oxfordjournals.pcp.a028949','Uemura 1997':'10.1006/bbrc.1997.6497'})
u['status']='Bibliographic overlay verified from publisher/PubMed; individual numeric values not checked against full original tables.'
u.to_csv('results/rubisco_Uemura_bibliographic_overlay.csv',index=False)
print('Study identifier records',len(result),'checked',sum(bool(r.get('identifier_checked')) for r in result))
