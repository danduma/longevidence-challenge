#!/usr/bin/env python3
"""C07 deterministic offline audit and mapping patch; no fitting or network calls."""
import hashlib,json,re,sys,platform
from pathlib import Path
from html.parser import HTMLParser
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
import sklearn
from fetch import norm,select
R=Path(__file__).resolve().parent
O=R/'results';O.mkdir(exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def load(p):return json.loads((R/p).read_text())
class Text(HTMLParser):
 def __init__(self):super().__init__();self.data=[]
 def handle_data(self,s):self.data.append(s)
def htmltext(path):
 p=Text();p.feed((R/path).read_text());return ' '.join(' '.join(p.data).split())
def fasta(path):
 out={};lab=None
 for line in (R/path).read_text().splitlines():
  if line.startswith('>'):lab=line[1:].split()[0];out.setdefault(lab,'')
  else:
   # repeated labels are separately checked using metadata lengths
   out[lab]+=line.strip()
 return out
# FASTA repeated records must not concatenate.
def fasta_records(path):
 out=[]
 for block in (R/path).read_text().split('>')[1:]:
  lines=block.splitlines();out.append((lines[0].split()[0],''.join(lines[1:])))
 return out

def run():
 for rec in load('inputs/manifest.json'):assert sha((R/rec['file']).read_bytes())==rec['sha256']
 m0=pd.read_csv(R/'inputs/meta.csv');m=m0.drop_duplicates().copy();m['key']=m.species.map(norm)
 records=fasta_records('inputs/sequences.fasta');seq={}
 assert len(records)==len(m0)
 for (label,s),(_,row) in zip(records,m0.iterrows()):
  assert label==row.label and len(s)==row.seq_len
  assert label not in seq or seq[label]==s
  seq[label]=s
 m['sequence']=m.label.map(seq);sus=select(m)
 assert len(sus)==60 and len(m)==207
 a=pd.read_csv(R/'inputs/assays.csv');a['key']=a.species.map(norm);a['source_row']=np.arange(len(a));wt=a[(a['primary']==1)&(a.mutant==0)].copy()
 crop=htmltext('evidence/crop2016.html');rp=htmltext('evidence/palustris2014.html')
 assert 'Avena sativa ‘Forridena’' in crop
 assert 'Avena sativa 10.8 ± 0.9 18.1 ± 2.0 2.3 ± 0.3 99.9 ± 3.0' in crop
 avrows=wt[wt.key==norm('Avena sativa')]
 assert len(avrows)==2 and avrows.KC.dropna().tolist()==[10.8] and avrows.vC.dropna().tolist()==[2.3] and avrows.S.dropna().tolist()==[99.9]
 rprow=wt[wt.key==norm('Rhodopseudomonas palustris')]
 assert len(rprow)==1 and tuple(rprow.iloc[0][['S','vC','KC','KO','KRuBP']])==(12.0,4.9,104.0,277.0,10.0)
 assert 'R. palustris strain CGA009' in rp
 assert 'Wild type 12 ± 1 4.9 ± 0.2 104 ± 22 277 ± 15' in rp
 old=load('evidence/uniprot/Q07N61.json');new=load('evidence/uniprot/Q6N0W9.json');pdb=load('evidence/pdb_4LF1.json')
 assert 'BisA53' in old['organism']['scientificName'] and 'CGA009' in new['organism']['scientificName']
 assert pdb['entity_src_gen'][0]['gene_src_strain']=='CGA009'
 assert pdb['entity_poly']['pdbx_seq_one_letter_code_can']=='MGSSHHHHHHSSGLVPRGSH'+new['sequence']['value']
 assert old['sequence']['value']==m.set_index('species').loc['Rhodopseudomonas palustris','sequence']
 assert old['sequence']['value']!=new['sequence']['value']
 av=(R/'evidence/unisave_A0A510B6V9_v17.txt').read_text();assert 'OS   Avena nuda' in av and 'NCBI_TaxID=4497' in av
 ledger=[];patch=[]
 for _,r in sus.iterrows():
  d=load('evidence/uniprot/'+r.uniprot_accession+'.json');inactive=d.get('entryType')=='Inactive';tax=d.get('organism',{}).get('scientificName','');taxid=d.get('organism',{}).get('taxonId','');ev=d.get('entryAudit',{}).get('entryVersion','');sv=d.get('entryAudit',{}).get('sequenceVersion','');snapshot='evidence/uniprot/'+r.uniprot_accession+'.json';hist=''
  if inactive:
   fs=sorted((R/'evidence/history').glob(r.uniprot_accession+'_v*.txt'))
   if fs:
    hist=str(fs[0].relative_to(R));s=fs[0].read_text();m_os=re.search(r'^OS   (.+)$',s,re.M);tax=' '.join(re.findall(r'^OS   (.+)$',s,re.M)).rstrip('.');ids=re.search(r'NCBI_TaxID=(\d+)',s);taxid=ids.group(1) if ids else '';v=re.search(r'entry version (\d+)',s);ev=v.group(1) if v else '';v=re.search(r'sequence version (\d+)',s);sv=v.group(1) if v else ''
  if inactive and hist:
   historical_seq=''.join(re.findall(r'[A-Z]+',s.split('SQ   SEQUENCE',1)[1].split('\n',1)[1].split('//',1)[0]))
   assert historical_seq==r.sequence, (r.species,'historical sequence differs')
  if not inactive:
   assert d['sequence']['value']==r.sequence, r.species
  refs=[c.get('value','') for ref in d.get('references',[]) for c in ref.get('referenceComments',[]) if c.get('type')=='STRAIN']
  reason=[]
  if r.seq_len<400:reason.append('length_lt400')
  if norm(r.species)!=norm(r.taxon_name):reason.append('literal_taxon_difference')
  if re.search(r'cv\.|cultivar|strain|sp\.|subsp\.|var\.|PCC|\d',r.species):reason.append('literal_qualifier')
  action='retain_unresolved';taxonomy='literal_difference_unresolved' if norm(r.species)!=norm(r.taxon_name) else 'name_agrees';assay='primary_assay_not_independently_resolved';why='No error inferred from name difference alone; retain for exploratory species linkage, not preparation-certified.';newacc='';newevidence='';assayev=''
  if r.seq_len<400:
   action='quarantine';why='Operational full-sequence threshold: deposited RbcL <400 aa. This is not a new fragment discovery.'
   if r.species in ['Capsicum annuum','Medicago sativa']:assay='crop_species_and_cultivar_confirmed; exact_sequence_bridge_unresolved';assayev='evidence/crop2016.html#MATERIALS-AND-METHODS'
  if r.species=='Avena sativa':
   action='quarantine';taxonomy='archived_accession_is_Avena_nuda';assay='Avena_sativa_Forridena_confirmed; accession_bridge_unresolved';why='Original crop assay names Avena sativa Forridena; UniSave v17 names Avena nuda taxid4497. No documented Forridena-to-accession bridge; synonymy alone cannot certify cultivar.';assayev='evidence/crop2016.html#MATERIALS-AND-METHODS';newevidence='evidence/unisave_A0A510B6V9_v17.txt'
  if r.species=='Rhodopseudomonas palustris':
   action='replace_native_core';taxonomy='species_agrees_strain_conflicts';assay='CGA009_cbbM_and_WT_Table2_confirmed';why='BisA53 is not the assayed CGA009 strain. Q6N0W9 native core matches 4LF1 after removing its documented 20-aa N-terminal tag. Repair applies to native encoded core, not an assertion of identical tag status for every assay.';newacc='Q6N0W9';assayev='evidence/palustris2014.html#EXPERIMENTAL-PROCEDURES;TABLE2';newevidence='evidence/pdb_4LF1.json; evidence/uniprot/Q6N0W9.json'
  if r.species=='Arabidopsis thaliana cv. Columbia':
   action='retain_reference_supported';taxonomy='cultivar_supported_in_accession_reference';assay='source_table_Columbia; exact_preparation_not_verified';why='O03042 reference 2 has cv. Columbia nucleotide-sequence support; do not quarantine merely because organism field omits cultivar.'
  if r.species.startswith('Nicotiana tabacum'):
   action='retain_reference_supported';taxonomy='cultivar_supported_in_accession_reference';assay='source_table_W38; exact_preparation_not_verified';why='P00876 reference 7 identifies cv. Wisconsin-38; released W38 label is reference-supported, not a proven distinct mismatch.'
  if r.species=='Glycine max cv. Wayne':
   action='quarantine';taxonomy='species_agrees_Wayne_not_supported_in_accession';assay='source_table_Wayne; primary_fulltext_identity_unresolved';why='P27066 sequence references do not identify Wayne (one names PI437654). This does not prove sequence difference. Quarantine strict cultivar link until an exact bridge is documented.';assayev='inputs/assays.csv#max_wayne_laing'
  if r.species=='Flaveria australasica':
   taxonomy='possible_synonym_or_sample_conflict_unresolved';why='E9N7A7 names trinervia and links Kapralov2011; Kubien2008 distinguishes australasica and trinervia assay labels. Taxonomic synonym reports conflict with assay sample distinction; no repair without specimen/sequence bridge.';assay='Kubien2008_abstract_distinguishes_labels; exact_sample_unresolved';assayev='https://pubmed.ncbi.nlm.nih.gov/18227079/'
  if r.species=='P. dilatatum':action='prior_known_excluded';why='Established Paspalum/Orthomnion mapping error; excluded from all new patch counts.'
  linked=wt[wt.key==r.key]
  z={'species':r.species,'old_accession':r.uniprot_accession,'selection_reasons':';'.join(reason),'released_taxon':r.taxon_name,'actual_accession_taxon':tax,'actual_taxid':taxid,'active_status':d['entryType'],'entry_version':ev,'sequence_version':sv,'sequence_length':int(r.seq_len),'accession_strain_reference_comments':';'.join(refs),'taxonomy_verdict':taxonomy,'assay_verdict':assay,'action':action,'rationale':why,'new_accession':newacc,'evidence_accession':snapshot,'evidence_history':hist,'evidence_primary_assay':assayev,'evidence_replacement':newevidence,'wt_source_records':len(linked),'source_identifiers':';'.join(linked.identifier.astype(str))}
  ledger.append(z)
  if action in ['quarantine','replace_native_core']:patch.append(z)
 pd.DataFrame(ledger).to_csv(O/'link_audit.csv',index=False)
 pd.DataFrame(patch).to_csv(O/'mapping_patch_v1.csv',index=False)
 (O/'mapping_patch_v1.json').write_text(json.dumps({'contribution':'C07','patch_version':'1.0.0','protocol':'protocol.md','operations':patch},indent=2))
 quarantine={norm(x['species']) for x in patch if x['action']=='quarantine'}
 corrected=m.copy();corrected['mapping_status']='not_audited'
 for z in ledger:corrected.loc[corrected.key==norm(z['species']),'mapping_status']=z['action']
 ix=corrected.species=='Rhodopseudomonas palustris';corrected.loc[ix,'uniprot_accession']='Q6N0W9';corrected.loc[ix,'sequence']=new['sequence']['value'];corrected.loc[ix,'seq_len']=len(new['sequence']['value']);corrected.loc[ix,'uniprot_length']=len(new['sequence']['value']);corrected.loc[ix,'taxon_name']=new['organism']['scientificName'];corrected.loc[ix,'taxon_id']=new['organism']['taxonId']
 # Derived sequence patch only; full corrected mapping is local-only because upstream data terms differ.
 corrected.drop(columns=['sequence']).to_csv(O/'mapping_with_status_v1.csv',index=False)
 (O/'replacement_native_core_v1.fasta').write_text('>Rhodopseudomonas_palustris_CGA009_Q6N0W9_native_core\n'+new['sequence']['value']+'\n')
 # Prior error is removed in both primary comparison arms, as required.
 prior={norm('P. dilatatum')};basekeys=set(m.key)-prior;afterkeys=basekeys-quarantine;additionalkeys=basekeys-{norm('Avena sativa'),norm('Glycine max cv. Wayne')}
 source=wt[wt.key.isin(set(m.key))].copy();source['prior_known']=source.key.isin(prior);source['c07_additional_only_eligible']=source.key.isin(additionalkeys);source['c07_eligible']=source.key.isin(afterkeys);source.to_csv(O/'source_record_eligibility.csv',index=False)
 summary={'contribution':'C07','selected_links':len(sus),'selected_accessions':len(set(sus.uniprot_accession)),'selected_inactive_accessions':len({x['old_accession'] for x in ledger if x['active_status']=='Inactive'}),'new_replaced_links':sum(x['action']=='replace_native_core' for x in patch),'patch_quarantined_links':len(quarantine),'fragment_quarantines_prior_recognized':7,'new_nonfragment_quarantines':len(quarantine)-7,'baseline_links_excluding_prior_error':len(basekeys),'additional_only_links':len(additionalkeys),'additional_only_wt_source_records':int(wt.key.isin(additionalkeys).sum()),'after_links':len(afterkeys),'baseline_wt_source_records':int(wt.key.isin(basekeys).sum()),'after_wt_source_records':int(wt.key.isin(afterkeys).sum()),'repaired_wt_source_records':int((wt.key==norm('Rhodopseudomonas palustris')).sum())}
 # Fixed folds from full released prediction order; never change folds after exclusions.
 impacts=[];exposures=[]
 def groups(mm,keys):
  tab=mm.set_index('key');parent={k:k for k in keys}
  def root(k):
   while parent[k]!=k:k=parent[k]
   return k
  for i,k in enumerate(sorted(keys)):
   for l in sorted(keys)[:i]:
    if ' '.join(k.split()[:2])==' '.join(l.split()[:2]) or tab.loc[k,'sequence']==tab.loc[l,'sequence']:parent[root(k)]=root(l)
  return {k:root(k) for k in keys}
 gb=groups(m,basekeys);ga=groups(corrected,afterkeys);gnew=groups(corrected,additionalkeys)
 for target in ['S','KC','vC','KO','vO','KRuBP']:
  p=pd.read_csv(R/'inputs'/('pred_'+target+'.csv'));p['key']=p.species.map(norm);p['released_row']=np.arange(len(p));fold=np.zeros(len(p),dtype=int)
  for f,(_,te) in enumerate(KFold(n_splits=5,shuffle=True,random_state=0).split(p)):fold[te]=f
  p['fold']=fold
  for arm,keys,mm,gg in [('baseline_excluding_prior',basekeys,m,gb),('C07_additional_only',additionalkeys,corrected,gnew),('C07_patch',afterkeys,corrected,ga)]:
   q=p[p.key.isin(keys)].copy();tab=mm.set_index('key');q['sequence_hash']=q.key.map(lambda k:sha(tab.loc[k,'sequence'].encode()));q['group']=q.key.map(gg)
   for col in ['key','sequence_hash','group']:
    count=0
    for _,row in q.iterrows():
     exposed=bool(((q[col]==row[col])&(q.fold!=row.fold)).any());count+=exposed
     exposures.append({'target':target,'arm':arm,'released_row':row.released_row,'species':row.species,'fold':int(row.fold),'convention':col,'exposed':exposed})
    impacts.append({'target':target,'arm':arm,'convention':col,'prediction_records':len(q),'exposed_records':count,'fraction':count/len(q) if len(q) else None})
  p['c07_eligible']=p.key.isin(afterkeys);p.to_csv(O/('fixed_fold_'+target+'.csv'),index=False)
 pd.DataFrame(impacts).to_csv(O/'exposure_impact.csv',index=False);pd.DataFrame(exposures).to_csv(O/'exposure_records.csv',index=False)
 (O/'summary.json').write_text(json.dumps(summary,indent=2));(O/'environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__},indent=2))
 print(json.dumps(summary,indent=2));return summary
if __name__=='__main__':run()
