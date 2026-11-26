"""Endpoint-specific provenance and control ledger; primary audits never overwrite raw inputs."""
from pathlib import Path
import pandas as pd,numpy as np,json,re,hashlib
R=Path(__file__).resolve().parents[1];S=R/'raw/shared_inputs'
m=pd.read_csv(R/'derived/historical_measurements.csv')
over=pd.read_csv(S/'results/rubisco_bibliographic_overlay.csv').drop_duplicates('identifier').set_index('identifier')
# Preserve identifiers, with a separate corrected DOI.
m['source_doi']=m.identifier.map(over.reviewed_correct_doi).fillna('unknown_or_unresolved_here')
m.loc[m.study.eq('Kubien 2008'),'source_doi']='10.1093/jxb/erm283'
m['origin_source']=m.study;m['origin_status']='compilation_asserted_primary_not_independently_verified'
m['strain_cultivar']='unknown';m['pH_endpoint']=m.pH_compiled;m['pH_evidence']='compilation_field_may_mix_endpoints'
m['gas_convention']='unknown';m['active_site_normalization']='unknown';m['n_range']='unknown'
m['primary_numeric_status']='not_rechecked';m['uncertainty_usable']='no_unverified_error_kind';m['paired_challenge_controls']='no_documented_pair'
m['control_role']='species_bridge_only';m['same_preparation_cross_study']='not_established'
# Detailed primary condition audit (values are kept as deposited).
a=m.study.eq('Kubien 2008');m.loc[a,'error_kind']='SE';m.loc[a,'n_range']=np.where(m.loc[a,'trait'].eq('S'),'3-5 repeat assays on purified sample','5-7 leaf samples')
m.loc[a,'pH_endpoint']=np.where(m.loc[a,'trait'].eq('S'),'8.3 method inheritance (Kane 1994); not restated in 2008 specificity paragraph','8.0 activation/extract; assay via Kubien2003')
m.loc[a,'pH_evidence']='primary methods plus compilation distinction';m.loc[a,'gas_convention']='S: 500 ppm CO2 in O2; dissolved concentrations; solubilities CO2=0.0334 and O2=0.00126 mol/(L bar) in Fig1 caption'
m.loc[a,'active_site_normalization']='vC: 14C-CABP catalytic sites, acid-stable 14CO2 fixation; S: 3H-RuBP product HPLC ratio'
m.loc[a,'primary_numeric_status']='primary text verified S means for maize/spinach; figure numeric vC remains compilation, PDF access incomplete'
m.loc[a & m.species_key.eq('Zea mays'),'strain_cultivar']='not stated in accessible 2008 methods'
m.loc[a,'origin_status']='primary methods and figure caption read; per-value table audit status retained'
# Whitney Figure 1a was visually inspected in the embedded image extracted from the primary PDF.
a=m.study.eq('Whitney 2001');m.loc[a,'primary_numeric_status']='verified against primary Fig1a image';m.loc[a,'error_kind']='unknown: Fig1a plus/minus kind not defined; do not inherit Fig6 SD'
m.loc[a,'pH_endpoint']=np.where(m.loc[a,'trait'].eq('S'),'8.3','unknown assay pH; extraction8.0; activation25C')
m.loc[a,'pH_evidence']='primary methods p544';m.loc[a,'gas_convention']='S: Kane1994 cited, exact gases/solubility inherited not rechecked here'
m.loc[a,'active_site_normalization']='vC: 14C-CABP stoichiometric active sites in extract; S on purified enzyme'
m.loc[a & m.species_key.eq('Griffithsia monilis'),'strain_cultivar']='wild intertidal Sydney Harbour collection; no culture identifier'
m.loc[a & m.species_key.eq('Galdieria sulfuraria'),'strain_cultivar']='UTEX2393; Allen strain107.79 (G. sulphuraria spelling in primary)'
m.loc[a & m.species_key.eq('Phaeodactylum tricornutum'),'strain_cultivar']='CS-29c'
a=m.study.eq('Whitney 2001') & m.species_key.eq('Nicotiana tabacum');m.loc[a,'origin_source']='Whitney 1999';m.loc[a,'origin_status']='Figure1 caption explicitly imports tobacco comparator from Whitney1999';m.loc[a,'control_role']='imported_comparator_not_contemporaneous_control'
# Uemura source resolved correctly, but metadata-only Elsevier API is not full text.
a=m.study.eq('Uemura 1997');m.loc[a,'primary_numeric_status']='abstract verifies rounded S=238 for G.partita only; detailed table not accessed';m.loc[a,'error_kind']='unknown (_SD compilation not certified)';m.loc[a,'origin_status']='original full numerical table inaccessible'
m.loc[a & m.species_key.eq('Galdieria partita Tokara'),'strain_cultivar']='Tokara from workbook; primary strain field unverified'
a=m.study.eq('Uemura 1997') & m.species_key.eq('Porphyridium purpureum') & m.trait.eq('vC');m.loc[a,'origin_source']='Read 1994';m.loc[a,'origin_status']='workbook note assigns turnover/Kc to Read1994; not a paired Uemura1997 turnover assay'
# Cultivar correction from primary 2016 text, applies where compilation omitted it.
cv={'Triticum aestivum':'Cajeme','Glycine max':'Ransom','Spinacia oleracea':'Butterfly','Phaseolus vulgaris':'Contender','Manihot esculenta':'unspecified'}
a=m.study.eq('Hermida-Carrera 2016');m.loc[a,'strain_cultivar']=m.loc[a,'species_key'].map(cv).fillna('see primary plant-material list');m.loc[a,'error_kind']='SE from primary TableI';m.loc[a,'n_range']='3-9 TableI; per-endpoint count not identified here'
# New table cohort, preserving imported cells and source error ambiguity.
z=pd.read_csv(R/'derived/cyano2026_all_parameters.csv',dtype={'value_string':str,'error_string':str})
identity={'C. aponinum':('Cyanobacterium aponinum','PCC10605'),'C. thermalis':('Chroococcidiopsis thermalis','KOMAREK1964/111; 2023 CCALA048'),'G. violaceus':('Gloeobacter violaceus','PCC7421'),'O. acuminata':('Oscillatoria acuminata','PCC6304'),'Synechococcus sp.':('Synechococcus sp.','PCC6301; 2023 UTEX625/CCALA188'),'T. aestivum':('Triticum aestivum','unknown cultivar')}
new=[]
for i,r in z.iterrows():
 if r['column'] not in ['Sc/o (mol·mol−1)','kcatc (s−1)']:continue
 trait='S' if r['column'].startswith('Sc/o') else 'vC';species,strain=identity[r.species_short]
 reused=bool(r.reused_italic)
 new.append({'measurement_id':f'kiag259:{r.species_short}:{trait}','excel_row':None,'identifier':r.species_short,'species_raw':r.species_short,'species_key':species,'study':'Aguilo-Nicolau 2026','trait':trait,'value':r.value,'unit':'mol/mol dissolved CO2/O2 convention' if trait=='S' else 's^-1 per catalytic site','compiled_error':r.error,'error_kind':'SE per Table1 caption; unresolved as reused conversion' if (r.species_short=='C. thermalis' and trait=='vC') else 'SE per Table1 caption','temp_C':25,'pH_compiled':'not applicable','heterologous_expression':'native','note':r.cell_text,'source_identifier':'10.1093/plphys/kiag259','cohort':'2026_table1_2023_reuse' if reused else '2026_table1_new','strict_shared_preparation':'not established across studies','source_doi':'10.1093/plphys/kiag259','origin_source':'Aguilo-Nicolau 2023' if reused else 'Aguilo-Nicolau 2026','origin_status':'italics plus methods specify imported kinetics' if reused else 'present-work measurement; wheat explicitly described as measured in present work in modeling section','strain_cultivar':strain,'pH_endpoint':'8.3' if trait=='S' else '8.1','pH_evidence':'2026 primary methods; for reused cells consult 2023 source','gas_convention':'99.95% O2/0.05% CO2; paper labels solubility ratio0.038 CO2/O2, orientation ambiguous' if trait=='S' else 'N2 for kcat; Kc also21% O2; pKa6.11; CO2 solubility0.034 mol/L/atm','active_site_normalization':'specificity via 3H-RuBP HPLC product partition' if trait=='S' else 'Vmax under N2 / 14C-CABP active sites; paper states8 sites/holoenzyme','n_range':'3-5 table2026; reused source3-6, per-trait count unknown' if reused else '3-5, exact per-trait n unknown','primary_numeric_status':'verified Table1 HTML incl typography','uncertainty_usable':'no_unresolved_SD_to_SE_transport' if (r.species_short=='C. thermalis' and trait=='vC') else 'reported_SE_preserved_no_covariance_or_exact_n','paired_challenge_controls':'no G.partita or maize in cohort','control_role':'contemporaneous wheat comparator' if r.species_short=='T. aestivum' else 'measured enzyme; reuse is not a new assay','same_preparation_cross_study':'not_established'})
m=pd.concat([m,pd.DataFrame(new)],ignore_index=True)
m['identity_warning']=m.species_raw.astype(str).str.contains(r'\b[A-Z]\d+[A-Z]\b|mutant',case=False,regex=True).map({True:'mutation_in_name_conflicts_with_WT_flag; quarantine_for_WT_use',False:'no_suffix_conflict_detected_not_full_identity_validation'})
m.to_csv(R/'derived/measurement_ledger.csv',index=False)
# Direct within-study comparator ratios are useful benchmark facts, not discoveries.
ctl=[]
for (study,cohort,trait),g in m.groupby(['study','cohort','trait']):
 control=None;status=None
 if study=='Aguilo-Nicolau 2026':control='Triticum aestivum';status='same_table_control; imported rows not contemporaneously assayed'
 elif study=='Kubien 2008':control='Zea mays';status='same_study distinct assays; cultivar unspecified'
 elif study=='Uemura 1997':control='Spinacia oleracea';status='compilation association; numerical primary source unavailable'
 elif study=='Whitney 2001':control='Nicotiana tabacum';status='IMPORTED comparator from Whitney1999; not a contemporaneous control'
 if control is None:continue
 # 2026 imported/new cohorts share table wheat control but retain original origin and no new-assay assertion.
 base=m[m.study.eq(study)&m.trait.eq(trait)&m.species_key.eq(control)]
 if len(base)!=1:continue
 c=base.iloc[0]
 for _,r in g.iterrows():
  if r.measurement_id==c.measurement_id:continue
  ctl.append({'target_id':r.measurement_id,'control_id':c.measurement_id,'trait':trait,'target':r.species_key,'control':control,'ratio':r.value/c.value,'source':study,'target_origin':r.origin_source,'control_origin':c.origin_source,'status':status,'ratio_error':'not calculated: covariance/identity/origin/n incomplete','certifies_challenge':'no'})
pd.DataFrame(ctl).to_csv(R/'derived/control_ledger.csv',index=False)
# Source inventory for every historical source, retaining unresolved sources and overlay explicitly.
src=pd.read_csv(S/'results/rubisco_source_pairs_linked.csv').drop_duplicates('short_ref')[['short_ref','verified_source_doi','verified_source_link','bibliographic_note']]
allrefs=pd.DataFrame({'short_ref':sorted(m.study.unique())}).merge(src,how='left',on='short_ref')
primary={'Aguilo-Nicolau 2026':('10.1093/plphys/kiag259','Table1 HTML; new3 cyano + new wheat, imported2 cyano'),'Aguilo-Nicolau 2023':('10.1007/s11120-023-01008-y','PDF + SI Table1 retrieved; 2026 correction affects photosynthesis not kinetic SI'),'dePins 2025':('10.1073/pnas.2501433122','XML, deposited library and code; 30C coupled cohort'),'Davidi correction2025':('10.1038/s44318-025-00419-y','Primary XML/table; dilution and active-site corrections distinct'),'cyano correction2026':('10.1007/s11120-025-01185-y','Primary correction concerns photosynthesis Fig2/text; no kinetic SI replacement specified')}
for label,(doi,note) in primary.items():
 row={'short_ref':label,'verified_source_doi':doi,'verified_source_link':'https://doi.org/'+doi,'bibliographic_note':note}
 if label in allrefs.short_ref.values:
  for k,v in row.items():allrefs.loc[allrefs.short_ref.eq(label),k]=v
 else:allrefs=pd.concat([allrefs,pd.DataFrame([row])],ignore_index=True)
allrefs.to_csv(R/'derived/primary_source_ledger.csv',index=False)
# Source/condition graph origin audit: don't silently treat paper grouping as source-independent.
a=m[m.study.isin(['Uemura 1997','Whitney 2001','Kubien 2008','Aguilo-Nicolau 2026'])]
a.to_csv(R/'derived/anchor_and_current_audit.csv',index=False)
print('Historical+2026 measurement rows',len(m),'control ratios',len(ctl),'audited rows',len(a))
