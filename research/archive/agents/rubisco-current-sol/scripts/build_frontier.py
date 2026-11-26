"""Within-cohort mean dominance only; no invented cross-source assay corrections."""
from pathlib import Path
import pandas as pd,numpy as np,json
R=Path(__file__).resolve().parents[1]
m=pd.read_csv(R/'derived/measurement_ledger.csv');hist=m[m.cohort.eq('historical2019')]
pairs=[]
for row,g in hist.groupby('excel_row'):
 if set(g.trait)!=set(['S','vC']):continue
 a=g[g.trait.eq('S')].iloc[0];b=g[g.trait.eq('vC')].iloc[0]
 pairs.append({'pair_id':f'2019:row{int(row)}','species':a.species_key,'strain_cultivar':a.strain_cultivar,'study':a.study,'S':a.value,'vC':b.value,'S_error':a.compiled_error,'vC_error':b.compiled_error,'S_error_kind':a.error_kind,'vC_error_kind':b.error_kind,'S_id':a.measurement_id,'vC_id':b.measurement_id,'S_origin':a.origin_source,'vC_origin':b.origin_source,'temp_C':25,'pH_S':a.pH_endpoint,'pH_vC':b.pH_endpoint,'source_evidence':a.primary_numeric_status+'; '+b.primary_numeric_status,'cohort_type':'historical_same_original_row','same_assay_origin':a.origin_source==b.origin_source})
new=m[m.study.eq('Aguilo-Nicolau 2026')]
for name,g in new.groupby('species_key'):
 a=g[g.trait.eq('S')].iloc[0];b=g[g.trait.eq('vC')].iloc[0]
 pairs.append({'pair_id':f'kiag259:{name}','species':name,'strain_cultivar':a.strain_cultivar,'study':a.study,'S':a.value,'vC':b.value,'S_error':a.compiled_error,'vC_error':b.compiled_error,'S_error_kind':a.error_kind,'vC_error_kind':b.error_kind,'S_id':a.measurement_id,'vC_id':b.measurement_id,'S_origin':a.origin_source,'vC_origin':b.origin_source,'temp_C':25,'pH_S':a.pH_endpoint,'pH_vC':b.pH_endpoint,'source_evidence':a.primary_numeric_status+'; '+b.primary_numeric_status,'cohort_type':'2023_imported' if 'reuse' in a.cohort else '2026_new','same_assay_origin':True})
p=pd.DataFrame(pairs)
p['frontier_cohort']=p.S_origin.astype(str)+' | vC_origin='+p.vC_origin.astype(str)+' | pH_S='+p.pH_S.astype(str)+' | pH_vC='+p.pH_vC.astype(str)
p['within_cohort_mean_frontier']=False;p['within_cohort_n']=0
for _,g in p.groupby('frontier_cohort'):
 a=g[['S','vC']].to_numpy()
 for j,i in enumerate(g.index):
  dom=(a>=a[j]).all(axis=1)&(a>a[j]).any(axis=1)
  p.loc[i,'within_cohort_mean_frontier']=not dom.any();p.loc[i,'within_cohort_n']=len(g)
p['historical_threshold_diagnostic']=p.S.ge(238.1)&p.vC.ge(4.05)
p['wt_identity_eligible']=~p.species.str.contains(r'\b[A-Z]\d+[A-Z]\b|mutant',case=False,regex=True)
p['challenge_label']='not_identified_without_paired_G.partita_and_maize_controls'
p.to_csv(R/'derived/assay_frontier_benchmark.csv',index=False)
s={'same_original_row_historical_pairs':int(p.cohort_type.eq('historical_same_original_row').sum()),'appended_new2026_pairs_including_wheat':int(p.cohort_type.eq('2026_new').sum()),'reused2023_pairs':int(p.cohort_type.eq('2023_imported').sum()),'total_pairs':len(p),'mutation_name_conflict_pairs':int((~p.wt_identity_eligible).sum()),'joint_historical_diagnostic_positives':int(p.historical_threshold_diagnostic.sum()),'certified_challenge_labels':0,'warning':'A zero historical-threshold diagnostic is not a measured challenge negative. Frontier membership is within source/endpoint condition only; singleton cohorts are uninformative.'}
(R/'results/frontier_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
