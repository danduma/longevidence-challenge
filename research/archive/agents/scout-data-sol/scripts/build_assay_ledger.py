"""Start an assay ledger without pretending unreviewed historical methods are known."""
import pandas as pd
pairs=pd.read_csv('results/rubisco_source_pairs_linked.csv')
front=pd.read_csv('results/rubisco_pareto_source_pairs.csv')[['species','short_ref']]
edge=pairs.merge(front,on=['species','short_ref'],how='inner')
anchors=pairs[pairs.species.str.contains('Zea mays',case=False)]
old=pd.concat([edge,anchors]).drop_duplicates(['species','short_ref'])
old=old[['species','short_ref','S','vC','S_SD','vC_SD','verified_source_doi','pair_type','pH_min','pH_max']].rename(columns={'verified_source_doi':'source_doi'})
old['temp_C']=25;old['S_error_kind']='SD as reported/compiled';old['vC_error_kind']='SD as reported/compiled'
old['active_site_quantification']='requires primary-method audit';old['specificity_method']='requires primary-method audit'
old['raw_progress_curves_available']='not verified';old['paired_G_partita_and_maize_controls']='not demonstrated in these compilation rows';old['cohort']='2019 compilation anchors/frontier'
new=pd.read_csv('data/rubisco_cyano2026_kinetics.csv').rename(columns={'Species':'species'})
new['short_ref']='Aguilo-Nicolau et al2026';new['pair_type']='same study, distinct specificity/carboxylation assays'
new['pH_carboxylation']=8.1;new['pH_specificity']=8.3;new['S_error_kind']='SE';new['vC_error_kind']='SE'
new['active_site_quantification']='14C-CABP binding; paper assumes8 sites/holoenzyme';new['specificity_method']='3H-RuBP product partition/HPLC under controlled gases'
new['raw_progress_curves_available']='not verified';new['paired_G_partita_and_maize_controls']='wheat control only';new['cohort']='2026 cyanobacterial primary cohort'
ledger=pd.concat([old,new],ignore_index=True);ledger.to_csv('results/rubisco_assay_ledger.csv',index=False)
print('Assay ledger rows:',len(ledger))
