"""Reviewed identifier corrections only; original measurements remain untouched."""
from pathlib import Path
import pandas as pd
fix={
'Uemura 1996':('10.1093/oxfordjournals.pcp.a028949','A Rapid and Sensitive Method for Determination of Relative Specificity of RuBisCO from Various Species by Anion-Exchange Chromatography','Source table cites the 1997 PMID/title for 1996 measurements.'),
'Uemura 1997':('10.1006/bbrc.1997.6497','Ribulose-1,5-bisphosphate carboxylase/oxygenase from thermophilic red algae with a strong specificity for CO2 fixation','Source DOI links to the 1996 assay-method paper.'),
'Whitney 2001':('10.1046/j.1365-313x.2001.01056.x','Form I Rubiscos from non-green algae are expressed abundantly but not assembled in tobacco chloroplasts','Source PMID26662726 resolves to Occhialini2015, not the cited2001 paper.'),
'Bainbridge 1995':('10.1093/jxb/46.special_issue.1269','Engineering Rubisco to change its catalytic properties','Source increments DOI page suffix across rows; paired wild-type row uses nonexistent1272; correct article begins1269.')}
audit=pd.read_csv('results/rubisco_citation_audit.csv')
for study,(doi,title,note) in fix.items():
 mask=audit.study.eq(study)
 audit.loc[mask,'reviewed_correct_doi']=doi;audit.loc[mask,'reviewed_correct_title']=title;audit.loc[mask,'audit_note']=note
# Do not call a link correct merely because it resolves; two links resolve to wrong titles.
audit['reviewed_link_status']=audit.apply(lambda r:'corrected_title_or_identifier_mismatch' if r.study in fix else ('resolved_title_checked_against_cited_title' if str(r.identifier_checked)=='True' else 'unresolved'),axis=1)
audit.to_csv('results/rubisco_citation_audit_reviewed.csv',index=False)
pairs=pd.read_csv('results/rubisco_source_pairs.csv')
meta=audit.set_index('study')
pairs['verified_source_doi']=pairs.short_ref.map(lambda s:fix[s][0] if s in fix else meta.loc[s,'resolved_doi'])
pairs['verified_source_link']=pairs.verified_source_doi.map(lambda s:'https://doi.org/'+str(s))
pairs['bibliographic_note']=pairs.short_ref.map(lambda s:fix[s][2] if s in fix else 'Identifier resolves to cited study; numerical source-table transcription not independently rechecked.')
pairs.to_csv('results/rubisco_source_pairs_linked.csv',index=False)
raw=pd.read_excel('data/rubisco2019_original.xlsx',sheet_name='kinetic_data')
u=raw[raw.short_ref.isin(fix)][['species','identifier','short_ref','pmid_or_doi','citation']].copy()
u['reviewed_correct_doi']=u.short_ref.map(lambda s:fix[s][0]);u['reviewed_correct_title']=u.short_ref.map(lambda s:fix[s][1]);u['audit_note']=u.short_ref.map(lambda s:fix[s][2]);u['numerical_validation']='No measurement value altered; original numerical source tables not completely rechecked.'
u.to_csv('results/rubisco_bibliographic_overlay.csv',index=False)
print('Linked source-pairs:',len(pairs),'reviewed reference records:',len(audit),'metadata overlay source rows:',len(u))
