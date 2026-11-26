"""Source-aware audit; all rates remain separated by study/temperature.
Run from the job directory: .venv/bin/python scripts/analyze_rubisco.py
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
OUT=Path('results');OUT.mkdir(exist_ok=True)
rng=np.random.default_rng(20261126)
raw=pd.read_excel('data/rubisco2019_original.xlsx',sheet_name='kinetic_data')
for c in ['S','vC','S_SD','vC_SD','temp_C','pH']:raw[c]=pd.to_numeric(raw[c],errors='coerce')
raw['species']=raw.species.str.strip()
f=raw[(raw['use']==True)&(raw.primary==1)&(raw.mutant==0)&(raw.temp_C==25)&(raw.isoform.astype(str)=='1')].copy()
strict=f[(f.S>0)&(f.vC>0)].copy()
# One source may report separate specificity/turnover rows at different pH.
# Keeping species + source preserves more provenance than a global median.
g=f.groupby(['species','short_ref'],dropna=False)
pairs=g.agg(S=('S','median'),vC=('vC','median'),S_SD=('S_SD','median'),vC_SD=('vC_SD','median'),pH_min=('pH','min'),pH_max=('pH','max'),taxonomy=('taxonomy','first'),pmid_or_doi=('pmid_or_doi','first'),citation=('citation','first')).reset_index()
pairs=pairs[(pairs.S>0)&(pairs.vC>0)].copy()
pairs['pair_type']=np.where(pd.MultiIndex.from_frame(pairs[['species','short_ref']]).isin(pd.MultiIndex.from_frame(strict[['species','short_ref']])),'same_original_row_available','same_study_separate_rows')
pairs.to_csv(OUT/'rubisco_source_pairs.csv',index=False)

def frontier(d):
 a=d[['S','vC']].to_numpy(float)
 dom=(a[:,None,:]<=a[None,:,:]).all(axis=2)&(a[:,None,:]<a[None,:,:]).any(axis=2)
 return ~dom.any(axis=1)

def corr(d):
 x=np.log(d.vC.to_numpy()); y=np.log(d.S.to_numpy())
 r,p=stats.pearsonr(x,y);rho,sp=stats.spearmanr(x,y)
 return {'n':len(d),'pearson_log_r':float(r),'pearson_p_naive':float(p),'spearman_rho':float(rho),'log_slope':float(np.polyfit(x,y,1)[0])}

def within_slope(d,block):
 x=np.log(d.vC);y=np.log(d.S)
 dx=x-x.groupby(d[block]).transform('mean');dy=y-y.groupby(d[block]).transform('mean')
 return float(np.dot(dx,dy)/np.dot(dx,dx)) if np.dot(dx,dx)>0 else float('nan')

def study_ci(d,B=3000):
 # Cluster bootstrap of centered source contributions; no phylogenetic correction.
 d=d.copy();d['x']=np.log(d.vC);d['y']=np.log(d.S)
 d['x']-=d.groupby('short_ref').x.transform('mean');d['y']-=d.groupby('short_ref').y.transform('mean')
 blocks=d.groupby('short_ref').apply(lambda z:pd.Series({'xy':np.dot(z.x,z.y),'xx':np.dot(z.x,z.x)}),include_groups=False)
 a=blocks[['xy','xx']].to_numpy(); idx=rng.integers(len(a),size=(B,len(a)));s=a[idx].sum(axis=1);v=s[:,0]/s[:,1]
 return [float(x) for x in np.quantile(v,[.025,.5,.975])]

Z=pd.get_dummies(pairs[['short_ref','taxonomy']],drop_first=True).to_numpy(float)
D=np.column_stack([np.ones(len(pairs)),Z,np.log(pairs.vC)])
combined_adjusted_slope=float(np.linalg.lstsq(D,np.log(pairs.S),rcond=None)[0][-1])

species=pairs.groupby('species',as_index=False).agg(S=('S','median'),vC=('vC','median'),taxonomy=('taxonomy','first'))
front=pairs[frontier(pairs)].sort_values('S',ascending=False)
front.to_csv(OUT/'rubisco_pareto_source_pairs.csv',index=False)
species[frontier(species)].to_csv(OUT/'rubisco_pareto_species.csv',index=False)
lo=[]
for ref in sorted(pairs.short_ref.unique()):
 d=pairs[pairs.short_ref!=ref];fp=d[frontier(d)]
 lo.append({'removed_study':ref,'n_remaining':len(d),'within_study_slope':within_slope(d,'short_ref'),'log_r':corr(d)['pearson_log_r'],'frontier_species_count':fp.species.nunique(),'new_frontier_species':'; '.join(sorted(set(fp.species)-set(front.species)))})
pd.DataFrame(lo).to_csv(OUT/'rubisco_leave_one_study.csv',index=False)
# Sensitivity: fixed +/- relative-error bounds, NOT a confidence interval.
unc=[]
a=pairs[['S','vC']].to_numpy(float)
for delta in [0,.05,.10,.20]:
 lower=a*(1-delta);upper=a*(1+delta)
 dom=(upper[:,None,:]<=lower[None,:,:]).all(axis=2)&(upper[:,None,:]<lower[None,:,:]).any(axis=2)
 possible=~dom.any(axis=1)
 unc.append({'relative_error_bound':delta,'not_certainly_dominated_source_pairs':int(possible.sum())})
# Monte Carlo for rows with both actually reported SDs, not SEM/CI estimates.
# Marginal measurement perturbation only; no assumed lab effects, no confidence claims.
u=pairs[(pairs.S_SD>0)&(pairs.vC_SD>0)].copy()
counts=np.zeros(len(u))
for _ in range(3000):
 d=u.copy()
 for col,sd in [('S','S_SD'),('vC','vC_SD')]:
  sigma=np.sqrt(np.log1p((u[sd]/u[col])**2));mu=np.log(u[col])-.5*sigma**2
  d[col]=rng.lognormal(mu,sigma)
 counts+=frontier(d)
u['frontier_fraction_independent_reported_SD_perturbation']=counts/3000
u.to_csv(OUT/'rubisco_reported_SD_sensitivity.csv',index=False)
# Original control observations; these do not replace mandated paired assays.
controls=raw[raw.species.str.contains('partita|mays',case=False,na=False)&(raw.mutant==0)]
controls.to_csv(OUT/'rubisco_control_observations.csv',index=False)
Sref=float(f[f.species.str.contains('Galdieria partita',case=False,na=False)].S.max())
Kref=float(strict[strict.species.str.fullmatch('Zea mays',case=False,na=False)].vC.median())
qual=pairs[(pairs.S>=Sref)&(pairs.vC>=Kref)]
# 2025 table has corrected coupled-assay rates at 30 C, no oxygen specificity.
lib=pd.read_csv('data/rubisco2025_supp3.csv');whole=pd.read_csv('data/rubisco2025_library.csv')
lib[~lib.identifier.str.contains('_ANC_',na=False)].sort_values('rate_lab_corrected',ascending=False).head(25).to_csv(OUT/'rubisco2025_fastest_missing_specificity.csv',index=False)
recent={'rows_in_supplementary_data3':len(lib),'rows_with_numeric_rates':int(lib.rate_lab_corrected.notna().sum()),'rows_flagged_below_0_5':int(lib['active<0.5s-1'].notna().sum()),'new_to_this_study':int((lib.new_to_this_study==1).sum()),'columns':lib.columns.to_list(),'max_corrected_rate_30C_including_reconstructed_ancestor':float(lib.rate_lab_corrected.max()),'max_natural_corrected_rate_30C':float(lib[~lib.identifier.str.contains('_ANC_',na=False)].rate_lab_corrected.max()),'ancestral_reconstruction_rows':int(lib.identifier.str.contains('_ANC_',na=False).sum()),'rate_correction_factor_median':float((whole.rate_lab_corrected/whole.rate_lab).median()),'specificity_column_present':False}
res={'input_rows':len(raw),'curation_use_true':int(raw['use'].sum()),'primary_wild_formI_25C_rows':len(f),'same_original_row_pairs':len(strict),'source_species_pairs':len(pairs),'source_count':pairs.short_ref.nunique(),'species_count':pairs.species.nunique(),'species_median_naive':corr(species),'source_paired_naive':corr(pairs),'strict_same_row_naive':corr(strict),'within_study_log_slope':within_slope(pairs,'short_ref'),'within_study_cluster_bootstrap_2_5_50_97_5':study_ci(pairs),'within_taxonomy_log_slope':within_slope(pairs,'taxonomy'),'study_plus_taxonomy_log_slope':combined_adjusted_slope,'same_pH8_to8_5_only':corr(pairs[(pairs.pH_min>=8)&(pairs.pH_max<=8.5)]),'frontier_source_pairs':len(front),'frontier_source_studies':front.short_ref.nunique(),'frontier_species':front.species.to_list(),'reported_SD_pairs':len(u),'uncertainty_bounds_sensitivity':unc,'leave_one_study_slope_range':[min(x['within_study_slope'] for x in lo),max(x['within_study_slope'] for x in lo)],'historical_diagnostic_thresholds':{'specificity_G_partita':Sref,'kcat_Z_mays':Kref,'qualifying_pairs':len(qual),'max_specificity_above_maize_turnover':float(pairs[pairs.vC>=Kref].S.max()),'max_turnover_above_G_partita_specificity':float(pairs[pairs.S>=Sref].vC.max())},'recent2025':recent,'caveats':['No phylogenetic independence assumed; study-block correction is not phylogenetic correction.','No harmonized 2026 dataset merged here.','Original workbook and public CSV have small numeric differences; original workbook used throughout.','Source-paired means study pairing, not simultaneous matched-control or same preparation assays.','SD perturbation fractions are model sensitivities, not posterior probabilities of true Pareto optimality.']}
cyano=pd.read_csv('data/rubisco_cyano2026_kinetics.csv')
cyano=cyano[~cyano.Species.str.contains('control')]
r,p=stats.pearsonr(np.log(cyano.vC),np.log(cyano.S));rho,sp=stats.spearmanr(cyano.vC,cyano.S)
res['primary_2026_cyanobacteria']={'source_doi':'10.1093/plphys/kiag259','n_strains':len(cyano),'temp_C':25,'errors_are':'standard errors','pearson_log_r':float(r),'pearson_p_exploratory':float(p),'spearman_rho':float(rho),'spearman_p_exploratory':float(sp),'historical_rectangle_qualifiers':int(((cyano.S>=Sref)&(cyano.vC>=Kref)).sum()),'max_S':float(cyano.S.max()),'max_vC':float(cyano.vC.max()),'note':'Separate same-study cohort; not pooled into 2019 data; small n and clade sampling prevent a global tradeoff conclusion.'}
(OUT/'rubisco_summary.json').write_text(json.dumps(res,indent=2,allow_nan=False))
print(json.dumps(res,indent=2))
