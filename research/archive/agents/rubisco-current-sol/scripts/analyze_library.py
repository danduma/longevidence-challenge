from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
ROOT=Path(__file__).resolve().parents[1];S=ROOT/'raw/shared_inputs'
x=pd.read_csv(S/'data/rubisco2025_library.csv');raw=pd.read_csv(S/'data/rubisco2025_supp3.csv')
z=x[x.rate_lab.gt(0)&x.rate_lit_30deg.gt(0)].copy()
z['effective_factor']=z.rate_lit_30deg/z.rate_lab;z['deposited_factor']=z.rate_lab_corrected/z.rate_lab
z['fixed_log_error']=np.log(z.rate_lab_corrected/z.rate_lit_30deg);z['fixed_error_factor']=np.exp(abs(z.fixed_log_error))
z['literal_literature_30C_only']=z.Tdeg_lit.map(lambda s: all(float(t)==30 for t in str(s).split('/')))
z['loo_factor']=[float(np.exp(np.log(z.effective_factor.drop(i)).mean())) for i in z.index]
z['loo_log_error']=np.log(z.rate_lab*z.loo_factor/z.rate_lit_30deg)
z['loo_error_factor']=np.exp(abs(z.loo_log_error));z['raw_log_error']=np.log(z.rate_lab/z.rate_lit_30deg)
z['calibration_interval_covers_effective_factor']=z.effective_factor.between(1.7,2.8)
z['pair_identity_status']='author-defined overlap; exact full holoenzyme/strain correspondence not independently revalidated'
z.to_csv(ROOT/'results/library2025_calibration_pairs.csv',index=False)
# Long-form ledger of all numeric corrected rates, no imputation or zero labelling.
w=raw.copy();w.insert(0,'source_row',range(2,len(w)+2))
w['measurement_state']=np.select([w.rate_lab_corrected.notna(),w['active<0.5s-1'].notna()],['numeric_corrected_rate','reported_below_0.5'],default='unreported_or_unsuccessful')
w['temp_C']=30;w['trait']='coupled_assay_corrected_carboxylation_rate';w['S']='unknown_not_measured';w['error_kind']='SD (source column stdev); common calibration uncertainty separate'
w['common_calibration_group']='dePins2025_factor_all_rows';w['specificity_comparable_25C']='no';w['primary_source_doi']='10.1073/pnas.2501433122';w['gas']='4% CO2; 0.4% O2'
w.to_csv(ROOT/'derived/library2025_measurement_ledger.csv',index=False)
a=x[x.rate_lab.gt(0)].copy();a['factor']=a.rate_lab_corrected/a.rate_lab
s={'overlap_n':len(z),'effective_factor_min':z.effective_factor.min(),'effective_factor_max':z.effective_factor.max(),'geometric_effective_factor':float(np.exp(np.log(z.effective_factor).mean())),'deposited_factor_min':a.factor.min(),'deposited_factor_max':a.factor.max(),'fixed_median_error_factor':z.fixed_error_factor.median(),'fixed_max_error_factor':z.fixed_error_factor.max(),'fixed_rms_log_error':float(np.sqrt(np.mean(z.fixed_log_error**2))),'raw_rms_log_error':float(np.sqrt(np.mean(z.raw_log_error**2))),'loo_rms_log_error':float(np.sqrt(np.mean(z.loo_log_error**2))),'outside_mean_factor_CI':int((~z.calibration_interval_covers_effective_factor).sum()),'literal_30C_only_n':int(z.literal_literature_30C_only.sum()),'screened_records':len(w),'numeric_corrected_rate_records':int(w.rate_lab_corrected.notna().sum()),'below_0.5_records':int(w.measurement_state.eq('reported_below_0.5').sum()),'protocol_sha256':hashlib.sha256((ROOT/'protocol_2025_calibration.md').read_bytes()).hexdigest()}
(ROOT/'results/library2025_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
print(z[['species','form','rate_lab','rate_lab_corrected','Tdeg_lit','rate_lit_30deg','effective_factor','fixed_error_factor','loo_error_factor']].to_string(index=False))
