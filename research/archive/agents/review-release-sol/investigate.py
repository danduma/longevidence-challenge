from audit import *
f,m,rec=load();wt=f[(f.primary==1)&(f.mutant==0)];raw=wt.groupby('species',as_index=False)[T].median();r=raw.merge(m[['species','sequence']],on='species');d=pd.read_csv(P/'manuscript_assets/flamholz/embedding_vs_kinetics/pred_S.csv');r=r[r.S.gt(0)].reset_index(drop=True)
print('length',len(r),len(d));print(pd.DataFrame({'expected':r.species,'deposit':d.species,'exp_y':np.log10(r.S),'deposit_y':d.y_true}).query('expected != deposit').head(30).to_string(index=False))
print('mismatch raw targets');p=pd.read_csv('outputs/target_provenance.csv');print(p[~p.raw_log_match].to_string(index=False));print('normalization changes',p.query('normalization_changes_target').groupby(['branch','target']).size())
print('meta species with >2 tokens',m[m.species.str.split().str.len()>2][['species','taxon_name','uniprot_accession']].to_string(index=False))
print('WT selected named mutants',wt[wt.species.str.contains('K332|F92|G91|I87|R88')][['species','mutant','primary','variant']].to_string(index=False))
