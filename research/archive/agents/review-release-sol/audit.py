#!/usr/bin/env python3
"""Audit pinned release without modifying inputs. Run with shared scientific Python."""
from pathlib import Path
from collections import defaultdict
import hashlib,json,re
import numpy as np,pandas as pd
from Bio import SeqIO, Align
from scipy.stats import spearmanr
from sklearn.metrics import r2_score,mean_squared_error
from sklearn.model_selection import KFold
P=Path('/path/to/longevidence-checkout/research/millennium/sources/rubisco-active-learning/supplementary/doubling_time_package_v2')
O=Path('outputs');O.mkdir(exist_ok=True)
T=['S','vC','KC','KO','vO','KRuBP']
def norm(s):return ' '.join(str(s).split())
def base(s):
 t=norm(s).split();return ' '.join(t[:2])
def metrics(y,p):return dict(spearman=float(spearmanr(y,p).statistic),r2=float(r2_score(y,p)),rmse=float(np.sqrt(mean_squared_error(y,p))))
def load():
 f=pd.read_csv(P/'inputs/flamholz_dataset_S1.csv');m=pd.read_csv(P/'manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv')
 rec=list(SeqIO.parse(P/'manuscript_assets/flamholz/flamholz_rubisco_wt.fasta','fasta'))
 by=defaultdict(set)
 for r in rec:by[r.id].add(str(r.seq).upper())
 assert all(len(v)==1 for v in by.values()),'ambiguous FASTA label'
 m['sequence']=m.label.map({k:next(iter(v)) for k,v in by.items()})
 assert m.sequence.notna().all();assert (m.sequence.str.len()==m.seq_len).all()
 return f,m,rec

def main():
 f,m,rec=load();wt=f[(f.primary==1)&(f.mutant==0)].copy()
 f[f.identifier.duplicated(keep=False)].to_csv(O/'reused_source_identifiers.csv',index=False)
 m[m.duplicated('label',keep=False)].to_csv(O/'duplicate_meta.csv',index=False)
 m[m.seq_len<400].to_csv(O/'short_sequences.csv',index=False)
 sm=m[['species','label','uniprot_accession','sequence']].drop_duplicates();sm['sha256']=sm.sequence.map(lambda x:hashlib.sha256(x.encode()).hexdigest());sm.to_csv(O/'sequence_map.csv',index=False)
 sm[sm.duplicated('sequence',keep=False)].sort_values('sha256').to_csv(O/'identical_sequence_species.csv',index=False)
 raw=wt.groupby('species')[T].median(); wn=wt.copy();wn.species=wn.species.map(norm);clean=wn.groupby('species')[T].median()
 provenance=[];overlap=[];recomputed=[];recon=[]
 for branch,folder in [('fixed',P/'manuscript_assets/flamholz/embedding_vs_kinetics'),('nested_archive',Path('archive_inspection/manuscript_assets/flamholz/tune_stepA'))]:
  for t in T:
   path=folder/f'pred_{t}.csv'
   if not path.exists():continue
   d=pd.read_csv(path);recomputed.append(dict(branch=branch,target=t,n=len(d),unique_species=d.species.nunique(),**metrics(d.y_true,d.y_pred)))
   # Reconstruct missing embedding IDs as every FASTA record in deposited order.
   emb_species=pd.DataFrame({'species':[dict(zip(m.label,m.species))[r.id] for r in rec],'sequence':[str(r.seq) for r in rec]})
   literal=raw.reset_index().merge(emb_species,on='species',how='inner');literal=literal[literal[t].gt(0)].reset_index(drop=True)
   # Alternative historical mapping assumption: raw last spelling wins after whitespace canonicalization.
   lookup={norm(s):s for s in raw.index}
   emb_species.species=emb_species.species.map(lambda s:lookup.get(norm(s),s))
   merged=raw.reset_index().merge(emb_species,on='species',how='inner');sub=merged[merged[t].gt(0)].reset_index(drop=True)
   same_species=len(sub)==len(d) and sub.species.equals(d.species)
   same_y=same_species and np.allclose(np.log10(sub[t]),d.y_true,atol=1e-12)
   recon.append(dict(branch=branch,target=t,n_literal=len(literal),literal_order_match=literal.species.equals(d.species),n_reconstructed=len(sub),species_order_match=same_species,y_true_match=same_y))
   for i,row in d.iterrows():
    sp=row.species;r=raw.loc[sp,t] if sp in raw.index else np.nan;c=clean.loc[norm(sp),t] if norm(sp) in clean.index else np.nan
    provenance.append(dict(branch=branch,target=t,row=i,species=sp,deposited_raw=10**row.y_true,raw_median=r,normalized_median=c,raw_log_match=bool(np.isclose(row.y_true,np.log10(r))),normalization_changes_target=bool(not np.isclose(r,c,equal_nan=True))))
   if not same_y:continue
   folds=np.empty(len(sub),int)
   for k,(_,te) in enumerate(KFold(5,shuffle=True,random_state=0).split(sub)):folds[te]=k
   sp=sub.species.to_numpy();sq=sub.sequence.to_numpy();bs=np.array([base(x) for x in sp])
   exposures={key:np.array([any((v==v[i])&(folds!=folds[i])) for i in range(len(v))]) for key,v in [('species',sp),('base_species',bs),('exact_sequence',sq)]}
   overlap.append(dict(branch=branch,target=t,n=len(sub),unique_species=len(set(sp)),unique_sequences=len(set(sq)),**{k+'_exposed_test_rows':int(v.sum()) for k,v in exposures.items()}))
   sub[['species','sequence']].assign(fold=folds,**{key+'_exposed':value for key,value in exposures.items()}).to_csv(O/f'reconstructed_{branch}_{t}_folds.csv',index=False)
 pd.DataFrame(provenance).to_csv(O/'target_provenance.csv',index=False);pd.DataFrame(recomputed).to_csv(O/'recomputed_metrics.csv',index=False);pd.DataFrame(overlap).to_csv(O/'reconstructed_overlap.csv',index=False);pd.DataFrame(recon).to_csv(O/'reconstruction_checks.csv',index=False)
 # Normalized target construction, no ambiguous sequence resolution.
 wn=wn.drop_duplicates(subset=[x for x in f.columns if x!=f.columns[0]])
 mm=m.copy();mm.species=mm.species.map(norm);mm=mm.drop_duplicates(['species','sequence'])
 amb=mm.groupby('species').sequence.nunique();amb=amb[amb>1].index.tolist();mm=mm[~mm.species.isin(amb)]
 d=wn.groupby('species',as_index=False)[T].median().merge(mm[['species','sequence','uniprot_accession']],on='species',validate='one_to_one').sort_values('species').reset_index(drop=True)
 d['base_species']=d.species.map(base);d['seq_sha256']=d.sequence.map(lambda s:hashlib.sha256(s.encode()).hexdigest());d.to_csv(O/'baseline_data.csv',index=False)
 # Full-length global alignments only; do not inflate fragment identity with short overlap denominators.
 aligner=Align.PairwiseAligner(mode='global',match_score=1,mismatch_score=-1,open_gap_score=-2,extend_gap_score=-0.5)
 pairs=[]
 for i in range(len(d)):
  for j in range(i):
   a,b=d.sequence[i],d.sequence[j]
   if min(len(a),len(b))/max(len(a),len(b))<.99:continue
   if a==b:ident=1.
   else:
    al=aligner.align(a,b)[0];ident=al.counts().identities/max(len(a),len(b))
   if ident>=.99:pairs.append(dict(i=i,j=j,species_i=d.species[i],species_j=d.species[j],identity=ident))
 pd.DataFrame(pairs,columns=['i','j','species_i','species_j','identity']).to_csv(O/'sequence_pairs_99.csv',index=False)
 info=dict(source_rows=len(f),wt_rows=len(wt),raw_wt_species=wt.species.nunique(),normalized_wt_species=wn.species.nunique(),metadata_rows=len(m),metadata_species=m.species.nunique(),fasta_records=len(rec),unique_sequences=m.sequence.nunique(),baseline_species=len(d),ambiguous_species_excluded=amb,source_whitespace_names=[repr(x) for x in f.species.unique() if norm(x)!=x],unmatched_wt_species=sorted(set(wn.species)-set(mm.species)),temp_counts=wt.temp_C.value_counts(dropna=False).to_dict(),pH_counts=wt.pH.value_counts(dropna=False).to_dict())
 (O/'audit_summary.json').write_text(json.dumps(info,indent=2,default=str));print(pd.DataFrame(recomputed).to_string(index=False));print(pd.DataFrame(overlap).to_string(index=False));print(info)
if __name__=='__main__':main()
