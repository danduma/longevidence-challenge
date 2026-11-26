#!/usr/bin/env python3
from audit import *
from baseline import groups
import sys,platform,importlib.metadata as md

def main():
 f,m,rec=load();b=pd.read_csv(O/'baseline_data.csv');g99=groups(b,pd.read_csv(O/'sequence_pairs_99.csv'),True);byseq=dict(zip(b.sequence,g99));rows=[]
 for branch,folder in [('fixed',P/'manuscript_assets/flamholz/embedding_vs_kinetics'),('nested_archive',Path('archive_inspection/manuscript_assets/flamholz/tune_stepA'))]:
  summary=pd.read_csv(folder/('summary_metrics.csv' if branch=='fixed' else 'summary_nestedcv.csv'))
  for t in summary.target:
   d=pd.read_csv(folder/f'pred_{t}.csv');computed=metrics(d.y_true,d.y_pred);shipped=summary.set_index('target').loc[t]
   assert int(shipped['n'])==len(d)
   for metric in ['spearman','r2','rmse']:assert np.isclose(computed[metric],shipped[metric],atol=1e-12)
   folds=pd.read_csv(O/f'reconstructed_{branch}_{t}_folds.csv');assert folds.species.equals(d.species)
   y=d.y_true.to_numpy();fold=folds.fold.to_numpy();mu=np.zeros(len(y));gg=folds.sequence.map(byseq).to_numpy()
   for k in range(5):mu[fold==k]=y[fold!=k].mean()
   exposed99=sum(any((gg==gg[i])&(fold!=fold[i])) for i in range(len(gg)))
   rows.append(dict(branch=branch,target=t,n=len(d),exposed99=exposed99,**{'mean_'+k:v for k,v in metrics(y,mu).items()}))
 pd.DataFrame(rows).to_csv(O/'original_fold_mean_control.csv',index=False)
 # Persist raw observation trace, including whitespace spelling and assay conditions.
 wt=f[(f.primary==1)&(f.mutant==0)].copy();wt['normalized_species']=wt.species.map(norm);wt.to_csv(O/'wt_observation_trace.csv',index=False)
 # Fresh integrity check against initial hashes.
 manifest=pd.read_csv(O/'release_manifest.csv');root=P.parents[1]
 for r in manifest.itertuples():assert hashlib.sha256((root/r.path).read_bytes()).hexdigest()==r.sha256
 env={'python':sys.version,'platform':platform.platform(),'packages':{x:md.version(x) for x in ['numpy','pandas','scipy','scikit-learn','biopython']}}
 (O/'environment.json').write_text(json.dumps(env,indent=2));print('All deposited metrics match summaries; all shared input hashes unchanged.');print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
