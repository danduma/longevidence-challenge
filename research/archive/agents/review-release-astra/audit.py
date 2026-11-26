"""Run from this job directory with ../../.venv/bin/python audit.py."""
from pathlib import Path
import hashlib, re, zipfile, io, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from Bio import SeqIO

root = Path('../../sources/rubisco-active-learning/supplementary/doubling_time_package_v2')
assets = root/'manuscript_assets/flamholz'
archive = zipfile.ZipFile('../../sources/rubisco-active-learning-Version3.zip')
def member(suffix):
    names = [n for n in archive.namelist() if n.endswith('/'+suffix)]
    assert len(names)==1
    return archive.read(names[0])
def save(rows,name):
    pd.DataFrame(rows).to_csv(name+'.csv',index=False)
def scores(y,p):
    return dict(r2=r2_score(y,p),spearman=spearmanr(y,p).statistic,rmse=np.mean((y-p)**2)**0.5)
def canonical(s):
    s=' '.join(s.split()).replace('H. vulgare','Hordeum vulgare').replace('S. cereale','Secale cereale')
    return re.split(r' (?:L\.|cv\.|var\.|subsp\.|ssp\.|Kotschy|MIT9313|PCC 7942)',s)[0]
manifest=[]
for path in sorted(root.rglob('*')):
    if path.is_file():
        b=path.read_bytes()
        manifest.append(dict(file=str(path.relative_to(root)),sha256=hashlib.sha256(b).hexdigest(),same=b==member(str(path.relative_to(root)))))
save(manifest,'manifest')
raw=pd.read_csv(root/'inputs/flamholz_dataset_S1.csv')
wt=raw.query('primary==1 and mutant==0').copy()
m=pd.read_csv(assets/'flamholz_rubisco_wt_meta.csv')
fa=list(SeqIO.parse(assets/'flamholz_rubisco_wt.fasta','fasta'))
assert [x.id for x in fa]==m.label.tolist()
m['seq']=[str(x.seq) for x in fa]
assert (m.seq.str.len()==m.seq_len).all()
m['hash']=m.seq.map(lambda s:hashlib.sha256(s.encode()).hexdigest())
m['canonical']=m.species.map(canonical)
wt['canonical']=wt.species.map(canonical)
save(m.drop(columns='seq'),'identities')
traits='S vC KC KO vO KRuBP'.split()
rawmed=wt.groupby('species')[traits].median()
med=wt.groupby('canonical')[traits].median()
counts=m.groupby('canonical').hash.nunique()
clean=m[m.canonical.isin(counts[counts==1].index)].drop_duplicates('canonical').merge(med,on='canonical')
save(wt[['species','canonical','identifier','short_ref',*traits]],'source_rows')
print('raw wt names normalized',len(raw),len(wt),wt.species.nunique(),wt.canonical.nunique(),'wt repeated IDs',wt.identifier.duplicated().sum())
print('meta',len(m),m.nunique().to_dict(),'ambiguous',counts[counts>1].to_dict())
print('control source',med.loc[med.index.str.contains('Galdieria|Zea mays'),['S','vC']])
print('clean control',clean[clean.canonical.str.contains('Galdieria|Zea mays')][['canonical','S','vC']])
metrics=[]; overlap=[]; provenance=[]; baseline=[]; oofs=[]
for t in traits:
    sub=m.merge(rawmed[t],on='species').sort_values('species',kind='stable')
    sub=sub[sub[t]>0].reset_index(drop=True)
    for variant,folder,summ in [('fixed','embedding_vs_kinetics','summary_metrics'),('nested','tune_stepA','summary_nestedcv')]:
        prefix='manuscript_assets/flamholz/'+folder+'/'
        if not any(n.endswith(prefix+'pred_'+t+'.csv') for n in archive.namelist()):continue
        p=pd.read_csv(io.BytesIO(member(prefix+'pred_'+t+'.csv')))
        reported=pd.read_csv(io.BytesIO(member(prefix+summ+'.csv'))).set_index('target').loc[t]
        expected=sub[['species',t]].copy()
        if variant=='fixed':
            expected.to_csv('literal_join_'+t+'.csv',index=False)
        by_name=m.drop_duplicates('species').set_index('species')
        sub=p[['species']].copy()
        for k in ['hash','canonical','uniprot_accession']:
            sub[k]=sub.species.str.strip().map(by_name[k])
        assert sub.hash.notna().all()
        sub[t]=sub.species.map(rawmed[t])
        assert sub[t].notna().all()
        calc=scores(p.y_true,p.y_pred)
        metrics.append(dict(target=t,analysis=variant,n=len(p),**calc,max_metric_delta=max(abs(calc[k]-reported[k]) for k in calc),max_y_delta=np.max(abs(np.log10(sub[t])-p.y_true))))
    for _,r in sub.iterrows():
        provenance.append(dict(target=t,species=r.species,canonical=r.canonical,original=r[t],normalized=med.loc[r.canonical,t],changed=not np.isclose(r[t],med.loc[r.canonical,t])))
    seen={k:0 for k in ('species','canonical','uniprot_accession','hash')}
    for tr,te in KFold(5,shuffle=True,random_state=0).split(sub):
        for k in seen:seen[k]+=sub.iloc[te][k].isin(sub.iloc[tr][k]).sum()
    overlap.append(dict(target=t,n=len(sub),unique_species=sub.species.nunique(),unique_sequences=sub.hash.nunique(),**seen))
    d=clean[clean[t]>0].reset_index(drop=True)
    if "--quality" in sys.argv:
        d=d[(d.species!="P. dilatatum") & d.seq_len.between(400,550)].reset_index(drop=True)
    X=np.array([[s.count(a)/len(s) for a in 'ACDEFGHIKLMNPQRSTVWY']+[len(s)] for s in d.seq])
    y=np.log10(d[t].to_numpy())
    groups=np.unique(d.hash)
    for seed in range(20):
        pred=np.zeros(len(y)); mean=pred.copy()
        for fold,(a,b) in enumerate(KFold(5,shuffle=True,random_state=seed).split(groups)):
            tr=d.hash.isin(groups[a]); te=d.hash.isin(groups[b])
            assert set(d[tr].canonical).isdisjoint(d[te].canonical) and set(d[tr].hash).isdisjoint(d[te].hash)
            pred[te]=make_pipeline(StandardScaler(),Ridge(alpha=10)).fit(X[tr],y[tr]).predict(X[te])
            mean[te]=y[tr].mean()
            for i in np.flatnonzero(te):
                oofs.append(dict(target=t,seed=seed,fold=fold,species=d.canonical[i],seq_hash=d.hash[i],y=y[i],ridge=pred[i],mean=mean[i]))
        for model,prediction in [('ridge',pred),('mean',mean)]:
            baseline.append(dict(target=t,seed=seed,n=len(y),groups=len(groups),model=model,**scores(y,prediction)))
for rows,name in [(metrics,'deposited'),(overlap,'overlap'),(provenance,'provenance'),(baseline,'baseline'),(oofs,'baseline_oof')]:save(rows,name+('_quality' if '--quality' in sys.argv and name.startswith('baseline') else ''))
print(pd.DataFrame(provenance).groupby('target').changed.sum())
print(pd.DataFrame(baseline).groupby(['target','model'])[['n','groups','r2','spearman','rmse']].median().round(3))
# Hypothesis: a whitespace-preserving metadata precursor was used when producing the outputs.
raw_names=sorted(wt.species.unique())
last={s.strip():s for s in raw_names}
for t in traits:
    alt=m.assign(species=m.species.map(last)).merge(rawmed[t],on='species').sort_values('species',kind='stable')
    alt=alt[alt[t]>0]
    p=pd.read_csv(assets/'embedding_vs_kinetics'/f'pred_{t}.csv')
    print('last raw alias ordered match',t,alt.species.tolist()==p.species.tolist())
print('ambiguous residues',sum(bool(set(s)-set('ACDEFGHIKLMNPQRSTVWY')) for s in m.seq))
print('identical metadata surplus',m.duplicated().sum(),'content duplicate WT surplus',wt.drop(columns=['Unnamed: 0']).duplicated().sum())
print('source examples',wt[wt.species.str.contains('Galdieria partita|Zea mays')][['species','S','vC','short_ref']].to_string(index=False))
