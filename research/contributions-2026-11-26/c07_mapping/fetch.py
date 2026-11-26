#!/usr/bin/env python3
import concurrent.futures,datetime,hashlib,json,re,shutil
from pathlib import Path
import urllib.request
import pandas as pd
ROOT=Path(__file__).resolve().parent
SHARED=Path('../../../millennium/sources/rubisco-active-learning')
BASE=SHARED/'supplementary/doubling_time_package_v2'
norm=lambda s:' '.join(s.lower().replace('−','-').split())
def select(m):
 return m[(m.seq_len<400)|m.apply(lambda x:norm(x.species)!=norm(x.taxon_name),axis=1)|m.species.str.contains(r'cv\.|cultivar|strain|sp\.|subsp\.|var\.|PCC|\d',regex=True)].copy()
def get(item):
 name,url=item;p=ROOT/'evidence'/name;p.parent.mkdir(parents=True,exist_ok=True)
 try:
  r=urllib.request.urlopen(url,timeout=55); body=r.read(); rec={'file':str(p.relative_to(ROOT)),'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'http_status':r.status,'headers':{k:v for k,v in r.headers.items() if k.lower() in ['x-uniprot-release','x-uniprot-release-date','content-type','last-modified']}}
  if r.status==200:p.write_bytes(body);rec['sha256']=hashlib.sha256(body).hexdigest()
  return rec
 except Exception as e:return {'url':url,'file':str(p.relative_to(ROOT)),'error':str(e),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
if __name__=='__main__':
 files={'meta.csv':'manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv','sequences.fasta':'manuscript_assets/flamholz/flamholz_rubisco_wt.fasta','assays.csv':'inputs/flamholz_dataset_S1.csv'}
 for t in ['S','KC','vC','KO','vO','KRuBP']:files['pred_'+t+'.csv']='manuscript_assets/flamholz/embedding_vs_kinetics/pred_'+t+'.csv'
 manifest=[]
 for dst,src in files.items():
  b=(BASE/src).read_bytes();(ROOT/'inputs'/dst).write_bytes(b);manifest.append({'file':'inputs/'+dst,'shared_path':str(BASE/src),'release_relative_path':'supplementary/doubling_time_package_v2/'+src,'sha256':hashlib.sha256(b).hexdigest(),'redistribution':'local verification snapshot; third-party terms not established; do not redistribute source CSV/predictions'})
 (ROOT/'inputs/manifest.json').write_text(json.dumps(manifest,indent=2))
 m=select(pd.read_csv(ROOT/'inputs/meta.csv').drop_duplicates());m.to_csv(ROOT/'results/suspects.csv',index=False)
 urls=[('uniprot/'+a+'.json','https://rest.uniprot.org/uniprotkb/'+a+'.json') for a in sorted(set(m.uniprot_accession))]
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:log=list(pool.map(get,urls))
 (ROOT/'evidence/retrieval.json').write_text(json.dumps(log,indent=2));print('requests',len(log),'success',sum(x.get('http_status')==200 for x in log))
