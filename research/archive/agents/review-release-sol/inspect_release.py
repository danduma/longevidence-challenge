from pathlib import Path
import zipfile,hashlib,json
import pandas as pd
from Bio import SeqIO
S=Path('/path/to/longevidence-checkout/research/millennium/sources/rubisco-active-learning')
P=S/'supplementary/doubling_time_package_v2'
Z=Path('/path/to/longevidence-checkout/research/millennium/sources/rubisco-active-learning-Version3.zip')
O=Path('outputs');O.mkdir(exist_ok=True)
z=zipfile.ZipFile(Z); entries=[]
for p in sorted(S.rglob('*')):
 if p.is_file():
  rel=str(p.relative_to(S)); candidates=[n for n in z.namelist() if n.endswith('/'+rel)]
  b=p.read_bytes(); entries.append(dict(path=rel,size=len(b),sha256=hashlib.sha256(b).hexdigest(),archive_match=len(candidates)==1 and b==z.read(candidates[0])))
pd.DataFrame(entries).to_csv(O/'release_manifest.csv',index=False)
(O/'archive_listing.txt').write_text('\n'.join(z.namelist()))
# Extract only relevant supplementary code/text/csv to this job, preserving sources.
for n in z.namelist():
 if '/doubling_time_package_v2/' in n and ('flamholz' in n.lower()) and n.endswith(('.py','.csv','.txt','.fasta')):
  dest=Path('archive_inspection')/n.split('/doubling_time_package_v2/')[1];dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(n))
f=pd.read_csv(P/'inputs/flamholz_dataset_S1.csv');m=pd.read_csv(P/'manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv');seq=list(SeqIO.parse(P/'manuscript_assets/flamholz/flamholz_rubisco_wt.fasta','fasta'))
print('source',f.shape,'meta',m.shape,'fasta',len(seq),'unique IDs',len(set(x.id for x in seq)),'unique sequences',len(set(str(x.seq) for x in seq)))
print('meta duplicate extra rows',m.duplicated().sum(),'species unique',m.species.nunique(),'label unique',m.label.nunique())
print(m[m.duplicated('label',keep=False)].to_string(index=False))
print('whitespace source',[(repr(x),repr(x.strip())) for x in f.species.unique() if x!=x.strip()])
print('whitespace meta',[(repr(x),repr(x.strip())) for x in m.species.unique() if x!=x.strip()])
print('duplicate observations',f.drop(columns=[f.columns[0]]).duplicated().sum(),'duplicate identifiers',f.identifier.duplicated().sum())
print('lengths',m.seq_len.describe().to_dict());print(m[m.seq_len<400][['species','label','uniprot_accession','seq_len']].to_string(index=False))
print('FASTA examples',[(x.id,x.description) for x in seq[:4]])
print('archive embeddings',[x for x in z.namelist() if 'flamholz' in x and x.endswith('.npy')])
print('archive SHA256',hashlib.sha256(Z.read_bytes()).hexdigest())
print('nonmatching',pd.DataFrame(entries).query('archive_match == False').path.tolist())
