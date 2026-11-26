from pathlib import Path
import urllib.request,urllib.parse,json,hashlib,csv
P=Path('sources'); tree=json.loads((P/'protrec_tree.json').read_text())['tree']; commit='8b9a1865108fa44c3e0a24ac20a59754222847ca'
rows=[]
for t in tree:
 path=t['path']
 if t['type']!='blob': continue
 if ('screen4_DMS/' in path and (path.endswith('.R') or any(path.endswith(x) for x in ['Data_DMS_unmapped.txt.zip','Data_DMS_controls_unmapped.txt','TEVs_BCs_suspicious.txt','Data_DMS_validation_PacBio.txt','Data_DMS_validation_manual.txt','Data_DMS_internal_controls.txt']))) or path.endswith('README.md') or path.endswith('README.txt'):
  url='https://raw.githubusercontent.com/JeschekLab/ProtRec/'+commit+'/'+urllib.parse.quote(path)
  dest=P/(path.replace('/', '__') if Path(path).name in ['README.md','README.txt'] else Path(path).name)
  if not dest.exists() or hashlib.sha1(b'blob '+str(dest.stat().st_size).encode()+b'\0'+dest.read_bytes()).hexdigest()!=t['sha']: dest.write_bytes(urllib.request.urlopen(url).read())
  b=dest.read_bytes();assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==t['sha'], path
  rows.append([str(dest),path,url,len(b),hashlib.sha256(b).hexdigest(),t['sha']])
with (P/'manifest.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['local','repository_path','url','bytes','sha256','git_blob']);w.writerows(rows)
print('verified',len(rows),'files')
