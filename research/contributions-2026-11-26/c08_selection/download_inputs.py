from pathlib import Path
import hashlib,json,urllib.request
root=Path(__file__).resolve().parent;out=root/'downloaded-inputs';out.mkdir(exist_ok=True)
for spec in json.loads((root/'sources.json').read_text())['inputs'].values():
 if 'download_url' not in spec:continue
 p=out/Path(spec['upstream_path']).name
 data=urllib.request.urlopen(spec['download_url'],timeout=60).read()
 assert hashlib.sha256(data).hexdigest()==spec['sha256'],p.name
 p.write_bytes(data);print('verified',p.name)
