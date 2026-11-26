from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
for name,sha in json.loads((root/'pinned_inputs.json').read_text()).items():
 assert hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name
print('C05 pinned input SHA-256 checks passed')
