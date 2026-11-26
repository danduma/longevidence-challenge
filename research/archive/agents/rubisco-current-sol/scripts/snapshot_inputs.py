"""Copy exact read-only shared inputs locally and record their cryptographic identity."""
from pathlib import Path
import shutil,hashlib,json
R=Path(__file__).resolve().parents[1];SH=R.parent/'scout-data-sol';D=R/'raw/shared_inputs';D.mkdir(exist_ok=True)
files=['data/rubisco2019_original.xlsx','data/rubisco2025_library.csv','data/rubisco2025_supp3.csv','sources/rubisco_cyano2026_pmc.html','sources/rubisco2025.xml','sources/rubisco_correction2025.xml','sources/rubisco2025_tree.json','results/rubisco_bibliographic_overlay.csv','results/rubisco_source_pairs_linked.csv','report.md','source_inventory.csv','scripts/analyze_rubisco.py','scripts/build_assay_ledger.py','scripts/extract_cyano2026.py']
manifest=[]
for f in files:
 p=SH/f;dest=D/f;dest.parent.mkdir(parents=True,exist_ok=True)
 if not dest.exists():shutil.copyfile(p,dest)
 b=p.read_bytes();assert b==dest.read_bytes()
 manifest.append({'shared_path':str(p),'snapshot':str(dest.relative_to(R)),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
for p in [R.parent/'review-release-sol/report.md',R.parent/'review-release-astra/report.md',R.parents[1]/'RANKING.md']:
 b=p.read_bytes();dest=D/'context'/p.parent.name/p.name;dest.parent.mkdir(parents=True,exist_ok=True)
 if not dest.exists():shutil.copyfile(p,dest)
 manifest.append({'shared_path':str(p),'snapshot':str(dest.relative_to(R)),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
(R/'results/shared_input_snapshots.json').write_text(json.dumps(manifest,indent=2))
print(len(manifest),'shared inputs copied unchanged')
