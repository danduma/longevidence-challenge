from pathlib import Path
import json,hashlib
R=Path(__file__).resolve().parents[1]
paths=[]
for folder in ['scripts','raw','derived','results']:
 for p in (R/folder).rglob('*'):
  if p.is_file() and '__pycache__' not in str(p) and p.suffix not in ['.log'] and p.name!='artifact_manifest.json':paths.append(p)
for name in ['run.sh','report.md','README.md','dataset_card.md','protocol.md','protocol_2025_calibration.md','refusals.md','access_failures.md','notes/amendment_01.md','notes/primary_reuse_audit.md']:paths.append(R/name)
rows=[{'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(set(paths))]
(R/'results/artifact_manifest.json').write_text(json.dumps({'release':'rubisco-control-benchmark-v1.0.0','date':'2026-11-26','artifacts':rows},indent=2));print(len(rows),'artifacts pinned')
