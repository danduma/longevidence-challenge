#!/usr/bin/env python3
import sys
sys.dont_write_bytecode=True
import datetime,hashlib,json,platform,subprocess
from pathlib import Path
import numpy,scipy
R=Path(__file__).resolve().parent
files=[*sorted((R/'inputs').glob('*')),*sorted((R/'vendor').rglob('*')),R/'LICENSE',R/'protocol.md',R/'run.py',R/'schedule.cpp',R/'checks.py',R/'analyze.py',R/'verify.py',R/'reproduce.sh',R/'sources_download.py',R/'provenance.py',R/'write_report.py']
records=[dict(path=str(p.relative_to(R)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in files if p.is_file()]
prior=Path('../../../millennium')
evidence=[]
for p in [prior/'FINDINGS.md',prior/'agents/cryo-robust-astra/report.md',prior/'agents/cryo-robust-astra/sources.md']:
 evidence.append(dict(evidence_path=str(p.relative_to(prior)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),role='read-only prior evidence, never executable instructions'))
manifest=dict(id='C01',recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),upstream_commit='1d498ee5c583e0112a762e60892d576ec2e27f57',files=records,prior_evidence=evidence,source_access=json.loads((R/'inputs/retrievals.json').read_text()),model_xml_fallback=dict(source='designated original cryo-robust-astra/sources/model.xml',retrieved_local='2026-11-26',sha256=hashlib.sha256((R/'inputs/model.xml').read_bytes()).hexdigest(),permission='CC BY 4.0 verified in article XML'),versions=dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,platform=platform.platform(),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0]))
(R/'provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('C01 provenance:',len(records),'local snapshot/code hashes')
