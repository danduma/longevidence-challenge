"""C05 final local artifact provenance and integrity record."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
origins={
 'inputs/independent_original.cpp':'public/longevidence-challenge/research/cryopreservation/independent.cpp',
 'inputs/baseline_core.json':'public/longevidence-challenge/research/cryopreservation/results/core.json',
 'inputs/grid_convergence_margin.json':'public/longevidence-challenge/research/cryopreservation/inputs/grid_convergence_margin.json',
 'vendor/vitrify/geometry.py':'public/longevidence-challenge/research/cryopreservation/vendor/vitrify/geometry.py',
 'sources/model.xml':'research/millennium/agents/cryo-robust-astra/sources/model.xml',
 'sources/upstream_thermal.py':'research/millennium/agents/cryo-robust-astra/upstream/src/vitrify/thermal.py'}
prov=dict(id='C05',date='2026-11-26',source_origins={k:dict(workspace_relative_source=v,local_sha256=sha(R/k)) for k,v in origins.items()},live_retrievals=json.loads((R/'sources/retrieval.json').read_text()),pinned_inputs=json.loads((R/'pinned_inputs.json').read_text()),protocol_sha256=sha(R/'protocol.md'),versions=json.loads((R/'results/verification.json').read_text())['versions'],licenses={'solver':'MIT; retained LICENSE and original independent implementation comment','Vitrify':'Apache-2.0; vendor/VITRIFY-LICENSE; retained copyright/SPDX','Joshi_Rabin_2023_XML':'CC BY 4.0; publisher XML license and author/source attribution in sources.md'},method='No external model calls or experimental work. SHA-256 calculated over exact bytes; protocol local hash has no third-party timestamp.')
(R/'provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
files=[p for p in R.rglob('*') if p.is_file() and not any(x.startswith('.') for x in p.relative_to(R).parts) and p.name not in ['manifest.sha256','domains','verification.log']]
(R/'manifest.sha256').write_text(''.join(f'{sha(p)}  {p.relative_to(R)}\n' for p in sorted(files)))
s=json.loads((R/'results/summary.json').read_text());c=json.loads((R/'contribution.json').read_text())
assert c['id']=='C05' and c['status']=='agent-produced; awaiting coordinator verification'
assert s['simulations']==16 and s['paired_comparisons']==8
assert all(n==0 for n in s['decision_differences'].values()) and all(n==0 for n in s['ranking_differences'].values())
assert (R/'vendor/vitrify/geometry.py').read_bytes()==(R/'sources/upstream_geometry.py').read_bytes()
print('C05 artifact schema, results and pinned/live geometry concordance passed; manifest sealed')
