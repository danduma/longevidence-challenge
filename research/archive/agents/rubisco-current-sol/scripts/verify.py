"""Independent numerical/provenance checks for the delivered scientific artifacts."""
from pathlib import Path
import json,hashlib,math,itertools,sys
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];checks=[]
def ok(name,cond):
 assert cond,name
 checks.append(name)
# Frozen protocols are not edited by subsequent curation.
ok('historical protocol frozen',hashlib.sha256((R/'protocol.md').read_bytes()).hexdigest()=='b70f20a097316dad007c184d22b2edf8183ee895d1dc57622730d1fc2ccf29ed')
ok('2025 extension frozen',hashlib.sha256((R/'protocol_2025_calibration.md').read_bytes()).hexdigest()=='3d631df0bc608149d286e170f4f2f003889c68bb2f4a2910e4881d777a630b36')
for r in json.loads((R/'results/shared_input_snapshots.json').read_text()):
 p=R/r['snapshot'];ok('snapshot '+p.name,hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'])
 original=Path(r['shared_path'])
 if original.exists():ok('shared unchanged '+p.name,hashlib.sha256(original.read_bytes()).hexdigest()==r['sha256'])
x=pd.read_excel(R/'raw/shared_inputs/data/rubisco2019_original.xlsx')
sel=x['use'].eq(True)&x.primary.eq(True)&x.mutant.eq(False)&x.temp_C.eq(25)
ok('historical per-endpoint measured counts',int((sel&x.S.gt(0)).sum())==337 and int((sel&x.vC.gt(0)).sum())==274)
# Enumerate rectangles with independent pivot/common-node loops.
ed=pd.read_csv(R/'results/graph_edges.csv');rr=pd.read_csv(R/'results/control_rectangles.csv')
for trait in ['S','vC']:
 e=ed[(ed['mode']=='condition')&(ed.trait==trait)];p=e.pivot(index='enzyme',columns='source',values='value');found={}
 for a,b in itertools.combinations(p.index,2):
  q=p.loc[[a,b]].dropna(axis=1)
  for s,t in itertools.combinations(q.columns,2):
   if s.split(' | pH=')[0]==t.split(' | pH=')[0]:continue
   d=math.log(q.loc[a,s])-math.log(q.loc[b,s])-math.log(q.loc[a,t])+math.log(q.loc[b,t])
   found[(a,b,s,t)]=math.exp(abs(d))
 r=rr[(rr['mode']=='condition')&(rr.trait==trait)]
 ok('rectangle enumeration '+trait,len(found)==len(r))
 ok('every cycle factor '+trait,all(abs(found[(v.a,v.b,v.s,v.t)]-v.factor)<1e-12 for v in r.itertuples()))
lib=pd.read_csv(R/'raw/shared_inputs/data/rubisco2025_library.csv');pair=lib[(lib.from_lab==1)&(lib.from_literature==1)]
ok('authors flags yield 11 calibration rows',len(pair)==11 and pair.rate_lab.gt(0).all() and pair.rate_lit_30deg.gt(0).all())
r=pair.rate_lit_30deg.to_numpy()/pair.rate_lab.to_numpy();c=np.prod(r)**(1/len(r))
ok('reproduce deposited factor',np.allclose(lib.loc[lib.rate_lab.notna(),'rate_lab_corrected'],lib.loc[lib.rate_lab.notna(),'rate_lab']*c,rtol=1e-12))
ok('SD columns scaled by same factor',np.allclose(lib.std_rate_lab_corrected,lib.std_rate_lab*c,equal_nan=True,rtol=1e-12))
s=json.loads((R/'results/calibration_identifiability.json').read_text());ok('two-column calibration design rank1',s['rank']==1 and s['n_original_30C']==0)
z=pd.read_csv(R/'results/reused2023_2026_audit.csv');ok('all reused common cells audited',len(z)==12 and z.new_marked_reused.all())
t=z[(z.species=='C. thermalis')&(z.trait=='vC')].iloc[0]
ok('turnover SD/SE incompatible under rounding n>=3',t.implied_n_max<3 and t.compatible_n_3_to_6=='none')
ok('exact rounding upper bound',abs(t.implied_n_max-25/9)<1e-12)
m=pd.read_csv(R/'derived/measurement_ledger.csv');ok('unique endpoint IDs',m.measurement_id.is_unique and len(m)==623)
ok('unresolved C. thermalis error quarantined',m.loc[m.measurement_id.eq('kiag259:C. thermalis:vC'),'uncertainty_usable'].iloc[0]=='no_unresolved_SD_to_SE_transport')
ok('2026 cohort is 4 new table entities plus2 reused',m[m.study.eq('Aguilo-Nicolau 2026')].cohort.value_counts().to_dict()=={'2026_table1_new':8,'2026_table1_2023_reuse':4})
f=pd.read_csv(R/'derived/assay_frontier_benchmark.csv')
ok('frontier paired rows and identity quarantine',len(f)==145 and int((~f.wt_identity_eligible).sum())==1)
ok('no historical diagnostic positives',not f.historical_threshold_diagnostic.any())
ok('historical mutation-name conflicts preserved',int(m.identity_warning.str.startswith('mutation_in_name').sum())==6)
orig=pd.read_csv(R/'results/origin_graph_summary.csv');s=orig[(orig.regime=='endpoint_origin')&(orig.trait=='S')].iloc[0]
ok('Griffithsia false contemporaneous bridge removed',not s.Griffithsia_maize_connected)
# Code blob hash checked against already downloaded GitLab tree, avoiding unpinned branch assumptions.
blob=(R/'raw/library2025_calibration.ipynb').read_bytes();gitsha=hashlib.sha1(b'blob '+str(len(blob)).encode()+b'\0'+blob).hexdigest()
tree=json.loads((R/'raw/shared_inputs/sources/rubisco2025_tree.json').read_text());want=next(v['id'] for v in tree if v['path']=='1_regression_measured_lit.ipynb')
ok('2025 notebook Git blob matches scout inventory',gitsha==want)
summary={'checks_passed':len(checks),'checks':checks,'notebook_git_blob':gitsha,'limitations':['Original Uemura1997 full numerical table inaccessible','Kubien2008 exact rate/error values retain compilation provenance; original figure text could not be fully transcribed','No independent biological replication or challenge validation performed']}
(R/'results/verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
