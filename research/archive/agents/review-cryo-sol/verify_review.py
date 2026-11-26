"""Final artifact and claim checks; independent simulations must already exist."""
from pathlib import Path
import json,re,hashlib,subprocess
ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'results/all.json').read_text())
assert len(rows)==22
def find(**criteria):
 found=[r for r in rows if all(r[k]==v for k,v in criteria.items())]
 assert len(found)==1,(criteria,len(found));return found[0]
for r in rows:
 assert r['missing']==0 and r['min_rate']>0
 assert abs(r['energy_rel_error'])<1e-6
 assert abs(r['muscle_mean']-r['mean'])<1e-12
 assert r['pass']==(r['min_rate']>=50 and r['cold']>=-35 and r['hot']<=0)
assert all(r['duty_s']==rows[0]['duty_s'] for r in rows)
for mm in [.75,.5,.3]:
 for d,expected in [('ramp',False),('flat',True)]:
  r=find(mm=mm,design=d,mean=1.83,amplitude=0,mode='clamped',dt_requested=.004 if mm==.3 else .01,h=0)
  assert r['pass']==expected
nominal={d:find(mm=.5,design=d,mean=1.47,amplitude=0,mode='extrapolated',dt_requested=.01,h=0) for d in ['ramp','flat']}
assert nominal['flat']['min_rate']>nominal['ramp']['min_rate']
assert nominal['flat']['window_rate']<nominal['ramp']['window_rate']
split=[find(mm=.3,design='ramp',mean=1.47,amplitude=.36,sign=s,mode='clamped') for s in [1,-1]]
assert split[0]['cold']< -35 <split[1]['cold']
assert abs(split[0]['iron_g']-split[1]['iron_g'])<1e-12
for amp,expected in [(.01,True),(.015,False)]:
 assert find(mm=.5,design='flat',mean=1.83,amplitude=amp,mode='clamped')['pass']==expected
cross=json.loads((ROOT/'results/crossing.json').read_text())
assert len(cross)==5
for r in cross[:4]:assert 0<=r['window_rate']-r['crossing_rate']<.005
assert cross[4]['cold']< -37 and cross[4]['min_rate']<50
base=json.loads((ROOT/'results/original_baselines.json').read_text())
assert len(base)==70 and max(abs(r['error']) for r in base)<1e-10
analytic=json.loads((ROOT/'results/analytic.json').read_text())
assert len(analytic['checks'])==3
scalar,symmetry,zero=analytic['checks']
assert abs(scalar['cold']-analytic['scalar_exact'])<.0002
assert scalar['cold']==symmetry['cold']==symmetry['hot']
assert zero['cold']==zero['hot']==-150 and zero['source_J']==0
assert not analytic['triplicate_conditional']['empirical_link_established']
assert abs(analytic['duty']['continuous_ramp_s']-analytic['duty']['continuous_flat_s'])<1e-12
for path,digest in analytic['readonly_input_hashes'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
assert subprocess.check_output(['git','-C',str(ROOT.parent/'cryo-robust-astra/upstream'),'rev-parse','HEAD'],text=True).strip()=='1d498ee5c583e0112a762e60892d576ec2e27f57'
for name in ['report.md','sources.md','access-failures.md','refusals.md']:
 assert (ROOT/name).exists()
 for link in re.findall(r'\]\(([^)]+)\)',(ROOT/name).read_text()):
  if not link.startswith('https://'):assert (ROOT/link).exists(),(name,link)
print('PASS: 22 main simulations, 5 crossing/boundary checks, 3 scalar fixtures, 70 original baseline metrics, document links and read-only input hashes.')
print('Maximum main-run relative energy error:',max(abs(r['energy_rel_error']) for r in rows))
print('Maximum original baseline absolute error:',max(abs(r['error']) for r in base))
print('0.30-mm spatial terminal difference:',split[1]['cold']-split[0]['cold'])
print('Maximum explicit crossing correction:',max(r['window_rate']-r['crossing_rate'] for r in cross[:4]))
