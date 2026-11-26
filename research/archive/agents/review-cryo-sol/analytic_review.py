"""Independent dimensional, finite-sample, scalar, symmetry and horizon checks."""
import json, subprocess, hashlib, sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
sys.dont_write_bytecode=True
root=Path(__file__).resolve().parent;target=root.parent/'cryo-robust-astra';up=target/'upstream'
p=json.loads((up/'results/grid_convergence_margin.json').read_text())['design_parameters'];flat=(p[3]+2*p[4]+p[5])/4
C=(1078+.43*150)*(2402.3-7.6*150)
onset={}
for mode,sar in [('extrapolated',1013),('clamped',691)]:
 for name,duty in [('ramp',p[3]),('flat',flat),('full',1)]:
  multiplier=60*1000*sar*duty/C
  onset[mode+'_'+name]=dict(rate_per_mg_ml=multiplier,necessary_loading=50/multiplier,nominal_rate=1.47*multiplier,high_rate=1.83*multiplier)
mu,s=1.47,.36;bound=mu-s/np.sqrt(3);rounded=mu+.005-(s-.005)/np.sqrt(3)
equality=np.array([bound,bound,mu+2*s/np.sqrt(3)])
rng=np.random.default_rng(756990);x=rng.normal(size=(100000,3));x-=x.mean(axis=1,keepdims=True);x*=s/x.std(axis=1,ddof=1,keepdims=True);x+=mu
assert x.min(axis=1).max()<=bound
assert abs(equality.mean()-mu)<1e-14 and abs(equality.std(ddof=1)-s)<1e-14
# Original schedule uses n linspace ordinates but n Euler updates, a rectangular integral.
dt_original=.012140447055757573;n=int(np.ceil(70/dt_original));u=np.interp(np.linspace(0,1,n),[0,.5,1],p[3:])
duty=dict(flat=flat,continuous_ramp_s=35*(p[3]+2*p[4]+p[5])/2,continuous_flat_s=70*flat,
          original_steps=n,original_horizon_s=n*dt_original,original_ramp_s=float(dt_original*u.sum()),original_flat_s=n*dt_original*flat)
def capacity(t):
 cp=2402.3+7.6*t if t< -120 else 66747+545.24*t if t<=-119 else 2678.8+11.73*t+.039*t*t
 return (1078-.43*t)*cp
def integral(t):
 return sum(quad(capacity,a,min(t,b),epsabs=1e-6)[0] for a,b in [(-150,-120),(-120,-119),(-119,0)] if t>a)
q=1000*1.83*691*flat;scalar_exact=brentq(lambda t:integral(t)-q*70,-150,-80)
checks=[]
for name,nn,fe in [('scalar',1,1.83),('symmetry',2,1.83),('zero',2,0)]:
 f=root/'results'/('analytic_'+name+'.txt');nf=nn-1
 f.write_text(f'{nn} {nf} .0005 .005 70 1 2 0 -150 {flat} {flat} {flat}\n'+''.join(f'2 {fe}\n' for _ in range(nn))+('0 1\n' if nf else ''))
 res=json.loads(subprocess.check_output([str(root/'independent'),str(f)],text=True));res['name']=name;checks.append(res)
 if fe:assert abs(res['cold']-scalar_exact)<.005 and res['hot']==res['cold']
 else:assert res['cold']==res['hot']==-150 and res['source_J']==0
files=['report.md','preregistration.md','thermal.cpp','investigate.py','verify.py','assay_bound.py','sources/heart.xml','sources/heart_fig1.jpg','upstream/src/vitrify/geometry.py','upstream/src/vitrify/properties.py','upstream/src/vitrify/thermal.py','upstream/coupler/workflow.py','upstream/results/grid_convergence_margin.json']
hashes={str(target/name):hashlib.sha256((target/name).read_bytes()).hexdigest() for name in files}
out=dict(capacity_J_m3K=C,onset=onset,triplicate_conditional=dict(mean=mu,sample_sd=s,sharp_bound=bound,rounded_bound=rounded,equality_triplet=equality.tolist(),random_max_min=float(x.min(axis=1).max()),rounded_ramp_rate=rounded*onset['extrapolated_ramp']['rate_per_mg_ml'],empirical_link_established=False),duty=duty,scalar_exact=scalar_exact,checks=checks,readonly_input_hashes=hashes)
(root/'results/analytic.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
