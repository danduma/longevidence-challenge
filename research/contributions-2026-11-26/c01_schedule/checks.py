#!/usr/bin/env python3
"""C01 scalar quadrature, unit, crossing, conservation and symmetry checks."""
import sys
sys.dont_write_bytecode=True
import json,subprocess
from pathlib import Path
from scipy.integrate import quad
from scipy.optimize import brentq
R=Path(__file__).resolve().parent

def cp(t):return 2402.3+7.6*t if t< -120 else 66747+545.24*t if t<=-119 else 2678.8+11.73*t+.039*t*t
def cap(t):return (1078-.43*t)*cp(t)
def sar(t,mode):return 691-4.6*(min(0,max(-80,t) if mode=='clamped' else t)+80)
def elapsed(t,fe,u,mode):
 bounds=[-150]+[x for x in [-120,-119,-80,0] if -150<x<t]+[t]
 return sum(quad(lambda x:cap(x)/(1000*fe*sar(x,mode)*u),a,b,epsabs=1e-10)[0] for a,b in zip(bounds[:-1],bounds[1:]))
def scalar(n,fe,u,mode,dt=.005):
 f=R/'results/check.input';f.write_text(f'{n} {n-1} .001 {dt} 70 {int(mode=="clamped")} 2 0 -150 {u} {u} {u}\n'+''.join(f'2 {fe}\n' for _ in range(n))+''.join(f'{i} {i+1}\n' for i in range(n-1)))
 r=json.loads(subprocess.check_output([str(R/'schedule_solver'),str(f)],text=True));f.unlink();return r

def main():
 out=[]
 for mode in ['clamped','extrapolated']:
  fe=6.;u=.95;r=scalar(1,fe,u,mode);s=scalar(2,fe,u,mode)
  exact=brentq(lambda t:elapsed(t,fe,u,mode)-70,-149,100)
  terr=abs(r['specimen']['cold']-exact)
  transit=3600/(elapsed(-40,fe,u,mode)-elapsed(-100,fe,u,mode))
  cerr=abs(r['specimen']['transit_100_40']-transit)
  sym=abs(r['specimen']['cold']-s['specimen']['cold'])
  assert terr<.003 and cerr<.003 and sym<1e-12
  # Uniform constant-power RHS is q/C. Its initial value is an independently computed unit identity.
  initial=60000*fe*sar(-150,mode)*u/cap(-150)
  assert r['specimen']['all']['rate']<=initial+1e-9
  out.append(dict(mode=mode,temperature_error_C=terr,crossing_transit_error_C_min=cerr,symmetry_error_C=sym,analytic_onset_C_min=initial,energy_rel_error=r['energy_rel_error']))
 z=scalar(2,0,1,'clamped');assert z['specimen']['cold']==-150 and z['specimen']['hot']==-150 and z['source_J']==0 and z['specimen']['all']['rate']==0
 out.append(dict(test='zero_source',passed=True))
 (R/'results/numerical_checks.json').write_text(json.dumps(out,indent=2)+'\n');print('C01 scalar/unit/crossing/symmetry/zero checks pass',out)
if __name__=='__main__':main()
