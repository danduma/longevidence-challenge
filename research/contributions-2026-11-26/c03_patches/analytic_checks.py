"""C03 one-cell quadrature/crossing check, independent of full-field comparison."""
import json,subprocess
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parent

def check():
 def capacity(t):
  cp=2402.3+7.6*t if t<-120 else 66747+545.24*t if t<=-119 else 2678.8+11.73*t+.039*t*t
  return (1078-.43*t)*cp
 def source(t):return 1000*5*(691-4.6*(max(-80,min(0,t))+80))*.9
 def time_to(t):
  cuts=sorted(set([-150]+[v for v in [-120,-119,-80,0] if -150<v<t]+[t]))
  return sum(quad(lambda v:capacity(v)/source(v),a,b,epsabs=1e-10)[0] for a,b in zip(cuts[:-1],cuts[1:]))
 exact=brentq(lambda t:time_to(t)-100,-100,50)
 inp=ROOT/'results/scalar.input';vox=ROOT/'results/scalar.voxels.csv'
 inp.write_text('1 0 0.0005 0.001 100 1 2 0 -150 .9 .9 .9\n2 5\n')
 r=json.loads(subprocess.check_output([str(ROOT/'solver'),str(inp),str(vox)],text=True))
 v=np.genfromtxt(vox,delimiter=',',names=True)
 errors={key:abs(float(v[key])-time_to(temp)) for key,temp in [('glass',-123),('endglass',-35),('entry',-100),('leave',-40)]}
 assert max(errors.values())<1e-4,errors
 assert abs(float(v['terminal'])-exact)<1e-4
 result=dict(id='C03',crossing_absolute_errors_s=errors,terminal_error_C=abs(float(v['terminal'])-exact),energy_relative_error=r['energy_rel_error'])
 (ROOT/'results/analytic_checks.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':print(json.dumps(check(),indent=2))
