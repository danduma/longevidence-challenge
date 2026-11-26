#!/usr/bin/env python3
"""C02 independent boundary-flux unit fixture for the compiled integrator."""
import json
import math
from pathlib import Path
import subprocess
import tempfile
ROOT=Path(__file__).resolve().parent

def check(executable):
    # A single source-free compartment cell has six exposed faces. Use one RK2 step.
    dx,dt,h,amb = .0005,.001,10.,22.
    def capacity(t):
        return (1078-.43*t)*(2402.3+7.6*t)
    def power(t):
        return 6*dx**2*h*(amb-t)
    start=-150.
    mid=start+.5*dt*power(start)/(dx**3*capacity(start))
    final=start+dt*power(mid)/(dx**3*capacity(mid))
    with tempfile.NamedTemporaryFile(mode='w',dir=ROOT) as f:
        f.write(f'1 0 {dx} {dt} {dt} 1 2 {h} {amb} 0 0 0\n2 0\n')
        f.flush()
        r=json.loads(subprocess.check_output([str(executable),f.name]))
    assert math.isclose(r['boundary_J'],dt*power(mid),rel_tol=1e-12,abs_tol=1e-18)
    assert math.isclose(r['cold'],final,rel_tol=0,abs_tol=1e-11)
    assert r['source_J']==0 and r['missing']==1
    print('C02 six-face unit fixture passed',r['boundary_J'],'J;',r['cold'],'deg C')
if __name__=='__main__':
    check(ROOT/'build/boundary')
