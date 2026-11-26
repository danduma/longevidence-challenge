"""Nonlinear adiabatic quadrature and null-source limits on an actual multi-cell grid."""
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from investigate import coefficients,library

rows=[]
# All 3 connected cells carry the same source and material; face flux must vanish.
for mode in ['clamped','extrapolated']:
    for iron in [0.,1.47,10.]:
        sar=lambda t:691.-4.6*(max(-80,t) + 80 if mode=='clamped' else min(t,0)+80)
        for dt in [.02,.01,.005]:
            steps=int(20/dt)
            lab=np.ones(3,dtype='i')
            faces=np.array([[0,1],[1,2]],dtype='i')
            out=np.zeros(15);end=np.zeros(3)
            library(3,2,faces,lab,np.full(3,iron),.001,dt,steps,np.ones(steps),int(mode=='clamped'),0.,100,out,end)
            def elapsed(T):
                return quad(lambda t: coefficients(np.array(t),False)[0]/(iron*1000*sar(t)),-150,T,
                            points=[x for x in [-120,-119,-80] if x<T],epsabs=1e-9)[0]
            target=brentq(lambda t:elapsed(t)-20,-150,-.01) if iron else -150.
            rows.append(dict(mode=mode,iron=iron,dt=dt,expected=target,calculated=end[0],error=end[0]-target,
                             spread=np.ptp(end),source_J=out[13]))
f=pd.DataFrame(rows)
f.to_csv('outputs/analytic.csv',index=False)
print(f.to_string(index=False))
assert f.spread.max()==0 and f.query('iron==0').error.abs().max()==0 and f.query('dt==.005').error.abs().max()<.02
