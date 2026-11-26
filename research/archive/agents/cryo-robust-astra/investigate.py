"""Run from the job directory. No dependency on JAX or container services."""
import ast
import ctypes
import hashlib
import json
import sys
from pathlib import Path
from functools import lru_cache

import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial as Poly
from scipy.ndimage import distance_transform_edt
import numpy.typing as npt

U = Path('upstream')
sys.path.insert(0,str(U/'src'))
from vitrify.geometry import rat_heart_cryobag, Phantom, CHAMBER, MYOCARDIUM, SURROUNDING_CPA

# Execute only the pinned pure NumPy design class, excluding Docker client imports.
tree = ast.parse((U/'coupler/workflow.py').read_text())
node = next(n for n in tree.body if getattr(n,'name',None)=='DesignParameterization')
CHAMBER_LOADING_MG_PER_ML, MYOCARDIUM_LOADING_MG_PER_ML = 10., 1.47
exec(compile(ast.Module(body=[node],type_ignores=[]),'upstream/coupler/workflow.py','exec'))

def coefficients(t,wall):
    if wall:
        return (931-.49*t-.00091*t*t)*(1031+3.57*t-.014*t*t-.000061710*t**3), .39+.00052*t
    cp=np.select([t< -120,t<=-119],[7.6*t+2402.3,545.24*t+66747],default=.039*t*t+11.73*t+2678.8)
    k=np.select([t<=-116,t< -96,t<=-77,t< -34],[.000816*t+.440,-.001*t+.234,.000459*t+.374,.0005*t+.374],default=.000459*t+.374)
    return (1078-.43*t)*cp,k

def enthalpy(t,wall):
    """Exact integral of piecewise volumetric capacity from storage to T."""
    if wall:
        p=(Poly([931,-.49,-.00091])*Poly([1031,3.57,-.014,-.000061710])).integ()
        return p(t)-p(-150)
    val=np.zeros_like(t)
    for lo,hi,cp in [(-150,-120,[2402.3,7.6]),(-120,-119,[66747,545.24]),(-119,np.inf,[2678.8,11.73,.039])]:
        p=(Poly([1078,-.43])*Poly(cp)).integ()
        val+=p(np.clip(t,lo,hi))-p(lo)
    return val

def read(name):
    return json.loads((U/'results'/f'{name}.json').read_text())

margin=read('grid_convergence_margin')['design_parameters']
equal=(margin[3]+2*margin[4]+margin[5])/4
designs={'caseII':[7.5]*3+[1.]*3,'original':read('grid_convergence')['design_parameters'],
         'margin':margin,'equal':margin[:3]+[equal]*3,'caseI':[10.]*3+[1.]*3,'caseIII':[5.]*3+[1.]*3}

@lru_cache(None)
def mesh(mm):
    p=rat_heart_cryobag(mm)
    active=p.active_mask
    number=np.full(active.shape,-1,dtype='i')
    number[active]=np.arange(active.sum())
    pairs=[]
    for axis in range(3):
        a=np.moveaxis(number,axis,0)
        left,right=a[:-1],a[1:]
        keep=(left>=0)&(right>=0)
        pairs.append(np.column_stack((left[keep],right[keep])))
    labels=p.labels[active].astype('i')
    dist=distance_transform_edt(~p.mask(CHAMBER,SURROUNDING_CPA))[active]
    muscle=np.flatnonzero(labels==MYOCARDIUM)
    order=muscle[np.argsort(dist[muscle],kind='stable')]
    offsets=np.zeros(len(labels))
    offsets[order[:len(order)//2]]=0.36
    offsets[order[-(len(order)//2):]]=-0.36
    return p,DesignParameterization(p),np.ascontiguousarray(np.concatenate(pairs)),labels,offsets

def setup(mm,design,scenario,half,duration):
    p,parameter,faces,labels,offsets=mesh(mm)
    fe,controls=parameter.fields(designs[design])
    fe=np.ascontiguousarray(fe[p.active_mask])
    muscle=labels==MYOCARDIUM
    fe[muscle]={'low':1.11,'nominal':1.47,'high':1.83,'deepLow':1.47,'deepHigh':1.47}[scenario]
    if scenario.startswith('deep'):
        fe+=offsets*(1 if scenario=='deepLow' else -1)
        assert abs(fe[muscle].mean()-1.47)<1e-12
    if design=='caseI':fe[muscle]=10
    t=np.linspace(-180,20,400)
    mats=[coefficients(t,w) for w in [True,False]]
    dt=.4*min(c.min() for c,k in mats)*p.spacing_m**2/(6*max(k.max() for c,k in mats))/half
    steps=int(np.ceil(duration/dt))
    duty=np.interp(np.linspace(0,1,steps),[0,.5,1],controls)
    return p,faces,labels,fe,dt,duty

library=ctypes.CDLL(str(Path('thermal.so').resolve())).solve
integer=ctypes.c_int
real=ctypes.c_double
ip=np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS')
dp=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
library.argtypes=[integer,integer,ip,ip,dp,real,real,integer,dp,integer,real,integer,dp,dp]
library.restype=None

def run(mm,design,scenario='nominal',sar='extrapolated',half=1,duration=70):
    p,faces,labels,fe,dt,duty=setup(mm,design,scenario,half,duration)
    end=np.zeros(len(fe))
    out=np.zeros(15)
    library(len(fe),len(faces),faces,labels,fe,p.spacing_m,dt,len(duty),duty,int(sar=='clamped'),15. if design=='caseIII' else 0.,
            int(np.ceil(len(duty)/120)),out,end)
    v,t,temp,i,w,wt,wtemp,j,cold,hot,mean,transit,missing,joules,glass=out
    H=sum(enthalpy(end[labels==0],True))+sum(enthalpy(end[labels>0],False))
    err=H*p.cell_volume_m3/joules-1
    idx=np.argwhere(p.active_mask)
    return dict(mm=mm,design=design,scenario=scenario,sar=sar,dt_divisor=half,dt_s=dt,actual_end_s=len(duty)*dt,
                min_rate=v,min_time=t,min_temp=temp,min_label=labels[int(i)],min_ijk=idx[int(i)].tolist(),
                window_rate=w,window_time=wt,window_temp=wtemp,window_label=labels[int(j)],window_ijk=idx[int(j)].tolist(),
                cold=cold,hot=hot,mean_rate=(mean+150)*60/(len(duty)*dt),transit_rate=transit,incomplete=int(missing),energy_J=joules,
                enthalpy_relative_error=err,glass_rate=glass,passes=v>=50 and cold>=-35 and hot<=0,
                window_passes=w>=50 and cold>=-35 and hot<=0)

def main():
    gates=[]
    for mode,values in read('thermal_validation')['readings'].items():
        for case,expected in values.items():
            r=run(.5,'case'+case,sar=mode,duration=30)
            gates.append((case,mode,.5,'glass_rate',r['glass_rate'],expected['mean_rate_c_per_min']))
    for file,design in [('grid_convergence','original'),('grid_convergence_margin','margin')]:
        for level in read(file)['levels']:
            for key,d in [('baseline','caseII'),('optimized',design)]:
                r=run(level['grid_spacing_mm'],d)
                for s,e in [('min_rate','min_warming_rate_c_per_min'),('cold','final_min_temperature_c'),('hot','final_max_temperature_c'),('mean_rate','mean_warming_rate_c_per_min')]:
                    gates.append((d,key,r['mm'],s,r[s],level['designs'][key][e]))
    g=pd.DataFrame(gates,columns=['design','comparison','mm','metric','calculated','committed'])
    g['error']=g.calculated-g.committed
    g.to_csv('outputs/baselines.csv',index=False)
    print('baseline maximum absolute discrepancy',g.error.abs().max(),flush=True)
    assert g.error.abs().max()<.02
    rows=[]
    for mm in [1.4,.9,.75,.5]:
        for d in ['caseII','original','margin','equal']:
            for s in ['low','nominal','high','deepLow','deepHigh']:
                for mode in ['extrapolated','clamped']:
                    rows.append(run(mm,d,s,mode))
        pd.DataFrame(rows).to_csv('outputs/sweep.csv',index=False)
        print('finished',mm,flush=True)
    for d in ['caseII','margin','equal','original']:
        for s in ['nominal','deepLow','deepHigh','low','high']:
            rows.append(run(.5,d,s,half=2))
    frame=pd.DataFrame(rows)
    frame.to_csv('outputs/sweep.csv',index=False)
    print('max energy discrepancy',frame.enthalpy_relative_error.abs().max(),flush=True)
    assert frame.enthalpy_relative_error.abs().max()<.005

if __name__=='__main__':main()
