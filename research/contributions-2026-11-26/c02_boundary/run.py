#!/usr/bin/env python3
"""C02 deterministic phase-gated envelope calculation. Run from any directory."""
import os
import sys
sys.dont_write_bytecode = True
import csv
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT/'pinned/vendor'))
from vitrify.geometry import rat_heart_cryobag
PAR = json.loads(Path('pinned/inputs/grid_convergence_margin.json').read_text())['design_parameters']
BUILD = Path('build')
BUILD.mkdir(exist_ok=True)
HS = [0, 2, 5, 10, 20, 40, 80]
AMES = [-150, -80, -20, 22]
FIELDS = ['rate_ok','cold_ok','hot_ok']
def geometry(mm):
    p = rat_heart_cryobag(mm)
    active = p.active_mask
    label = p.labels[active]
    coord = np.meshgrid(*[np.arange(n)+.5 for n in p.shape], indexing='ij')
    radius = np.sqrt(sum(((coord[j]-p.shape[j]/2)/(p.spec.heart_dims_mm[j]/(2*mm)))**2 for j in range(3)))
    w = np.stack([np.exp(-.5*((radius-c)/s)**2) for c,s in zip([1.05,1.45,2],[.32,.42,.65])])
    fe = np.sum((w/w.sum(axis=0))*np.array(PAR[:3])[:,None,None,None],axis=0)[active]
    for lab,value in [(0,0),(2,1.83),(3,10)]:
        fe[label==lab] = value
    ids = np.full(p.shape,-1)
    ids[active] = np.arange(len(label))
    faces = []
    for axis in range(3):
        arr = ids.swapaxes(0,axis)
        left,right = arr[:-1],arr[1:]
        valid = (left>=0)&(right>=0)
        faces.append(np.column_stack((left[valid],right[valid])))
    faces = np.concatenate(faces)
    exposed = 6-np.bincount(faces.ravel(),minlength=len(label))
    body = ''.join(f'{l} {f:.17g}\n' for l,f in zip(label,fe))+''.join(f'{a} {b}\n' for a,b in faces)
    return (len(label),len(faces),int(exposed.sum()),int(exposed[label!=0].sum()),body)

def calculate(mm,design,h,amb,dt):
    key = (mm,design,h,amb,dt)
    if key in cache:
        return cache[key]
    n,nf,ext,bare,body = geometries[mm]
    u = PAR[3:] if design=='ramp' else [(PAR[3]+2*PAR[4]+PAR[5])/4]*3
    path = BUILD / ('_'.join(map(str,key))+'.txt')
    path.write_text(f'{n} {nf} {mm*.001} {dt} 70 1 2 {h} {amb} '+ ' '.join(map(str,u))+'\n'+body)
    r = json.loads(subprocess.check_output(['build/boundary',str(path)]))
    path.unlink()
    r.update(mm=mm,design=design,h=h,ambient=amb,exposed_faces=ext,exposed_nonwall=bare,
             rate_ok=r['min_rate']>=50,cold_ok=r['cold']>=-35,hot_ok=r['hot']<=0)
    assert r['bad_cross']==0 and abs(r['energy_rel_error'])<1e-5
    assert abs(r['duty_s']-66.9762264567)<1e-8
    assert h or r['boundary_J']==0
    cache[key] = r
    return r

def envelope(task):
    mm,amb,design = task
    dt = .02 if mm==.75 else .01
    def at(h): return calculate(mm,design,h,amb,dt)
    brackets = []
    for lo0,hi0 in zip(HS,HS[1:]):
        for field in FIELDS+['pass']:
            lo,hi=lo0,hi0
            if at(lo)[field]==at(hi)[field]: continue
            while hi-lo>2:
                mid=(lo+hi)/2
                if at(mid)[field]==at(lo)[field]: lo=mid
                else: hi=mid
            brackets.append(dict(mm=mm,ambient=amb,design=design,criterion=field,
                                 lo=lo,hi=hi,lo_pass=at(lo)[field],hi_pass=at(hi)[field]))
    if design=='flat' and amb in [-150,22]:
        joint=[b for b in brackets if b['criterion']=='pass']
        for h in ([joint[0]['lo'],joint[0]['hi']] if joint else [80]):
            calculate(mm,design,h,amb,dt/2)
    return brackets

def table(path,rows):
    with Path(path).open('w') as out:
        writer=csv.DictWriter(out,fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
if __name__=='__main__':
    subprocess.run(['shasum','-a','256','-c','protocol.sha256'],check=True)
    subprocess.run(['c++','-O3','-std=c++17','boundary.cpp','-o','build/boundary'],check=True)
    cache={}
    geometries={mm:geometry(mm) for mm in [.5,.75]}
    for design,expected in [('ramp',[44.43576,-32.99651,-12.84565]),('flat',[50.31942,-32.66986,-12.99979])]:
        r=calculate(.5,design,0,-150,.01)
        assert np.max(np.abs(np.array([r[k] for k in ['min_rate','cold','hot']])-expected))<.015
        assert r['pass']==(design=='flat')
        if design=='ramp':
            rate=60000*1.83*691*PAR[3]/1442177.75
            assert abs(r['min_rate']-rate)<1e-6
        print('C02 baseline',design,r['min_rate'],r['cold'],r['hot'],flush=True)
    tasks=[(mm,amb,d) for mm in [.75,.5] for amb in (AMES if mm==.75 else [-150,22]) for d in ['flat','ramp']]
    with ThreadPoolExecutor(max_workers=int(os.environ.get('C02_WORKERS','4'))) as pool:
        brackets=sum(pool.map(envelope,tasks),[])
    table('results/envelope.csv',[cache[k] for k in sorted(cache)])
    table('results/brackets.csv',brackets)
    print('C02',len(cache),'scenarios including diagnostics;',len(brackets),'switch intervals')
