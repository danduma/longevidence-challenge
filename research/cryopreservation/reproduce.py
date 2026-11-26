"""Review calculations: own loading basis, topology and solver; reuse geometry only.
Run with PYTHONDONTWRITEBYTECODE=1 ../../.venv/bin/python run_review.py
"""
import sys, json, subprocess, hashlib, platform, csv
from pathlib import Path
import numpy as np
from scipy.ndimage import distance_transform_edt
ROOT=Path(__file__).resolve().parent
UP=ROOT/'vendor'
sys.dont_write_bytecode=True
sys.path.insert(0,str(UP))
from vitrify.geometry import rat_heart_cryobag
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
params=json.loads((ROOT/'inputs/grid_convergence_margin.json').read_text())['design_parameters']
flat=(params[3]+2*params[4]+params[5])/4
def run(mm,design,mean=1.47,amplitude=0,sign=1,mode='extrapolated',dt=.01,h=0,env=-150,scheme=2,solver='independent'):
    p=rat_heart_cryobag(mm); active=p.active_mask
    labels=p.labels[active]; coords=np.meshgrid(*[np.arange(n)+.5 for n in p.shape],indexing='ij')
    rad=np.sqrt(sum(((coords[j]-p.shape[j]/2)/(p.spec.heart_dims_mm[j]/(2*mm)))**2 for j in range(3)))
    weights=np.stack([np.exp(-.5*((rad-c)/w)**2) for c,w in zip([1.05,1.45,2],[.32,.42,.65])]);weights/=weights.sum(axis=0)
    iron=np.sum(weights*np.array(params[:3])[:,None,None,None],axis=0)[active]
    iron[labels==0]=0;iron[labels==3]=10;iron[labels==2]=mean
    muscle=np.flatnonzero(labels==2)
    distance=distance_transform_edt(~((p.labels==1)|(p.labels==3)))[active]
    order=muscle[np.lexsort((muscle,distance[muscle]))]; nh=len(order)//2
    if amplitude:iron[order[:nh]]+=sign*amplitude;iron[order[-nh:]]-=sign*amplitude
    index=np.full(p.shape,-1,dtype=int);index[active]=np.arange(active.sum());faces=[]
    for axis in range(3):
        left=[slice(None)]*3;right=left.copy();left[axis]=slice(None,-1);right[axis]=slice(1,None)
        a=index[tuple(left)];b=index[tuple(right)];ok=(a>=0)&(b>=0);faces.extend(zip(a[ok],b[ok]))
    controls=params[3:] if design=='ramp' else [flat]*3
    name=f'{mm}_{design}_{mean}_{amplitude}_{sign}_{mode}_{dt}_{h}_{env}_{scheme}'+('' if solver=='independent' else '_'+solver)
    f=OUT/'input.txt'
    with f.open('w') as stream:
        stream.write(f'{len(labels)} {len(faces)} {mm*.001:.17g} {dt} 70 {int(mode=="clamped")} {scheme} {h} {env} '+' '.join(map(str,controls))+'\n')
        for lab,fe in zip(labels,iron):stream.write(f'{lab} {fe:.17g}\n')
        for a,b in faces:stream.write(f'{a} {b}\n')
    digest=hashlib.sha256(f.read_bytes()).hexdigest()
    proc=subprocess.run([str(ROOT/solver),str(f)],text=True,capture_output=True,check=True)
    data=json.loads(proc.stdout);data.update(mm=mm,design=design,mean=mean,amplitude=amplitude,sign=sign,mode=mode,dt_requested=dt,h=h,env=env,scheme=scheme,input_sha256=digest,muscle_cells=len(muscle),muscle_mean=float(iron[muscle].mean()),volumes=p.report())
    (OUT/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n')
    print(name,proc.stdout.strip(),flush=True)
    return data
if __name__=='__main__':
    subprocess.run(['c++','-O3','-std=c++17',str(ROOT/'independent.cpp'),'-o',str(ROOT/'independent')],check=True)
    rows=[]
    for mode,mean in [('clamped',1.83),('extrapolated',1.47)]:
        for design in ['ramp','flat']: rows.append(run(.5,design,mean=mean,mode=mode))
    assert not rows[0]['pass'] and rows[1]['pass']
    assert rows[1]['window_rate']<rows[0]['window_rate']
    assert rows[3]['min_rate']>rows[2]['min_rate'] and rows[3]['window_rate']<rows[2]['window_rate']
    assert all(r['missing']==0 and abs(r['energy_rel_error'])<1e-6 for r in rows)
    (OUT/'core.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('PASS: four independently integrated thermal cases. No biological validation claimed.')
