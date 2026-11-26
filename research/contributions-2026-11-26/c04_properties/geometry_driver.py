"""C04 fixed design input builder adapted from Longevidence reproduce.py (MIT)."""
import sys, json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
UP=ROOT/'vendor'
sys.dont_write_bytecode=True
sys.path.insert(0,str(UP))
from vitrify.geometry import rat_heart_cryobag

params=json.loads((ROOT/'inputs/grid_convergence_margin.json').read_text())['design_parameters']
flat=(params[3]+2*params[4]+params[5])/4
def write_input(f, mm, design, dt):
    p=rat_heart_cryobag(mm); active=p.active_mask
    labels=p.labels[active]; coords=np.meshgrid(*[np.arange(n)+.5 for n in p.shape],indexing='ij')
    rad=np.sqrt(sum(((coords[j]-p.shape[j]/2)/(p.spec.heart_dims_mm[j]/(2*mm)))**2 for j in range(3)))
    weights=np.stack([np.exp(-.5*((rad-c)/w)**2) for c,w in zip([1.05,1.45,2],[.32,.42,.65])]);weights/=weights.sum(axis=0)
    iron=np.sum(weights*np.array(params[:3])[:,None,None,None],axis=0)[active]
    iron[labels==0]=0;iron[labels==3]=10;iron[labels==2]=1.83
    index=np.full(p.shape,-1,dtype=int);index[active]=np.arange(active.sum());faces=[]
    for axis in range(3):
        left=[slice(None)]*3;right=left.copy();left[axis]=slice(None,-1);right[axis]=slice(1,None)
        a=index[tuple(left)];b=index[tuple(right)];ok=(a>=0)&(b>=0);faces.extend(zip(a[ok],b[ok]))
    controls=params[3:] if design=='ramp' else [flat]*3
    with f.open('w') as stream:
        stream.write(f'{len(labels)} {len(faces)} {mm*.001:.17g} {dt} 70 1 2 0 -150 '+' '.join(map(str,controls))+'\n')
        for lab,fe in zip(labels,iron):stream.write(f'{lab} {fe:.17g}\n')
        for a,b in faces:stream.write(f'{a} {b}\n')
