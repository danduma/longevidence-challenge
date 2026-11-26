#!/usr/bin/env python3
"""C01 local exact-horizon equal-command family benchmark. MIT; see LICENSE."""
import sys
sys.dont_write_bytecode=True
import argparse, csv, hashlib, json, platform, subprocess, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
from vitrify.geometry import rat_heart_cryobag
P=json.loads((ROOT/'inputs/grid_convergence_margin.json').read_text())['design_parameters']
D=4-(P[3]+2*P[4]+P[5]); FLAT=(4-D)/4
FAMILY={'original':P[3:],'flat':[FLAT]*3,
 'deficit_start':[1-D,1,1],'deficit_mid':[1,1-D/2,1], 'deficit_end':[1,1,1-D],
 'edge_start_mid':[1-D/2,1-D/4,1], 'edge_start_end':[1-D/2,1,1-D/2],
 'edge_mid_end':[1,1-D/4,1-D/2], 'barycenter':[1-D/3,1-D/6,1-D/3]}
SUBSET=['original','flat','deficit_start','deficit_mid','deficit_end']
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare(mm,mean):
 p=rat_heart_cryobag(mm); active=p.active_mask; labels=p.labels[active]
 coords=np.meshgrid(*[np.arange(n)+.5 for n in p.shape],indexing='ij')
 rad=np.sqrt(sum(((coords[j]-p.shape[j]/2)/(p.spec.heart_dims_mm[j]/(2*mm)))**2 for j in range(3)))
 weights=np.stack([np.exp(-.5*((rad-c)/w)**2) for c,w in zip([1.05,1.45,2],[.32,.42,.65])]);weights/=weights.sum(axis=0)
 iron=np.sum(weights*np.array(P[:3])[:,None,None,None],axis=0)[active]
 iron[labels==0]=0;iron[labels==3]=10;iron[labels==2]=mean
 index=np.full(p.shape,-1,dtype=int);index[active]=np.arange(active.sum());faces=[]
 for axis in range(3):
  left=[slice(None)]*3;right=left.copy();left[axis]=slice(None,-1);right[axis]=slice(1,None)
  a=index[tuple(left)];b=index[tuple(right)];ok=(a>=0)&(b>=0);faces.extend(zip(a[ok],b[ok]))
 body=''.join(f'{lab} {fe:.17g}\n' for lab,fe in zip(labels,iron))+''.join(f'{a} {b}\n' for a,b in faces)
 ijk=np.argwhere(active)
 return p,labels,iron,faces,body,ijk
GEOM={}
def run(mm,name,mean,mode,dt,phase,solver='schedule_solver'):
 key=(mm,mean)
 if key not in GEOM:GEOM[key]=prepare(mm,mean)
 p,labels,iron,faces,body,ijk=GEOM[key]
 tag=f'{phase}_{mm}_{mean}_{mode}_{name}_{dt}'
 f=ROOT/'results'/f'{tag}.input'; controls=FAMILY[name]
 f.write_text(f'{len(labels)} {len(faces)} {mm*.001:.17g} {dt} 70 {int(mode=="clamped")} 2 0 -150 '+' '.join(f'{u:.17g}' for u in controls)+'\n'+body)
 start=time.monotonic();r=json.loads(subprocess.run([str(ROOT/solver),str(f)],text=True,capture_output=True,check=True).stdout)
 r.update(id='C01',phase=phase,mm=mm,mean=mean,mode=mode,schedule=name,controls=controls,input_sha256=digest(f),elapsed_s=time.monotonic()-start,volumes=p.report())
 if solver=='schedule_solver':
  for domain in ['specimen','myocardium']:
   for w in ['all','glass','glass35','window']:
    m=r[domain][w];m['ijk']=ijk[m['cell']].tolist();m['xyz_mm']=((ijk[m['cell']]+.5)*mm).tolist()
  initial=60000*mean*(691 if mode=='clamped' else 1013)*controls[0]/((1078+.43*150)*(2402.3-7.6*150))
  r['analytic_myocardium_onset_rate']=initial
 (ROOT/'results'/f'{tag}.json').write_text(json.dumps(r,indent=2)+'\n');f.unlink()
 print(f'C01 {tag}: '+(f"all={r['specimen']['all']['rate']:.5f} glass={r['specimen']['glass']['rate']:.5f}" if solver=='schedule_solver' else str(r['min_rate'])),flush=True)
 return r

def compile_all():
 for src,exe in [('schedule.cpp','schedule_solver'),('inputs/original_independent.cpp','baseline_solver')]:subprocess.run(['c++','-O3','-std=c++17',str(ROOT/src),'-o',str(ROOT/exe)],check=True)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=2);args=ap.parse_args()
 compile_all()
 assert digest(ROOT/'protocol.md')==(ROOT/'inputs/protocol.sha256').read_text().split()[0]
 for u in FAMILY.values():assert all(0<=x<=1 for x in u) and abs(17.5*(u[0]+2*u[1]+u[2])-17.5*(4-D))<1e-12
 (ROOT/'results/family.json').write_text(json.dumps(FAMILY,indent=2)+'\n')
 baselines=[run(.5,n,1.83,'clamped',.01,'baseline',solver='baseline_solver') for n in ['original','flat']]
 assert abs(baselines[0]['min_rate']-44.4358)<.03 and abs(baselines[1]['min_rate']-50.3155)<.03
 assert all(abs(r['energy_rel_error'])<1e-5 and r['missing']==0 for r in baselines)
 (ROOT/'results/baseline_verified.json').write_text(json.dumps(baselines,indent=2)+'\n')
 tasks=[]
 for mean in [1.47,1.83]:
  for mode in ['clamped','extrapolated']:
   tasks += [(.75,n,mean,mode,.02,'discovery') for n in FAMILY]
   tasks += [(.5,n,mean,mode,.01,'spatial') for n in SUBSET]
 for mode in ['clamped','extrapolated']:tasks += [(.5,n,1.83,mode,.005,'time') for n in SUBSET]
 # Precompute shared geometry before threads; workers only write unique local paths.
 for mm,_,mean,_,_,_ in tasks:
  if (mm,mean) not in GEOM:GEOM[(mm,mean)]=prepare(mm,mean)
 with ThreadPoolExecutor(max_workers=args.jobs) as pool: rows=list(pool.map(lambda t:run(*t),tasks))
 (ROOT/'results/all_runs.json').write_text(json.dumps(rows,indent=2)+'\n')
 versions={'id':'C01','python':sys.version,'numpy':np.__version__,'platform':platform.platform(),'compiler':subprocess.check_output(['c++','--version'],text=True).splitlines()[0]}
 (ROOT/'results/versions.json').write_text(json.dumps(versions,indent=2)+'\n')
 subprocess.run([sys.executable,str(ROOT/'verify.py')],check=True)
if __name__=='__main__':main()
