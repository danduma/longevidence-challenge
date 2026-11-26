#!python3
"""C03 executable full reproduction using only pinned local snapshots."""
import os, sys
os.environ['PYTHONDONTWRITEBYTECODE']='1'
sys.dont_write_bytecode=True
from pathlib import Path
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
import json, hashlib, subprocess, csv, platform
import numpy as np
import scipy
from scipy.ndimage import distance_transform_edt, label
sys.path.insert(0,str(ROOT/'snapshots/vendor'))
from vitrify.geometry import rat_heart_cryobag
P=json.loads(Path('snapshots/grid_convergence_margin.json').read_text())['design_parameters']
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def setup(mm):
 p=rat_heart_cryobag(mm); active=p.active_mask; labs=p.labels[active]; muscle=np.flatnonzero(labs==2)
 xyz=np.stack(np.meshgrid(*[(np.arange(n)+.5)*mm for n in p.shape],indexing='ij'),axis=-1)[active]
 # Preserve released radial shell convention, based on grid midpoint, as in portable runner.
 coords=np.meshgrid(*[np.arange(n)+.5 for n in p.shape],indexing='ij')
 rad=np.sqrt(sum(((coords[j]-p.shape[j]/2)/(p.spec.heart_dims_mm[j]/(2*mm)))**2 for j in range(3)))
 w=np.stack([np.exp(-.5*((rad-c)/v)**2) for c,v in zip([1.05,1.45,2],[.32,.42,.65])]);w/=w.sum(axis=0)
 iron=np.sum(w*np.array(P[:3])[:,None,None,None],axis=0)[active];iron[labs==0]=0;iron[labs==3]=10;iron[labs==2]=1.47
 idx=np.full(p.shape,-1,int);idx[active]=np.arange(len(labs));faces=[]
 for ax in range(3):
  a=[slice(None)]*3;b=a.copy();a[ax]=slice(None,-1);b[ax]=slice(1,None);aa=idx[tuple(a)];bb=idx[tuple(b)];ok=(aa>=0)&(bb>=0);faces.extend(zip(aa[ok],bb[ok]))
 dc=distance_transform_edt(p.labels!=3,sampling=mm)[active];de=distance_transform_edt(p.labels!=1,sampling=mm)[active]
 center=np.array(p.spec.container_dims_mm)/2; eligible=muscle[(xyz[muscle,1]>=center[1])&(abs(xyz[muscle,0]-center[0])<=1)&(abs(xyz[muscle,2]-center[2])<=1)]
 tie=np.sum((xyz[eligible]-(center+[0,3.5,0]))**2,axis=1)
 seeds={}
 for name,criterion in [('chamber',dc),('deep',-np.minimum(dc,de)),('exterior',de)]:
  seeds[name]=int(eligible[np.lexsort((eligible,tie,criterion[eligible]))[0]])
 return p,labs,muscle,xyz,iron,faces,seeds,dc,de

def run(mm,name,volume=0,dt=.01,mean=1.47,zero=False):
 tag=f'{mm}_{name}_{volume}_{dt}_{mean}'
 if '--resume' in sys.argv and (OUT/(tag+'.json')).exists():
  cached=json.loads((OUT/(tag+'.json')).read_text()); assert digest(OUT/(tag+'.voxels.csv'))==cached['voxel_sha256']; return cached
 p,labs,m,xyz,iron,faces,seeds,dc,de=setup(mm); iron[m]=mean
 n=round(volume*1000/mm**3);patch=np.array([],int)
 meta={}
 if name in seeds:
  seed=seeds[name];order=m[np.lexsort((m,np.sum((xyz[m]-xyz[seed])**2,axis=1)))];patch=order[:n]
  iron[m]+= .36*n/(len(m)-n);iron[patch]=mean-.36
  mask=np.zeros(p.shape,bool);tmp=np.zeros(len(labs),bool);tmp[patch]=True;mask[p.active_mask]=tmp
  meta=dict(seed_xyz_mm=xyz[seed].tolist(),components=int(label(mask)[1]),patch_max_radius_mm=float(np.sqrt(np.max(np.sum((xyz[patch]-xyz[seed])**2,axis=1)))),patch_chamber_distance_mean_mm=float(dc[patch].mean()),patch_exterior_distance_mean_mm=float(de[patch].mean()))
  assert abs(iron[m].mean()-mean)<1e-12
  assert label(mask)[1]==1
 if name=='lower_mean':iron[m]=mean-.36*n/len(m)
 if zero:iron[:]=0
 tag=f'{mm}_{name}_{volume}_{dt}_{mean}'; inp=OUT/(tag+'.input');vox=OUT/(tag+'.voxels.csv')
 with inp.open('w') as f:
  f.write(f'{len(labs)} {len(faces)} {mm*.001:.17g} {dt} 70 1 2 0 -150 '+ ' '.join(map(str,P[3:]))+'\n')
  for l,fe in zip(labs,iron):f.write(f'{l} {fe:.17g}\n')
  for a,b in faces:f.write(f'{a} {b}\n')
 proc=subprocess.run([str(ROOT/'solver'),str(inp),str(vox)],capture_output=True,text=True,check=True)
 r=json.loads(proc.stdout); v=np.genfromtxt(vox,delimiter=',',names=True)
 r.update(mm=mm,name=name,target_patch_ml=volume,requested_dt=dt,mean_input=mean,patch_n=len(patch),target_n=n,muscle_n=len(m),patch_ml=len(patch)*mm**3/1000,deficit_mg=len(patch)*mm**3/1000*.36,muscle_mean=float(iron[m].mean()),input_sha256=digest(inp),voxel_sha256=digest(vox),**meta)
 for group,mask in [('myo',labs==2),('whole',labs>0)]:
  r[group+'_cold']=float(v['terminal'][mask].min());r[group+'_hot']=float(v['terminal'][mask].max())
  for key,entry,leave,span in [('post','glass','endglass',88),('window','entry','leave',60)]:
   ok=(v[entry]>=0)&(v[leave]>v[entry]); r[group+'_'+key+'_missing']=int(np.sum(mask&~ok))
   r[group+'_'+key+'_longest_s']=float(np.max(v[leave][mask&ok]-v[entry][mask&ok])) if np.any(mask&ok) else None
   r[group+'_'+key+'_rate']=float(np.min(v['postrate' if key=='post' else 'winrate'][mask]))
   assert np.all(v[leave][mask&ok]<=70) and np.all(v[entry][mask&ok]>=0)
 if len(patch):
  r['patch_terminal_min']=float(v['terminal'][patch].min());np.savez_compressed(OUT/(tag+'.patch.npz'),indices=patch,xyz=xyz[patch],fe=iron[patch])
 assert abs(r['energy_rel_error'])<1e-5
 if zero:assert np.all(v['terminal']==-150)
 (OUT/(tag+'.json')).write_text(json.dumps(r,indent=2)+'\n')
 inp.unlink() # deterministically regenerable; retain hash, voxel fields and compact maps
 print(tag,'cold',r['myo_cold'],'post_s',r['myo_post_longest_s'],'energy',r['energy_rel_error'],flush=True)
 return r

def main():
 for p,h in json.loads(Path('pinned_hashes.json').read_text()).items():assert digest(p)==h,p
 subprocess.run([sys.executable,'build_solver.py'],check=True)
 subprocess.run(['c++','-O3','-std=c++17','solver.cpp','-o','solver'],check=True)
 # Central reproduction must succeed before scoring any new cases.
 baseline=run(.5,'baseline',mean=1.83);old=json.loads(Path('snapshots/baseline.json').read_text())
 for key in ['min_rate','window_rate','cold','hot','transit','source_J']:
  assert abs(baseline[key]-old[key])<1e-8,(key,baseline[key],old[key])
 # Explicit source-unit check: 1 mg/mL =1000 g/m3; W/g * g/m3 gives W/m3.
 cap=(1078-.43*(-150))*(2402.3+7.6*(-150)); expected=60000*1.83*691*P[3]/cap
 assert abs(baseline['min_rate']-expected)<1e-10
 zero=run(.75,'zero',zero=True)
 rows=[]
 rows.append(run(.75,'uniform'))
 for volume in [.010,.040]:
  rows.append(run(.75,'lower_mean',volume))
  for name in ['chamber','deep','exterior']:rows.append(run(.75,name,volume))
 rows.append(run(.5,'uniform'))
 for dt in [.01,.005]:
  for name in ['chamber','deep','exterior']:rows.append(run(.5,name,.040,dt))
 # Incomplete crossings remain censored and are reported for both domains.
 for mm in [.75,.5]:
  trio=[r for r in rows if r['mm']==mm and r['target_patch_ml']==.04 and r['requested_dt']==.01 and r['name'] in ['chamber','deep','exterior']]
  assert len(set(r['patch_n'] for r in trio))==1
  assert max(r['iron_g'] for r in trio)-min(r['iron_g'] for r in trio)<1e-12
 keys=sorted(set().union(*(r.keys() for r in rows)))
 with (OUT/'summary.csv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)
 (OUT/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
 from summarize import summarize
 summarize()
 from analytic_checks import check
 analytic=check()
 checks={'analytic_checks':analytic,'id':'C03','baseline':'six metrics matched within 1e-8','source_unit_check':True,'zero_source_check':True,'max_energy_relative_error':max(abs(r['energy_rel_error']) for r in rows),'protocol_sha256':digest('protocol.md'),'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'compiler':subprocess.check_output(['c++','--version'],text=True).splitlines()[0]}
 (OUT/'verification.json').write_text(json.dumps(checks,indent=2)+'\n');print('C03 PASS',json.dumps(checks),flush=True)
if __name__=='__main__':main()
