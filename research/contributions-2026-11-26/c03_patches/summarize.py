"""C03 deterministic contrasts, geometry audit, and local perturbation summaries."""
import csv,json
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parent/'results'
def write(name,rows):
 with (R/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=sorted(set().union(*(x.keys() for x in rows))));w.writeheader();w.writerows(rows)
def summarize():
 rows=json.loads((R/'summary.json').read_text());contrasts=[];maps=[];local=[]
 metrics=['myo_cold','whole_cold','myo_post_longest_s','whole_post_longest_s','myo_window_longest_s','whole_window_longest_s','myo_window_rate','whole_window_rate']
 for mm in [.75,.5]:
  for vol in [.01,.04]:
   trio=[r for r in rows if r['mm']==mm and r['target_patch_ml']==vol and r['name'] in ['chamber','deep','exterior'] and r['requested_dt']==.01]
   if not trio:continue
   by={r['name']:r for r in trio}
   for a,b in [('deep','chamber'),('deep','exterior'),('exterior','chamber')]:
    x,y=by[a],by[b];contrasts.append(dict(mm=mm,volume=vol,pair=a+' minus '+b,censored_post=(x['myo_post_missing']+y['myo_post_missing'])>0,**{k:x[k]-y[k] for k in metrics}))
   base=np.genfromtxt(R/f'{mm}_uniform_0_0.01_1.47.voxels.csv',delimiter=',',names=True)
   for name,r in by.items():
    tag=f'{mm}_{name}_{vol}_0.01_1.47';idx=np.load(R/(tag+'.patch.npz'))['indices'];v=np.genfromtxt(R/(tag+'.voxels.csv'),delimiter=',',names=True)
    local.append(dict(mm=mm,volume=vol,name=name,patch_mean_terminal_delta=float(np.mean(v['terminal'][idx]-base['terminal'][idx])),patch_max_terminal_depression=float(np.max(base['terminal'][idx]-v['terminal'][idx])),patch_terminal_min=r['patch_terminal_min']))
    for other in by:
     idx2=np.load(R/f'{mm}_{other}_{vol}_0.01_1.47.patch.npz')['indices']
     maps.append(dict(mm=mm,volume=vol,name=name,other=other,patch_n=len(idx),overlap_fraction=len(np.intersect1d(idx,idx2))/len(idx),seed_xyz_mm=r['seed_xyz_mm'],patch_chamber_distance_mean_mm=r['patch_chamber_distance_mean_mm'],patch_exterior_distance_mean_mm=r['patch_exterior_distance_mean_mm'],components=r['components'],patch_ml=r['patch_ml'],deficit_mg=r['deficit_mg'],muscle_mean=r['muscle_mean'],iron_g=r['iron_g']))
 temporal=[]
 for name in ['chamber','deep','exterior']:
  a=next(r for r in rows if r['mm']==.5 and r['name']==name and r['requested_dt']==.01)
  b=next(r for r in rows if r['mm']==.5 and r['name']==name and r['requested_dt']==.005)
  temporal.append(dict(name=name,**{k:b[k]-a[k] for k in metrics}))
 write('contrasts.csv',contrasts);write('maps.csv',maps);write('local_effects.csv',local);write('time_refinement.csv',temporal)
 return dict(id='C03',placement_contrasts=contrasts,local_effects=local,time_refinement=temporal)
if __name__=='__main__':print(json.dumps(summarize(),indent=2))
