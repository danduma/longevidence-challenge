#!/usr/bin/env python3
"""C05 deterministic local verifier; all paths anchored to this job directory."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json, subprocess, hashlib, csv, platform, os
import numpy as np
import scipy
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
sys.path.insert(0,str(ROOT/'vendor'))
from vitrify.geometry import rat_heart_cryobag
OUT=ROOT/'results'
P=json.loads((ROOT/'inputs/grid_convergence_margin.json').read_text())['design_parameters']
DOMAINS={'wall':[0],'surrounding_CPA':[1],'myocardium':[2],'chamber_CPA':[3],'heart_contents':[2,3],'whole_specimen':[1,2,3]}
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def mesh(mm):
 p=rat_heart_cryobag(mm); active=p.active_mask; labels=p.labels[active]
 coords=np.meshgrid(*[np.arange(n)+.5 for n in p.shape],indexing='ij')
 rad=np.sqrt(sum(((coords[j]-p.shape[j]/2)/(p.spec.heart_dims_mm[j]/(2*mm)))**2 for j in range(3)))
 weights=np.stack([np.exp(-.5*((rad-c)/w)**2) for c,w in zip([1.05,1.45,2],[.32,.42,.65])]);weights/=weights.sum(axis=0)
 iron=np.sum(weights*np.array(P[:3])[:,None,None,None],axis=0)[active]
 iron[labels==0]=0;iron[labels==3]=10
 index=np.full(p.shape,-1,dtype=int);index[active]=np.arange(active.sum());faces=[]
 for axis in range(3):
  left=[slice(None)]*3;right=left.copy();left[axis]=slice(None,-1);right[axis]=slice(1,None)
  a=index[tuple(left)];b=index[tuple(right)];ok=(a>=0)&(b>=0);faces.extend(zip(a[ok],b[ok]))
 return p,labels,iron,np.array(faces),np.argwhere(active)
def write_input(path,labels,iron,faces,mm,dt,mode,controls,duration=70):
 with path.open('w') as f:
  f.write(f'{len(labels)} {len(faces)} {mm*.001:.17g} {dt} {duration} {int(mode=="clamped")} 2 0 -150 '+' '.join(map(str,controls))+'\n')
  for lab,fe in zip(labels,iron): f.write(f'{lab} {fe:.17g}\n')
  for a,b in faces:f.write(f'{a} {b}\n')
def solve(inp,cells):
 return json.loads(subprocess.run([str(ROOT/'domains'),str(inp),str(cells)],check=True,text=True,capture_output=True).stdout)
def summarize(cells,idx,mm,metadata):
 rows=[];vol=mm**3/1000
 def loc(i):
  return dict(index=int(i),label=int(cells['label'][i]),ijk=idx[i].tolist(),xyz_mm=((idx[i]+.5)*mm).tolist())
 for domain,labs in DOMAINS.items():
  selected=np.flatnonzero(np.isin(cells['label'],labs)); c=cells[selected];n=len(c)
  r=dict(metadata,domain=domain,n=n,volume_ml=n*vol)
  for key in ['all','window','glass']:
   j=int(np.argmin(c[key+'_rate']));i=selected[j]
   r[key+'_rate']=float(c[key+'_rate'][j]);r[key+'_location']=json.dumps(dict(loc(i),time_s=float(c[key+'_time'][j]),temperature_C=float(c[key+'_temp'][j])))
   bad=c[key+'_rate']<50
   r[key+'_rate_fail_ml']=int(bad.sum())*vol
   r[key+'_rate_fail_fraction']=float(bad.mean())
   ties=np.abs(c[key+'_rate']-r[key+'_rate'])<1e-8
   r[key+'_min_ties']=int(ties.sum());r[key+'_min_tie_ijk_min']=json.dumps(idx[selected[ties]].min(axis=0).tolist());r[key+'_min_tie_ijk_max']=json.dumps(idx[selected[ties]].max(axis=0).tolist())
   r[key+'_failure_temp_min']=float(c[key+'_temp'][bad].min()) if bad.any() else None
   r[key+'_failure_temp_max']=float(c[key+'_temp'][bad].max()) if bad.any() else None
  for key,fun in [('cold',np.argmin),('hot',np.argmax)]:
   j=int(fun(c['terminal']));r[key]=float(c['terminal'][j]);r[key+'_location']=json.dumps(dict(loc(selected[j]),time_s=70))
   r[key+'_ties']=int((np.abs(c['terminal']-r[key])<1e-8).sum())
  cold=c['terminal']< -35;hot=c['terminal']>0
  r['cold_fail_ml']=int(cold.sum())*vol;r['hot_fail_ml']=int(hot.sum())*vol
  for key,lo,hi in [('window','cross100','cross40'),('glass','cross123','cross35')]:
   valid=(c[lo]>=0)&(c[hi]>=0);missing=~valid
   r[key+'_missing_entry']=int((c[lo]<0).sum());r[key+'_missing_exit']=int((c[hi]<0).sum());r[key+'_complete']=int(valid.sum());r[key+'_missing_entry_ml']=int((c[lo]<0).sum())*vol;r[key+'_missing_exit_ml']=int((c[hi]<0).sum())*vol
   if valid.any():
    transit=c[hi][valid]-c[lo][valid];assert np.all(transit>0)
    delta=60 if key=='window' else 88
    r[key+'_transit_seconds_min']=float(transit.min());r[key+'_transit_seconds_max']=float(transit.max())
    r[key+'_transit_rate_min']=float(delta*60/transit.max());r[key+'_transit_rate_max']=float(delta*60/transit.min())
    j=np.flatnonzero(valid)[int(np.argmax(transit))]
    r[key+'_slow_transit_location']=json.dumps(dict(loc(selected[j]),entry_s=float(c[lo][j]),exit_s=float(c[hi][j])))
    j=np.flatnonzero(valid)[int(np.argmin(transit))]
    r[key+'_fast_transit_location']=json.dumps(dict(loc(selected[j]),entry_s=float(c[lo][j]),exit_s=float(c[hi][j])))
   else:
    for k in ['transit_seconds_min','transit_seconds_max','transit_rate_min','transit_rate_max','slow_transit_location']:r[key+'_'+k]=None
   fail=(c[key+'_rate']<50)|cold|hot|missing
   r[key+'_union_fail_ml']=int(fail.sum())*vol;r[key+'_union_fail_fraction']=float(fail.mean());r[key+'_pass']=not fail.any()
  fail=(c['all_rate']<50)|cold|hot
  r['all_union_fail_ml']=int(fail.sum())*vol;r['all_union_fail_fraction']=float(fail.mean());r['all_pass']=not fail.any()
  rows.append(r)
 return rows

def run(mm,design,mean,mode,dt):
 name=f'{mm}_{design}_{mean}_{mode}_{dt}';p,labels,iron,faces,idx=mesh(mm);iron[labels==2]=mean
 controls=P[3:] if design=='original' else [(P[3]+2*P[4]+P[5])/4]*3
 inp=OUT/(name+'.input');cp=OUT/(name+'.cells.csv')
 write_input(inp,labels,iron,faces,mm,dt,mode,controls)
 energy=solve(inp,cp);cells=np.genfromtxt(cp,delimiter=',',names=True)
 assert np.array_equal(cells['label'],labels)
 meta=dict(case=name,mm=mm,design=design,mean=mean,mode=mode,dt_requested=dt,**energy)
 (OUT/(name+'.json')).write_text(json.dumps(dict(meta,geometry=p.report(),input_sha256=digest(inp),cells_sha256=digest(cp)),indent=2)+'\n')
 rows=summarize(cells,idx,mm,meta)
 print(name,'specimen',[(r['all_rate'],r['cold'],r['all_pass']) for r in rows if r['domain']=='whole_specimen'],flush=True)
 return rows

def tests():
 # One isolated CPA voxel: independently integrate C(T)/q(T) using scalar quadrature.
 inp=OUT/'scalar.input';cp=OUT/'scalar.cells.csv';write_input(inp,[2],[1.83],[],1,.001,'clamped',[.9568032350959227]*3)
 e=solve(inp,cp);a=np.genfromtxt(cp,delimiter=',',names=True);T=float(a['terminal'])
 def capacity(t):
  cp=2402.3+7.6*t if t< -120 else 66747+545.24*t if t<=-119 else 2678.8+11.73*t+.039*t*t
  return (1078-.43*t)*cp
 u=.9568032350959227
 time=quad(lambda t:capacity(t)/(1000*1.83*(691-4.6*(max(-80,min(0,t))+80))*u),-150,T,points=[-120,-119,-80],epsabs=1e-8)[0]
 assert abs(time-70)<.0002,(time,T)
 assert abs(e['energy_rel_error'])<1e-6
 # Unit conversion and independent initial RHS.
 initial=60000*1.83*691*u/capacity(-150)
 assert abs(float(a['initial_rate'])-initial)<1e-9
 crossing_exact=quad(lambda t:capacity(t)/(1000*1.83*691*u),-150,-123,epsabs=1e-10)[0]
 crossing_error=float(a['cross123'])-crossing_exact
 assert abs(crossing_error)<1e-6,crossing_error
 # Two connected cells, zero source: preserve initial state, equal/opposite internal flux.
 inp=OUT/'zero.input';cp=OUT/'zero.cells.csv'
 write_input(inp,[2,2],[0,0],[(0,1)],1,.01,'clamped',[1]*3)
 e0=solve(inp,cp);zero=np.genfromtxt(cp,delimiter=',',names=True)
 assert np.all(zero['terminal']==-150) and e0['source_J']==0 and e0['enthalpy_J']==0
 # Two identical heated cells stay identical.
 inp=OUT/'symmetric.input';cp=OUT/'symmetric.cells.csv'
 write_input(inp,[2,2],[1.83,1.83],[(0,1)],1,.001,'clamped',[u]*3)
 es=solve(inp,cp);sym=np.genfromtxt(cp,delimiter=',',names=True)
 assert np.ptp(sym['terminal'])==0 and abs(sym['terminal'][0]-T)<1e-10
 assert abs(es['energy_rel_error'])<1e-6
 return dict(scalar_quadrature_time_s=time,scalar_time_error_s=time-70,scalar_crossing_error_s=crossing_error,initial_unit_rate=initial,zero_source=True,symmetry=True)

def csvwrite(path,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def main():
 assert digest('protocol.md')==Path('protocol.sha256').read_text().split()[0]
 subprocess.run(['c++','-O3','-std=c++17','domains.cpp','-o','domains'],check=True)
 checks=tests()
 # Execute baseline before discovery or scoring new comparisons.
 base=[]
 for d in ['original','flat']:base+=run(.5,d,1.83,'clamped',.01)
 expected=json.loads(Path('inputs/baseline_core.json').read_text())
 baseline_errors=[]
 for d in ['original','flat']:
  r=next(r for r in base if r['design']==d and r['domain']=='whole_specimen')
  b=next(b for b in expected if b['mode']=='clamped' and b['design']==('ramp' if d=='original' else 'flat'))
  for k,bk in [('all_rate','min_rate'),('window_rate','window_rate'),('cold','cold'),('hot','hot')]:
   error=r[k]-b[bk];assert abs(error)<.02,(d,k,error);baseline_errors.append(dict(design=d,metric=k,error=error))
 csvwrite(OUT/'baseline_check.csv',baseline_errors)
 checks['baseline_max_abs_error']=max(abs(r['error']) for r in baseline_errors)
 rows=[]
 for mean in [1.47,1.83]:
  for mode in ['extrapolated','clamped']:
   for d in ['original','flat']:rows+=run(.75,d,mean,mode,.01)
 rows+=base
 for mean,mode in [(1.83,'clamped'),(1.47,'extrapolated')]:
  for d in ['original','flat']:
   if mean!=1.83:rows+=run(.5,d,mean,mode,.01)
   rows+=run(.5,d,mean,mode,.005)
 csvwrite(OUT/'domains.csv',rows)
 assert all(abs(r['energy_rel_error'])<1e-6 for r in rows)
 assert all(r['minimum_observed_rate']>=0 for r in rows)
 assert all(r['window_complete']==r['n'] for r in rows if r['domain']!='wall')
 assert all(r[key+'_complete']+r[key+'_missing_exit']==r['n'] and r[key+'_missing_entry']<=r[key+'_missing_exit'] for r in rows for key in ['window','glass'])
 comparisons=[]
 for r in rows:
  if r['design']!='original':continue
  f=next(f for f in rows if all(f[k]==r[k] for k in ['mm','mean','mode','dt_requested','domain']) and f['design']=='flat')
  comparisons.append({k:r[k] for k in ['mm','mean','mode','dt_requested','domain']}|{key+'_flat_minus_original':f[key]-r[key] for key in ['all_rate','window_rate','glass_rate','cold','hot']}|{key+'_original_pass':r[key+'_pass'] for key in ['all','window','glass']}|{key+'_flat_pass':f[key+'_pass'] for key in ['all','window','glass']})
 csvwrite(OUT/'paired.csv',comparisons)
 refine=[]
 for r in rows:
  if r['mm']!=.5 or r['dt_requested']!=.005:continue
  b=next(b for b in rows if all(b[k]==r[k] for k in ['mm','mean','mode','domain','design']) and b['dt_requested']==.01)
  refine.append({k:r[k] for k in ['mean','mode','domain','design']}|{k+'_delta':r[k]-b[k] for k in ['all_rate','window_rate','glass_rate','cold','hot','window_transit_rate_min','glass_transit_rate_min']})
 csvwrite(OUT/'time_refinement.csv',refine)
 checks.update(max_energy_rel_error=max(abs(r['energy_rel_error']) for r in rows),simulations=len(rows)//6,domain_rows=len(rows),max_time_rate_delta=max(abs(r[k+'_delta']) for r in refine for k in ['all_rate','window_rate','glass_rate']),max_time_endpoint_delta=max(abs(r[k+'_delta']) for r in refine for k in ['cold','hot']),nonwall_window_complete=True,crossing_accounting=True,monotone=True,versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0]))
 (OUT/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
 print('C05 VERIFIED',json.dumps(checks),flush=True)
if __name__=='__main__':main()
