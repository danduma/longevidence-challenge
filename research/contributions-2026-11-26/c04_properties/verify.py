#!/usr/bin/env python3
"""C04. Run from any directory. No network calls."""
import sys
sys.dont_write_bytecode = True
import csv
import hashlib
import json
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from itertools import product
from scipy.integrate import quad
from scipy.optimize import brentq
from geometry_driver import ROOT, write_input, flat

OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)
designs = ['ramp', 'flat']
scales = [.95, 1., 1.05]
def table(name, records):
    with (OUT/name).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(records)

def main():
    for line in (ROOT/'inputs.sha256').read_text().splitlines():
        digest, name = line.split(None, 1)
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest
    with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
        work = Path(tmp)
        exe = work/'solver'
        subprocess.run(['c++','-O3','-std=c++17',ROOT/'solver.cpp','-o',exe], check=True)
        def solve(f, s, c):
            return json.loads(subprocess.check_output([exe, f, str(s), str(c)]))
        # Constant clamped SAR below -80, exact energy and nonlinear T from quadrature.
        scalar = work/'scalar'
        scalar.write_text('1 0 0.001 0.005 15 1 2 0 -150 1 1 1\n2 1.83\n')
        q = .95*1000*1.83*691
        def capacity(t):
            return 1.05*(1078-.43*t)*(2402.3+7.6*t)
        sc = solve(scalar,.95,1.05)
        exact = brentq(lambda t: quad(capacity,-150,t)[0]-q*15,-150,-120)
        assert abs(sc['cold']-exact)<1e-7 and abs(sc['source_J']-q*15*1e-9)<1e-12
        # Independent scalar crossing check through both capacity joins and SAR knee.
        def seconds_per_degree(t):
            cp = (2402.3+7.6*t if t < -120 else
                  66747+545.24*t if t <= -119 else 2678.8+11.73*t+.039*t*t)
            sar = 691-4.6*(max(-80.,t)+80.)
            return 1.05*(1078-.43*t)*cp/(.95*1000*1.83*sar)
        joins = [-120,-119,-80]
        duration = quad(seconds_per_degree,-150,-35,points=joins)[0]
        exact_transit = 83*60/quad(seconds_per_degree,-123,-40,points=joins)[0]
        scalar.write_text(f'1 0 0.001 0.005 {duration} 1 2 0 -150 1 1 1\n2 1.83\n')
        crossing = solve(scalar,.95,1.05)
        assert crossing['missing']==0 and abs(crossing['transit']-exact_transit)<.002
        inputs = {}
        for mm,dt in [(.5,.01),(.5,.02),(.75,.02)]:
            for design in designs:
                f=work/f'{mm}_{dt}_{design}'
                write_input(f,mm,design,dt)
                inputs[mm,dt,design]=f
        def case(key):
            mm,dt,d,s,c=key
            r=solve(inputs[mm,dt,d],s,c)
            assert abs(r['energy_rel_error'])<1e-5 and r['boundary_J']==0
            assert abs(r['duty_s']-flat*70)<1e-10
            assert r['post']<1e90 and r['mus_post']<1e90 and r['transit']>0
            r['mus_pass']=r['mus_all']>=50 and r['mus_cold']>=-35 and r['mus_hot']<=0
            return dict(mm=mm,dt_s=dt,design=d,s=s,c=c,**r)
        baseline=[]
        for d in designs:
            r=case((.5,.01,d,1.,1.))
            with (ROOT/'inputs/baselines.csv').open() as stream:
                ref=next(row for row in csv.DictReader(stream) if row['design']==d)
            for k in ['min_rate','cold','hot','window_rate','source_J']:
                assert abs(r[k]-float(ref[k]))<1e-7, (d,k,r[k],ref[k])
            assert r['pass']==(ref['pass']=='True')
            baseline.append(r)
        print('C04 central baselines recovered', flush=True)
        keys=list(product([.75],[.02],designs,scales,scales))
        keys += [(.75,.02,'flat',round(c*ratio,7),c) for c in scales for ratio in [.992,.996]]
        for dt in [.02,.01]:
            keys.extend((.5,dt,d,s,c) for d in designs for s,c in [(.95,1.05),(1.05,.95),(1.,1.)] if dt!=.01 or s!=1.)
            keys.extend((.5,dt,'flat',s,1.) for s in [.992,.996])
        with ThreadPoolExecutor(max_workers=3) as pool:
            rows=baseline+list(pool.map(case, keys))
    table('cases.csv',rows)
    pairs=[]
    for ramp in rows:
        if ramp['design']!='ramp': continue
        key={k:ramp[k] for k in ['mm','dt_s','s','c']}
        constant=next(r for r in rows if r['design']=='flat' and all(r[k]==v for k,v in key.items()))
        pairs.append(dict(**key,delta_post=ramp['post']-constant['post'],delta_mus_post=ramp['mus_post']-constant['mus_post'],delta_all=ramp['min_rate']-constant['min_rate'],ramp_pass=ramp['pass'],flat_pass=constant['pass']))
    table('pairs.csv',pairs)
    bounds=[]
    for mm,dt,c in sorted({(r['mm'],r['dt_s'],r['c']) for r in rows}):
        subset=sorted((r for r in rows if (r['mm'],r['dt_s'],r['c'],r['design'])==(mm,dt,c,'flat')),key=lambda r:r['s'])
        for lo,hi in zip(subset,subset[1:]):
            if lo['pass']!=hi['pass']:
                bounds.append(dict(mm=mm,dt_s=dt,c=c,s_lo=lo['s'],lo_pass=lo['pass'],s_hi=hi['s'],hi_pass=hi['pass'],ratio_lo=lo['s']/c,ratio_hi=hi['s']/c))
    table('flat_brackets.csv',bounds)
    (OUT/'checks.json').write_text(json.dumps(dict(scalar_transit_exact=exact_transit,scalar_transit_RK2=crossing['transit'],scalar_exact_C=exact,scalar_solver_C=sc['cold'],scalar_source_J=sc['source_J'],max_abs_relative_energy_error=max(abs(r['energy_rel_error']) for r in rows),max_missing_crossings=max(r['missing'] for r in rows),onset_required_ratio=50*(1078+.43*150)*(2402.3-7.6*150)/(60000*flat*1.83*691)),indent=2)+'\n')
    # Regression checks on the observed C04 conclusions, after baseline and physics checks.
    discovery=[p for p in pairs if p['mm']==.75]
    assert len(rows)==40 and len(discovery)==9
    assert sum(p['flat_pass'] for p in discovery)==6
    assert all(not p['ramp_pass'] and p['delta_post']>0 and p['delta_mus_post']>0 for p in pairs)
    assert abs(min(p['delta_post'] for p in discovery)-2.554295)<1e-5
    assert abs(max(p['delta_post'] for p in discovery)-2.647500)<1e-5
    assert len(bounds)==5 and all(abs(b['ratio_lo']-.992)<1e-10 and
        abs(b['ratio_hi']-.996)<1e-10 and not b['lo_pass'] and b['hi_pass'] for b in bounds)
    print(f'C04 completed {len(rows)} deterministic calculations and {len(pairs)} paired contrasts.')
if __name__=='__main__': main()
