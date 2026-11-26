#!/usr/bin/env python3
"""C02 offline verifier. Default repeats every integration in a private job subdirectory."""
import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))
def audit(path):
    rows=read(path/'envelope.csv')
    lookup={(float(r['mm']),r['design'],float(r['h']),int(r['ambient']),float(r['dt'])):r for r in rows}
    assert len(rows)==len(lookup)
    for mm,ambients in [(.75,[-150,-80,-20,22]),(.5,[-150,22])]:
        dt=.02 if mm==.75 else .01
        for amb in ambients:
            for design in ['flat','ramp']:
                for h in [0,2,5,10,20,40,80]:
                    assert (mm,design,h,amb,dt) in lookup
    for r in rows:
        assert int(r['bad_cross'])==0 and int(r['exposed_nonwall'])==0
        assert abs(float(r['energy_rel_error']))<1e-5
        assert all(math.isfinite(float(r[k])) for k in ['min_rate','cold','hot','boundary_J','source_J','enthalpy_J'])
        flags=[float(r['min_rate'])>=50,float(r['cold'])>=-35,float(r['hot'])<=0]
        assert (r['pass']=='True')==all(flags)
        assert [r[k]=='True' for k in ['rate_ok','cold_ok','hot_ok']]==flags
        assert float(r['myo_cold'])>=float(r['cold']) and float(r['myo_hot'])<=float(r['hot'])
        assert float(r['myo_rate'])>=float(r['min_rate'])
    for b in read(path/'brackets.csv'):
        mm=float(b['mm'])
        lo,hi=float(b['lo']),float(b['hi'])
        assert 0<hi-lo<=2 and b['lo_pass']!=b['hi_pass']
        for h,val in [(lo,b['lo_pass']),(hi,b['hi_pass'])]:
            r=lookup[mm,b['design'],h,int(b['ambient']),.01 if mm==.5 else .02]
            assert r[b['criterion']]==val
    return len(rows)
def compare(first,second):
    a,b=read(first),read(second)
    assert len(a)==len(b)
    for left,right in zip(a,b):
        assert left.keys()==right.keys()
        for k,v in left.items():
            if v!=right[k]:
                assert math.isclose(float(v),float(right[k]),abs_tol=1e-7,rel_tol=1e-8),(k,v,right[k])
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--tables-only',action='store_true')
    options=parser.parse_args()
    os.chdir(ROOT)
    subprocess.run(['shasum','-a','256','-c','manifest.sha256'],check=True)
    count=audit(ROOT/'results')
    if not options.tables_only:
        with tempfile.TemporaryDirectory(dir=ROOT,prefix='verify-') as tmp:
            work=Path(tmp)
            for f in ['run.py','boundary.cpp','protocol.md','protocol.sha256']:
                shutil.copy2(ROOT/f,work/f)
            shutil.copytree(ROOT/'pinned',work/'pinned')
            (work/'results').mkdir()
            subprocess.run([sys.executable,'-B',str(work/'run.py')],check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
            audit(work/'results')
            for name in ['envelope.csv','brackets.csv']:
                compare(ROOT/'results'/name,work/'results'/name)
    print('C02 verified',count,'rows;', 'table audit' if options.tables_only else 'full offline integration and table comparison')
