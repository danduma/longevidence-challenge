"""Verify numerical decisions, independent checks, and the provenance of source bytes."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import pandas as pd

def table(name):
    return pd.read_csv(f'outputs/{name}.csv')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

pin=subprocess.check_output(['git','-C','upstream','rev-parse','HEAD'],text=True).strip()
assert pin=='1d498ee5c583e0112a762e60892d576ec2e27f57'
assert sha('preregistration.md')==Path('outputs/preregistration.sha256').read_text().split()[0]
assert sha('sources/liter.xlsx')==sha('sources/publisher_check.xlsx')=='c059e1beef528efd5a3aa04fc51d039d64e1c8583376a3596160e837069d2898'
base=table('baselines')
assert len(base)==70 and base.error.abs().max()<.02
s=table('sweep')
assert len(s)==180 and s.enthalpy_relative_error.abs().max()<.005
for _,g in s.query("mm==.5 and dt_divisor==1").groupby(['sar','scenario']):
    g=g.set_index('design')
    assert g.loc['margin','window_rate']>g.loc['equal','window_rate']
r=s.query("mm==.5 and dt_divisor==1 and sar=='extrapolated' and design=='margin'").set_index('scenario')
assert r.passes.to_dict()=={'low':False,'high':True,'nominal':True,'deepLow':False,'deepHigh':False}
assert r.window_passes.all() and r.incomplete.max()==0
fine=pd.concat([s,table('refinement')]).query("design=='equal' and scenario=='high' and dt_divisor==1 and sar=='clamped'")
assert len(fine)==8 and (fine.passes==(fine.mm<1.4)).all()
bound=table('loading_boundary');assert len(bound)==12 and (bound.passes==(bound.fe>=1.82)).all()
tol=table('local_tolerance');assert len(tol)==18 and (tol.passes==(tol.amplitude<=.01)).all()
split=table('split_fine');assert len(split)==4 and (split.window_passes==(split.scenario=='deepHigh')).all()
check=table('independent')
assert len(check)==12 and check.difference.abs().max()<.032
raw=table('raw_measurement_check')
assert len(raw)==6 and raw.raw_min_C_min.max()==0 and raw.transit_to_zero.min()>85
assert table('analytic').query('dt==.005').error.abs().max()<.006
b=json.loads(Path('outputs/assay_bound.json').read_text())
assert b['tests']['extrapolated_margin']['rounding_conservative_slowest_start_rate']<45.21
for doc in ['report.md','sources.md']:
    for link in re.findall(r'\]\(([^)]+)\)',Path(doc).read_text()):
        if not link.startswith('https:'):
            assert Path(link).exists(),link
print('Pinned upstream, preregistration hash and original publisher bytes verified.')
print('70 baseline metrics verified; max abs discrepancy',base.error.abs().max())
print('180 thermal runs; max enthalpy discrepancy',s.enthalpy_relative_error.abs().max())
print('Verified eight-grid clamped comparison, 12 loading-bracket runs, 18 equal-mean tolerance runs and four sub-0.35-mm split runs.')
print('Independent RK2, analytic quadrature, assay bound, six raw traces and document links verified.')
