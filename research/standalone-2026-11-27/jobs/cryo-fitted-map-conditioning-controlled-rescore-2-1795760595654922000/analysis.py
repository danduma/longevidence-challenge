"""Paired mechanics with common physical points and a conditioning mask."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json

import numpy as np
from scipy.special import logsumexp

from fitted import SHAPES, Mesh, Solver, overlap, physical
from heat_cache import DT, H, Q, compute
from inherited import SCALE, Sample, cloud
from regular import RegularMesh
from thermal_graph import cut_graph

ROOT = Path(__file__).resolve().parent
CONFIGS = [(family, level, q) for q in (12,20) for level in (0,1)
           for family in ('cube','regular')]


def support():
    ref, xyz, weight = cloud(3)
    differences = np.eye(3)*1e-5
    jac = np.stack([(physical(ref+h)-physical(ref-h))/(2e-5) for h in differences],axis=-1)
    condition = np.linalg.cond(jac)
    return ref, xyz, weight, condition <= 10


def prepare():
    cache = ROOT/'thermal-cache.npz'
    if cache.exists():
        graph = cut_graph(H,Q)
        with np.load(cache,allow_pickle=False) as data:
            return graph, dict(data), json.loads((ROOT/'events.json').read_text())
    graph, histories, events = compute()
    stride = round(.3/DT)
    histories = {name: np.concatenate((history[stride:-1:stride], history[-1:]))
                 for name, history in histories.items()}
    np.savez_compressed(cache,**histories)
    (ROOT/'events.json').write_text(json.dumps(events,indent=2)+'\n')
    return graph,histories,events


def measure(p, w):
    return np.array([p[:, -1].max(), SCALE*logsumexp(p/SCALE+np.log(w[:,None]))])/1000


def run(family, level, q, inputs, scorer, sampling, cloud_data):
    graph,histories,events = inputs
    mesh = Mesh(SHAPES[level]) if family == 'cube' else RegularMesh(4*2**level)
    solver = Solver(mesh)
    matrix, diagnostic = overlap(mesh,graph,q)
    _, xyz, weights, keep = cloud_data
    sample = sampling(mesh, cloud_data)
    masks = {'all':np.ones(len(xyz),bool),'interior':keep}
    measures = {key:weights[mask]/weights[mask].sum() for key,mask in masks.items()}
    scores = {}
    for name, history in histories.items():
        cadence = np.arange(1,len(history))*round(.3/DT)*DT
        times = np.r_[cadence, events[name]['time_s']]
        assert 0 < times[-1]-times[-2] <= .3+1e-10
        field = (matrix@history.T).T
        traces = {key:[] for key in masks}
        def observe(time, u, memory, shear, eigenstrain):
            p = sample.evaluate(u,memory,shear,eigenstrain,solver.bulk)
            for key, mask in masks.items():
                traces[key].append(scorer(p[mask],measures[key]))
        solver.solve(times,field,observe)
        scores[name] = dict(times=times.tolist(),
                            **{key:np.asarray(trace).tolist() for key,trace in traces.items()})
        print(family,level,q,name,flush=True)
    conditioning = None
    if family == 'regular':
        points = np.stack(np.meshgrid(*([np.linspace(0,1,5)]*3),indexing='ij'),axis=-1).reshape(-1,3)
        _, j, _ = mesh.quadrature(points)
        s = np.linalg.svd(j*mesh.width[:,None,None,:],compute_uv=False)
        conditioning = dict(max_condition=float((s[...,0]/s[...,2]).max()),
                            min_singular_mm=float(s[...,2].min()))
    return dict(family=family,level=level,q=q,cells=len(mesh.cell),
                volume_mm3=float(solver.volume.sum()),conditioning=conditioning,
                transfer=diagnostic,scores=scores)


class RegularSample(Sample):
    def __init__(self, mesh, xyz):
        self.owner, self.gradient = mesh.locate(xyz)


def make_sample(mesh, cloud_data):
    if isinstance(mesh, RegularMesh):
        return RegularSample(mesh, cloud_data[1])
    return Sample(mesh, 3)


def regenerate(scorer=measure, sampling=make_sample, cloud_source=support):
    inputs = prepare()
    cloud_data = cloud_source()
    _, _, w, keep = cloud_data
    with ThreadPoolExecutor(max_workers=4) as pool:
        tasks = [pool.submit(run,*config,inputs,scorer,sampling,cloud_data) for config in CONFIGS]
        rows = [task.result() for task in tasks]
    return dict(points=len(keep),retained_points=int(keep.sum()),
                retained_weight=float(w[keep].sum()/w.sum()),runs=rows)


if __name__ == '__main__':
    (ROOT/'results.json').write_text(json.dumps(regenerate(),allow_nan=False)+'\n')
