"""Cross transfer and mechanical spacing on the same extrapolated SAR histories."""
import json
from pathlib import Path

import numpy as np

from cut_geometry import CENTRE, NATIVE, classify
from heat_cache import DT, H, Q, compute
from mechanics import Mechanics, rat_heart_cryobag, summarize
from overlap import patch_tests, transfers

ROOT = Path(__file__).resolve().parent
SETTINGS = {'native': (1.4, 6, 0), 'coarse': (1.4, 6, 0),
            'medium': (1.1, 6, 0), 'fine': (.9, 6, 0),
            'fine_quadrature': (.9, 9, 2)}
CADENCE = .3


def target_grid(name, spacing):
    if name == 'native':
        phantom = rat_heart_cryobag(spacing)
        return phantom, [np.arange(n+1)*spacing for n in phantom.shape]
    axes = [np.arange(int(np.ceil(d/spacing))+1)*spacing for d in 2*CENTRE]
    centers = [(a[1:]+a[:-1])/2 for a in axes]
    labels = classify(np.stack(np.meshgrid(*centers, indexing='ij'), axis=-1))
    return type(NATIVE)(labels, spacing/1000, NATIVE.spec, NATIVE.pillow_exponent), axes


def compare(name, graph, histories, events):
    spacing, q, axis = SETTINGS[name]
    phantom, axes = target_grid(name, spacing)
    maps, centroid, diagnostics = transfers(graph,H,axes,phantom.labels,q,axis)
    tests = patch_tests(maps,graph,centroid,phantom.labels)
    row = dict(spacing_mm=spacing, shape=phantom.shape,
               label_counts=[int((phantom.labels==p).sum()) for p in range(4)],
               patches=tests, overlap=diagnostics, schedules={})
    if name == 'fine_quadrature':
        maps = {'overlap_affine': maps['overlap_affine']}
    solver = Mechanics(phantom)
    for schedule, history in histories.items():
        event = events[schedule]
        stride = round(CADENCE/DT)
        regular = np.arange(stride,len(history)-1,stride)
        temperatures = history[np.r_[regular,len(history)-1]]
        times = np.r_[regular*DT,event['time_s']]
        row['schedules'][schedule] = {}
        for mapping, matrix in maps.items():
            scores = summarize(times, solver.solve(times, (matrix@temperatures.T).T))
            row['schedules'][schedule][mapping] = scores
            print(name,schedule,mapping,scores,flush=True)
    return row


def main():
    from concurrent.futures import ThreadPoolExecutor
    from evaluate import record

    graph, histories, events = compute()
    with ThreadPoolExecutor(max_workers=len(SETTINGS)) as pool:
        jobs = {name: pool.submit(compare, name, graph, histories, events) for name in SETTINGS}
        runs = {name: job.result() for name, job in jobs.items()}
    output = record(runs, events)
    (ROOT/'results.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
