"""Independently regenerate the inputs and check every matched contrast."""
import json
import math
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

from analysis import SETTINGS, compare, target_grid
from evaluate import METRICS, OVERLAP, PAIR, REFERENCE, REFINED, record
from heat_cache import DT, compute
from cut_geometry import CENTRES, POWERS, RADII, bounds
from overlap import moments

ROOT = Path(__file__).resolve().parent


def close_tree(a, b, path=''):
    if isinstance(b, dict):
        assert isinstance(a, dict) and a.keys() == b.keys(), path
        for k in b:
            close_tree(a[k], b[k], f'{path}/{k}')
    elif isinstance(b, (tuple, list)):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a,b)):
            close_tree(x, y, f'{path}/{i}')
    elif isinstance(b, float):
        if path.endswith(('/residual','/heat_balance')):
            assert 0 <= a < 2e-7 and 0 <= b < 2e-7, path
        else:
            assert math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-8), (path,a,b)
    else:
        assert a == b, (path,a,b)


def audit_overlap_boxes():
    # This rectangular patch lies inside the inner bag and entirely outside
    # the heart. The minimum and maximum powers give exact tests over a box.
    lo, hi = np.array([2.6,8.3,5.9]), np.array([3.4,9.2,6.7])
    for shape in (0, 1, 2):
        c, r, p = CENTRES[shape], RADII[shape], POWERS[shape]
        distance = (np.maximum(abs(lo-c),abs(hi-c)) if shape < 2
                    else np.maximum(np.maximum(lo-c,c-hi),0))
        inside = np.sum((distance/r)**p) < 1
        assert inside == (shape < 2)

    donor = bounds(.75)
    axes = [np.array([a,b]) for a,b in zip(lo,hi)]
    labels = np.ones((1,1,1),int)
    uid = np.arange(4*math.prod(len(a)-1 for a in donor)).reshape(
        *(len(a)-1 for a in donor),4)
    intersecting = [np.flatnonzero((a[:-1] < hi[d]) & (a[1:] > lo[d]))
                    for d,a in enumerate(donor)]
    expected = {}
    for cell in product(*intersecting):
        a = np.maximum(lo,[donor[d][i] for d,i in enumerate(cell)])
        b = np.minimum(hi,[donor[d][i+1] for d,i in enumerate(cell)])
        volume = np.prod(b-a)
        expected[uid[(*cell,1)]] = volume*np.r_[1.,(a+b)/2]

    for axis in range(3):
        row, column, integrals = moments(donor,axes,uid,labels,3,axis)
        assert (row == 0).all() and set(column) == set(expected)
        for code, exact in expected.items():
            np.testing.assert_allclose(integrals[column == code].sum(axis=0), exact,
                                       rtol=1e-12,atol=1e-12)


def audit_history(graph, history, event):
    phase = graph[0]
    heart = history[:,phase >= 2].min(axis=1)
    assert np.all(heart[:-1] < -35) and abs(heart[-1]+35) < 1e-10
    assert np.all(history[0] == -150.)
    specimen = history[:,phase > 0]
    elapsed = event['time_s']-(len(specimen)-2)*DT
    assert 0 < elapsed <= DT+1e-10
    crossed = np.zeros(specimen.shape[1],bool)
    floor = np.inf
    for i in range(1,len(specimen)):
        valid = crossed | (specimen[i] > -123.)
        if valid.any():
            delta = (specimen[i,valid]-specimen[i-1,valid]).min()
            floor = min(floor, delta*60/(elapsed if i == len(specimen)-1 else DT))
        crossed |= specimen[i] >= -123.
    assert abs(floor-event['min_c_min']) < 1e-8
    assert abs(history[-1].max()-event['hot_any_c']) < 1e-10
    assert event['passed'] and floor >= 50 and history[-1].max() < 0


def audit_result(runs, archived):
    previous = json.loads((ROOT/'sources/predecessor-cut-fine.json').read_text())['extrapolated']
    for schedule, entry in previous.items():
        for metric in METRICS:
            current = runs['native']['schedules'][schedule]['idw8'][metric]
            assert abs(current-entry['stress']['idw8'][metric]) < 1e-5
    for name, run in runs.items():
        phantom, _ = target_grid(name,run['spacing_mm'])
        assert sum(run['label_counts']) == phantom.active_mask.sum()
        assert run['overlap']['total_missing_fraction'] < 1e-5
        assert run['patches'][OVERLAP]['affine'] < 1e-9
        for mapping in run['schedules'][REFERENCE]:
            for metric in METRICS:
                a, b = (run['schedules'][s][mapping][metric] for s in PAIR)
                ref = run['schedules'][REFERENCE][mapping][metric]
                pair = archived['contrasts'][name][mapping][metric]
                assert abs((a-ref)-(b-ref)-pair['graded_minus_strict']) < 2e-8
                assert abs(a-ref-pair['graded_minus_reference']) < 2e-8
                assert abs(b-ref-pair['strict_minus_reference']) < 2e-8

    def delta(grid, mapping, metric):
        values = runs[grid]['schedules']
        return values[PAIR[0]][mapping][metric]-values[PAIR[1]][mapping][metric]

    for metric in METRICS:
        sample = [delta(grid,OVERLAP,metric) for grid in REFINED]
        fine = sample[-1]
        checks = archived['sensitivity'][metric]
        assert abs(abs(fine)-checks['fine_abs_contrast']) < 2e-8
        assert abs(np.ptp(sample)-checks['sensitivity']['affine_three_grid_range']) < 2e-8
        transfer = [delta('fine', m, metric) for m in runs['fine']['schedules'][REFERENCE]]
        assert abs(np.ptp(transfer)-checks['sensitivity']['fine_transfer_range']) < 2e-8


def audit_decision(runs):
    def score(grid, schedule, mapping, metric):
        return runs[grid]['schedules'][schedule][mapping][metric]

    gaps = [score(grid, PAIR[0], OVERLAP, 'objective_mpa') -
            score(grid, PAIR[1], OVERLAP, 'objective_mpa') for grid in REFINED]
    assert gaps[0] < 0 < gaps[1] and gaps[2] < 0, gaps
    assert np.ptp(gaps) >= abs(gaps[-1]), gaps
    for grid, run in runs.items():
        for mapping in run['schedules'][REFERENCE]:
            for metric in METRICS:
                reference = score(grid, REFERENCE, mapping, metric)
                assert all(reference < score(grid, schedule, mapping, metric)
                           for schedule in PAIR), (grid, mapping, metric)


def main():
    audit_overlap_boxes()
    graph, histories, events = compute()
    for schedule, history in histories.items():
        audit_history(graph, history, events[schedule])
    with ThreadPoolExecutor(max_workers=len(SETTINGS)) as pool:
        jobs = {name: pool.submit(compare, name, graph, histories, events) for name in SETTINGS}
        runs = {name: task.result() for name, task in jobs.items()}
    archived = json.loads((ROOT/'results.json').read_text())
    close_tree(record(runs, events), archived)
    audit_result(runs, archived)
    audit_decision(runs)
    print('Verified fresh trajectories, overlap patches, coupled stresses and signed contrasts.')


if __name__ == '__main__':
    main()
