"""Rebuild geometry, integrate both fields and assert the archived coupled result."""
import json
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from math import isclose
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from analysis import run, ROOT, FLOOR
from baseline import baseline
from mechanics import UPSTREAM, read

METHODS = ('voxel', 'cut')
MODES = ('extrapolated', 'clamped')
MAPS = ('nearest', 'idw8')
SCORES = ('objective_mpa', 'peak_mpa')
RUNS = {'main': (.1,.6,16), 'refined': (.05,.3,24)}


def same(actual, expected, path=''):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), path
        for name, value in expected.items():
            same(actual[name], value, f'{path}/{name}')
    elif isinstance(expected, list):
        assert len(actual) == len(expected), path
        for index, (left, right) in enumerate(zip(actual, expected)):
            same(left, right, f'{path}/{index}')
    elif type(expected) is float:
        # Residuals are diagnostics and need only satisfy their own solver bound.
        if path.endswith(('/residual', '/heat_balance')):
            assert actual < 2e-7, (path, actual)
        else:
            assert isclose(actual,expected,rel_tol=2e-8,abs_tol=2e-8), (path,actual,expected)
    else:
        assert actual == expected, (path, actual, expected)


def conclusions(records):
    output = {}
    for key, resolution in records.items():
        output[key] = {}
        for mode in MODES:
            detail = {'thermal': {}, 'ranks': {}, 'pairs': {}}
            output[key][mode] = detail
            for method in METHODS:
                records_by_schedule = resolution[method][mode]
                feasible = []
                for name, row in records_by_schedule.items():
                    event = row['event']
                    assert abs(event['cold_heart_c'] + 35) < 1e-10
                    assert event['min_location']['time_s'] <= event['time_s']
                    passed = event['min_c_min'] >= FLOOR and event['hot_any_c'] < 0
                    assert passed == event['passed']
                    if passed:
                        feasible.append(name)
                detail['thermal'][method] = sorted(feasible)
                detail['ranks'][method] = {
                    mapping: {score: sorted(feasible, key=lambda n: records_by_schedule[n]['stress'][mapping][score])
                              for score in SCORES}
                    for mapping in MAPS}
            candidates = sorted(set.intersection(*(set(s) for s in detail['thermal'].values())))
            for a,b in combinations(candidates,2):
                detail['pairs'][f'{a} minus {b}'] = {
                    mapping: {score: {
                        method: resolution[method][mode][a]['stress'][mapping][score] -
                                resolution[method][mode][b]['stress'][mapping][score]
                        for method in METHODS}
                        for score in SCORES}
                    for mapping in MAPS}
    return output


def verify():
    reproduced = baseline()
    reference = read('objective_horizon.json')['designs']
    for name, row in reproduced.items():
        assert abs(row['peak_mpa']-reference[name]['12']['peak_tension_mpa']) < 1e-8
        assert abs(row['peak_time_s']-reference[name]['12']['peak_time_s']) < 1e-8
    assert abs(reproduced['optimized']['objective_mpa'] -
               read('optimization_extrapolated.json')['optimizer']['objective_mpa']) < 1e-8
    with TemporaryDirectory(dir=ROOT) as tmp:
        with ProcessPoolExecutor(max_workers=2) as pool:
            jobs = {resolution: {
                method: pool.submit(run,*settings,method,Path(tmp)/f'{resolution}-{method}.json',audit_thermal)
                for method in METHODS}
                for resolution,settings in RUNS.items()}
            records = {resolution: {method: future.result() for method,future in group.items()}
                       for resolution,group in jobs.items()}
    expected = json.loads((ROOT/'results.json').read_text())
    same(reproduced,expected['baseline'])
    same(records,expected['runs'])
    comparison = conclusions(records)
    same(comparison,expected['comparison'])
    assert_result(comparison)
    print('Verified released objective, fresh coupled trajectories and schedule comparison.')


def audit_thermal(graph, history, event, dt):
    """Calculate the screen directly from temperatures, without solver rates."""
    phase = graph[0]
    heart, specimen = phase >= 2, phase > 0
    assert (history[:-1, heart].min(axis=1) < -35).all()
    assert abs(history[-1, heart].min()+35.) < 1e-10
    last_step = event['time_s'] - (len(history)-2)*dt
    assert 0 < last_step <= dt + 1e-10
    # First output sample at or above Tg identifies the eligible crossing
    # interval independently for each cell.
    sampled = history[:, specimen]
    above = sampled >= -123.
    first = above.argmax(axis=0)
    assert above.any(axis=0).all() and (first > 0).all()
    rate_min = np.inf
    for frame in range(1,len(sampled)):
        eligible = (first < frame) | ((first == frame) & (sampled[frame] > -123.))
        if eligible.any():
            delta = sampled[frame,eligible] - sampled[frame-1,eligible]
            rate_min = min(rate_min, delta.min() * 60 / (last_step if frame==len(sampled)-1 else dt))
    assert isclose(rate_min,event['min_c_min'],rel_tol=1e-8,abs_tol=1e-8)
    assert isclose(history[-1,specimen].max(),event['hot_c'],abs_tol=1e-10)
    assert isclose(history[-1].max(),event['hot_any_c'],abs_tol=1e-10)


def assert_result(comparison):
    winner = 'uniform_released_policy'
    all_schedules = {'graded', winner, 'uniform_strict'}
    for settings in comparison.values():
        for details in settings.values():
            for method in METHODS:
                assert set(details['thermal'][method]) == all_schedules
                for mapping in MAPS:
                    for score in SCORES:
                        assert details['ranks'][method][mapping][score][0] == winner
        extrapolated = settings['extrapolated']['pairs']['graded minus uniform_strict']
        for mapping in MAPS:
            objectives = extrapolated[mapping]['objective_mpa']
            assert objectives['voxel'] > 0 > objectives['cut']


if __name__ == '__main__':
    verify()
