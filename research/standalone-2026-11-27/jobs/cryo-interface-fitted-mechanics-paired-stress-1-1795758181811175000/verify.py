"""Recompute thermal paths, both fitted mechanical levels and overlap sensitivity."""
import json
from math import gamma, prod
from pathlib import Path

import numpy as np

from analysis import CONFIGURATIONS, SHAPES, physical, regenerate
from cut_geometry import CENTRE, POWERS, RADII
from heat_cache import DT

ROOT = Path(__file__).resolve().parent
METRICS = ('objective_mpa', 'peak_mpa')
PAIR = ('graded', 'uniform_strict')
REF = 'uniform_released_policy'


def comparisons(record):
    """Evaluate signed contrasts directly from the absolute stress scores."""
    rows = record['runs']
    gap = np.array([[1000*(row['schedules'][PAIR[0]][key] -
                            row['schedules'][PAIR[1]][key])
                     for key in METRICS] for row in rows])
    lead = np.array([[[1000*(row['schedules'][name][key] -
                              row['schedules'][REF][key])
                       for key in METRICS] for name in PAIR] for row in rows])
    return gap, lead


def geometry_check():
    coarse, fine = SHAPES
    for axis, (a,b) in enumerate(zip(coarse,fine)):
        coarse_axis = np.linspace(-1,1,a+1)
        fine_axis = np.linspace(-1,1,b+1)
        np.testing.assert_array_equal(coarse_axis, fine_axis[::2])
        uv = np.linspace(-1,1,37)
        u = np.zeros((len(uv)**2,3))
        u[:,[d for d in range(3) if d != axis]] = np.array(
            np.meshgrid(uv,uv,indexing='ij')).reshape(2,-1).T
        for sign in (-1,1):
            u[:,axis] = sign
            level = (abs((physical(u)-CENTRE)/RADII[1])**POWERS[1]).sum(axis=1)
            np.testing.assert_allclose(level,1,atol=2e-14,rtol=0)


def audit(actual, expected):
    geometry_check()
    power = POWERS[1]
    exact = 8*prod(RADII[1])*gamma(1+1/power)**3/gamma(1+3/power)
    for idx, (row, wanted) in enumerate(zip(actual['runs'],expected['runs'],strict=True)):
        shape, order = CONFIGURATIONS[idx]
        assert row['shape'] == list(shape) == wanted['shape']
        assert row['q'] == order == wanted['q']
        assert abs(row['volume_mm3']/exact-1) < 4e-4
        assert row['uniform_patch_kpa'] < 1e-8
        assert row['overlap']['patch_error'] < 1e-9
        assert row['overlap']['missing_fraction'] < 1e-4
        for name in (*PAIR,REF):
            event = actual['events'][name]
            baseline = expected['events'][name]
            assert event['passed'] and event['min_c_min'] >= 50 and event['hot_any_c'] < 0
            np.testing.assert_allclose([event[k] for k in ('time_s','cold_heart_c','min_c_min')],
                                       [baseline[k] for k in ('time_s','cold_heart_c','min_c_min')],
                                       rtol=0,atol=2e-8)
            fresh = row['schedules'][name]
            prior = wanted['schedules'][name]
            assert fresh['steps'] == int(event['time_s']/(6*DT))+1 == prior['steps']
            np.testing.assert_allclose([fresh[k] for k in (*METRICS,'peak_time_s')],
                                       [prior[k] for k in (*METRICS,'peak_time_s')],
                                       rtol=0,atol=2e-8)
            assert fresh['release_12s']['steps'] == prior['release_12s']['steps'] == 40
            np.testing.assert_allclose(
                [fresh['release_12s'][k] for k in (*METRICS,'peak_time_s')],
                [prior['release_12s'][k] for k in (*METRICS,'peak_time_s')],
                rtol=0,atol=2e-8)
    old = json.loads((ROOT/'sources/staircase-results.json').read_text())
    old_gap = np.array([[1000*(old['runs'][grid]['schedules'][PAIR[0]]['overlap_affine'][key]-
                               old['runs'][grid]['schedules'][PAIR[1]]['overlap_affine'][key])
                        for key in METRICS] for grid in ('coarse','medium','fine')])
    assert list(np.sign(old_gap[:,0])) == [-1,1,-1]
    for name in (*PAIR,REF):
        assert abs(old['events'][name]['time_s']-actual['events'][name]['time_s']) < 2e-8
    gap, lead = comparisons(actual)
    assert (gap[:2] > 0).all() and (gap[2:] < 0).all()
    assert (lead > 0).all()
    for a,b in ((0,1),(2,3)):
        assert (abs(gap[a]-gap[b]) < np.minimum(abs(gap[a]),abs(gap[b]))).all()
    for index, row in enumerate(actual['runs']):
        early = row['schedules']
        for metric in METRICS:
            delta = early[PAIR[0]]['release_12s'][metric] - early[PAIR[1]]['release_12s'][metric]
            assert np.sign(delta) == np.sign(gap[index, METRICS.index(metric)])
    archived_gap, archived_lead = comparisons(expected)
    np.testing.assert_allclose(gap, archived_gap, rtol=0, atol=4e-5)
    np.testing.assert_allclose(lead, archived_lead, rtol=0, atol=4e-5)
    print('Staircase graded minus strict, kPa\n',old_gap)
    print('Fitted graded minus strict, kPa\n',gap)
    print('Fitted leads over tuned reference, kPa\n',lead)


if __name__ == '__main__':
    actual = regenerate()
    expected = json.loads((ROOT/'results.json').read_text())
    audit(actual,expected)
