"""Re-solve the benchmark and independently evaluate the spatial stress scores."""
import json
from itertools import product
from math import gamma

import numpy as np

from analysis import ROOT, regenerate
from cut_geometry import CENTRE, POWERS, RADII


def map_and_derivative(points):
    """Analytic Jacobian of the inherited spherical and radial maps."""
    sphere = np.empty_like(points)
    jac = np.empty((len(points), 3, 3))
    for i in range(3):
        j, k = (i+1) % 3, (i+2) % 3
        a, b, c = points[:, [i, j, k]].T
        f = np.sqrt(1-(b*b+c*c)/2+b*b*c*c/3)
        sphere[:, i] = a*f
        jac[:, i, i] = f
        jac[:, i, j] = a*b*(-.5+c*c/3)/f
        jac[:, i, k] = a*c*(-.5+b*b/3)/f
    power = POWERS[1]
    l2_squared = np.sum(sphere**2, axis=1)
    lp_power = np.sum(abs(sphere)**power, axis=1)
    radial = np.sqrt(l2_squared)/lp_power**(1/power)
    dlog = (sphere/l2_squared[:, None]
            - np.sign(sphere)*abs(sphere)**(power-1)/lp_power[:, None])
    radial_jac = radial[:, None, None]*(np.eye(3)+sphere[:, :, None]*dlog[:, None])
    return CENTRE+sphere*radial[:, None]*RADII[1], RADII[1][:, None]*(radial_jac@jac)


class IndependentSample:
    def __init__(self, mesh, order):
        fine = np.array([16, 12, 8])
        if order == 2:
            roots, weights = [-1/np.sqrt(3), 1/np.sqrt(3)], [1, 1]
        else:
            r = np.sqrt(3/5)
            roots, weights = [-r, 0, r], [5/9, 8/9, 5/9]
        boxes = np.array(list(product(*(range(n) for n in fine))))
        combos = np.array(list(product(range(order), repeat=3)))
        x = -1+(2*boxes[:, None]+1+np.asarray(roots)[combos])/fine
        x = x.reshape(-1, 3)
        self.xyz, jac = map_and_derivative(x)
        self.volume = np.linalg.det(jac)*np.tile(np.asarray(weights)[combos].prod(1),len(boxes))/fine.prod()
        assert np.all(self.volume > 0)
        exact_volume = 8*np.prod(RADII[1])*gamma(1+1/POWERS[1])**3/gamma(1+3/POWERS[1])
        assert abs(sum(self.volume)/exact_volume-1) < 1e-3
        self.log_weight = np.log(self.volume/sum(self.volume))
        grid = (x+1)*mesh.shape/2
        i, j, k = owner = np.floor(grid).astype(int).T
        _, ny, nz = mesh.shape
        self.owner = (i*ny+j)*nz+k
        assert np.all((owner >= 0)&(owner < mesh.shape[:, None]))
        assert np.all(np.bincount(self.owner) == order**3*np.prod(fine//mesh.shape))
        local = grid-owner.T
        assert np.all((local > 0)&(local < 1))
        signs = np.array(list(product((-1, 1), repeat=3)))

        def basis(at):
            return np.prod(.5+(at[:, None]-.5)*signs, axis=2)

        # Centered differences give the exact directional derivative of a
        # trilinear basis in reference space, up to floating-point roundoff.
        delta = .00001
        ref_gradient = np.stack([
            (basis(local+d)-basis(local-d))/(2*delta)
            for d in np.eye(3)*delta
        ], axis=2)*mesh.shape/2
        self.gradient = np.linalg.solve(jac.transpose(0, 2, 1),
                                        ref_gradient.transpose(0, 2, 1))

    def evaluate(self, u, history, shear, eigenstrain, bulk):
        idx = self.owner
        derivative = np.einsum('pan,pnb->pab', self.gradient, u[idx].reshape(-1,8,3))
        strain = .5*(derivative+derivative.transpose(0, 2, 1))
        trace = np.einsum('pii->p', strain)
        strain[:, range(3), range(3)] -= trace[:, None]/3
        stress = 2*shear[idx, None, None]*strain+history[idx]
        stress[:, range(3), range(3)] += bulk*(trace-eigenstrain[idx])[:, None]
        return np.linalg.eigvalsh(stress)


def measures(record, order, early):
    if order == 'native':
        value = record['native_12s' if early else 'native']
        return np.array([value['objective_mpa'], value['peak_mpa']])
    data = record['shared'][order]
    trace = [row for row in data['time_history'] if not early or row['time_s'] <= 12]
    logs = np.array([row['weighted_logsum'] for row in trace])
    smooth = .02*(logs.max()+np.log(np.exp(logs-logs.max()).sum()))
    stored_key = 'weighted_objective_12s_mpa' if early else 'weighted_objective_mpa'
    np.testing.assert_allclose(smooth, data[stored_key], atol=1e-12, rtol=0)
    return np.array([smooth, max(row['peak_mpa'] for row in trace)])


def audit(actual, expected, baseline):
    assert len(actual['runs']) == len(expected['runs']) == len(baseline['runs']) == 4
    assert actual['events'].keys() == expected['events'].keys()
    for name in actual['events']:
        np.testing.assert_allclose(actual['events'][name]['time_s'],
                                   expected['events'][name]['time_s'], atol=1e-8, rtol=0)
    for row, archived, native in zip(actual['runs'], expected['runs'], baseline['runs'], strict=True):
        assert row['shape'] == archived['shape'] == native['shape']
        assert row['q'] == archived['q'] == native['q']
        assert row['schedules'].keys() == archived['schedules'].keys() == native['schedules'].keys()
        for name, rec in row['schedules'].items():
            target = archived['schedules'][name]
            for key in ('native', 'native_12s'):
                old = native['schedules'][name]
                if key.endswith('12s'):
                    old = old['release_12s']
                for metric in ('objective_mpa', 'peak_mpa', 'peak_time_s', 'steps'):
                    np.testing.assert_allclose(rec[key][metric], [old[metric], target[key][metric]],
                                               atol=1e-8, rtol=0)
            for order in ('2', '3'):
                a, b = rec['shared'][order], target['shared'][order]
                assert a['points'] == b['points'] == 1536*int(order)**3
                np.testing.assert_allclose(a['volume_mm3'], b['volume_mm3'], atol=1e-4, rtol=0)
                assert len(a['time_history']) == len(b['time_history']) == 220
                # The production map Jacobian uses finite differences. Allow
                # 0.5 Pa for the analytic comparison near ill-conditioned corners.
                for metric, tol in (('time_s', 1e-8), ('peak_mpa', 5e-7), ('weighted_logsum', 2.5e-5)):
                    np.testing.assert_allclose([step[metric] for step in a['time_history']],
                                               [step[metric] for step in b['time_history']],
                                               atol=tol, rtol=0)
                for trace in (a, b):
                    peak = max(trace['time_history'], key=lambda step: step['peak_mpa'])
                    assert 0 < peak['time_s'] <= 12
                np.testing.assert_allclose(
                    max(a['time_history'], key=lambda step: step['peak_mpa'])['peak_xyz_mm'],
                    max(b['time_history'], key=lambda step: step['peak_mpa'])['peak_xyz_mm'],
                    atol=1e-8, rtol=0)
            for order, early in product(('native', '2', '3'), (False, True)):
                np.testing.assert_allclose(measures(rec, order, early),
                                           measures(target, order, early), atol=5e-7, rtol=0)


def contrasts(result):
    for row in result['runs']:
        for order, early in product(('native', '2', '3'), (False, True)):
            s = row['schedules']
            graded, strict, tuned = (measures(s[name], order, early)
                                     for name in ('graded', 'uniform_strict', 'uniform_released_policy'))
            yield row['shape'], row['q'], order, early, (graded-strict)*1000, (np.minimum(graded, strict)-tuned)*1000


if __name__ == '__main__':
    baseline = json.loads((ROOT/'sources/fitted-results.json').read_text())
    actual = regenerate(IndependentSample)
    expected = json.loads((ROOT/'results.json').read_text())
    audit(actual, expected, baseline)
    for shape, q, order, early, gap, lead in contrasts(actual):
        print(shape, q, order, early, 'gap', gap, 'reference lead', lead)
