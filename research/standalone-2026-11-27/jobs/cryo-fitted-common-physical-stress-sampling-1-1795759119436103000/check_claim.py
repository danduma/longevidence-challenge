"""Check the interpretation against the archive after verify.py has reproduced it."""
import json
from itertools import product

import numpy as np

from analysis import ROOT, SHAPES, Mesh, cloud
from verify import IndependentSample, contrasts, map_and_derivative


def check():
    result = json.loads((ROOT/'results.json').read_text())
    pattern = {
        ('native', False): (1, 1), ('native', True): (-1, -1),
        ('2', False): (-1, 1), ('2', True): (-1, -1),
        ('3', False): (1, 1), ('3', True): (1, 1),
    }
    for shape, q, order, early, gap, lead in contrasts(result):
        fine = shape == list(SHAPES[-1])
        np.testing.assert_array_equal(np.sign(gap), pattern[order, fine])
        expected_lead = (-1, -1) if order == '3' else (-1, 1) if order == '2' and fine else (1, 1)
        np.testing.assert_array_equal(np.sign(lead), expected_lead)

    corners = np.array(list(product((-1., 1.), repeat=3)))
    singular_values = np.linalg.svd(map_and_derivative(corners)[1], compute_uv=False)
    assert np.all(singular_values[:, 1] > 6) and np.all(singular_values[:, 2] < 1e-12)
    for order in (2, 3):
        reference, xyz, weights = cloud(order)
        _, jac = map_and_derivative(reference)
        cond = np.linalg.cond(jac)
        # All points and weights agree with independently generated clouds.
        for shape in SHAPES:
            independent = IndependentSample(Mesh(shape), order)
            np.testing.assert_allclose(independent.xyz, xyz, atol=1e-12, rtol=0)
            np.testing.assert_allclose(independent.volume, weights, atol=3e-9, rtol=0)
        for row in result['runs']:
            for item in row['schedules'].values():
                data = item['shared'][str(order)]
                peak = max(data['time_history'], key=lambda step: step['peak_mpa'])
                index = np.linalg.norm(xyz-peak['peak_xyz_mm'], axis=1).argmin()
                np.testing.assert_allclose(xyz[index], peak['peak_xyz_mm'], atol=1e-10, rtol=0)
                weight = weights[index]/sum(weights)
                fraction = weight*np.exp((peak['peak_mpa']-data['weighted_objective_mpa'])/.02)
                if order == 3:
                    assert 460 < cond[index] < 465
                    assert .14 < fraction < .38
                    assert np.all(abs(reference[index]) > .97)
                    assert peak['peak_mpa'] > 2*item['native']['peak_mpa']
                else:
                    assert cond[index] < 3
                    assert fraction < .025
    print('Both horizons: ranking assertions and corner Jacobian diagnostics passed.')


if __name__ == '__main__':
    check()
