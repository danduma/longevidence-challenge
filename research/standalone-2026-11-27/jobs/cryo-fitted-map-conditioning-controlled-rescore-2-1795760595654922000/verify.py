"""Re-solve both maps and independently score pointwise principal stresses."""
import json
from itertools import product
from math import gamma
from pathlib import Path

import numpy as np

from analysis import regenerate
from cut_geometry import CENTRE, POWERS, RADII
from regular import INNER, P, R, RegularMesh, map_block


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


def independent_cloud():
    fine = np.array([16,12,8])
    boxes = np.array(list(product(*(range(n) for n in fine))))
    combinations = np.array(list(product(range(3),repeat=3)))
    root = np.sqrt(3/5)
    nodes = np.array([-root,0,root])
    weights = np.array([5/9,8/9,5/9])
    reference = (-1+(2*boxes[:,None]+1+nodes[combinations])/fine).reshape(-1,3)
    xyz,jac = map_and_derivative(reference)
    weight = np.linalg.det(jac)*np.tile(weights[combinations].prod(axis=1),len(boxes))/fine.prod()
    singular = np.linalg.svd(jac,compute_uv=False)
    assert np.all(weight > 0)
    return reference,xyz,weight,singular[:,0] <= 10*singular[:,2]


def differential(local,jac):
    signs = np.array(list(product((-1,1),repeat=3)))
    def basis(points):
        return np.prod(.5+(points[:,None]-.5)*signs,axis=2)
    delta = 1e-5
    reference = np.stack([(basis(local+h)-basis(local-h))/(2*delta)
                          for h in np.eye(3)*delta],axis=1)
    return np.linalg.solve(jac.transpose(0,2,1),reference).transpose(0,2,1)


class IndependentSample:
    def __init__(self,mesh,data):
        ref,xyz,_,_ = data
        if isinstance(mesh,RegularMesh):
            self.owner,self.gradient = mesh.locate(xyz,differential)
        else:
            grid = (ref+1)*mesh.shape/2
            owner = np.floor(grid).astype(int)
            self.owner = np.ravel_multi_index(owner.T,mesh.shape)
            jac = map_and_derivative(ref)[1]
            self.gradient = differential(grid-owner,jac*mesh.width)

    def evaluate(self,u,history,shear,eigenstrain,bulk):
        i = self.owner
        du = np.einsum('pnk,pnj->pjk',self.gradient,u[i].reshape(-1,8,3))
        strain = .5*(du+du.transpose(0,2,1))
        trace = np.trace(strain,axis1=1,axis2=2)
        strain[:,range(3),range(3)] -= trace[:,None]/3
        sigma = 2*shear[i,None,None]*strain+history[i]
        sigma[:,range(3),range(3)] += bulk*(trace-eigenstrain[i])[:,None]
        return np.linalg.eigvalsh(sigma)


def independent_score(p,w):
    peak = p.max()
    return np.array([peak,peak+20_000*np.log(w@np.exp((p-peak)/20_000).sum(axis=1))])/1000


def aggregate(trace):
    values = np.array(trace)
    center = values[:,1].max()
    return np.array([values[:,0].max(), center+20*np.log(np.exp((values[:,1]-center)/20).sum())])


def validate_geometry():
    # Distinguish the element-axis scaling from the unit block Jacobian.
    q = np.array(list(product(np.linspace(0,1,5),repeat=3)))
    for block in range(7):
        x,j = map_block(block,q)
        numerical = np.stack([(map_block(block,q+h)[0]-map_block(block,q-h)[0])/(2e-6)
                              for h in np.eye(3)*1e-6],axis=-1)
        np.testing.assert_allclose(j,numerical,atol=2e-8,rtol=0)
        if block:
            outer = q[:,2] == 1
            np.testing.assert_allclose(np.sum(abs((x[outer]-CENTRE)/R)**P,axis=1),1)
            analytic_lower = 4*INNER**2*(3**(-1/P)-INNER)*R.prod()
            assert np.min(abs(np.linalg.det(j))) >= analytic_lower


def compare(actual,expected,baseline):
    assert actual['points'] == expected['points'] == 41472
    assert actual['retained_points'] == expected['retained_points'] < actual['points']
    np.testing.assert_allclose(actual['retained_weight'],expected['retained_weight'],atol=1e-9,rtol=0)
    assert .99349 < actual['retained_weight'] < .99351
    exact = 8*R.prod()*gamma(1+1/P)**3/gamma(1+3/P)
    assert len(actual['runs']) == len(expected['runs']) == 8
    for a,b in zip(actual['runs'],expected['runs'],strict=True):
        f,l,q = [a[k] for k in ('family','level','q')]
        for k in ('family','level','q','cells'):
            assert a[k] == b[k]
        assert abs(a['volume_mm3']/exact-1) < .001
        if f == 'regular':
            assert a['conditioning']['max_condition'] < 10
            np.testing.assert_allclose(list(a['conditioning'].values()),list(b['conditioning'].values()))
        incumbent = next(row for row in baseline['runs'] if row['q'] == q
                         and np.prod(row['shape']) == (192,1536)[l])
        assert a['scores'].keys() == b['scores'].keys() == incumbent['schedules'].keys()
        aggregates = {}
        early = {}
        for name,scores in a['scores'].items():
            ref = incumbent['schedules'][name]['shared']['3']['time_history']
            np.testing.assert_allclose(scores['times'],[r['time_s'] for r in ref],atol=1e-9,rtol=0)
            aggregates[name] = {}
            early[name] = {}
            first_12s = np.asarray(scores['times']) <= 12
            for region in ('all','interior'):
                np.testing.assert_allclose(scores[region],b['scores'][name][region],atol=.003,rtol=0)
                aggregates[name][region] = aggregate(scores[region])
                early[name][region] = aggregate(np.asarray(scores[region])[first_12s])
            if f == 'cube':
                np.testing.assert_allclose(scores['all'],[[r['peak_mpa']*1000,20*r['weighted_logsum']]
                                                          for r in ref],atol=.003,rtol=0)
        for region in ('all','interior'):
            pair = aggregates['graded'][region]-aggregates['uniform_strict'][region]
            early_pair = early['graded'][region] - early['uniform_strict'][region]
            assert np.all(np.sign(early_pair) == np.sign(pair))
            tuned = aggregates['uniform_released_policy'][region]
            peak = [round(by_region[region][0],2) for by_region in aggregates.values()]
            if f == 'cube' and region == 'all':
                assert np.all(pair > 70)
                assert min(peak) > 2800
                for regions in aggregates.values():
                    retained_fraction = actual['retained_weight'] * np.exp(
                        (regions['interior'][1] - regions['all'][1]) / 20
                    )
                    assert retained_fraction < 5e-28
            elif f == 'cube':
                assert np.all(pair < -1)
                assert max(peak) < 1600
            else:
                for regions in aggregates.values():
                    difference = regions['interior'] - regions['all']
                    expected_shift = -20*np.log(actual['retained_weight'])
                    np.testing.assert_allclose(
                        difference, [0, expected_shift], atol=5e-9, rtol=0
                    )
                if l == 0:
                    assert np.all(pair > 6)
                else:
                    assert np.all(pair < -2)
                assert 1300 < min(peak) <= max(peak) < 1460
                assert np.all(tuned < aggregates['graded'][region])
                assert np.all(tuned < aggregates['uniform_strict'][region])
            print(f,l,q,region,'graded minus strict',pair,'tuned',tuned,'peaks',peak,flush=True)


if __name__ == '__main__':
    root = Path(__file__).resolve().parent
    validate_geometry()
    fresh = regenerate(independent_score,IndependentSample,independent_cloud)
    compare(fresh,json.loads((root/'results.json').read_text()),
            json.loads((root/'sources/common-results.json').read_text()))
