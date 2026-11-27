"""Rebuild geometry and independently march the decisive temperature trajectories."""
import json
from functools import lru_cache

import numpy as np
from scipy.special import gamma
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu

from analysis import CASES, ROOT, cut_graph, voxel_graph
from cut_geometry import CAVITY, CENTRE, HEART, POWERS, RADII
from lattice import VECTORS


def exact_volumes():
    p = POWERS[0]

    def super_volume(radii):
        return 8 * radii.prod() * gamma(1+1/p)**3 / gamma(1+3/p)

    inner, outer = super_volume(RADII[1]), super_volume(CENTRE)
    # On the ellipsoid, the inner-bag level is a convex function of
    # squared normalized coordinates. Its maximum occurs on an axis.
    ratio = HEART/RADII[1]
    assert p >= 2 and ratio.max() < 1
    # The inner superellipsoid contains its smallest-axis ball. Convexity
    # puts that ball, scaled by (1-ratio.max()), around every heart point.
    # The bound is attained on z, the axis extremizing both radii tests.
    assert np.argmax(ratio) == np.argmin(RADII[1]) == 2
    assert .80 < RADII[1][2]-HEART[2] < .82
    heart = 4*np.pi*HEART.prod()/3
    separation = .8*HEART[0]/CAVITY[0]
    overlap = np.pi*(4+separation)*(2-separation)**2/12
    cavity = CAVITY.prod()*(8*np.pi/3-overlap)
    return np.array([outer-inner,inner-heart,heart-cavity,cavity])


def exact_heart_area():
    # Integrate the surface Jacobian of the ellipsoid's map from the sphere.
    z, w = np.polynomial.legendre.leggauss(100)
    phi = 2*np.pi*np.arange(200)/200
    r = np.sqrt(1-z*z)
    normal = np.broadcast_arrays(r[:,None]*np.cos(phi),
                                 r[:,None]*np.sin(phi),z[:,None])
    jacobian = HEART.prod()*np.sqrt(sum((n/a)**2 for n,a in zip(normal,HEART)))
    return 2*np.pi*np.dot(w,jacobian.mean(axis=1))


def independent(mesh, vector, dt, theta):
    """Sparse incidence fluxes and increment solve, with independent material evaluation."""
    labels, vol, xyz, (i,j,area,dl,dr), _ = mesh
    wall = labels == 0
    position = (xyz-np.array([11.55,8.8,6.05]))/np.array([8.75,6.5,5.5])
    r = np.sqrt(np.sum(position**2, axis=1))
    weights = np.vstack([np.exp(-((r-c)/s)**2/2)
                         for c,s in ((1.05,.32),(1.45,.42),(2.,.65))])
    shell = (np.array(vector[:3])[:,None]*weights).sum(axis=0)/weights.sum(axis=0)
    loading = np.choose(labels,[0.,shell,1.47,10.])
    cols = np.arange(len(i))
    incidence = coo_matrix((np.r_[np.ones(len(i)), -np.ones(len(i))],
                            (np.r_[i,j],np.tile(cols,2))),
                           shape=(len(labels),len(i))).tocsr()

    incidence_t = incidence.T.tocsr()

    def capacity(t):
        wax = ((-6.171e-5*t-.014)*t+3.57)*t+1031
        return vol*np.where(wall,wax*(931-.49*t-.00091*t*t),
                            (2402.3+7.6*t)*(1078-.43*t))

    def conductance(t):
        k = np.where(wall,.39+.00052*t,.440+.000816*t)
        return 1e6*area/(dl/k[i]+dr/k[j])

    t = np.full(len(labels),-150.)
    mass0 = capacity(t)
    g0 = conductance(t)
    if theta:
        ref = diags(mass0)+theta*dt*(incidence.multiply(g0)@incidence.T)
        factor = splu(ref.tocsc())
    minimum = np.inf
    for n in range(round(1/dt)):
        mass = capacity(t)
        g = conductance(t)
        assert -151 < t.min() <= t.max() < -120
        duty = vector[3]+(vector[4]-vector[3])*n*dt/35
        source = 1000*vol*loading*(323-4.6*t)*duty
        rhs = source-incidence@(g*(incidence_t@t))
        if theta:
            v = factor.solve(rhs)
            for _ in range(6):
                residual = rhs-mass*v-theta*dt*(incidence@(g*(incidence_t@v)))
                if np.max(abs(residual))/np.max(abs(rhs)) < 3e-9:
                    break
                v += factor.solve(residual)
            else:
                raise AssertionError("Implicit flux update did not converge")
        else:
            v = rhs/mass
        assert abs((mass*v).sum()-source.sum())/source.sum() < 1e-7
        minimum = min(minimum,60*v[~wall].min())
        t += dt*v
    return minimum


@lru_cache
def build(name):
    h,q = CASES[name]
    return cut_graph(h,q) if q > 0 else voxel_graph(h,q == -1)


def main():
    saved = json.loads((ROOT/"results.json").read_text())
    assert set(saved) == set(CASES)
    volumes = exact_volumes()
    area = exact_heart_area()
    rates = {}
    for name, (_, q) in CASES.items():
        mesh = build(name)
        phase,volume,_,_,stats = mesh
        row = saved[name]
        assert len(phase) == row["states"]
        totals = np.bincount(phase,volume)
        np.testing.assert_allclose(totals,row["phase_volumes_mm3"],rtol=1e-10)
        if q > 0:
            np.testing.assert_allclose(totals,volumes,rtol=.005,atol=.1)
            assert abs(stats["interface_mm2"]["1-2"]-area) < .005*area
            assert stats["interface_mm2"]["0-2"] == 0
            np.testing.assert_allclose(list(stats["interface_mm2"].values()),
                                       list(row["interface_mm2"].values()),atol=1e-8)
        for design,vector in VECTORS.items():
            for dt in (.002,.001):
                stored = row["designs"][design][str(dt)]
                val = independent(mesh,vector,dt,float(q>0))
                assert abs(val-stored["rate"]) < 2e-5, (name,design,dt,val,stored)
            assert abs(row["designs"][design]["0.002"]["rate"]-
                       row["designs"][design]["0.001"]["rate"]) < .025
        rates[name] = row["designs"]["uniform_strict"]["0.001"]["rate"]
    assert rates["voxel_1.1"] < 51.5
    assert all(value > 51.5 for name,value in rates.items() if name != "voxel_1.1")
    for name, mesh in [("off-sequence quadrature",cut_graph(1.1,24)),
                       ("refined cut",build("cut_0.75_q16"))]:
        rate = independent(mesh,VECTORS["uniform_strict"],.0005,.5)
        assert rate > 51.5 and abs(rate-rates["cut_1.1_q32"]) < .1, (name,rate)
        print(name,round(rate,6))
    print("Verified startup rates in °C/min",rates)


if __name__ == "__main__":
    main()
