"""Recompute the full comparison and independently integrate startup."""
import json
import math

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix

from analysis import CASES, calculate, face_geometry, material
from lattice import ROOT, VECTORS, geometry, loading_field, sample


def independent_startup(name, mode, dt=.002, design="uniform_strict", original_only=False):
    """Use predicate roots, sparse incidence and a common explicit step."""
    spacing, phase = CASES[name]
    native = geometry.rat_heart_cryobag(spacing)
    phantom, axes = sample(spacing, native, phase)
    labels = phantom.labels
    coords = np.argwhere(labels >= 0)
    index = {tuple(x): i for i, x in enumerate(coords)}
    kinds = labels[tuple(coords.T)].tolist()
    vector = VECTORS[design]
    load = loading_field(phantom, axes, native.shape, vector)[labels >= 0].tolist()
    original_count = len(load)
    volumes = [1.] * original_count
    centre = np.asarray(native.spec.container_dims_mm)/2
    heart = np.asarray(native.spec.heart_dims_mm)/2
    heart *= (native.spec.heart_volume_ml*1000/(4*np.pi/3*np.prod(heart)))**(1/3)
    inner = centre-native.spec.wall_thickness_mm
    edges = []
    intervals = []
    for i, xyz in enumerate(coords):
        for direction in np.eye(3,dtype=int):
            j = index.get(tuple(xyz+direction))
            if j is None:
                continue
            if {kinds[i], kinds[j]} != {0,2}:
                edges.append((i,j,.5,.5))
                continue
            h,w = (i,j) if kinds[i] == 2 else (j,i)
            start = np.array([ax[a] for ax,a in zip(axes,coords[h])])
            end = np.array([ax[a] for ax,a in zip(axes,coords[w])])

            def root(radius, power):
                def level(s):
                    return np.sum(abs((start+s*(end-start)-centre)/radius)**power)-1
                return brentq(level, 0, 1) if level(1) > 0 else 1.

            b = root(inner, native.pillow_exponent)
            a = min(b,root(heart,2))
            intervals.append((b-a)*spacing)
            if mode == "voxel":
                edges.append((h,w,.5,.5))
            elif mode == "series" or b-a <= 1e-9:
                edges.append((h,w,b,1-b))
            else:
                g = len(volumes)
                volumes[h] -= max(0,min(.5,b)-a)/3
                volumes[w] -= max(0,b-max(.5,a))/3
                volumes.append((b-a)/3)
                kinds.append(1)
                midpoint = start + (a+b)/2 * (end-start)
                origin = np.asarray(native.shape) * spacing / 2
                radius = np.linalg.norm((midpoint-origin)/(np.asarray(native.spec.heart_dims_mm)/2))
                weights = np.exp(-.5*((radius-np.array([1.05,1.45,2.]))/[.32,.42,.65])**2)
                gap_loading = np.dot(weights, vector[:3])/weights.sum()
                load.append(0. if mode == "unheated_gap" else gap_loading)
                edges.extend([(h,g,a,(b-a)/2),(g,w,(b-a)/2,1-b)])
    analytic = face_geometry(phantom,axes)[2]
    np.testing.assert_allclose(sorted(intervals),sorted((b-a)*spacing for _,_,_,a,b,_ in analytic),atol=1e-10)
    left,right,dl,dr = np.asarray(edges).T
    left,right = left.astype(int),right.astype(int)
    n = len(volumes)
    columns = np.arange(len(left))
    divergence = coo_matrix((np.r_[np.ones(len(left)), -np.ones(len(left))],
                             (np.r_[left,right],np.tile(columns,2))),shape=(n,len(left))).tocsr()
    t = np.full(n,-150.)
    wall = np.array(kinds) == 0
    measured = ~wall
    if original_only:
        measured[original_count:] = False
    volumes, load = np.asarray(volumes), np.asarray(load)
    assert np.all(volumes > 0) and abs(volumes.sum()-original_count) < 1e-8
    minimum = math.inf
    for step in range(round(1/dt)):
        cp,k = material(t,wall)
        duty = np.interp(2*step*dt/70,[0,1,2],vector[3:])
        power = volumes*load*1000*(691-4.6*(t+80))*duty
        flux = (t[right]-t[left])/((spacing/1000)**2*(dl/k[left]+dr/k[right]))
        rate = (power+divergence@flux)/(cp*volumes)
        minimum = min(minimum,60*rate[measured].min())
        t += dt*rate
    return minimum


def main():
    cached = json.loads((ROOT/"results.json").read_text())
    prior = json.loads((ROOT/"sources/prior-thermal.json").read_text())
    fresh = calculate()
    assert set(cached) == set(fresh)
    for name, designs in fresh.items():
        spacing, phase = CASES[name]
        anchor = next(row for row in prior if row["spacing_mm"] == spacing and tuple(row["phase"]) == phase)
        assert set(cached[name]) == set(designs)
        for design, modes in designs.items():
            assert set(cached[name][design]) == set(modes)
            for mode, row in modes.items():
                assert set(row) == set(cached[name][design][mode])
                for key, value in row.items():
                    assert math.isclose(value,cached[name][design][mode][key],rel_tol=1e-9,abs_tol=1e-9), (name,design,mode,key)
                assert row["flux_residual"] < 1e-12
                assert row["event_s"] < 70 and -35 <= row["cold70"] <= row["hot70"] < 0
                assert row["time_s"] < 1 and row["temp_c"] < -149.6
            assert abs(modes["voxel"]["rate"]-anchor["designs"][design]["rate_70"]) < 1e-8
            if name in ("z","fine_y"):
                assert modes["voxel"] == modes["gap"] == modes["series"]
            else:
                assert modes["gap"]["gap_faces"] == 10
                assert modes["series"]["rate"] <= modes["voxel"]["rate"]
                assert modes["gap"]["full_rate"]-modes["voxel"]["rate"] > 1
                assert abs(modes["gap"]["full_rate"]-fresh["fine_y"][design]["voxel"]["rate"]) < .01
        strict = designs["uniform_strict"]
        assert (strict["voxel"]["rate"] >= 51.5) == (name not in ("x","y"))
        assert strict["gap"]["full_rate"] >= 51.5

    ablation = json.loads((ROOT/"startup-original-nodes.json").read_text())
    refined_cache = json.loads((ROOT/"startup-refined.json").read_text())
    for name in CASES:
        rates = {mode: independent_startup(name,mode) for mode in ("voxel","series","gap")}
        half_step = {mode: independent_startup(name,mode,.001) for mode in rates}
        for mode, value in rates.items():
            assert abs(value-refined_cache[name][mode]) < 1e-9
            assert abs(value-half_step[mode]) < .01
        original_node_rates = {}
        for dt, baseline in ((.002, rates), (.001, half_step)):
            measured = dict(baseline)
            for mode in ("gap", "unheated_gap"):
                measured[mode] = independent_startup(name, mode, dt, original_only=True)
            for mode, value in measured.items():
                assert abs(value-ablation[name][str(dt)][mode]) < 1e-9
            original_node_rates[dt] = measured
            if name in ("z", "fine_y"):
                assert measured["gap"] == measured["unheated_gap"]
            else:
                assert measured["unheated_gap"] < measured["voxel"] < measured["gap"]
            if name == "y":
                assert measured["unheated_gap"] < 50 < measured["series"] < 51.5 < measured["gap"]
        ablation_delta = abs(original_node_rates[.002]["unheated_gap"]
                             - original_node_rates[.001]["unheated_gap"])
        assert ablation_delta < .01
        strict = fresh[name]["uniform_strict"]
        assert abs(rates["gap"]-strict["gap"]["full_rate"]) < .015
        assert (rates["voxel"] >= 51.5) == (name != "y")
        assert (rates["series"] >= 51.5) == (name != "y")
        assert rates["gap"] >= 51.5
        if name in ("native","x","y"):
            assert rates["gap"]-rates["voxel"] > .8
            assert rates["series"] < rates["voxel"]
        else:
            assert rates["voxel"] == rates["series"] == rates["gap"]
        print(name, {mode: round(float(value), 6) for mode, value in rates.items()},
              "unheated gap on original nodes", round(float(original_node_rates[.002]["unheated_gap"]), 6))
    print("Verified full-run records, matched controls, explicit startup thresholds and the heat-source ablation.")


if __name__ == "__main__":
    main()
