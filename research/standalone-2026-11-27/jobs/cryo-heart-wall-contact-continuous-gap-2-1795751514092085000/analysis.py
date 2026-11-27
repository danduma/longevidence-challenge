"""Matched Vitrify phantom face-flux experiment. See report.md."""
import json
import math
import numpy as np
from lattice import ROOT, VECTORS, geometry, loading_field, sample

def material(t, wall):
    cp = np.select([t < -120, t <= -119],
                   [2402.3 + 7.6*t, 66747 + 545.24*t],
                   2678.8 + 11.73*t + 0.039*t*t)
    k = np.select([t <= -116, t < -96, t <= -77, t < -34],
                  [0.440 + 0.000816*t, 0.234 - 0.001*t,
                   0.374 + 0.000459*t, 0.374 + 0.0005*t],
                  0.374 + 0.000459*t)
    wall_cp = 1031 + 3.57*t - 0.014*t*t - 6.171e-5*t**3
    wall_rho = 931 - 0.49*t - 0.00091*t*t
    return np.where(wall, wall_cp*wall_rho, cp*(1078-0.43*t)), np.where(wall, 0.39+0.00052*t, k)


def face_geometry(phantom, points):
    """Return center-to-center links and the two analytic interval endpoints."""
    labels = phantom.labels
    active = labels >= 0
    index = np.full(labels.shape, -1)
    index[active] = np.arange(active.sum())
    spec = phantom.spec
    centre = np.asarray(spec.container_dims_mm)/2
    inner = centre - spec.wall_thickness_mm
    heart = np.asarray(spec.heart_dims_mm)/2
    heart *= np.cbrt(spec.heart_volume_ml*1000/(4*np.pi/3*np.prod(heart)))
    spacing = phantom.spacing_m*1000
    left, right, links = [], [], []
    for axis in range(3):
        lo = np.take(index, range(labels.shape[axis]-1), axis)
        hi = np.take(index, range(1, labels.shape[axis]), axis)
        valid = (lo >= 0) & (hi >= 0)
        offset = len(left)
        left.extend(lo[valid])
        right.extend(hi[valid])
        types = labels[active]
        for row, (i, j) in enumerate(zip(lo[valid], hi[valid])):
            if {types[i], types[j]} != {0, 2}:
                continue
            h, w = (i, j) if types[i] == 2 else (j, i)
            grid = np.unravel_index(np.flatnonzero(valid)[row], valid.shape)
            xyz = np.array([p[g] for p, g in zip(points, grid)])
            sign = 1 if h == i else -1
            if sign < 0:
                xyz[axis] += spacing
            transverse = np.arange(3) != axis

            def crossing(axes, exponent):
                rem = 1 - np.sum(abs((xyz-centre)[transverse]/axes[transverse])**exponent)
                plane = centre[axis] + sign*axes[axis]*max(0, rem)**(1/exponent)
                return np.clip(sign*(plane-xyz[axis])/spacing, 0, 1)

            b = crossing(inner, phantom.pillow_exponent)
            a = min(b, crossing(heart, 2))
            midpoint = xyz.copy()
            midpoint[axis] += sign * spacing * (a+b)/2
            links.append((offset+row, h, w, a, b, midpoint))
    return np.array(left), np.array(right), links


def run(phantom, points, native_shape, vector, mode, refinement=1):
    active = phantom.active_mask
    labels = phantom.labels[active]
    i, j, links = face_geometry(phantom, points)
    dx = phantom.spacing_m
    loading = loading_field(phantom, points, native_shape, vector)[active]
    volumes = np.ones(len(labels))
    fraction = np.full(len(i), .5)
    ids = np.array([row for row, *_ in links], dtype=int)
    gap = [link for link in links if link[4]-link[3] > 1e-9]
    if mode != "voxel":
        for row, h, w, a, b, _ in links:
            fraction[row] = b if i[row] == h else 1-b
    if mode == "gap" and gap:
        remove = np.array([row for row, *_ in gap])
        keep = ~np.isin(np.arange(len(i)), remove)
        i, j, fraction = i[keep], j[keep], fraction[keep]
        extra_load = []
        for _, h, w, a, b, midpoint in gap:
            # Allocate one third of the axial gap length as voxel volume,
            # corresponding to the two directional shares out of six faces.
            volumes[h] -= max(0, min(b,.5)-a)/3
            volumes[w] -= max(0, b-max(a,.5))/3
            node = len(labels) + len(extra_load)
            volumes = np.append(volumes, (b-a)/3)
            i, j = np.r_[i, h, node], np.r_[j, node, w]
            fraction = np.r_[fraction, a/((a+b)/2),
                              ((b-a)/2)/(1-(a+b)/2)]
            radial = np.linalg.norm((midpoint-np.asarray(native_shape)*dx*500) /
                                    (np.asarray(phantom.spec.heart_dims_mm)/2))
            weights = np.exp(-.5*((radial-np.array([1.05,1.45,2.]))/[.32,.42,.65])**2)
            extra_load.append(np.dot(weights, vector[:3])/weights.sum())
        lengths = np.r_[np.ones(keep.sum()), np.ravel(
            [[(a+b)/2, 1-(a+b)/2] for _,_,_,a,b,_ in gap])]
        loading = np.r_[loading, extra_load]
        labels = np.r_[labels, np.ones(len(gap), dtype=int)]
    else:
        lengths = np.ones(len(i))
    assert np.all(volumes > 0) and abs(volumes.sum()-active.sum()) < 1e-8
    sweep = np.linspace(-180,20,400)
    cap, k = material(np.tile(sweep,2), np.repeat([False,True],400))
    base_dt = .4*cap.min()*dx**2/(6*k.max())
    degree = np.bincount(np.r_[i,j], np.tile(1/lengths,2), minlength=len(labels))
    # Gap nodes may be stiff. Integrate each independently by eliminating its
    # incident flux implicitly within the otherwise explicit conservative step.
    wall = labels == 0
    specimen = labels > 0
    original = specimen & (np.arange(len(labels)) < active.sum())
    heart = labels >= 2
    dt = min(base_dt, .4*cap.min()*dx**2/k.max() * min(volumes[original | wall]/degree[original | wall]))/refinement
    n = math.ceil(70/dt)
    duty = np.interp(np.linspace(0,1,n), [0,.5,1], vector[3:])
    t = np.full(len(labels), -150.)
    minimum, full_minimum, startup = np.inf, np.inf, np.inf
    maximum_residual = 0.
    result = {}
    if mode == "gap" and gap:
        cut = len(i)-2*len(gap)
        h = i[cut::2]
        g = j[cut::2]
        w = j[cut+1::2]
    else:
        cut = len(i)
    for step in range(n):
        cap, k = material(t,wall)
        capacity = volumes*cap
        q = volumes*1000*loading*(691-4.6*(np.minimum(t,0)+80))*duty[step]
        conductance = 1/(dx**2*lengths*(fraction/k[i]+(1-fraction)/k[j]))
        flux = conductance*(t[j]-t[i])
        if cut != len(i):
            gh, gw = conductance[cut::2], conductance[cut+1::2]
            dgap = (q[g]+gh*(t[h]-t[g])+gw*(t[w]-t[g]))/(capacity[g]+dt*(gh+gw))
            flux[cut::2] += gh*dt*dgap
            flux[cut+1::2] -= gw*dt*dgap
        power = q + np.bincount(i, flux, minlength=len(t)) - np.bincount(j, flux, minlength=len(t))
        maximum_residual = max(maximum_residual, abs(power.sum()-q.sum())/q.sum())
        rate = 60*power/capacity
        limiting = np.argmin(np.where(original,rate,np.inf))
        full_minimum = min(full_minimum, rate[specimen].min())
        if step*dt < 1:
            startup = min(startup, rate[limiting])
        if rate[limiting] < minimum:
            minimum = rate[limiting]
            result.update(time_s=step*dt, temp_c=t[limiting], label=int(labels[limiting]))
        after = t + dt*rate/60
        if "event_s" not in result and after[heart].min() >= -35:
            cold = heart & (t < -35)
            share = np.max((-35-t[cold])/(after[cold]-t[cold]))
            result.update(event_s=(step+share)*dt,
                          hot_event_c=(t+share*(after-t))[specimen].max())
        t = after
    result.update(rate=minimum, full_rate=full_minimum, startup=startup,
                  cold70=t[specimen].min(), hot70=t[specimen].max(),
                  gap_faces=len(gap), contact_faces=len(ids)-len(gap), dt=dt,
                  minimum_gap_mm=min(((b-a)*dx*1000 for _,_,_,a,b,_ in gap), default=0.),
                  maximum_gap_mm=max(((b-a)*dx*1000 for _,_,_,a,b,_ in gap), default=0.),
                  flux_residual=maximum_residual)
    return result


CASES = {"native": (1.1, (0,0,0)), "x": (1.1,(.5,0,0)),
         "y": (1.1,(0,.5,0)), "z": (1.1,(0,0,.5)), "fine_y": (.75,(0,.5,0))}


def calculate(refinement=1):
    results = {}
    for name, (spacing, phase) in CASES.items():
        native = geometry.rat_heart_cryobag(spacing)
        phantom, points = sample(spacing, native, phase)
        results[name] = {
            design: {mode: run(phantom, points, native.shape, vector, mode, refinement)
                     for mode in ("voxel", "series", "gap")}
            for design, vector in VECTORS.items()
        }
    return results


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2, allow_nan=False))
