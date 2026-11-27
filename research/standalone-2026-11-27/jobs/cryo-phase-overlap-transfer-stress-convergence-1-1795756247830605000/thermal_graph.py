"""Matched one-second heat-balance comparison. All geometry inputs are local."""
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu

from cut_geometry import CENTRE, NATIVE, ORIGIN, PAIRS, bounds, classify, integrate
from gap_model import face_geometry, material
from lattice import VECTORS

ROOT = Path(__file__).resolve().parent


def centre_grid(axes):
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)


def cut_graph(h, q):
    cells, moments, interfaces, faces = integrate(h, q)
    vol = moments[..., 0]
    present = vol > 0
    uid = np.full(vol.shape, -1)
    uid[present] = np.arange(present.sum())
    cell, phase = np.where(present)
    volume = vol[present]
    xyz = moments[..., 1:][present] / volume[:, None]
    edges = []
    lost = 0.
    for axis, aperture in enumerate(faces):
        pair = np.moveaxis(cells, axis, -1)
        low, high = pair[..., :-1].ravel(), pair[..., 1:].ravel()
        for p in range(4):
            i, j = uid[low, p], uid[high, p]
            area = aperture[low, p]
            valid = (i >= 0) & (j >= 0)
            lost += area[~valid].sum()
            ok = valid & (area > 0)
            i, j, area = i[ok], j[ok], area[ok]
            # Use the normal projection of each centroid distance on grid faces.
            plane = bounds(h)[axis][1:-1]
            plane = np.broadcast_to(plane, pair[..., :-1].shape).ravel()[ok]
            edges.append((i, j, area, plane - xyz[i, axis], xyz[j, axis] - plane))
    for n, (p, r) in enumerate(PAIRS):
        area = interfaces[:, n, 0]
        i, j = uid[:, p], uid[:, r]
        valid = (i >= 0) & (j >= 0)
        lost += area[~valid].sum()
        ok = valid & (area > 0)
        i, j = i[ok], j[ok]
        surface = interfaces[ok, n, 1:] / area[ok, None]
        edges.append((i, j, area[ok], np.linalg.norm(xyz[i] - surface, axis=1),
                      np.linalg.norm(xyz[j] - surface, axis=1)))
    assert lost < 1e-10, lost
    stats = dict(interface_mm2={f"{p}-{r}": float(interfaces[:, k, 0].sum())
                               for k, (p, r) in enumerate(PAIRS)})
    return phase, volume, xyz, tuple(map(np.concatenate, zip(*edges))), stats


def voxel_graph(h, gap=False, phase_shift=(0,.5,0)):
    if phase_shift == (0,.5,0):
        axes = [(a[1:] + a[:-1]) / 2 for a in bounds(h)]
    else:
        axes = [(np.arange(-1, int(np.ceil(d/h))+1) + .5 + s)*h
                for d,s in zip(2*CENTRE, phase_shift)]
    points = centre_grid(axes)
    labels = classify(points)
    phantom = type(NATIVE)(labels, h/1000, NATIVE.spec, NATIVE.pillow_exponent)
    active = labels >= 0
    xyz = points[active]
    phase = labels[active]
    # Refined voxel runs omit the centerline gap state.
    i, j, links = face_geometry(phantom, axes)
    dl = np.full(len(i), h/2)
    dr = dl.copy()
    volume = np.full(len(phase), h**3)
    if gap:
        keep = np.ones(len(i), bool)
        extra = []
        for row, heart, wall, a, b, midpoint in links:
            if b-a <= 1e-9:
                dl[row], dr[row] = (b*h, (1-b)*h) if i[row] == heart else ((1-b)*h,b*h)
                continue
            k = len(phase)+len(extra)
            keep[row] = False
            volume[heart] -= max(0,min(b,.5)-a)*h**3/3
            volume[wall] -= max(0,b-max(a,.5))*h**3/3
            extra.append((heart, wall, k, a, b, midpoint))
        for heart, wall, k, a, b, midpoint in extra:
            i,j = np.r_[i,heart,k],np.r_[j,k,wall]
            dl,dr = np.r_[dl,a*h,(b-a)*h/2],np.r_[dr,(b-a)*h/2,(1-b)*h]
            xyz = np.vstack([xyz, midpoint])
            volume = np.r_[volume,(b-a)*h**3/3]
        keep = np.r_[keep,np.ones(len(extra)*2,bool)]
        i,j,dl,dr = i[keep],j[keep],dl[keep],dr[keep]
        phase = np.r_[phase,np.ones(len(extra),int)]
        assert np.isclose(volume.sum(), active.sum()*h**3)
    return phase, volume, xyz, (i, j, np.full(len(i), h**2), dl, dr), {}


def loads(phase, xyz, vector):
    radii = np.array(NATIVE.spec.heart_dims_mm)/2
    radial = np.linalg.norm((xyz-ORIGIN)/radii,axis=1)
    weights = np.exp(-.5*((radial[:,None]-[1.05,1.45,2.])/[.32,.42,.65])**2)
    shell = (weights @ vector[:3])/weights.sum(axis=1)
    return np.select([phase==1,phase==2,phase==3],[shell,1.47,10.])


def startup(graph, vector, dt, implicit=False, theta=1.):
    phase, volume, xyz, (i,j,area,dl,dr), _ = graph
    wall = phase == 0
    specimen = ~wall
    t = np.full(len(phase), -150.)
    load = loads(phase, xyz, vector)
    heat = volume * load * 1000
    _, k0 = material(t, wall)
    area = area * 1e6  # mm conversion, with volumes stored in mm³
    g0 = area / (dl/k0[i] + dr/k0[j])
    if implicit:
        lap = coo_matrix((np.r_[-g0,-g0,g0,g0],
                          (np.r_[i,j,i,j],np.r_[i,j,j,i])),
                         shape=(len(t),len(t))).tocsc()
        cap0, _ = material(t, wall)
        mass0 = cap0*volume
        stepper = splu(diags(mass0,format="csc")-theta*dt*lap,
                       permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.)

    def flux_sum(flux):
        return (np.bincount(i,flux,minlength=len(t))-
                np.bincount(j,flux,minlength=len(t)))

    minimum = float('inf')
    imbalance = 0.
    for step in range(round(1/dt)):
        cap,k = material(t, wall)
        mass = cap*volume
        duty = np.interp(step*dt,[0,35,70],vector[3:])
        power = heat*(691-4.6*(t+80))*duty
        g = area / (dl/k[i]+dr/k[j])
        rhs = power+flux_sum(g*(t[j]-t[i]))
        if implicit:
            rate = stepper.solve(rhs)
            # Five fixed-point corrections keep conductivity and heat capacity
            # at the current temperature while reusing the reference factor.
            for _ in range(5):
                residual = rhs-mass*rate+theta*dt*flux_sum(g*(rate[j]-rate[i]))
                rate += stepper.solve(residual)
        else:
            rate = rhs/mass
        imbalance = max(imbalance, abs((mass*rate).sum()-power.sum())/power.sum())
        limit = np.argmin(np.where(specimen,rate,np.inf))
        value = float(60*rate[limit])
        if value < minimum:
            minimum = value
            result = dict(rate=value, time_s=step*dt, temp_c=float(t[limit]),
                          phase=int(phase[limit]), xyz_mm=xyz[limit].tolist(),
                          volume_mm3=float(volume[limit]))
        t += dt*rate
    assert np.isfinite(t).all() and imbalance < 2e-7, imbalance
    return result


CASES = {
    "voxel_1.1": (1.1, 0), "gap_1.1": (1.1, -1),
    "voxel_0.75": (.75, 0), "voxel_0.55": (.55, 0), "voxel_0.375": (.375, 0),
    "cut_1.1_q8": (1.1, 8), "cut_1.1_q16": (1.1, 16),
    "cut_1.1_q32": (1.1, 32), "cut_0.75_q16": (.75, 16),
}


def compare(names=CASES):
    output = {}
    for name in names:
        h,q = CASES[name]
        mesh = cut_graph(h,q) if q > 0 else voxel_graph(h,q == -1)
        phase,volume,_,_,stats = mesh
        record = {**stats, "phase_volumes_mm3": np.bincount(phase,volume).tolist(),
                  "states":len(phase), "designs":{}}
        for design,vector in VECTORS.items():
            record["designs"][design] = {str(dt): startup(mesh,vector,dt,q>0)
                                         for dt in (.002,.001)}
        output[name] = record
        print(name, {k: v["0.001"]["rate"] for k,v in record["designs"].items()},flush=True)
    return output


if __name__ == "__main__":
    import sys
    names = sys.argv[1:] or CASES
    records = compare(names)
    path = ROOT / "results.json"
    saved = json.loads(path.read_text()) if path.exists() else {}
    path.write_text(json.dumps({**saved, **records},indent=2,allow_nan=False)+"\n")
