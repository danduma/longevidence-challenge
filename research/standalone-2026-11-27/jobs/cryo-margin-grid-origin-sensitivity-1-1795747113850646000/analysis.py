"""Prespecified translation of the sampling lattice, with no fitting across phases."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

from concurrent.futures import ProcessPoolExecutor
import json
import math
import sys

import numpy as np

from forward import ROOT, Mechanics, read, summarize, thermal
from vitrify import geometry

CONFIG = json.loads((ROOT / "protocol.json").read_text())
GRADED = read("grid_convergence_margin.json")["design_parameters"]
CANDIDATES = json.loads((ROOT / "candidates.json").read_text())
VECTORS = {"graded": GRADED}
for name, (loading, *duty) in CANDIDATES.items():
    VECTORS[name] = [10 * loading] * 3 + duty


def sample(spacing, native, phase):
    """Shift queries to source predicates and freeze the native pillow exponent."""
    shift = np.asarray(phase) * spacing
    # Use one guard cell on both sides to detect support beyond the native array.
    coordinates = [(np.arange(-1, n + 1) + .5) * spacing + delta
                   for n, delta in zip(native.shape, shift)]
    grids = np.meshgrid(*coordinates, indexing="ij", sparse=True)
    spec = native.spec
    centre = np.array(spec.container_dims_mm) / 2
    exponent = native.pillow_exponent
    mask = geometry._superellipsoid
    outer = mask(grids, centre, centre, exponent)
    inner = mask(grids, centre, centre - spec.wall_thickness_mm, exponent)
    heart_axes = np.array(spec.heart_dims_mm) / 2
    heart_axes *= np.cbrt(spec.heart_volume_ml * 1000 / (4*np.pi/3*np.prod(heart_axes)))
    heart = mask(grids, centre, heart_axes, 2) & inner
    cav_axes = heart_axes * np.cbrt(spec.chamber_fraction / 2)
    separation = np.array([.4 * heart_axes[0], 0., 0.])
    cavity = heart & (mask(grids, centre - separation, cav_axes, 2)
                      | mask(grids, centre + separation, cav_axes, 2))
    labels = np.select([cavity, heart, inner, outer], [3, 2, 1, 0], -1).astype(np.int8)
    occupied = np.argwhere(outer)
    start = np.minimum(1, occupied.min(axis=0))
    end = np.maximum(np.asarray(native.shape)+1, occupied.max(axis=0)+1)
    slices = tuple(slice(a, b) for a, b in zip(start, end))
    phantom = geometry.Phantom(labels[slices], native.spacing_m, spec, exponent)
    points = tuple(coords[sl] for coords, sl in zip(coordinates, slices))
    return phantom, points


def loading_field(phantom, points, native_shape, vector):
    # The source loading profile is centred on the native array, which may
    # differ slightly from the analytic bag centre. Freeze that physical origin.
    spacing = phantom.spacing_m * 1000
    axes = np.asarray(phantom.spec.heart_dims_mm) / 2
    centred = [(a - n*spacing/2)/r for a,n,r in zip(points, native_shape, axes)]
    radius = np.sqrt(sum(x*x for x in np.meshgrid(*centred, indexing="ij", sparse=True)))
    centres = np.array([1.05, 1.45, 2.0])[:,None,None,None]
    widths = np.array([.32, .42, .65])[:,None,None,None]
    gaussian = np.exp(-.5*((radius-centres)/widths)**2)
    shell = np.tensordot(vector[:3], gaussian/gaussian.sum(axis=0), axes=1)
    labels = phantom.labels
    return np.select([labels == 1, labels == 2, labels == 3], [shell, 1.47, 10.])


def metrics(phantom, points, native_shape, vector, mechanics):
    field = loading_field(phantom, points, native_shape, vector)
    dt, n, states, rates, _, _, events = thermal(phantom, field, vector[3:])
    active = phantom.active_mask
    labels = phantom.labels[active]
    specimen = labels > 0
    step = int(rates[:n].argmin())
    cell = int(np.where(specimen, states[step+1]-states[step], np.inf).argmin())
    xyz = np.argwhere(active)[cell]
    reconstructed = (states[step+1,cell]-states[step,cell])*60/dt
    assert abs(reconstructed-rates[step]) < 1e-8
    stop, temp = events["heart"]
    result = dict(rate_70=float(rates[step]), rate_event=float(rates[:math.ceil(stop/dt)].min()),
                  time_s=step*dt, location_mm=[float(p[i]) for p,i in zip(points, xyz)],
                  limiting_label=int(labels[cell]), limiting_temp_c=float(states[step,cell]),
                  heart_stop_s=float(stop), hot_at_event_c=float(temp[specimen].max()),
                  cold_70=float(states[n,specimen].min()), hot_70=float(states[n,specimen].max()))
    if mechanics is not None:
        early_stride = round(CONFIG["early_cadence_s"]/dt)
        late_stride = round(CONFIG["late_cadence_s"]/dt)
        early = np.arange(early_stride, math.floor(CONFIG["cadence_change_s"]/dt)+1, early_stride)
        late = np.arange(early[-1]+late_stride, math.ceil(stop/dt), late_stride)
        indices = np.r_[early, late]
        times = np.r_[indices*dt, stop]
        stress = mechanics.solve(times, np.vstack((states[indices], temp)))
        result.update(summarize(times, stress))
    return result


def evaluate(job):
    spacing, phase, solve_mechanics = job
    native = geometry.rat_heart_cryobag(spacing)
    phantom, points = sample(spacing, native, phase)
    if tuple(phase) == (0, 0, 0):
        np.testing.assert_array_equal(phantom.labels, native.labels)
    mechanics = Mechanics(phantom) if solve_mechanics else None
    return dict(spacing_mm=spacing, phase=phase, voxel_counts=np.bincount(
                phantom.labels.ravel()+1, minlength=5)[1:].tolist(),
                pillow_exponent=native.pillow_exponent,
                designs={name: metrics(phantom, points, native.shape, VECTORS[name], mechanics)
                         for name in CONFIG["schedules"]})


def recompute(solve_mechanics=True):
    jobs = [(spacing, phase, solve_mechanics) for spacing in CONFIG["spacings_mm"]
            for phase in CONFIG["phase_fractions"]]
    workers = int(os.environ.get("CRYO_WORKERS", 2))
    if workers == 1:
        return list(map(evaluate, jobs))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(evaluate, jobs))


if __name__ == "__main__":
    thermal_only = "--thermal" in sys.argv
    print(json.dumps(recompute(not thermal_only), indent=2, allow_nan=False))
