import json
import math
import sys

import numpy as np

from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "sources/upstream/src"))

def read(name):
    return json.loads((ROOT / "sources/upstream/results" / name).read_text())
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


