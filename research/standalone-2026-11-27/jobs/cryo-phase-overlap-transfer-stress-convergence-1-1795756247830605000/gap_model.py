"""Matched Vitrify phantom face-flux experiment. See report.md."""
import numpy as np

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

