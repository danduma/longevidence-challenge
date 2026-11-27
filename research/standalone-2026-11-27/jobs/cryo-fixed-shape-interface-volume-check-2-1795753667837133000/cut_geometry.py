"""Ray quadrature of phase volumes, centroids, face apertures and interfaces."""
import numpy as np
from lattice import geometry

NATIVE = geometry.rat_heart_cryobag(1.1)
CENTRE = np.array(NATIVE.spec.container_dims_mm) / 2
ORIGIN = np.array(NATIVE.shape) * 1.1 / 2
HEART = np.array(NATIVE.spec.heart_dims_mm) / 2
HEART *= np.cbrt(690 / (4 * np.pi / 3 * HEART.prod()))
CAVITY = HEART * np.cbrt(NATIVE.spec.chamber_fraction / 2)
POWERS = [NATIVE.pillow_exponent] * 2 + [2.] * 3
RADII = [CENTRE, CENTRE - NATIVE.spec.wall_thickness_mm, HEART, CAVITY, CAVITY]
CENTRES = np.tile(CENTRE, (5, 1))
CENTRES[3:, 0] += np.array([-1, 1]) * .4 * HEART[0]
PAIRS = [(a, b) for a in range(4) for b in range(a + 1, 4)]


def phases(inside):
    outer, inner, heart, c1, c2 = inside
    return np.select([~outer, ~inner, ~heart, c1 | c2], [-1, 0, 1, 3], 2)


def classify(points):
    return phases([np.sum(abs((points - c) / r)**p, axis=-1) <= 1
                   for c, r, p in zip(CENTRES, RADII, POWERS)])


def bounds(h):
    return [np.arange(-1, int(np.ceil(d / h)) + 2) * h + shift
            for d, shift in zip(2 * CENTRE, [0, h / 2, 0])]


def integrate(h, q):
    axes = bounds(h)
    shape = tuple(len(a) - 1 for a in axes)
    cells = np.arange(np.prod(shape)).reshape(shape)
    moments = np.zeros((cells.size, 4, 4))
    interfaces = np.zeros((cells.size, len(PAIRS), 4))
    faces = []
    patch = (h / q)**2
    for axis in range(3):
        dims = [i for i in range(3) if i != axis]
        transverse = [axes[d][0] + (np.arange(shape[d] * q) + .5) * h / q
                      for d in dims]
        uv = np.stack(np.meshgrid(*transverse, indexing='ij'), -1).reshape(-1, 2)
        uv_idx = np.indices(tuple(len(t) for t in transverse)).reshape(2, -1).T // q
        face_area = np.zeros((cells.size, 4))
        ax = axes[axis]
        local_cells = np.moveaxis(cells, axis, -1)
        for start in range(0, len(uv), 512):
            sl = slice(start, start + 512)
            points = uv[sl]
            a, b = uv_idx[sl].T
            ids = local_cells[a, b]
            intervals = []
            for c, r, p in zip(CENTRES, RADII, POWERS):
                rest = 1 - np.sum(abs((points - c[dims]) / r[dims])**p, axis=1)
                half = r[axis] * np.maximum(0, rest)**(1/p)
                intervals.append((c[axis] - half, c[axis] + half))

            def at(x):
                return phases([(lo[:, None] < x) & (x < hi[:, None])
                               for lo, hi in intervals])

            roots = np.sort(np.column_stack([v for pair in intervals for v in pair]), axis=1)
            cuts = np.sort(np.concatenate([np.broadcast_to(ax, (len(points), len(ax))),
                                           roots], axis=1), axis=1)
            length = np.diff(cuts, axis=1)
            middle = (cuts[:, :-1] + cuts[:, 1:]) / 2
            index = np.clip(np.searchsorted(ax, middle) - 1, 0, len(ax) - 2)
            phase = at(middle)
            r0, s0 = np.where((phase >= 0) & (length > 0))
            phi = phase[r0, s0]
            cell = ids[r0, index[r0, s0]]
            weight = length[r0, s0] * patch / 3
            np.add.at(moments[..., 0], (cell, phi), weight)
            np.add.at(moments[..., axis + 1], (cell, phi), weight * middle[r0, s0])
            for d, coords in zip(dims, points.T):
                np.add.at(moments[..., d + 1], (cell, phi), weight * coords[r0])

            phase = at(ax[None, 1:-1])
            rr, ff = np.where(phase >= 0)
            np.add.at(face_area, (ids[rr, ff], phase[rr, ff]), patch)

            for (lo, hi), c, radii, p in zip(intervals, CENTRES, RADII, POWERS):
                root = np.stack([lo, hi], axis=1)
                left, right = at(root - 1e-9), at(root + 1e-9)
                pair_lo, pair_hi = np.minimum(left, right), np.maximum(left, right)
                for number, (p0, p1) in enumerate(PAIRS):
                    rr, kk = np.where((pair_lo == p0) & (pair_hi == p1) & (root[:, :1] < root[:, 1:]))
                    if not len(rr):
                        continue
                    xyz = np.empty((len(rr), 3))
                    xyz[:, dims] = points[rr]
                    xyz[:, axis] = root[rr, kk]
                    grad = abs((xyz - c) / radii)**(p - 1) / radii
                    # A projected patch contributes |n_axis| dA_projected.
                    # Summing all directions integrates sum(n_axis**2) dA.
                    area = patch * grad[:, axis] / np.linalg.norm(grad, axis=1)
                    cell = ids[rr, np.searchsorted(ax, xyz[:, axis]) - 1]
                    # A buried endpoint can coincide with another shape only on a
                    # measure-zero curve. Check that this root changes its own predicate.
                    np.add.at(interfaces[:, number, 0], cell, area)
                    for d in range(3):
                        np.add.at(interfaces[:, number, d + 1], cell, area * xyz[:, d])
        faces.append(face_area)
    return cells, moments, interfaces, faces
