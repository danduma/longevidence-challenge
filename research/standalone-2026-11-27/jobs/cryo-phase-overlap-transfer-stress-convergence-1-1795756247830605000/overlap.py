"""Integrate a same-phase donor reconstruction over each target element."""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree

from cut_geometry import bounds, CENTRES, POWERS, RADII, phases


def common_axis(donor, target, q):
    edges = np.unique(np.round(np.r_[donor, target], 12))
    lower, upper = edges[:-1], edges[1:]
    mask = (lower >= target[0] - 1e-10) & (upper <= target[-1] + 1e-10)
    lower, upper = lower[mask], upper[mask]
    points = (lower[:, None] + (upper-lower)[:, None] *
              ((np.arange(q)+.5)/q)).ravel()
    return points, np.repeat((upper-lower)/q, q)


def moments(donor, target, donor_uid, labels, q, axis):
    """Analytic segments along one axis, composite midpoint across the other two."""
    dims = [d for d in range(3) if d != axis]
    transverse = [common_axis(donor[d], target[d], q) for d in dims]
    u, v = np.meshgrid(*(p for p, _ in transverse), indexing='ij')
    points = np.column_stack((u.ravel(), v.ravel()))
    areas = np.multiply.outer(*(w for _, w in transverse)).ravel()
    cuts = np.unique(np.round(np.r_[donor[axis], target[axis]], 12))
    source_shape = donor_uid.shape[:3]
    target_shape = labels.shape
    rows, cols, values = [], [], []
    for start in range(0, len(points), 256):
        sl = slice(start, start+256)
        uv = points[sl]
        intervals = []
        for c, r, power in zip(CENTRES, RADII, POWERS):
            rest = 1 - np.sum(abs((uv-c[dims])/r[dims])**power, axis=1)
            half = r[axis]*np.maximum(rest, 0)**(1/power)
            intervals.extend((c[axis]-half, c[axis]+half))
        edges = np.sort(np.column_stack([np.broadcast_to(cuts, (len(uv), len(cuts))),
                                         *intervals]), axis=1)
        mid = (edges[:, 1:] + edges[:, :-1])/2
        phase = phases([(lo[:, None] < mid) & (mid < hi[:, None])
                        for lo, hi in zip(intervals[::2], intervals[1::2])])
        i, j = np.where((phase >= 0) & (np.diff(edges, axis=1) > 1e-12))
        xyz = np.empty((len(i), 3))
        xyz[:, axis] = mid[i, j]
        xyz[:, dims] = uv[i]
        src = np.column_stack([np.searchsorted(a, x)-1 for a,x in zip(donor, xyz.T)])
        dst = np.column_stack([np.searchsorted(a, x)-1 for a,x in zip(target, xyz.T)])
        inside = ((src >= 0) & (src < source_shape) &
                  (dst >= 0) & (dst < target_shape)).all(axis=1)
        src, dst, xyz, i, j = (a[inside] for a in (src,dst,xyz,i,j))
        p = phase[i,j]
        same = labels[tuple(dst.T)] == p
        src, dst, xyz, i, j, p = (a[same] for a in (src,dst,xyz,i,j,p))
        rows.append(np.ravel_multi_index(dst.T, target_shape))
        cols.append(donor_uid[(*src.T, p)])
        volume = (edges[i,j+1]-edges[i,j])*areas[sl][i]
        values.append(volume[:,None]*np.column_stack((np.ones(len(i)), xyz)))
    return np.concatenate(rows), np.concatenate(cols), np.concatenate(values)


def transfers(graph, thermal_h, target, labels, q, axis=0):
    phase, _, centroids, _, _ = graph
    donor = bounds(thermal_h)
    source_shape = tuple(len(a)-1 for a in donor)
    cell = np.column_stack([np.searchsorted(a, x)-1 for a,x in zip(donor,centroids.T)])
    uid = np.full((*source_shape, 4), -1)
    uid[(*cell.T, phase)] = np.arange(len(phase))
    assert np.count_nonzero(uid >= 0) == len(phase)

    active = labels.ravel() >= 0
    target_index = np.cumsum(active)-1
    centers = np.stack(np.meshgrid(*[(a[1:]+a[:-1])/2 for a in target],
                                   indexing='ij'), axis=-1)[labels >= 0]
    row, col, integrals = moments(donor, target, uid, labels, q, axis)
    row = target_index[row]
    volume = np.bincount(row, integrals[:,0], minlength=len(centers))
    assert (volume > 0).all()
    target_centroid = np.column_stack([
        np.bincount(row, v, minlength=len(centers))/volume for v in integrals[:,1:].T])

    nearest = np.empty(len(centers), int)
    donor_neighbors = np.empty((len(phase),16), int)
    gradient = np.empty((len(phase),3,16))
    idw_neighbors = np.empty((len(centers),8), int)
    idw_weights = np.empty_like(idw_neighbors, float)
    for p in range(4):
        group = np.flatnonzero(phase == p)
        tree = cKDTree(centroids[group])
        mask = labels.ravel()[active] == p
        distance, idx = tree.query(centers[mask], k=8)
        idw_neighbors[mask] = group[idx]
        nearest[mask] = group[idx[:,0]]
        weights = 1/np.maximum(distance, 1e-8)**2
        idw_weights[mask] = weights/weights.sum(axis=1, keepdims=True)
        distance, idx = tree.query(centroids[group], k=16)
        donor_neighbors[group] = group[idx]
        offset = centroids[group[idx]] - centroids[group,None]
        scale = 1/np.maximum(distance, .2)
        basis = np.concatenate([np.ones((*offset.shape[:2],1)), offset], axis=2)
        gradient[group] = np.linalg.pinv(basis * scale[:,:,None])[:,1:,:]*scale[:,None,:]

    def sparse(r, c, w):
        return coo_matrix((w.ravel(), (r.ravel(), c.ravel())),
                          shape=(len(centers),len(phase))).tocsr()

    missing = col < 0
    # A sub-quadrature sliver can be absent from the thermal graph. Use the
    # element's nearest same-phase donor for that sliver and record its weight.
    missing_fraction = np.bincount(row[missing], integrals[missing,0],
                                   minlength=len(centers))/volume
    col[missing] = nearest[row[missing]]
    weight = integrals[:,0]/volume[row]
    constant = sparse(row, col, weight)
    displacement = (integrals[:,1:] - integrals[:,:1]*centroids[col])/volume[row,None]
    # Aggregate first moments before expanding each donor's gradient stencil.
    offset = [sparse(row,col,d).tocoo() for d in displacement.T]
    linear = constant.copy()
    for axis_moment, component in enumerate(offset):
        r, c, data = component.row, component.col, component.data
        linear += sparse(np.broadcast_to(r[:,None], (len(r),16)), donor_neighbors[c],
                         data[:,None]*gradient[c,axis_moment])
    idw = sparse(np.broadcast_to(np.arange(len(centers))[:,None], idw_neighbors.shape),
                 idw_neighbors, idw_weights)
    diagnostics = dict(max_missing_fraction=float(missing_fraction.max()),
                       total_missing_fraction=float((missing_fraction@volume)/volume.sum()),
                       minimum_element_fraction=float(volume.min()/(target[0][1]-target[0][0])**3))
    return {'idw8': idw, 'overlap_constant': constant, 'overlap_affine': linear}, target_centroid, diagnostics


def patch_tests(maps, graph, target_centroid, labels):
    phase, _, xyz, _, _ = graph
    target_phase = labels[labels >= 0]
    constant = 10. + 31*phase
    target_constant = 10. + 31*target_phase
    # Independent basis fields with different jumps across phases.
    slopes = np.array([.11,.19,.23,.31])
    affine = constant[:,None] + slopes[phase,None]*xyz
    expected = target_constant[:,None] + slopes[target_phase,None]*target_centroid
    test = {}
    for name, matrix in maps.items():
        r, c = matrix.nonzero()
        assert np.array_equal(target_phase[r], phase[c])
        test[name] = dict(constant=float(abs(matrix@constant-target_constant).max()),
                          affine=float(abs(matrix@affine-expected).max()))
        assert test[name]['constant'] < 1e-9
    assert test['overlap_affine']['affine'] < 1e-9
    return test
