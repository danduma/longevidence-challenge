"""Nested curved Q1 specimen mechanics on the predecessor's thermal paths."""
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree

from cut_geometry import CENTRE, POWERS, RADII, bounds, classify
from heat_cache import DT, H, compute
from mechanics import summarize

ROOT = Path(__file__).resolve().parent
CORNERS = np.array(list(product((0, 1), repeat=3)))
SHAPES = ((8, 6, 4), (16, 12, 8))
CONFIGURATIONS = tuple((shape, q) for shape in SHAPES for q in (12, 20))


def physical(u):
    """Spherified cube followed by a radial map to the analytic superellipsoid."""
    x, y, z = np.moveaxis(u, -1, 0)
    a, b, c = x*x, y*y, z*z
    sphere = u * np.sqrt(np.stack((1-(b+c)/2+b*c/3,
                                  1-(a+c)/2+a*c/3,
                                  1-(a+b)/2+a*b/3), axis=-1))
    l2 = np.linalg.norm(sphere, axis=-1)
    lp = np.sum(abs(sphere)**POWERS[1], axis=-1)**(1/POWERS[1])
    ratio = np.divide(l2, lp, out=np.ones_like(l2), where=lp > 0)
    return CENTRE + RADII[1] * sphere * ratio[..., None]


class Mesh:
    def __init__(self, shape):
        self.shape = np.asarray(shape)
        self.cell = np.indices(shape).reshape(3, -1).T
        self.width = 2/self.shape
        self.origin = self.cell*self.width - 1
        self.label = classify(physical(self.origin+self.width/2))
        assert np.all(self.label > 0)
        nodes = np.ravel_multi_index((self.cell[:, None]+CORNERS).transpose(2, 0, 1),
                                     tuple(self.shape+1))
        self.dofs = (3*nodes[..., None]+np.arange(3)).reshape(-1, 24)
        self.ndof = 3*np.prod(self.shape+1)
        self.free = np.ones(self.ndof, dtype=bool)
        for point, comps in [((0,0,0), (0,1,2)), ((shape[0],0,0),(1,2)),
                             ((0,shape[1],0),(2,))]:
            node = np.ravel_multi_index(point, tuple(self.shape+1))
            self.free[3*node+np.array(comps)] = False

    def quadrature(self, points):
        u = self.origin[:, None] + points*self.width
        xyz = physical(u)
        jac = np.stack([(physical(u+d)-physical(u-d))/(2e-5)
                        for d in np.eye(3)*1e-5], axis=-1)
        det = np.linalg.det(jac)
        assert np.all(det > 0)
        return xyz, jac, det*np.prod(self.width)

    def strain(self, point):
        _, jac, volume = self.quadrature(point[None])
        factors = np.where(CORNERS, point, 1-point)
        reference = np.column_stack([
            (2*CORNERS[:, a]-1)*np.prod(np.delete(factors, a, axis=1), axis=1)
            / self.width[a] for a in range(3)
        ])
        grad = reference@np.linalg.inv(jac[:, 0])
        trace = grad.reshape(-1, 24)
        gradient = np.einsum('ena,bc->ebanc', grad, np.eye(3)).reshape(-1,3,3,24)
        sym = (gradient+gradient.swapaxes(1,2))/2
        dev = sym - np.eye(3)[None,:,:,None]*trace[:,None,None,:]/3
        return dev.reshape(-1,9,24), trace, volume[:,0]


def overlap(mesh, graph, q):
    """Composite quadrature of same-phase donor/curved-target intersections."""
    phase, _, xyz, _, _ = graph
    donor = bounds(H)
    cell = tuple(np.searchsorted(a, x)-1 for a,x in zip(donor, xyz.T))
    uid = np.full(tuple(len(a)-1 for a in donor)+(4,), -1)
    uid[(*cell, phase)] = np.arange(len(phase))
    neighbors = np.empty((len(phase),16), int)
    gradient = np.empty((len(phase),3,16))
    trees = {}
    for p in (1,2,3):
        group = np.flatnonzero(phase == p)
        tree = cKDTree(xyz[group])
        trees[p] = (tree, group)
        dist, ids = tree.query(xyz[group], k=16)
        neighbors[group] = group[ids]
        delta = xyz[group[ids]]-xyz[group,None]
        scale = 1/np.maximum(dist,.2)
        basis = np.concatenate([np.ones((*delta.shape[:2],1)),delta], axis=-1)
        gradient[group] = np.linalg.pinv(basis*scale[...,None])[:,1:,:]*scale[:,None,:]

    points = np.array(list(product((np.arange(q)+.5)/q, repeat=3)))
    shape = (len(mesh.cell), len(phase))
    integrals = [csr_matrix(shape) for _ in range(4)]
    missing = total = 0.
    for block in np.array_split(points, max(1, len(points)//48)):
        coords, _, measure = mesh.quadrature(block)
        r, k = np.where(classify(coords) == mesh.label[:,None])
        pos = coords[r,k]
        weight = measure[r,k]/len(points)
        ijk = tuple(np.searchsorted(a,x)-1 for a,x in zip(donor,pos.T))
        col = uid[(*ijk,mesh.label[r])]
        total += weight.sum()
        missing += weight[col<0].sum()
        for p, (tree, group) in trees.items():
            mask = (col < 0) & (mesh.label[r] == p)
            col[mask] = group[tree.query(pos[mask])[1]]
        values = np.column_stack((weight, weight[:,None]*pos))
        for d, value in enumerate(values.T):
            integrals[d] += coo_matrix((value,(r,col)),shape=shape).tocsr()
    volume = np.asarray(integrals[0].sum(axis=1)).ravel()
    assert np.all(volume > 0)
    average = integrals[0].multiply(1/volume[:,None]).tocsr()
    matrix = average.copy()
    centroid = []
    for d in range(3):
        moment = integrals[d+1].multiply(1/volume[:,None]).tocsr()
        centroid.append(np.asarray(moment.sum(axis=1)).ravel())
        offset = (moment-average.multiply(xyz[:,d])).tocoo()
        data = offset.data[:,None]*gradient[offset.col,d]
        matrix += coo_matrix((data.ravel(), (np.repeat(offset.row,16),
                                             neighbors[offset.col].ravel())), shape=shape).tocsr()
    r,c = matrix.nonzero()
    assert np.array_equal(mesh.label[r],phase[c])
    source = np.column_stack((10+phase, xyz+phase[:,None]))
    target = np.column_stack((10+mesh.label,np.array(centroid).T+mesh.label[:,None]))
    error = float(np.max(abs(matrix@source-target)))
    assert error < 1e-9
    return matrix, dict(patch_error=error, missing_fraction=float(missing/total))


class Solver:
    def __init__(self, mesh):
        self.mesh = mesh
        self.bulk = 800e6/1.5
        count = len(mesh.cell)
        self.kdev = np.zeros((count,24,24))
        self.kbulk = np.zeros_like(self.kdev)
        self.int_dev = np.zeros((count,9,24))
        self.int_trace = np.zeros((count,24))
        self.volume = np.zeros(count)
        gauss, weights = np.polynomial.legendre.leggauss(3)
        for index in product(range(3),repeat=3):
            index = list(index)
            dev, trace, volume = mesh.strain((gauss[index]+1)/2)
            volume *= weights[index].prod()/8
            self.volume += volume
            self.kdev += 2*(dev.transpose(0,2,1)@dev)*volume[:,None,None]
            self.kbulk += self.bulk*trace[:,:,None]*trace[:,None,:]*volume[:,None,None]
            self.int_dev += dev*volume[:,None,None]
            self.int_trace += trace*volume[:,None]
        self.dev, self.trace, _ = mesh.strain(np.full(3,.5))
        row = np.repeat(mesh.dofs,24,axis=1).ravel()
        col = np.tile(mesh.dofs,(1,24)).ravel()
        self.keep = mesh.free[row]&mesh.free[col]
        idx = np.cumsum(mesh.free)-1
        self.row, self.col = idx[row[self.keep]],idx[col[self.keep]]
        n = mesh.free.sum()
        ids = self.row*n+self.col
        unique, self.scatter = np.unique(ids, return_inverse=True)
        self.row,self.col = np.divmod(unique,n)
        self.size = (n,n)

    def solve(self, times, field, observe=None):
        mesh = self.mesh
        prev_t, prev_time = -150., 0.
        center_s = np.zeros((len(mesh.cell),9))
        center_e = np.zeros_like(center_s)
        principal = []
        disp = np.zeros(mesh.ndof)
        last_shear = None
        for time,t in zip(times,field,strict=True):
            mid = (t+prev_t)/2
            eta = np.select([mid >= -100,mid >= -140],
                            [1.21e4,4.2783e-23*np.exp(-.6091*mid)],4.63e14)
            ratio = (time-prev_time)*320e6/eta
            decay = np.exp(-ratio)
            shear = 320e6*(-np.expm1(-ratio))/ratio
            history = decay[:,None]*center_s-2*shear[:,None]*center_e
            # Represent free uniform dilation analytically on the curved domain.
            # Physical temperature still sets the local viscosity.
            eigenstrain = 3*1.1e-4*(t-np.average(t, weights=self.volume))
            load = self.bulk*eigenstrain[:,None]*self.int_trace
            load -= np.einsum('ci,cij->cj', history, self.int_dev)
            rhs = np.bincount(mesh.dofs.ravel(),load.ravel(),minlength=mesh.ndof)
            if last_shear is None or not np.allclose(shear,last_shear,rtol=1e-13,atol=0):
                k = self.kbulk + shear[:,None,None]*self.kdev
                values = np.bincount(self.scatter,k.ravel()[self.keep])
                matrix = coo_matrix((values,(self.row,self.col)),shape=self.size).tocsc()
                factor = splu(matrix, permc_spec='MMD_AT_PLUS_A',
                              diag_pivot_thresh=0, options={'SymmetricMode': True})
                last_shear = shear
            disp[mesh.free] = factor.solve(rhs[mesh.free])
            u = disp[mesh.dofs]
            center_now = np.einsum('cij,cj->ci',self.dev,u)
            center_s = history+2*shear[:,None]*center_now
            center_e = center_now
            pressure = self.bulk*(np.einsum('ci,ci->c',self.trace,u)-eigenstrain)
            stress = center_s.reshape(-1,3,3)+pressure[:,None,None]*np.eye(3)
            principal.append(np.linalg.eigvalsh(stress))
            if observe is not None:
                observe(time, u, history.reshape(-1, 3, 3), shear, eigenstrain)
            prev_time,prev_t = time,t
        return np.array(principal)


def run(shape, graph, histories, events, q):
    mesh = Mesh(shape)
    solver = Solver(mesh)
    matrix, diagnostics = overlap(mesh,graph,q)
    records = {}
    # Each schedule uses the predecessor's identical regular and event samples.
    for name, history in histories.items():
        stride = round(.3/DT)
        regular = np.arange(stride,len(history)-1,stride)
        times = np.r_[regular*DT,events[name]['time_s']]
        field = (matrix@history[np.r_[regular,len(history)-1]].T).T
        principal = solver.solve(times, field)
        records[name] = summarize(times, principal)
        early = times <= 12
        records[name]['release_12s'] = summarize(times[early], principal[early])
        print(shape,q,name,records[name],flush=True)
    # The analytic dilation must remove spurious stress from uniform warming.
    uniform = solver.solve([.3],np.full((1,len(mesh.cell)),-149.))
    return dict(shape=list(shape), q=q, volume_mm3=float(solver.volume.sum()),
                uniform_patch_kpa=float(abs(uniform).max()/1e3),
                overlap=diagnostics, schedules=records)


def regenerate(inputs=None):
    graph,histories,events = compute() if inputs is None else inputs
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=len(CONFIGURATIONS)) as pool:
        jobs = [pool.submit(run,shape,graph,histories,events,q)
                for shape,q in CONFIGURATIONS]
        rows = [job.result() for job in jobs]
    return dict(events=events,runs=rows)


if __name__ == '__main__':
    import sys
    inputs = None
    if len(sys.argv) == 3:
        from heat_cache import Q
        from thermal_graph import cut_graph
        with np.load(sys.argv[1], allow_pickle=False) as cache:
            histories = dict(cache)
        inputs = cut_graph(H,Q), histories, json.loads(Path(sys.argv[2]).read_text())
    (ROOT/'results.json').write_text(json.dumps(regenerate(inputs),indent=2,allow_nan=False)+'\n')
