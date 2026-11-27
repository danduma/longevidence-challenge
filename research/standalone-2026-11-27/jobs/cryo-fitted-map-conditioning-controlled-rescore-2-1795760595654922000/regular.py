"""Seven conforming blocks on the same analytic superellipsoid."""
from itertools import product

import numpy as np

from fitted import CENTRE, CORNERS, Mesh, POWERS, RADII, classify

INNER = .45
P = POWERS[1]
R = RADII[1]
FACES = list(product(range(3), (-1, 1)))


def map_block(block, q):
    """q is in [0,1]^3, radial coordinate last on the shell blocks."""
    if block == 0:
        xyz = CENTRE + R*INNER*(2*q-1)
        return xyz, np.broadcast_to(np.diag(2*INNER*R), (*q.shape[:-1], 3, 3))
    axis, sign = FACES[block-1]
    tang = [i for i in range(3) if i != axis]
    v = np.empty_like(q)
    v[..., axis] = sign
    v[..., tang] = 2*q[..., :2]-1
    norm = np.sum(abs(v)**P, axis=-1)**(1/P)
    radius = 1/norm
    t = q[..., 2]
    s = INNER+t*(radius-INNER)
    xyz = CENTRE+R*s[..., None]*v
    dr = -2*np.sign(v[..., tang])*abs(v[..., tang])**(P-1)/norm[..., None]**(P+1)
    jac = np.empty((*q.shape[:-1], 3, 3))
    jac[..., :, :2] = 2*s[..., None, None]*np.eye(3)[:, tang]
    jac[..., :, :2] += v[..., :, None]*t[..., None, None]*dr[..., None, :]
    jac[..., :, 2] = v*(radius-INNER)[..., None]
    return xyz, R[:, None]*jac


class RegularMesh(Mesh):
    def __init__(self, n):
        self.shapes = [(n, n, n)] + [(n, n, n//2)]*6
        self.offset = np.r_[0, np.cumsum([np.prod(s) for s in self.shapes])]
        self.cell = np.concatenate([np.indices(s).reshape(3,-1).T for s in self.shapes])
        self.width = np.repeat(1/np.asarray(self.shapes), np.diff(self.offset), axis=0)
        self.origin = self.cell*self.width
        xyz, _, _ = self.quadrature(np.array([[.5]*3]))
        self.label = classify(xyz[:, 0])
        assert np.all(self.label > 0)
        vertices, _, _ = self.quadrature(CORNERS)
        flat = vertices.reshape(-1, 3)
        nodes, first, ids = np.unique(flat.round(11), axis=0, return_index=True, return_inverse=True)
        assert np.max(abs(flat-flat[first][ids])) < 1e-10
        self.dofs = (3*ids.reshape(-1,8,1)+np.arange(3)).reshape(-1,24)
        self.ndof = 3*len(nodes)
        self.free = np.ones(self.ndof, dtype=bool)
        for signs, comps in [((-1,-1,-1),(0,1,2)), ((1,-1,-1),(1,2)), ((-1,1,-1),(2,))]:
            point = CENTRE + R*np.asarray(signs)*3**(-1/P)
            node = np.linalg.norm(nodes-point, axis=1).argmin()
            assert np.linalg.norm(nodes[node]-point) < 1e-10
            self.free[3*node+np.array(comps)] = False

    def quadrature(self, points):
        q = self.origin[:,None]+points*self.width[:,None]
        xyz, jac = self.map_cells(q)
        det = abs(np.linalg.det(jac))
        assert np.all(det > 0)
        return xyz, jac, det*self.width.prod(axis=1)[:,None]

    def map_cells(self, q):
        blocks = [map_block(i, q[lo:hi])
                  for i, (lo, hi) in enumerate(zip(self.offset, self.offset[1:]))]
        return tuple(np.concatenate(parts) for parts in zip(*blocks))

    def strain(self, point):
        _, jac, volume = self.quadrature(point[None])
        grad = gradients(point[None], jac[:,0]*self.width[:,None,:])
        trace = grad.reshape(-1, 24)
        gradient = np.einsum('ena,bc->ebanc',grad,np.eye(3)).reshape(-1,3,3,24)
        sym = (gradient+gradient.swapaxes(1,2))/2
        dev = sym-np.eye(3)[None,:,:,None]*trace[:,None,None,:]/3
        return dev.reshape(-1,9,24), trace, volume[:,0]

    def locate(self, xyz, differentiate=None):
        if differentiate is None:
            differentiate = gradients
        y = (xyz-CENTRE)/R
        s = np.max(abs(y),axis=1)
        dom = abs(y).argmax(axis=1)
        block_id = np.where(s <= INNER, 0, 1+2*dom+(y[np.arange(len(y)),dom]>0))
        owner = np.empty(len(y), int)
        grad = np.empty((len(y),8,3))
        for block, shape in enumerate(self.shapes):
            use = block_id == block
            if block == 0:
                q = (y[use]/INNER+1)/2
            else:
                axis, _ = FACES[block-1]
                tang = [i for i in range(3) if i != axis]
                v = y[use]/s[use,None]
                radius = np.sum(abs(v)**P,axis=1)**(-1/P)
                q = np.column_stack(((v[:,tang]+1)/2,(s[use]-INNER)/(radius-INNER)))
            x, jac = map_block(block,q)
            assert np.max(abs(x-xyz[use]), initial=0) < 1e-10
            grid = q*shape
            cell = np.clip(np.floor(grid).astype(int),0,np.asarray(shape)-1)
            local = grid-cell
            assert np.min(local,initial=1) > 1e-10 and np.max(local,initial=0) < 1-1e-10
            owner[use] = self.offset[block]+np.ravel_multi_index(cell.T,shape)
            grad[use] = differentiate(local,jac/np.asarray(shape))
        return owner, grad


def gradients(local, jac):
    factors = np.where(CORNERS, local[:,None], 1-local[:,None])
    reference = np.stack([
        (2*CORNERS[:, a]-1)*np.prod(np.delete(factors,a,axis=2),axis=2)
        for a in range(3)
    ],axis=-1)
    return reference@np.linalg.inv(jac)
