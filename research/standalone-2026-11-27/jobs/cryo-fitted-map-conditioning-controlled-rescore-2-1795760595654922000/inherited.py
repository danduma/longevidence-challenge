"""Common physical stress cloud from the immediately preceding audit."""
from itertools import product

import numpy as np

from fitted import CORNERS, SHAPES, Mesh, physical

SCALE = 20_000.


def cloud(order):
    """Composite Gauss rule on the common refinement of the two partitions."""
    fine = Mesh(SHAPES[-1])
    nodes, weights = np.polynomial.legendre.leggauss(order)
    indices = np.array(list(product(range(order), repeat=3)))
    local = (1 + nodes[indices])/2
    xyz, _, volume = fine.quadrature(local)
    reference = fine.origin[:, None] + local*fine.width
    volume *= weights[indices].prod(axis=1)/8
    return reference.reshape(-1, 3), xyz.reshape(-1, 3), volume.ravel()


class Sample:
    def __init__(self, mesh, order):
        reference, self.xyz, self.volume = cloud(order)
        grid = (reference+1)/mesh.width
        cell = np.floor(grid).astype(int)
        self.owner = np.ravel_multi_index(cell.T, mesh.shape)
        local = grid-cell
        assert np.min(local) > 0 and np.max(local) < 1
        jac = np.stack([(physical(reference+d)-physical(reference-d))/(2e-5)
                        for d in np.eye(3)*1e-5], axis=-1)
        factors = np.where(CORNERS, local[:, None], 1-local[:, None])
        gradient = np.stack([
            (2*CORNERS[:, axis]-1)*np.prod(np.delete(factors, axis, axis=2), axis=2)
            / mesh.width[axis] for axis in range(3)
        ], axis=-1)
        self.gradient = gradient@np.linalg.inv(jac)
        self.log_weight = np.log(self.volume/self.volume.sum())

    def evaluate(self, u, history, shear, eigenstrain, bulk):
        cell = self.owner
        gradient = u[cell].reshape(-1, 8, 3).transpose(0, 2, 1)@self.gradient
        trace = np.trace(gradient, axis1=1, axis2=2)
        isotropic = bulk*(trace-eigenstrain[cell])-2*shear[cell]*trace/3
        stress = history[cell]+shear[cell, None, None]*(gradient+gradient.transpose(0, 2, 1))
        stress += isotropic[:, None, None]*np.eye(3)
        return np.linalg.eigvalsh(stress)


