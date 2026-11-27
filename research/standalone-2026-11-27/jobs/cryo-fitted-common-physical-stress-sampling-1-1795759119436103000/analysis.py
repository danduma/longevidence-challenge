"""Rescore the inherited fitted solves with a shared physical quadrature."""
import json
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from fitted import CONFIGURATIONS, CORNERS, SHAPES, Mesh, Solver, overlap, physical
from heat_cache import DT, compute
from mechanics import summarize

ROOT = Path(__file__).resolve().parent
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


def run(shape, graph, histories, events, transfer_order, sample_type):
    mesh = Mesh(shape)
    solver = Solver(mesh)
    matrix, diagnostic = overlap(mesh, graph, transfer_order)
    samples = {str(order): sample_type(mesh, order) for order in (2, 3)}
    records = {}
    for name, history in histories.items():
        stride = round(.3/DT)
        regular = np.arange(stride, len(history)-1, stride)
        times = np.r_[regular*DT, events[name]['time_s']]
        field = (matrix@history[np.r_[regular, len(history)-1]].T).T
        scores = {order: [] for order in samples}

        def observe(time, u, memory, shear, eigenstrain):
            for order, sample in samples.items():
                principal = sample.evaluate(u, memory, shear, eigenstrain, solver.bulk)
                peak = principal[:, -1].argmax()
                scores[order].append(dict(
                    peak_mpa=float(principal[peak, -1]/1e6),
                    peak_xyz_mm=sample.xyz[peak].tolist(),
                    weighted_logsum=float(logsumexp(principal/SCALE+sample.log_weight[:, None]))
                ))

        centers = solver.solve(times, field, observe=observe)
        early = times <= 12
        records[name] = {
            'native': summarize(times, centers),
            'native_12s': summarize(times[early], centers[early]),
            'shared': {}
        }
        for order, sample in samples.items():
            records[name]['shared'][order] = {
                'points': len(sample.xyz),
                'volume_mm3': float(sample.volume.sum()),
                'weighted_objective_mpa': float(SCALE/1e6*logsumexp(
                    [step['weighted_logsum'] for step in scores[order]])),
                'weighted_objective_12s_mpa': float(SCALE/1e6*logsumexp(
                    [step['weighted_logsum'] for step, use in zip(scores[order], early) if use])),
                'time_history': [dict(time_s=float(t), **step)
                                 for t, step in zip(times, scores[order])]
            }
        print(shape, transfer_order, name, flush=True)
    return dict(shape=list(shape), q=transfer_order, overlap=diagnostic, schedules=records)


def regenerate(sample_type=Sample):
    graph, histories, events = compute()
    with ThreadPoolExecutor(max_workers=len(CONFIGURATIONS)) as pool:
        jobs = [pool.submit(run, shape, graph, histories, events, q, sample_type)
                for shape, q in CONFIGURATIONS]
        rows = [job.result() for job in jobs]
    return dict(events=events, runs=rows)


if __name__ == '__main__':
    result = regenerate()
    (ROOT/'results.json').write_text(json.dumps(result, separators=(',', ':'),
                                              allow_nan=False)+'\n')
