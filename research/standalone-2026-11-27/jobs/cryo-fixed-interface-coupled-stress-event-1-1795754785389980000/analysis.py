"""Couple matched thermal discretizations to one fixed mechanical grid."""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import cg
from scipy.spatial import cKDTree

from mechanics import Mechanics, rat_heart_cryobag, summarize
from thermal_graph import cut_graph, voxel_graph, loads, material
from lattice import VECTORS

ROOT = Path(__file__).resolve().parent
GLASS = -123.
CUTOFF = -35.
FLOOR = 50.
SPACING = .75


def heat_history(graph, vector, mode, dt):
    phase, volume, xyz, (i, j, area, dl, dr), _ = graph
    heart, specimen, wall = phase >= 2, phase > 0, phase == 0
    area = area * 1e6
    load = 1000 * volume * loads(phase, xyz, vector)
    t = np.full(len(phase), -150.)
    history = [t.copy()]
    crossed = np.zeros(len(t), bool)
    minimum = (np.inf, None)
    max_balance = 0.
    worst_residual = 0.
    rows, cols = np.r_[i,j,i,j], np.r_[i,j,j,i]

    def flux(edge):
        return np.bincount(i, edge, minlength=len(t)) - np.bincount(j, edge, minlength=len(t))

    for step in range(round(110/dt)):
        time = step*dt
        capacity, conductivity = material(t, wall)
        mass = volume * capacity
        g = area / (dl/conductivity[i] + dr/conductivity[j])
        effective_t = np.clip(t, -80., 0.) if mode == 'clamped' else np.minimum(t, 0.)
        power = load * (691 - 4.6*(effective_t+80))
        # Midpoint integration of the prescribed piecewise linear duty.
        power *= np.interp(time+dt/2, [0.,35.,70.], vector[3:])
        rhs = power + flux(g*(t[j]-t[i]))
        lap = coo_matrix((np.r_[g,g,-g,-g], (rows,cols)),
                         shape=(len(t),len(t))).tocsr()
        system = diags(mass,format='csr') + dt*lap
        scale = 1/np.sqrt(system.diagonal())
        diagonal = diags(scale,format='csr')
        normalized = diagonal @ system @ diagonal
        solution, flag = cg(normalized, scale*rhs, rtol=2e-13, maxiter=300)
        if flag:
            raise RuntimeError(f'Thermal linear solve failed at {time:g} s: {flag}')
        rate = scale*solution
        delta = rhs - mass*rate + dt*flux(g*(rate[j]-rate[i]))
        residual = abs(delta).max()/max(abs(rhs).max(), 1.)
        assert residual < 2e-7, (time, residual)
        worst_residual = max(worst_residual, residual)
        max_balance = max(max_balance, abs((mass*rate).sum()-power.sum())/power.sum())
        nxt = t + dt*rate
        assert np.isfinite(nxt).all()
        fraction = 1.
        reached = nxt[heart].min() >= CUTOFF
        if reached:
            cold = heart & (t < CUTOFF)
            fraction = float(np.max((CUTOFF-t[cold])/(nxt[cold]-t[cold])))
            assert 0 < fraction <= 1
        end = t + fraction*dt*rate
        # Each linear step contributes only its portion after this state first
        # crosses Tg and before the coldest-heart event.
        eligible = specimen & (crossed | (end > GLASS))
        if eligible.any():
            idx = int(np.argmin(np.where(eligible, rate, np.inf)))
            value = float(60*rate[idx])
            if value < minimum[0]:
                onset = 0. if crossed[idx] else (GLASS-t[idx]) / rate[idx]
                minimum = (value, dict(time_s=time+onset, temp_c=max(GLASS,float(t[idx]))
                                       if not crossed[idx] else float(t[idx]),
                                       phase=int(phase[idx]), xyz_mm=xyz[idx].tolist()))
        crossed |= end >= GLASS
        history.append(end)
        if reached:
            hot = int(np.argmax(np.where(specimen,end,-np.inf)))
            hot_all = int(np.argmax(end))
            event = dict(time_s=time+fraction*dt, cold_heart_c=float(end[heart].min()),
                         min_c_min=minimum[0], min_location=minimum[1],
                         hot_c=float(end[hot]), hot_phase=int(phase[hot]),
                         hot_xyz_mm=xyz[hot].tolist(), hot_any_c=float(end[hot_all]),
                         passed=bool(minimum[0]>=FLOOR and end[hot_all]<0.),
                         residual=worst_residual, heat_balance=max_balance)
            return np.array(history), event
        t = nxt
    raise RuntimeError('Heart event not reached')


def transfers(graph, phantom):
    phase, _, xyz, _, _ = graph
    points = (np.indices(phantom.shape).reshape(3,-1).T + .5) * phantom.spacing_m * 1000
    active = phantom.labels.ravel() >= 0
    labels = phantom.labels.ravel()[active]
    points = points[active]
    indices = np.empty((len(points),8),int)
    weights = np.empty_like(indices,dtype=float)
    for p in np.unique(labels):
        mask = labels == p
        subset = np.flatnonzero(phase == p)
        distance, idx = cKDTree(xyz[subset]).query(points[mask], k=8)
        indices[mask] = subset[idx]
        w = 1/np.maximum(distance, 1e-8)**2
        weights[mask] = w / w.sum(axis=1,keepdims=True)
    return indices, weights


def mechanical_scores(solver, graph, phantom, history, event, dt, cadence):
    indices, weights = transfers(graph, phantom)
    stride = round(cadence/dt)
    assert np.isclose(stride*dt,cadence)
    regular = np.arange(stride,len(history)-1,stride)
    temperatures = history[np.r_[regular,len(history)-1]]
    times = np.r_[regular*dt,event['time_s']]
    result = {}
    for mapping in ('nearest','idw8'):
        if mapping == 'nearest':
            path = temperatures[:,indices[:,0]]
        else:
            path = np.einsum('tij,ij->ti',temperatures[:,indices],weights)
        result[mapping] = summarize(times, solver.solve(times,path))
    return result


def run(dt, cadence, q, method, output, audit=None):
    phantom = rat_heart_cryobag(1.4)
    solver = Mechanics(phantom)
    graph = cut_graph(SPACING,q) if method == 'cut' else voxel_graph(SPACING)
    records = {}
    for mode in ('extrapolated','clamped'):
        records[mode] = {}
        for name, vector in VECTORS.items():
            history, event = heat_history(graph, vector, mode, dt)
            if audit is not None:
                audit(graph, history, event, dt)
            print(method, dt, mode, name, 'thermal', event, flush=True)
            stress = mechanical_scores(solver, graph, phantom, history,event,dt,cadence)
            print(method, dt, mode, name, 'stress', stress, flush=True)
            records[mode][name] = dict(event=event, stress=stress)
            Path(output).write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    return records


if __name__ == '__main__':
    dt, cadence, q, method, output = sys.argv[1:]
    run(float(dt),float(cadence),int(q),method,output)
