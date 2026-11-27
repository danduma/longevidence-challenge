"""Constrained uniform-shell search and released-mesh comparison."""
from functools import lru_cache
import json
import math
import sys

import numpy as np
from scipy.optimize import minimize
from analysis import ROOT, Mechanics, design_fields, rat_heart_cryobag, read, summarize, thermal

LEVELS = tuple(x["grid_spacing_mm"] for x in read("grid_convergence_margin.json")["levels"])
GRADED = np.array(read("optimization_margin.json")["optimizer"]["parameters"])

@lru_cache(None)
def mesh(spacing):
    return rat_heart_cryobag(spacing)

@lru_cache(None)
def mechanics(spacing):
    return Mechanics(mesh(spacing))

@lru_cache(maxsize=2)
def histories(spacing, parameters):
    phantom = mesh(spacing)
    x = np.array(parameters)
    return thermal(phantom, design_fields(phantom, x), x[3:])

def thermal_metrics(spacing, parameters):
    dt, n, states, rates, _, _, endpoints = histories(spacing, tuple(parameters))
    labels = mesh(spacing).labels[mesh(spacing).active_mask]
    stop, temperature = endpoints["heart"]
    return dict(stop_s=float(stop), event_min_rate=float(rates[:math.ceil(stop/dt)].min()),
                event_hot_c=float(temperature[labels>0].max()),
                event_cold_heart_c=float(temperature[labels>=2].min()),
                rate_70=float(rates[:n].min()),
                cold_70=float(states[n,labels>0].min()),
                hot_70=float(states[n,labels>0].max()))

def stress(spacing, parameters, full=True, cadence=.24, late_cadence=None):
    dt, _, states, _, _, _, events = histories(spacing, tuple(parameters))
    stop, temperature = events["heart"]
    stride = round(cadence/dt)
    last = stop if full else 16.
    indices = np.arange(stride, math.ceil(last/dt), stride)
    if full and late_cadence is not None:
        early = indices[indices*dt <= 16]
        late_stride = round(late_cadence/dt)
        late = np.arange(early[-1]+late_stride, math.ceil(stop/dt), late_stride)
        indices = np.r_[early, late]
    times = indices*dt
    path = states[indices]
    if full:
        times = np.r_[times, stop]
        path = np.vstack((path, temperature))
    return summarize(times, mechanics(spacing).solve(times, path))

def uniform(x):
    return [x[0]*10]*3 + list(x[1:])

def margins(spacing, x, floor):
    r = thermal_metrics(spacing, uniform(x))
    return np.array([r["rate_70"]-floor, r["cold_70"]+35,
                     -r["hot_70"], -r["event_hot_c"]])

def optimize(strict=False, budget=65):
    floor = 51.5 if strict else 50.
    start = np.array([.94, .85, .995, 1.])
    bounds = [(0, 1), (.2, 1), (.2, 1), (.2, 1)]
    lower, upper = np.array(bounds).T
    def constraints(x):
        return np.r_[margins(1.4, x, 51.5), margins(1.1, x, floor)]
    best = [None, float("inf")]
    def objective(x):
        value = stress(1.4, uniform(x), full=False)["objective_mpa"]
        if (value < best[1] and np.all(lower <= x) and np.all(x <= upper)
                and constraints(x).min() >= 0):
            best[:] = [x.tolist(), value]
            print("best", strict, best, flush=True)
            (ROOT / f"bounded-search-{strict}.json").write_text(json.dumps(best)+"\n")
        return value
    minimize(objective, start, method="COBYLA", bounds=bounds,
             constraints={"type":"ineq","fun":constraints},
             options={"maxiter":budget,"rhobeg":.017,"tol":.00012,"catol":1e-8})
    assert best[0] is not None
    return best

if __name__ == "__main__":
    optimize(sys.argv[1] == "strict", int(sys.argv[2]))
