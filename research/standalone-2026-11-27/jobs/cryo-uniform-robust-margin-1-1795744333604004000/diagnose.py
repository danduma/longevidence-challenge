"""Locate the rate-limiting voxel from consecutive thermal states."""
import json
import math

import numpy as np

from search import GRADED, LEVELS, ROOT, histories, mesh, uniform


def limiting_voxel(spacing, parameters):
    dt, n, states, rates, _, _, endpoints = histories(spacing, tuple(parameters))
    stop = endpoints["heart"][0]
    use = max(n, math.ceil(stop/dt))
    step = int(np.argmin(rates[:use]))
    phantom = mesh(spacing)
    coords = np.argwhere(phantom.active_mask)
    labels = phantom.labels[phantom.active_mask]
    cell = int(np.argmin(np.where(labels > 0, states[step+1]-states[step], np.inf)))
    return dict(time_s=step*dt, cell=coords[cell].tolist(),
                label=int(labels[cell]), temperature_c=float(states[step,cell]),
                rate_c_min=float((states[step+1,cell]-states[step,cell])*60/dt))


if __name__ == "__main__":
    candidates = json.loads((ROOT/"candidates.json").read_text())
    cases = {"graded": GRADED, **{k: uniform(v) for k,v in candidates.items()}}
    result = {str(s): {k: limiting_voxel(s, x) for k,x in cases.items()} for s in LEVELS}
    print(json.dumps(result, indent=2))
