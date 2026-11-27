"""Recompute the original coarse voxel objective against the release."""
import json
import math
import numpy as np
from mechanics import Mechanics, design_fields, rat_heart_cryobag, thermal, summarize, UPSTREAM

def baseline():
    phantom = rat_heart_cryobag(1.4)
    solver = Mechanics(phantom)
    records = {}
    source = json.loads((UPSTREAM / 'results/objective_horizon.json').read_text())
    with np.load(UPSTREAM / 'results/render_fields.npz', allow_pickle=False) as fields:
        np.testing.assert_array_equal(phantom.labels, fields['labels'])
        for name in ('baseline', 'optimized'):
            parameters = fields[name + '_parameters']
            loading = design_fields(phantom, parameters)
            np.testing.assert_allclose(loading, fields[name + '_loading'], atol=1e-12)
            dt, steps, history, *_ = thermal(phantom, loading, parameters[3:])
            stride = math.ceil(steps / 40)
            indices = np.arange(stride, math.floor(12 / dt) + 1, stride)
            times = dt * indices
            records[name] = summarize(times, solver.solve(times, history[indices]))
    return records

if __name__ == '__main__':
    print(json.dumps(baseline(), indent=2))
