"""Recompute event-stopped comparisons from the pinned physical inputs."""
from concurrent.futures import ProcessPoolExecutor
import json
import math
import os
import sys

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

from search import (LEVELS, GRADED, ROOT, histories, mechanics, read, stress,
                    summarize, thermal_metrics, uniform)


def calculate_case(task):
    spacing, name, x = task
    assert np.all(np.r_[np.zeros(3), [.2]*3] <= x)
    assert np.all(x <= np.r_[[10]*3, np.ones(3)])
    thermal = thermal_metrics(spacing, x)
    assert abs(thermal["event_cold_heart_c"] + 35) < 1e-9
    assert thermal["event_hot_c"] < 0
    assert thermal["cold_70"] >= -35 and thermal["hot_70"] <= 0
    floor = 51.5 if spacing == 1.4 or name == "uniform_strict" else 50.
    assert min(thermal["event_min_rate"], thermal["rate_70"]) >= floor

    if name == "graded":
        reference = read("grid_convergence_margin.json")
        level = next(x for x in reference["levels"] if x["grid_spacing_mm"] == spacing)
        source = level["designs"]["optimized"]
        for field, key in (("rate_70", "min_warming_rate_c_per_min"),
                           ("cold_70", "final_min_temperature_c"),
                           ("hot_70", "final_max_temperature_c")):
            assert abs(thermal[field] - source[key]) < 1e-8

    late_cadence = 1.8 if spacing == .75 else None
    values = {**thermal, **stress(spacing, x, late_cadence=late_cadence)}
    print(spacing, name, values, flush=True)
    return str(spacing), name, values


def calculate(candidates):
    cases = {"graded": GRADED, **{name: uniform(x) for name, x in candidates.items()}}
    output = {str(spacing): {} for spacing in LEVELS}
    tasks = [(s, name, x) for s in reversed(LEVELS) for name, x in cases.items()]
    workers = max(1, int(os.environ.get("VITRIFY_VERIFY_WORKERS", "3")))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for spacing, name, values in executor.map(calculate_case, tasks):
            output[spacing][name] = values
    return output


def check_source_objective():
    dt, n, states, _, _, _, _ = histories(1.4, tuple(GRADED))
    stride = math.ceil(n/40)
    indices = np.arange(stride, math.floor(12/dt)+1, stride)
    times = indices*dt
    replay = summarize(times, mechanics(1.4).solve(times, states[indices]))
    reference = read("optimization_margin.json")["optimum"]
    assert abs(replay["objective_mpa"] - reference["stress_objective_mpa"]) < 1e-8
    assert abs(replay["peak_mpa"] - reference["peak_tension_mpa"]) < 1e-8


def assert_close(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            assert_close(actual[key], saved[key])
    else:
        assert math.isclose(actual, saved, rel_tol=2e-8, abs_tol=2e-9), (actual, saved)


def main():
    check_source_objective()
    candidates = json.loads((ROOT / "candidates.json").read_text())
    actual = calculate(candidates)
    for spacing, cases in actual.items():
        graded = cases["graded"]
        uniform_policy = cases["uniform_released_policy"]
        assert uniform_policy["objective_mpa"] < graded["objective_mpa"], spacing
        assert uniform_policy["peak_mpa"] < graded["peak_mpa"], spacing
        assert abs(cases["uniform_strict"]["objective_mpa"] - graded["objective_mpa"]) < .0013
        assert (graded["rate_70"] >= 51.5) == (spacing != "1.1")
    path = ROOT / "results.json"
    if "--write" in sys.argv:
        path.write_text(json.dumps(actual, indent=2) + "\n")
    else:
        assert_close(actual, json.loads(path.read_text()))
        print("Verified the event-stopped stress values and all released-mesh feasibility constraints.")


if __name__ == "__main__":
    main()
