"""Event-stopped matched-iron ablation of the released 1.5 C/min margin design."""
import argparse
import json
import math
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

from analysis import ROOT, Mechanics, design_fields, rat_heart_cryobag, read, summarize, thermal


def calculate():
    phantom = rat_heart_cryobag(1.1)
    parameters = np.array(read("optimization_margin.json")["optimizer"]["parameters"])
    loading = design_fields(phantom, parameters)
    surrounding = phantom.labels == 1
    solver = Mechanics(phantom)
    result = {}
    for name, field in (("graded", loading), ("uniform_matched", loading.copy())):
        if name == "uniform_matched":
            field[surrounding] = loading[surrounding].mean()
        dt, _, states, rates, _, glass_rates, endpoints = thermal(phantom, field, parameters[3:])
        stop, temperature = endpoints["heart"]
        labels = phantom.labels[phantom.active_mask]
        heart = labels >= 2
        after = np.flatnonzero(states[:, heart].min(axis=1) >= -35)[0]
        assert (after - 1) * dt < stop <= after * dt
        assert abs(temperature[heart].min() + 35) < 1e-10
        cadence = 4
        indices = np.arange(cadence, after, cadence)
        times = np.r_[indices * dt, stop]
        path = np.vstack((states[indices], temperature))
        result[name] = {
            "mean_surrounding_mg_ml": float(field[surrounding].mean()),
            "stop_s": stop,
            "max_heart_c": float(temperature[heart].max()),
            "all_specimen_min_rate_c_min": float(rates[:after].min()),
            "heart_min_rate_glass_to_stop_c_min": float(glass_rates[:after].min()),
            "mechanics_step_s": cadence * dt,
            **summarize(times, solver.solve(times, path)),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    actual = calculate()
    path = ROOT / "margin-results.json"
    if args.write:
        path.write_text(json.dumps(actual, indent=2) + "\n")
        return
    saved = json.loads(path.read_text())
    assert actual.keys() == saved.keys()
    for name, values in actual.items():
        assert values.keys() == saved[name].keys()
        for key, value in values.items():
            assert math.isclose(value, saved[name][key], rel_tol=2e-8, abs_tol=2e-9), (name, key)
    graded, uniform = actual["graded"], actual["uniform_matched"]
    assert abs(graded["mean_surrounding_mg_ml"] - uniform["mean_surrounding_mg_ml"]) < 1e-12
    convergence = read("grid_convergence_margin.json")
    reference = next(level for level in convergence["levels"] if level["grid_spacing_mm"] == 1.1)
    published_rate = reference["designs"]["optimized"]["min_warming_rate_c_per_min"]
    assert abs(graded["all_specimen_min_rate_c_min"] - published_rate) < 1e-8
    for design in (graded, uniform):
        assert design["all_specimen_min_rate_c_min"] > 51
        assert design["max_heart_c"] < 0
        assert design["peak_time_s"] < 12
    for key in ("peak_mpa", "objective_mpa"):
        assert abs(uniform[key] - graded[key]) < .002, key
    print("Verified matched-iron mechanical and temperature comparison for the released margin design.")


if __name__ == "__main__":
    main()
