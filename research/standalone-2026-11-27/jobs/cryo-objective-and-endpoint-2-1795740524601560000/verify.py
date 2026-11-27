"""Recompute the thermal histories, finite-element states and audit conclusions."""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import json
import math

import numpy as np

from analysis import ROOT, UPSTREAM, design_fields, rat_heart_cryobag, read, run, thermal


def check_tree(actual, saved, path=""):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys(), path
        for key in actual:
            check_tree(actual[key], saved[key], f"{path}/{key}")
    elif isinstance(actual, int):
        assert type(saved) is int and actual == saved, (path, actual, saved)
    else:
        assert math.isclose(actual, saved, rel_tol=2e-8, abs_tol=2e-9), (path, actual, saved)


def check_upstream(coarse):
    original = read("optimization_extrapolated.json")
    horizon = read("objective_horizon.json")
    fields = np.load(UPSTREAM / "results/render_fields.npz", allow_pickle=False)
    assert str(fields["sar_mode"]) == original["sar_mode"] == "extrapolated"
    np.testing.assert_array_equal(fields["baseline_parameters"], [7.5, 7.5, 7.5, 1, 1, 1])
    np.testing.assert_allclose(fields["optimized_parameters"], original["optimizer"]["parameters"],
                               rtol=0, atol=1e-13)
    for name in ("baseline", "optimized"):
        design = coarse[name]
        wrapped = original["baseline"] if name == "baseline" else original["optimum"]
        for horizon_seconds in ("12", "30"):
            source = horizon["designs"][name][horizon_seconds]
            reproduction = design["upstream_"+horizon_seconds]
            assert abs(reproduction["peak_mpa"] - source["peak_tension_mpa"]) < 1e-8
            assert abs(reproduction["peak_time_s"] - source["peak_time_s"]) < 1e-10
        expected_objective = wrapped["stress_objective_mpa"]
        assert abs(design["upstream_12"]["objective_mpa"] - expected_objective) < 1e-8
        assert abs(design["cold_specimen_70_c"] - wrapped["final_min_temperature_c"]) < 1e-8
        assert abs(design["min_rate_70_c_min"] - wrapped["min_warming_rate_c_per_min"]) < 1e-8


def check_endpoint_brackets(coarse):
    fields = np.load(UPSTREAM / "results/render_fields.npz", allow_pickle=False)
    phantom = rat_heart_cryobag(1.4)
    labels = phantom.labels[phantom.active_mask]
    for name in ("baseline", "optimized"):
        params = fields[name+"_parameters"]
        dt, _, history, _, _, _, events = thermal(phantom, design_fields(phantom, params), params[3:])
        for endpoint, mask in (("specimen", labels > 0), ("heart", labels >= 2)):
            time, state = events[endpoint]
            minimum = history[:, mask].min(axis=1)
            after = np.flatnonzero(minimum >= -35)[0]
            assert (after-1)*dt < time <= after*dt
            assert abs(state[mask].min()+35) < 1e-10
            assert abs(coarse[name][endpoint+"_event"]["time_s"] - time) < 1e-10


def check_conclusions(coarse, fine):
    def event(design):
        return design["heart_event"]

    baseline, optimized, uniform, full = (coarse[key] for key in
        ("baseline", "optimized", "uniform_matched", "unthrottled"))
    assert event(baseline)["time_s"] > 75 > event(optimized)["time_s"]
    assert event(baseline)["min_heart_rate_glass_to_event"] < 50
    for design in (optimized, uniform, full):
        assert event(design)["max_heart_c"] < 0
        assert event(design)["min_heart_rate_glass_to_event"] > 50
        assert event(design)["min_rate"] > 50
    for stride in ("10", "5", "2", "1"):
        original, staged, simple, unthrottled = [
            item["cadences"][stride]["heart"] for item in (baseline, optimized, uniform, full)]
        for key in ("peak_mpa", "objective_mpa"):
            assert original[key] - staged[key] > .04
            assert unthrottled[key] - staged[key] > .02
        assert max(item["peak_time_s"] for item in (original, staged, simple, unthrottled)) < 12
    for grid in (coarse, fine):
        staged, simple = grid["optimized"], grid["uniform_matched"]
        assert abs(staged["mean_surrounding_mg_ml"] - simple["mean_surrounding_mg_ml"]) < 1e-12
        assert 0 < event(staged)["time_s"] - event(simple)["time_s"] < .07
        for design in (staged, simple):
            assert event(design)["min_rate"] > 50
            assert event(design)["max_specimen_c"] < 0
        for stride in staged["cadences"]:
            staged_score = staged["cadences"][stride]["heart"]
            simple_score = simple["cadences"][stride]["heart"]
            for key in ("peak_mpa", "objective_mpa"):
                assert 0 < simple_score[key] - staged_score[key] < .001


def main():
    coarse = run()
    fine = run((5, 2), ("optimized", "uniform_matched"), spacing=1.2)
    actual = {"grid_1.4_mm": coarse, "grid_1.2_mm": fine}
    check_upstream(coarse)
    check_endpoint_brackets(coarse)
    check_conclusions(coarse, fine)
    check_tree(actual, json.loads((ROOT / "results.json").read_text()))
    print("Verified source reproduction, end-temperature brackets, and the event-stopped mechanical comparison.")


if __name__ == "__main__":
    main()
