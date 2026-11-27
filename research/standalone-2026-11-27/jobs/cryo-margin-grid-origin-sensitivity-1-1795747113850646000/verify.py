"""Recompute all temperature and stress trajectories and test the reported comparisons."""
import json
import math
import sys

from analysis import ROOT, geometry, np, read, recompute, sample
from summarize_results import paired_rows


def compare(expected, actual, path="results"):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), path
        for key in expected:
            compare(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert len(expected) == len(actual), path
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{path}[{i}]")
    elif isinstance(expected, float):
        assert math.isclose(expected, actual, abs_tol=2e-7, rel_tol=2e-9), (path, expected, actual)
    else:
        assert expected == actual, (path, expected, actual)


def check_thermal(rows):
    source = {level["grid_spacing_mm"]: level["designs"]["optimized"]
              for level in read("grid_convergence_margin.json")["levels"]}
    assert len(rows) == 10
    for row in rows:
        designs = row["designs"]
        for design in designs.values():
            assert design["rate_70"] > 50
            assert design["cold_70"] > -35 and design["hot_70"] < 0
            assert design["hot_at_event_c"] < 0
            assert design["limiting_label"] == 2
            assert design["time_s"] < .31 and design["limiting_temp_c"] < -149.6
            assert math.isclose(design["rate_70"], design["rate_event"], abs_tol=1e-8)
        if row["phase"] == [0, 0, 0]:
            graded = designs["graded"]
            reference = source[row["spacing_mm"]]
            for metric, key in (("rate_70", "min_warming_rate_c_per_min"),
                                ("cold_70", "final_min_temperature_c"),
                                ("hot_70", "final_max_temperature_c")):
                assert math.isclose(graded[metric], reference[key], abs_tol=1e-8)
    coarse = [row for row in rows if row["spacing_mm"] == 1.1]
    fine = [row for row in rows if row["spacing_mm"] == .75]
    assert [row["designs"]["uniform_strict"]["rate_70"] >= 51.5 for row in coarse] == [
        True, False, False, True, True]
    assert all(design["rate_70"] >= 51.5 for row in fine for design in row["designs"].values())
    for name in fine[0]["designs"]:
        rates = [row["designs"][name]["rate_70"] for row in fine]
        assert max(rates) - min(rates) < .003
    assert len({tuple(row["designs"]["uniform_strict"]["location_mm"]) for row in coarse}) == 5
    for row in rows:
        spacing = row["spacing_mm"]
        phantom, axes = sample(spacing, geometry.rat_heart_cryobag(spacing), row["phase"])
        location = row["designs"]["uniform_strict"]["location_mm"]
        cell = np.array([np.abs(axis-coordinate).argmin()
                         for axis, coordinate in zip(axes, location)])
        offsets = np.eye(3, dtype=int)
        neighbors = [phantom.labels[tuple(cell+offset)] for offset in (*offsets, *-offsets)]
        expected = [0, 2, 2, 2, 2, 2] if row in coarse[:3] else [2]*6
        assert sorted(neighbors) == expected
        if row in coarse[:3]:
            # Closed-form axial surface intersections at the limiting cell's x, y.
            centre = np.array(phantom.spec.container_dims_mm) / 2
            inner = centre - phantom.spec.wall_thickness_mm
            heart = np.array(phantom.spec.heart_dims_mm) / 2
            heart *= np.cbrt(1000 * phantom.spec.heart_volume_ml /
                             (4 * np.pi / 3 * np.prod(heart)))
            lateral = abs(np.array(location[:2]) - centre[:2])
            exponent = phantom.pillow_exponent
            bag_z = centre[2] + inner[2] * (1 - sum((lateral/inner[:2])**exponent))**(1/exponent)
            heart_z = centre[2] + heart[2] * (1 - sum((lateral/heart[:2])**2))**.5
            assert .809 < bag_z - heart_z < .834
            assert location[2] < heart_z < bag_z < location[2] + spacing


def main():
    thermal_only = "--thermal" in sys.argv
    rows = recompute(not thermal_only)
    check_thermal(rows)
    if thermal_only:
        compare(json.loads((ROOT / "thermal.json").read_text()), rows)
    else:
        compare(json.loads((ROOT / "results.json").read_text()), rows)
        check_mechanics(rows)
    print(f"Verified {len(rows)} lattices with freshly recomputed solutions.", flush=True)


def check_mechanics(rows):
    # Independent anchor from the previous full source-cadence comparison.
    native = next(row for row in rows if row["spacing_mm"] == .75 and row["phase"] == [0, 0, 0])
    assert math.isclose(native["designs"]["graded"]["objective_mpa"], 1.1853555698417106, abs_tol=1e-8)
    assert math.isclose(native["designs"]["graded"]["peak_mpa"], 1.1111467840967504, abs_tol=1e-8)


    differences = list(paired_rows(rows))
    assert all(-.0017 < objective < -.0007 for _, _, objective, _, _, _ in differences)
    assert [peak < 0 for _, _, _, peak, _, _ in differences] == [
        True, True, True, False, True, True, True, True, True, True]
    assert all(objective > 0 for _, _, _, _, objective, _ in differences)
    assert [peak < 0 for _, _, _, _, _, peak in differences[5:]] == [
        True, False, False, False, False]


if __name__ == "__main__":
    main()
