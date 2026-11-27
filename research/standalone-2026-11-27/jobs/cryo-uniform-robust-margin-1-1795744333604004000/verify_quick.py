"""Check the ranking reversal on 1.1 mm and both rate policies on all meshes."""
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import json

from search import GRADED, LEVELS, ROOT, mesh, stress, thermal_metrics, uniform
from analysis import design_fields
from verify import check_source_objective


def main():
    check_source_objective()
    candidates = json.loads((ROOT / "candidates.json").read_text())
    replay = {"graded": GRADED, **{k: uniform(x) for k, x in candidates.items()}}
    phen = mesh(1.1)
    mask = phen.labels == 1
    loading = float(design_fields(phen, GRADED)[mask].mean())
    fixed_duty = [loading]*3 + list(GRADED[3:])
    scored = {k: stress(1.1, x)["objective_mpa"] for k, x in replay.items()}
    matched = stress(1.1, fixed_duty)["objective_mpa"]
    assert matched - scored["graded"] > .001
    assert scored["graded"] - scored["uniform_released_policy"] > .001
    for s in LEVELS:
        original = thermal_metrics(s, GRADED)["rate_70"]
        assert (original >= 51.5) == (s != 1.1)
        for name, x in replay.items():
            m = thermal_metrics(s, x)
            floor = 51.5 if name == "uniform_strict" or s == 1.4 else 50.
            assert m["cold_70"] >= -35 and m["hot_70"] < 0 and m["event_hot_c"] < 0
            assert min(m["rate_70"], m["event_min_rate"]) >= floor
    print("Verified fixed-duty penalty, tuned reversal, and the all-mesh rates.")


if __name__ == "__main__":
    main()
