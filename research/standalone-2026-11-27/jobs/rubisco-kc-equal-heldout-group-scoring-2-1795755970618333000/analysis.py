"""Rescore frozen KC predictions under two empirical distributions."""
import json
from collections import Counter
from itertools import groupby
from math import fsum, sqrt
from pathlib import Path
from statistics import fmean

ROOT = Path(__file__).resolve().parent


def midranks(values, weights):
    result = [0.0] * len(values)
    below = 0.0
    ordered = sorted(range(len(values)), key=values.__getitem__)
    for _, members in groupby(ordered, key=values.__getitem__):
        tied = list(members)
        mass = fsum(weights[i] for i in tied)
        for i in tied:
            result[i] = below + mass / 2
        below += mass
    return result


def score(y, pred, weights):
    total = fsum(weights)
    a, b = (midranks(values, weights) for values in (y, pred))
    # Each marginal weighted midrank has mean equal to half the total mass.
    a = [rank - total / 2 for rank in a]
    b = [rank - total / 2 for rank in b]
    variance_a = fsum(w*x*x for w, x in zip(weights, a))
    variance_b = fsum(w*x*x for w, x in zip(weights, b))
    denominator = sqrt(variance_a * variance_b)
    rho = fsum(w*x*z for w, x, z in zip(weights, a, b)) / denominator if denominator else None
    rmse = sqrt(fsum(w*(x-z)**2 for w, x, z in zip(weights, y, pred)) / total)
    return dict(rho=rho, rmse=rmse)


def analyze(data):
    rows = data["rows"]
    counts = Counter(row["group"] for row in rows)
    folds = data["tests"]
    result = {
        "n_rows": len(rows),
        "n_groups": len(counts),
        "group_sizes": dict(sorted(Counter(counts.values()).items())),
        "models": {},
    }
    for name, pred in data["predictions"].items():
        model = {}
        for unit in ("row", "group"):
            weights = [1 / counts[row["group"]] if unit == "group" else 1 for row in rows]
            pooled = score([row["y"] for row in rows], pred, weights)
            within = [
                score([rows[i]["y"] for i in test], [pred[i] for i in test],
                      [weights[i] for i in test])["rho"]
                for test in folds
            ]
            defined = [value for value in within if value is not None]
            model[unit] = dict(**pooled, fold_rho=within, mean_fold_rho=fmean(defined) if defined else None)
        result["models"][name] = model
    return result


if __name__ == "__main__":
    data = json.loads((ROOT / "sources/frozen.json").read_text())
    comparators = json.loads((ROOT / "sources/comparators.json").read_text())
    data["predictions"].update(comparators["predictions"])
    (ROOT / "results.json").write_text(json.dumps(analyze(data), indent=2, allow_nan=False) + "\n")
