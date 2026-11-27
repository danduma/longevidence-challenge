"""Independently select and refit each arm, then check the saved results."""
import csv
import json
from collections import Counter
from itertools import chain
from math import isclose, lcm, sqrt
from statistics import correlation, fmean, median

from analysis import ROOT, SOURCES, kernel, load


def tied_ranks(values):
    counts = Counter(values)
    ranks, below = {}, 0
    for value in sorted(counts):
        ranks[value] = below + (counts[value] + 1) / 2
        below += counts[value]
    return [ranks[v] for v in values]


def rho(y, pred, repeats):
    if len(set(y)) == 1 or len(set(pred)) == 1:
        return None
    a, b = (tied_ranks(list(chain.from_iterable([v] * n for v, n in zip(x, repeats))))
            for x in (y, pred))
    return correlation(a, b)


def solve(matrix, response):
    """Gaussian elimination on a positive-definite system."""
    n = len(response)
    for k in range(n - 1):
        pivot = matrix[k]
        for i in range(k + 1, n):
            row = matrix[i]
            factor = row[k] / pivot[k]
            for j in range(k + 1, n):
                row[j] -= factor * pivot[j]
            response[i] -= factor * response[k]
    for i in range(n - 1, -1, -1):
        response[i] -= sum(matrix[i][j] * response[j] for j in range(i + 1, n))
        response[i] /= matrix[i][i]
    return response


def close(actual, expected, tolerance=2e-9):
    if actual is None or expected is None:
        assert actual is expected, (actual, expected)
    else:
        assert isclose(actual, expected, rel_tol=0, abs_tol=tolerance), (actual, expected)


def verify():
    data, sequences = load()
    gram = kernel(sequences)
    rows = data["rows"]
    groups = [row["group"] for row in rows]
    y = [row["y"] for row in rows]
    identity = data["identity"]
    component = data["component_id"]
    genus = [r["species"].split()[0].casefold() for r in rows]
    results = json.loads((ROOT / "results.json").read_text())
    reference = json.loads((SOURCES / "outer-reference.json").read_text())
    with (SOURCES / "pred_KC.csv").open() as stream:
        released = list(csv.DictReader(stream))
    assert len(released) == len(rows) == results["n_rows"]
    assert len(set(groups)) == results["n_groups"]
    for a, b in zip(rows, released):
        assert a["species"] == b["species"]
        close(a["y"], float(b["y_true"]), 1e-12)

    def multiplicity(indices):
        counts = Counter(groups[i] for i in indices)
        scale = lcm(*counts.values())
        return [scale // counts[groups[i]] for i in indices]

    def metrics(predictions):
        mse = {}
        for g in set(groups):
            mse[g] = fmean((y[i] - predictions[i])**2 for i in indices if groups[i] == g)
        return {"rho": rho(y, predictions, multiplicity(indices)), "rmse": sqrt(fmean(mse.values()))}

    def check_metrics(predictions, recorded):
        for name, value in metrics(predictions).items():
            close(value, recorded[name])

    def predict(model, setting, train, valid):
        if model in ("ridge", "ridge_balanced"):
            counts = Counter(groups[j] for j in train)
            denominator = [counts[groups[j]] if model == "ridge_balanced" else 1 for j in train]
            center = sum(y[j] / d for j, d in zip(train, denominator)) / sum(1 / d for d in denominator)
            matrix = [[gram[i][j] for j in train] for i in train]
            for i, d in enumerate(denominator):
                matrix[i][i] += setting * d
            coef = solve(matrix, [y[i] - center for i in train])
            return [round(center + sum(gram[i][j] * a for j, a in zip(train, coef)), 10) for i in valid]
        values = []
        donor = {i: i if model == "knn" else groups[i] for i in train}
        members = {g: [j for j in train if donor[j] == g] for g in set(donor.values())}
        for i in valid:
            best_by_group = {}
            for j in sorted(train, key=lambda j: (-identity[i][j], j)):
                best_by_group.setdefault(donor[j], j)
            nearest = list(best_by_group.items())
            cutoff = identity[i][nearest[min(setting, len(nearest)) - 1][1]]
            labels = [median(y[h] for h in members[g]) if model == "knn_median" else y[j]
                      for g, j in nearest if identity[i][j] >= cutoff]
            values.append(round(fmean(labels), 10))
        return values

    indices = list(range(len(y)))
    for pool, models in reference["baseline"].items():
        folds = data[pool]
        assert sorted(chain.from_iterable(f["test"] for f in folds)) == indices
        for model, baseline in models.items():
            record = results["pools"][pool][model]
            grid = (.001, .01, .1, 1, 10, 100) if model.startswith("ridge") else (1, 3, 5, 10)
            assert list(grid) == record["grid"]
            outer = {unit: {} for unit in ("row", "group")}
            fixed = {setting: {} for setting in grid}
            for fold, details in zip(folds, record["folds"], strict=True):
                train, test = fold["train"], fold["test"]
                def eligible(i, held):
                    return all(genus[i] != genus[j] and identity[i][j] < .95 for j in held) and (
                        pool == "direct" or all(component[i] != component[j] for j in held))
                assert train == [i for i in indices if eligible(i, test)]
                assert set(groups[i] for i in train).isdisjoint(groups[i] for i in test)
                held_sets = fold["validation"]
                assert sorted(chain.from_iterable(held_sets)) == train
                inner = {unit: [] for unit in outer}
                for setting in grid:
                    predicted = {}
                    for held in held_sets:
                        donors = [i for i in train if eligible(i, held)]
                        assert donors and set(groups[i] for i in donors).isdisjoint(groups[i] for i in held)
                        predicted.update(zip(held, predict(model, setting, donors, held)))
                    for unit, repeats in (("row", [1] * len(train)), ("group", multiplicity(train))):
                        inner[unit].append(rho([y[i] for i in train], [predicted[i] for i in train], repeats))
                    fixed[setting].update(zip(test, predict(model, setting, train, test)))
                for unit, scores in inner.items():
                    best = max(range(len(grid)), key=lambda j: scores[j] if scores[j] is not None else -2)
                    assert grid[best] == details[unit]["choice"], (pool, model, scores, details)
                    for value, saved in zip(scores, details[unit]["inner_rho"], strict=True):
                        close(value, saved)
                    outer[unit].update((i, fixed[grid[best]][i]) for i in test)
                    if unit == "row" and model in ("ridge", "knn"):
                        key = "trimer" if model == "ridge" else "identity_knn"
                        assert grid[best] == fold["chosen"][key][0]
            for unit, predictions in outer.items():
                vector = [predictions[i] for i in indices]
                for value, saved in zip(vector, record["predictions"][unit], strict=True):
                    close(value, saved)
                if unit == "row":
                    for value, saved in zip(vector, baseline, strict=True):
                        close(value, saved)
                check_metrics(vector, record["outer"][unit])
            for name in ("rho", "rmse"):
                row, group = (metrics([outer[u][i] for i in indices])[name] for u in ("row", "group"))
                close(group - row, record["delta"][name])
            for setting, pred in fixed.items():
                check_metrics([pred[i] for i in indices], record["fixed_outer"][str(setting)])
    for name, pred in reference["fixed"].items():
        check_metrics(pred, results["fixed_references"][name])
    print("Verified deposited targets, preserved nested eligibility, independently refitted all grids, "
          "recovered original choices, and recomputed paired outer results.")


if __name__ == "__main__":
    verify()
