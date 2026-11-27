"""Independently rebuild eligibility, nested predictions and the rank comparison."""
import csv
import json
from collections import Counter, defaultdict
from math import isclose, sqrt
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / "sources/kc-inputs.json").read_text())
saved = json.loads((ROOT / "results.json").read_text())
rows, identity = data["rows"], data["identity"]
archived = list(csv.DictReader((ROOT / "sources/pred_KC.csv").open()))
assert len(rows) == len(archived) == 197
assert [r["species"] for r in rows] == [r["species"] for r in archived]
y = [float(r["y_true"]) for r in archived]
assert all(isclose(a, row["y"], abs_tol=1e-12, rel_tol=0) for a, row in zip(y, rows))
# Use the predecessor's pre-rounding targets after verifying release agreement.
y = [row["y"] for row in rows]
groups = [r["group"] for r in rows]
genus = [r["species"].lower().split()[0] for r in rows]
all_rows = set(range(len(y)))
grid = [1, 3, 5, 10]


def direct_edge(i, j):
    return genus[i] == genus[j] or identity[i][j] >= .95


components = {}
for seed in sorted(all_rows):
    if seed in components:
        continue
    pending = [seed]
    components[seed] = seed
    while pending:
        i = pending.pop()
        for j in all_rows - components.keys():
            if groups[i] == groups[j] or direct_edge(i, j):
                components[j] = seed
                pending.append(j)
assert all((components[i] == components[j]) ==
           (data["component_id"][i] == data["component_id"][j])
           for i in all_rows for j in all_rows)


def split(indices, n):
    sizes = Counter(groups[i] for i in indices)
    loads = [0] * n
    allocation = {}
    for group in sorted(sizes, key=lambda g: (sizes[g], g), reverse=True):
        fold = loads.index(min(loads))
        loads[fold] += sizes[group]
        allocation[group] = fold
    return [[i for i in sorted(indices) if allocation[groups[i]] == f] for f in range(n)]


def purged(indices, held, pool):
    if pool == "direct":
        return [i for i in sorted(set(indices) - set(held)) if
                not any(direct_edge(i, j) for j in held)]
    blocked = {components[i] for i in held}
    return [i for i in sorted(indices) if components[i] not in blocked]


def predict(eligible, query, k, mode):
    ordered = sorted(eligible, key=lambda i: (identity[query][i], -i), reverse=True)
    if mode != "rows":
        group_labels = defaultdict(list)
        for i in ordered:
            group_labels[groups[i]].append(y[i])
        first = {}
        for i in ordered:
            first.setdefault(groups[i], i)
        ordered = list(first.values())
    threshold = identity[query][ordered[min(k, len(ordered)) - 1]]
    chosen = [i for i in ordered if identity[query][i] >= threshold]
    if mode == "median_group":
        labels = [median(group_labels[groups[i]]) for i in chosen]
    else:
        labels = [y[i] for i in chosen]
    return round(sum(labels) / len(labels), 10)


def centered_ranks(values):
    counts = Counter(values)
    center = (len(values) - 1) / 2
    mapping, previous = {}, 0
    for v, count in sorted(counts.items()):
        mapping[v] = previous + (count - 1) / 2 - center
        previous += count
    return [mapping[v] for v in values]


def spearman(obs, pred):
    a, b = centered_ranks(obs), centered_ranks(pred)
    norm = sqrt(sum(x*x for x in a) * sum(x*x for x in b))
    return sum(x*z for x, z in zip(a, b)) / norm if norm else None


def close(actual, expected):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys()
        for k, v in actual.items():
            close(v, expected[k])
    elif isinstance(actual, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            close(a, b)
    elif actual is None:
        assert expected is None
    else:
        assert isclose(actual, expected, abs_tol=1e-12, rel_tol=0), (actual, expected)


tests = split(all_rows, 5)


def metrics(vector):
    return {
        "rho": spearman(y, vector),
        "rmse": sqrt(sum((y[i] - vector[i])**2 for i in all_rows) / len(y)),
        "fold_rho": [spearman([y[i] for i in test], [vector[i] for i in test]) for test in tests],
    }


computed = {name: metrics(data[name]) for name in ("trimer", "embedding")}
for pool in ("direct", "component"):
    computed[pool] = {}
    for mode in ("rows", "nearest_group", "median_group"):
        nested, row_k = [None] * len(y), [None] * len(y)
        fixed = {k: [None] * len(y) for k in grid}
        choices, scores = [], []
        assert len(data[pool]) == len(tests)
        for test, fold in zip(tests, data[pool]):
            train = purged(all_rows, test, pool)
            validation = split(train, 4)
            assert (test, train, validation) == (fold["test"], fold["train"], fold["validation"])
            inner = [(purged(train, held, pool), held) for held in validation]
            order = [i for _, held in inner for i in held]
            inner_rho = [spearman([y[i] for i in order],
                                  [predict(fit, i, k, mode) for fit, held in inner for i in held])
                         for k in grid]
            best = max(range(len(grid)), key=lambda j: float("-inf") if inner_rho[j] is None else inner_rho[j])
            k = grid[best]
            choices.append(k)
            scores.append(inner_rho)
            row_choice = fold["chosen"]["identity_knn"][0]
            if mode == "rows":
                assert k == row_choice
            for i in test:
                for candidate in grid:
                    fixed[candidate][i] = predict(train, i, candidate, mode)
                nested[i] = fixed[k][i]
                row_k[i] = fixed[row_choice][i]
        if mode == "rows":
            close(nested, data["prior_knn"][pool])
        computed[pool][mode] = {
            "k": choices, "inner_rho": scores, "predictions": nested, "nested": metrics(nested),
            "at_row_selected_k": metrics(row_k),
            "fixed": {str(k): metrics(p) for k, p in fixed.items()},
        }

close(computed, saved)
row = computed["direct"]["rows"]["nested"]
for mode in ("nearest_group", "median_group"):
    result = computed["direct"][mode]
    held_k = result["at_row_selected_k"]
    assert held_k["rho"] < row["rho"] < result["nested"]["rho"] < computed["trimer"]["rho"]
    assert computed["trimer"]["rmse"] < result["nested"]["rmse"] < held_k["rmse"] < row["rmse"]
    assert computed["component"][mode]["predictions"] == computed["component"]["rows"]["predictions"]
    for metric in [*result["fixed"].values(), result["nested"]]:
        assert metric["rho"] < computed["trimer"]["rho"]
print("Verified nested KC predictions, fixed-k contrasts and preserved trimer lead on 197 rows.")
