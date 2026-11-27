"""Independent fold reconstruction, Gaussian solver and KNN replay."""
import csv
import json
from collections import Counter, defaultdict
from math import isclose, sqrt
from pathlib import Path
from statistics import fmean, median

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


def split(indices, n):
    sizes = Counter(groups[i] for i in indices)
    loads = [0] * n
    allocation = {}
    for group in sorted(sizes, key=lambda g: (sizes[g], g), reverse=True):
        fold = loads.index(min(loads))
        loads[fold] += sizes[group]
        allocation[group] = fold
    return [[i for i in sorted(indices) if allocation[groups[i]] == f] for f in range(n)]


def purged(indices, held):
    return [i for i in sorted(set(indices) - set(held))
            if not any(direct_edge(i, j) for j in held)]


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
        assert isclose(actual, expected, abs_tol=1.1e-10, rel_tol=0), (actual, expected)


tests = split(all_rows, 5)


def metrics(vector):
    return {
        "rho": spearman(y, vector),
        "rmse": sqrt(sum((y[i] - vector[i])**2 for i in all_rows) / len(y)),
        "fold_rho": [spearman([y[i] for i in test], [vector[i] for i in test]) for test in tests],
    }


def features():
    source = ROOT / "sources"
    sequences = defaultdict(str)
    for line in (source / "released.fasta").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = ""
        else:
            sequences[name] += line.strip()
    sequences["Rhodopseudomonas_palustris"] = json.loads(
        (source / "Q6N0W9.json").read_text())["sequence"]["value"]
    result = []
    for row in rows:
        sequence = sequences[row["label"]]
        counts = Counter(map("".join, zip(sequence, sequence[1:], sequence[2:])))
        norm = sqrt(sum(count**2 for count in counts.values()))
        result.append({word: count/norm for word, count in counts.items()})
    return result


vectors = features()
gram = [[sum(a[word] * b[word] for word in sorted(a.keys() & b.keys()))
         for b in vectors] for a in vectors]
alphas = [.001, .01, .1, 1, 10, 100]


def eliminate(system):
    n = len(system)
    for j in range(n):
        pivot = max(range(j, n), key=lambda i: abs(system[i][j]))
        system[j], system[pivot] = system[pivot], system[j]
        lead = system[j]
        for i in range(j+1, n):
            row = system[i]
            ratio = row[j] / lead[j]
            for k in range(j+1, n+1):
                row[k] -= ratio*lead[k]
    x = [0.] * n
    for i in reversed(range(n)):
        row = system[i]
        x[i] = (row[-1] - sum(row[j]*x[j] for j in range(i+1, n))) / row[i]
    return x


def ridge(train, held, alpha, balanced):
    grouped = defaultdict(list)
    for i in train:
        grouped[groups[i]].append(y[i])
    offset = fmean(map(fmean, grouped.values())) if balanced else fmean(y[i] for i in train)
    scale = [sqrt(1/len(grouped[groups[i]])) if balanced else 1. for i in train]
    # Solve the symmetrically rescaled dual, then recover the original coefficients.
    system = [[scale[a]*scale[b]*gram[i][j] + (alpha if a == b else 0.)
               for b, j in enumerate(train)] + [scale[a]*(y[i]-offset)]
              for a, i in enumerate(train)]
    coeff = [s*c for s, c in zip(scale, eliminate(system))]
    return [round(offset + sum(gram[i][j]*c for j, c in zip(train, coeff)), 10)
            for i in held]


ridge_modes = ("unweighted", "balanced", "balanced_fixed_alpha")
knn_modes = ("nearest_group", "median_group")
outputs = {name: {} for name in (*ridge_modes, *knn_modes)}
choices = {name: [] for name in knn_modes}
details = []
for index, test in enumerate(tests):
    train = purged(all_rows, test)
    validation = split(train, 4)
    cached = data["direct"][index]
    assert (test, train, validation) == (cached["test"], cached["train"], cached["validation"])
    inner = [(purged(train, v), v) for v in validation]
    pools = [train, *[a for a, _ in inner]]
    detail = {"train_sizes": list(map(len, pools)),
              "group_counts": [len({groups[i] for i in a}) for a in pools]}
    order = [i for _, v in inner for i in v]
    for name in ridge_modes[:2]:
        scores = [spearman([y[i] for i in order], [p for a, v in inner
                   for p in ridge(a, v, alpha, name == "balanced")]) for alpha in alphas]
        best = max(range(len(alphas)), key=scores.__getitem__)
        detail[name] = {"alpha": alphas[best], "inner_rho": scores}
    assert detail["unweighted"]["alpha"] == cached["chosen"]["trimer"][0]
    for name in ridge_modes:
        tuned = "balanced" if name == "balanced" else "unweighted"
        outputs[name].update(zip(test, ridge(train, test, detail[tuned]["alpha"], name != "unweighted")))
    for mode in knn_modes:
        scores = [spearman([y[i] for i in order],
                           [predict(a, i, k, mode) for a, v in inner for i in v]) for k in grid]
        k = max(grid, key=lambda k: scores[grid.index(k)])
        choices[mode].append(k)
        outputs[mode].update({i: predict(train, i, k, mode) for i in test})
    details.append(detail)

computed = {"folds": details}
for name, output in outputs.items():
    assert set(output) == all_rows
    vector = [output[i] for i in range(len(y))]
    computed[name] = {**metrics(vector), "predictions": vector}
    if name in choices:
        computed[name]["k"] = choices[name]
close(computed, saved)
close(computed["unweighted"]["predictions"], data["trimer"])
u, b, f = (computed[name] for name in ridge_modes)
for name in knn_modes:
    k = computed[name]
    assert all(model["rho"] > k["rho"] and model["rmse"] < k["rmse"] for model in (u, b, f))
assert b["rho"] < f["rho"] < u["rho"]
assert f["rmse"] < u["rmse"] < b["rmse"]
assert sum(d["balanced"]["alpha"] != d["unweighted"]["alpha"] for d in details) == 1
assert sum(a > c for a, c in zip(b["fold_rho"], computed["median_group"]["fold_rho"])) == 2
print("Verified nested ridge and capped KNN predictions, both pooled leads, and the fixed-alpha contrast.")
