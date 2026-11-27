"""Paired nested KC selection with fixed pools, grids and outer scoring."""
import json
from collections import Counter, defaultdict
from itertools import groupby
from math import fsum, sqrt
from pathlib import Path
from statistics import fmean, median

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"
ALPHAS = (.001, .01, .1, 1, 10, 100)
NEIGHBORS = (1, 3, 5, 10)
MODELS = {
    "direct": ("ridge", "ridge_balanced", "knn", "knn_nearest", "knn_median"),
    "component": ("ridge", "knn"),
}
def load():
    data = json.loads((SOURCES / "kc-inputs.json").read_text())
    fasta = {}
    for record in (SOURCES / "released.fasta").read_text().split(">")[1:]:
        header, *lines = record.splitlines()
        fasta[header.split()[0]] = "".join(lines)
    fasta["Rhodopseudomonas_palustris"] = json.loads(
        (SOURCES / "Q6N0W9.json").read_text())["sequence"]["value"]
    return data, [fasta[row["label"]] for row in data["rows"]]


def kernel(sequences):
    counts = [Counter(s[i:i + 3] for i in range(len(s) - 2)) for s in sequences]
    norms = [sqrt(sum(v*v for v in c.values())) for c in counts]
    return [[sum(v*b[word] for word, v in a.items()) / (na*nb)
             for b, nb in zip(counts, norms)] for a, na in zip(counts, norms)]


def cholesky_solve(matrix, rhs):
    lower = []
    for i, row in enumerate(matrix):
        factors = []
        for j in range(i):
            residual = row[j] - sum(factors[k]*lower[j][k] for k in range(j))
            factors.append(residual / lower[j][j])
        factors.append(sqrt(row[i] - sum(v*v for v in factors)))
        lower.append(factors)
    forward = []
    for row, value in zip(lower, rhs):
        i = len(forward)
        forward.append((value - sum(row[j]*forward[j] for j in range(i))) / row[i])
    solution = [0.] * len(rhs)
    for i in reversed(range(len(rhs))):
        solution[i] = (forward[i] - sum(lower[j][i]*solution[j]
                       for j in range(i+1, len(rhs)))) / lower[i][i]
    return solution


def fit(gram, y, groups, train, held, alpha, balanced):
    counts = Counter(groups[i] for i in train)
    weights = [1/counts[groups[i]] if balanced else 1. for i in train]
    center = sum(w*y[i] for w, i in zip(weights, train)) / sum(weights)
    system = [[gram[i][j] + (alpha/w if i == j else 0.)
               for j in train] for i, w in zip(train, weights)]
    coeff = cholesky_solve(system, [y[i]-center for i in train])
    return [round(center + sum(gram[i][j]*v for j, v in zip(train, coeff)), 10)
            for i in held]


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


def knn_predictions(y, groups, identity, train, held, model):
    donors = defaultdict(list)
    for index in train:
        donors[index if model == "knn" else groups[index]].append(index)
    output = {k: [] for k in NEIGHBORS}
    for query in held:
        votes = []
        for members in donors.values():
            nearest = min(members, key=lambda i: (-identity[query][i], i))
            label = median(y[i] for i in members) if model == "knn_median" else y[nearest]
            votes.append((identity[query][nearest], nearest, label))
        votes.sort(key=lambda v: (-v[0], v[1]))
        for k in NEIGHBORS:
            cutoff = votes[min(k, len(votes)) - 1][0]
            output[k].append(round(fmean(label for sim, _, label in votes if sim >= cutoff), 10))
    return output


def run():
    data, sequences = load()
    gram = kernel(sequences)
    y = [row["y"] for row in data["rows"]]
    groups = [row["group"] for row in data["rows"]]
    genus = [row["species"].casefold().split()[0] for row in data["rows"]]
    identity = data["identity"]
    components = data["component_id"]
    all_rows = list(range(len(y)))

    def weights(indices):
        counts = Counter(groups[i] for i in indices)
        return [1 / counts[groups[i]] for i in indices]

    def retain(pool, candidate, held):
        if pool == "direct":
            return [i for i in candidate if all(genus[i] != genus[j] and
                                               identity[i][j] < .95 for j in held)]
        blocked = {components[i] for i in held}
        return [i for i in candidate if components[i] not in blocked]

    def candidates(model, train, held):
        if model.startswith("ridge"):
            return {alpha: fit(gram, y, groups, train, held, alpha, model == "ridge_balanced")
                    for alpha in ALPHAS}
        return knn_predictions(y, groups, identity, train, held, model)

    result = {"n_rows": len(y), "n_groups": len(set(groups)), "pools": {}}
    for pool, models in MODELS.items():
        folds = data[pool]
        output = result["pools"][pool] = {}
        for model in models:
            grid = ALPHAS if model.startswith("ridge") else NEIGHBORS
            oof = {unit: [None] * len(y) for unit in ("row", "group", "error")}
            fixed = {setting: [None] * len(y) for setting in grid}
            history = []
            for fold in folds:
                train, test = fold["train"], fold["test"]
                assert train == retain(pool, all_rows, test)
                inner = {setting: {} for setting in grid}
                for held in fold["validation"]:
                    for setting, vector in candidates(model, retain(pool, train, held), held).items():
                        inner[setting].update(zip(held, vector))
                final = candidates(model, train, test)
                observed = [y[i] for i in train]
                detail = {}
                for unit in oof:
                    mass = [1] * len(train) if unit == "row" else weights(train)
                    metric = "rmse" if unit == "error" else "rho"
                    scores = [score(observed, [inner[s][i] for i in train], mass)[metric]
                              for s in grid]
                    if unit == "error":
                        winner = min(range(len(grid)), key=scores.__getitem__)
                    else:
                        winner = max(range(len(grid)),
                                     key=lambda j: scores[j] if scores[j] is not None else -2)
                    chosen = grid[winner]
                    detail[unit] = {"choice": chosen, "inner_" + metric: scores}
                    for i, prediction in zip(test, final[chosen]):
                        oof[unit][i] = prediction
                for setting, vector in final.items():
                    for i, prediction in zip(test, vector):
                        fixed[setting][i] = prediction
                history.append(detail)
            mass = weights(all_rows)
            scores = {unit: score(y, vector, mass) for unit, vector in oof.items()}
            output[model] = {
                "grid": grid,
                "folds": history,
                "outer": scores,
                "delta": {metric: scores["error"][metric] - scores["group"][metric]
                          for metric in ("rho", "rmse")},
                "predictions": oof,
                "outer_by_fold": [
                    {unit: score([y[i] for i in fold["test"]],
                                 [vector[i] for i in fold["test"]], weights(fold["test"]))
                     for unit, vector in {**oof, "fixed": fixed[3 if model.startswith("knn") else .01]}.items()}
                    for fold in folds
                ],
                "fixed_outer": {str(setting): score(y, vector, mass)
                                for setting, vector in fixed.items()},
            }
    references = json.loads((SOURCES / "outer-reference.json").read_text())
    result["fixed_references"] = {name: score(y, pred, weights(all_rows))
                                  for name, pred in references["fixed"].items()}
    return result


if __name__ == "__main__":
    (ROOT / "results.json").write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
