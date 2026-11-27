"""Paired nested direct-pool KC refit, using standard-library Python."""
import json
from collections import Counter
from math import sqrt
from pathlib import Path
from statistics import fmean

from donor_knn import evaluate, rho

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"
ALPHAS = (.001, .01, .1, 1, 10, 100)


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


def run():
    data, sequences = load()
    gram = kernel(sequences)
    y = [r["y"] for r in data["rows"]]
    groups = [r["group"] for r in data["rows"]]
    genus = [r["species"].casefold().split()[0] for r in data["rows"]]
    tests = [fold["test"] for fold in data["direct"]]

    def purge(pool, held):
        return [i for i in pool if all(genus[i] != genus[j] and
                data["identity"][i][j] < .95 for j in held)]

    def summary(pred):
        return {"rho": rho(y, pred),
                "rmse": sqrt(fmean((a-b)**2 for a, b in zip(y, pred))),
                "fold_rho": [rho([y[i] for i in t], [pred[i] for i in t]) for t in tests],
                "predictions": pred}

    names = ("unweighted", "balanced", "balanced_fixed_alpha")
    predictions = {name: [None]*len(y) for name in names}
    details = []
    for fold in data["direct"]:
        train, test = fold["train"], fold["test"]
        assert train == purge(range(len(y)), test)
        inner = [(purge(train, held), held) for held in fold["validation"]]
        detail = {"train_sizes": [len(train), *[len(a) for a, _ in inner]],
                  "group_counts": [len({groups[i] for i in a})
                                   for a in [train, *[a for a, _ in inner]]]}
        for name, balanced in zip(names, (False, True)):
            scores = []
            for alpha in ALPHAS:
                pred = {}
                for a, b in inner:
                    pred.update(zip(b, fit(gram, y, groups, a, b, alpha, balanced)))
                scores.append(rho([y[i] for i in train], [pred[i] for i in train]))
            index = max(range(len(scores)), key=scores.__getitem__)
            detail[name] = {"alpha": ALPHAS[index], "inner_rho": scores}
        for name, source in zip(names, ("unweighted", "balanced", "unweighted")):
            values = fit(gram, y, groups, train, test, detail[source]["alpha"], name != "unweighted")
            for i, value in zip(test, values):
                predictions[name][i] = value
        details.append(detail)
    knn = evaluate(data)["direct"]
    result = {name: summary(vector) for name, vector in predictions.items()}
    for mode in ("nearest_group", "median_group"):
        result[mode] = {**knn[mode]["nested"], "k": knn[mode]["k"],
                        "predictions": knn[mode]["predictions"]}
    result["folds"] = details
    return result


if __name__ == "__main__":
    (ROOT / "results.json").write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
