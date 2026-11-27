"""Nested identity KNN with one vote per eligible donor group."""
import json
import math
from pathlib import Path
from statistics import correlation, fmean, median

ROOT = Path(__file__).resolve().parent
GRID = (1, 3, 5, 10)
MODES = ("rows", "nearest_group", "median_group")


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0] * len(order)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[start]] == values[order[end]]:
            end += 1
        for index in order[start:end]:
            result[index] = (start + end - 1) / 2
        start = end
    return result


def rho(y, p):
    if len(set(p)) == 1 or len(set(y)) == 1:
        return None
    return correlation(ranks(y), ranks(p))


def evaluate(data):
    rows, identity = data["rows"], data["identity"]
    y = [row["y"] for row in rows]
    group = [row["group"] for row in rows]
    genus = [row["species"].casefold().split()[0] for row in rows]
    tests = [fold["test"] for fold in data["direct"]]

    def metrics(p):
        return {
            "rho": rho(y, p),
            "rmse": math.sqrt(fmean((a - b)**2 for a, b in zip(y, p))),
            "fold_rho": [rho([y[i] for i in test], [p[i] for i in test]) for test in tests],
        }

    def votes(query, eligible, mode):
        members = {}
        for i in eligible:
            members.setdefault(i if mode == "rows" else group[i], []).append(i)
        votes = []
        for members_of_group in members.values():
            best = min(members_of_group, key=lambda i: (-identity[query][i], i))
            label = median(y[i] for i in members_of_group) if mode == "median_group" else y[best]
            votes.append((identity[query][best], best, label))
        return sorted(votes, key=lambda vote: (-vote[0], vote[1]))

    def fits(eligible, held, mode):
        predictions = {k: {} for k in GRID}
        for query in held:
            ordered = votes(query, eligible, mode)
            for k in GRID:
                cutoff = ordered[min(k, len(ordered)) - 1][0]
                predictions[k][query] = round(fmean(v[2] for v in ordered if v[0] >= cutoff), 10)
        return predictions

    result = {"trimer": metrics(data["trimer"]), "embedding": metrics(data["embedding"])}
    for pool in ("direct", "component"):
        pool_results = {}
        for mode in MODES:
            choices, inner_scores = [], []
            fixed = {k: {} for k in GRID}
            nested, row_k = {}, {}
            for fold in data[pool]:
                train = fold["train"]
                valid_predictions = {k: {} for k in GRID}
                for held in fold["validation"]:
                    if pool == "direct":
                        eligible = [i for i in train if all(
                            genus[i] != genus[j] and identity[i][j] < .95 for j in held)]
                    else:
                        component = data["component_id"]
                        blocked = {component[j] for j in held}
                        eligible = [i for i in train if component[i] not in blocked]
                    for k, vector in fits(eligible, held, mode).items():
                        valid_predictions[k].update(vector)
                scores = [rho([y[i] for i in train], [valid_predictions[k][i] for i in train])
                          for k in GRID]
                chosen = max(zip(GRID, scores),
                             key=lambda pair: -math.inf if pair[1] is None else pair[1])[0]
                choices.append(chosen)
                inner_scores.append(scores)
                final = fits(train, fold["test"], mode)
                nested.update(final[chosen])
                row_k.update(final[fold["chosen"]["identity_knn"][0]])
                for k in GRID:
                    fixed[k].update(final[k])
            vector = [nested[i] for i in range(len(y))]
            pool_results[mode] = {
                "k": choices, "inner_rho": inner_scores, "predictions": vector,
                "nested": metrics(vector),
                "at_row_selected_k": metrics([row_k[i] for i in range(len(y))]),
                "fixed": {str(k): metrics([fixed[k][i] for i in range(len(y))]) for k in GRID},
            }
        result[pool] = pool_results
    return result


if __name__ == "__main__":
    data = json.loads((ROOT / "sources/kc-inputs.json").read_text())
    (ROOT / "results.json").write_text(json.dumps(evaluate(data), indent=2, allow_nan=False) + "\n")
