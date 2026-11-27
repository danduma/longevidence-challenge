"""Verify by reranking integer-replicated rows and averaging group MSEs."""
import csv
import hashlib
import json
import math
from collections import defaultdict, Counter
from pathlib import Path
from statistics import correlation, fmean

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def same(actual, expected):
    require(math.isclose(actual, expected, rel_tol=0, abs_tol=2e-12),
            f"{actual!r} differs from {expected!r}")


def ordinary_ranks(values):
    counts = Counter(values)
    rank_by_value = {}
    first = 1
    for value in sorted(counts):
        count = counts[value]
        rank_by_value[value] = first + (count - 1) / 2
        first += count
    return [rank_by_value[value] for value in values]


def independently_score(indices, rows, predictions, equal_groups):
    blocks = defaultdict(list)
    for i in indices:
        blocks[rows[i]["group"] if equal_groups else i].append(i)
    copies = math.lcm(*(len(members) for members in blocks.values()))
    replicated = [
        i for members in blocks.values() for i in members
        for _ in range(copies // len(members))
    ]
    ranks_y = ordinary_ranks([rows[i]["y"] for i in replicated])
    ranks_p = ordinary_ranks([predictions[i] for i in replicated])
    rho = correlation(ranks_y, ranks_p) if min(len(set(ranks_y)), len(set(ranks_p))) > 1 else None
    rmse = math.sqrt(fmean(
        fmean((rows[i]["y"] - predictions[i])**2 for i in members)
        for members in blocks.values()
    ))
    return rho, rmse


def main():
    payload = (ROOT / "sources/frozen.json").read_bytes()
    require(hashlib.sha256(payload).hexdigest() ==
            "a166d4dcd44e91f5db1b888478590f20d20afbb663d7da3ace65c1e5bd81e44d",
            "Frozen predecessor extract changed")
    inputs = json.loads(payload)
    extras = (ROOT / "sources/comparators.json").read_bytes()
    require(hashlib.sha256(extras).hexdigest() ==
            "579f481a8ef3c7606392892b47d0f834dbd6d28946b74606d4af0f5e2a9b1fab",
            "Frozen comparator extract changed")
    comparators = json.loads(extras)
    require(not (inputs["predictions"].keys() & comparators["predictions"].keys()),
            "Comparator overwrites an existing vector")
    inputs["predictions"].update(comparators["predictions"])
    inputs["prior_pooled"].update(comparators["prior_pooled"])
    report = json.loads((ROOT / "results.json").read_text())
    rows, folds = inputs["rows"], inputs["tests"]
    with (ROOT / "sources/pred_KC.csv").open(newline="", encoding="utf-8") as stream:
        archive = list(csv.DictReader(stream))
    require(len(rows) == len(archive) == report["n_rows"], "Row count mismatch")
    for row, reference in zip(rows, archive):
        require(row["species"] == reference["species"], "Archive row order mismatch")
        same(row["y"], float(reference["y_true"]))
    all_indices = list(range(len(rows)))
    require(sorted(i for fold in folds for i in fold) == all_indices,
            "Folds do not partition the panel")
    pools = comparators["training_pools"]
    require(len(pools) == 66, "Missing eligible training pools")
    for name, pool in pools.items():
        require(len(pool) == len(folds), "Missing comparator fold")
        for train, direct, component, test in zip(pool, pools["direct"], pools["component"], folds):
            require(set(train) <= set(direct) and not (set(train) & set(test)),
                    f"Training exclusion failed for {name}")
            if name.startswith("thin_"):
                require(len(train) == len(component), "Matched pool has wrong size")
    require(all(set(c) < set(d) for c, d in zip(pools["component"], pools["direct"])),
            "Component pool should be strictly smaller in each fold")
    groups = Counter(row["group"] for row in rows)
    require(sum(len({rows[i]["group"] for i in f}) for f in folds) ==
            len(groups) == report["n_groups"], "Groups split across held-out folds")
    require({str(k): n for k, n in Counter(groups.values()).items()} ==
            report["group_sizes"], "Group histogram mismatch")

    scores = {}
    require(report["models"].keys() == inputs["predictions"].keys(), "Model roster differs")
    for name, vector in inputs["predictions"].items():
        require(len(vector) == len(rows) and all(map(math.isfinite, vector)),
                f"Invalid vector: {name}")
        scores[name] = {}
        for unit in ("row", "group"):
            weighting = unit == "group"
            shown = report["models"][name][unit]
            rho, rmse = independently_score(all_indices, rows, vector, weighting)
            scores[name][unit] = rho, rmse
            same(rho, shown["rho"])
            same(rmse, shown["rmse"])
            by_fold = [independently_score(f, rows, vector, weighting)[0] for f in folds]
            require(len(shown["fold_rho"]) == len(folds), "Missing fold score")
            for expected, actual in zip(shown["fold_rho"], by_fold):
                if actual is None:
                    require(expected is None, "Undefined fold score must be null")
                else:
                    same(actual, expected)
            defined = [value for value in by_fold if value is not None]
            if defined:
                same(fmean(defined), shown["mean_fold_rho"])
            else:
                require(shown["mean_fold_rho"] is None, "Undefined fold mean must be null")
        old = inputs["prior_pooled"][name]
        same(scores[name]["row"][0], old["rho"])
        if "rmse" in old:
            same(scores[name]["row"][1], old["rmse"])

    alternatives = ("median_group", "nearest_group", "row_knn", "embedding")
    gaps = {}
    for unit in ("row", "group"):
        ridge, baseline = (scores[n][unit] for n in ("unweighted", "median_group"))
        for name in alternatives:
            require(scores[name][unit][0] <= baseline[0] and
                    scores[name][unit][1] >= baseline[1], "Unexpected strongest alternative")
        gaps[unit] = ridge[0] - baseline[0], baseline[1] - ridge[1]
    require(0 < gaps["group"][0] < gaps["row"][0], "Rank lead did not survive and shrink")
    require(0 < gaps["row"][1] < gaps["group"][1], "Error advantage did not grow")
    fixed_knn = [name for name in scores if name.startswith("fixed_")]
    for unit, expected_rank_winner in (("row", "fixed_rows_1"),
                                       ("group", "fixed_nearest_group_1")):
        best = max(fixed_knn, key=lambda name: scores[name][unit][0])
        require(best == expected_rank_winner, "Unexpected fixed-k rank envelope")
        require(scores["unweighted"][unit][0] > scores[best][unit][0] and
                scores["unweighted"][unit][1] < min(scores[n][unit][1] for n in fixed_knn),
                "A fixed-k comparator overtook ordinary ridge")
    ordinary = scores["unweighted"]
    fixed = scores["balanced_fixed_alpha"]
    require(fixed["row"][0] < ordinary["row"][0] < ordinary["group"][0] < fixed["group"][0],
            "Weighted rank ordering did not flip between trimer fits")
    require(fixed["group"][1] < ordinary["group"][1], "Fixed-alpha RMSE comparison differs")
    primary = [name for name in scores if not name.startswith("thin_") and
               ("fixed_" not in name or name == "balanced_fixed_alpha")]
    component_fixed = [name for name in scores if name.startswith("component_fixed_")]
    matched = [name for name in scores if name.startswith("thin_")]
    for unit in ("row", "group"):
        nested = scores["component_identity_knn"][unit]
        ridge = scores["unweighted"][unit]
        require(nested[0] > ridge[0] and nested[1] > ridge[1],
                "Nested rank and RMSE trade-off failed")
        require(max(primary, key=lambda n: scores[n][unit][0]) == "component_identity_knn",
                "Unexpected rank leader among full-pool models")
        require(min(primary, key=lambda n: scores[n][unit][1]) == "balanced_fixed_alpha",
                "Unexpected RMSE leader among full-pool models")
        require(nested[0] > max(scores[name][unit][0] for name in matched),
                "A matched draw exceeded component nested KNN rank")
        require(scores["component_fixed_exact_1"][unit][0] ==
                max(scores[name][unit][0] for name in component_fixed),
                "Unexpected fixed component rank leader")
        require(scores["component_fixed_ties_3"][unit][1] ==
                min(scores[name][unit][1] for name in component_fixed) <
                min(scores[name][unit][1] for name in primary + matched + fixed_knn),
                "Fixed component RMSE sensitivity failed")
    print("Verified all scores, frozen row-weighted replay, and the paired ranking comparisons.")


if __name__ == "__main__":
    main()
