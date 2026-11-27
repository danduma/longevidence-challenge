"""Reconstruct signed targets and independently check all headline comparisons."""
import csv
import json
from collections import Counter, defaultdict
from itertools import pairwise
from math import isclose, sqrt
from statistics import linear_regression
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "sources"
REFERENCE = dict(zip(("P6", "P5", "P4", "P3", "P2", "P1", "P1'"), "ENLYFQS"))


def read(name, delimiter=","):
    with (SOURCE / name).open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def average(values):
    values = list(values)
    return sum(values) / len(values)


def check(actual, expected):
    if isinstance(expected, (int, float)):
        assert isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12), (actual, expected)
    else:
        assert actual == expected, (actual, expected)


def samples():
    raw = {(r["TEVp"], r["pos"], r["TEVs_AA"]): float(r["viniabs"])
           for r in read("data_viniabs_Luzi.txt", "\t")}
    relative = {(r["TEVp"], r["pos"], r["TEVs_AA"]): float(r["vinirel"])
                for r in read("protease_vitro_relative.tsv", "\t")}
    assert raw.keys() == relative.keys() and len(raw) == 840
    normalized = {key: slope / raw[key[0], key[1], REFERENCE[key[1]]]
                  for key, slope in raw.items()}
    for key in normalized:
        check(relative[key], normalized[key])
    data = defaultdict(list)
    keys = set()
    canon = "".join(REFERENCE.values())
    for row in read("biochemical_common_keys.csv"):
        enzyme, motif = row["TEVp"], row["TEVs"]
        assert (enzyme, motif) not in keys and len(motif) == 7
        keys.add((enzyme, motif))
        assert min(int(row[k]) for k in ("min_raw", "min_corrected")) >= 100
        if motif == canon:
            check(float(row["vinirel"]), 1.0)
            continue
        (position, residue), = [(position, aa) for position, aa in zip(REFERENCE, motif)
                                if aa != REFERENCE[position]]
        y = normalized[enzyme, position, residue]
        check(float(row["vinirel"]), y)
        data[position].append((enzyme, residue, float(row["AUC_corrected"]), y))
    assert len(keys) == 786 and sum(map(len, data.values())) == 780
    assert set(data) == set(REFERENCE) and len({p for p, _ in keys}) == 6
    return data


def increasing_curve(training):
    """Slopes of the lower convex hull of cumulative count and response."""
    grouped = defaultdict(list)
    for _, _, auc, y in training:
        grouped[auc].append(y)
    aucs = sorted(grouped)
    hull = [(0, 0, 0.0)]
    count, total = 0, 0.0
    for index, auc in enumerate(aucs, 1):
        count += len(grouped[auc])
        total += sum(grouped[auc])
        while len(hull) > 1:
            _, a, b = hull[-2]
            _, c, d = hull[-1]
            if (d-b)*(count-c) < (total-d)*(c-a):
                break
            hull.pop()
        hull.append((index, count, total))
    levels = [(d-b)/(c-a) for (start, a, b), (stop, c, d) in pairwise(hull)
              for _ in range(start, stop)]

    def predict(x):
        if x <= aucs[0]:
            return levels[0]
        for i, (low, high) in enumerate(pairwise(aucs)):
            if x <= high:
                return ((high-x)*levels[i] + (x-low)*levels[i+1])/(high-low)
        return levels[-1]
    return predict


RESIDUE_CLASSES = dict(zip("DEAILMVFWYKRHNQSTCGP", "00111112223334444567"))
FEATURES = ("mean", "enzyme", "residue", "chemical", "chemical_group")


def anchor(position, residue):
    if position == "P5":
        return 1.0
    allowed = {"P1": "", "P4": "VMAIL", "P6": "ED"}
    if position in allowed and residue not in allowed[position]:
        return 0.0
    return None


def feature(name, position, row):
    enzyme, residue, _, _ = row
    residue_class = RESIDUE_CLASSES[residue]
    return {
        "mean": None, "enzyme": enzyme, "residue": residue,
        "chemical": residue_class == RESIDUE_CLASSES[REFERENCE[position]],
        "chemical_group": residue_class,
    }[name]


def predictions(data, fit_positions, holdout, subset, name, penalty=0, rules="four"):
    train = [(pos, row) for pos in fit_positions for row in data[pos]
             if subset == "all" or anchor(pos, row[1]) is None]
    ybar = average(row[-1] for _, row in train)
    if name == "linear":
        slope, intercept = linear_regression([row[2] for _, row in train],
                                             [row[-1] for _, row in train])
        evaluate = lambda row: slope * row[2] + intercept
    elif name == "isotonic":
        curve = increasing_curve([row for _, row in train])
        evaluate = lambda row: curve(row[2])
    else:
        groups = defaultdict(list)
        for pos, row in train:
            groups[feature(name, pos, row)].append(row[-1])

        def evaluate(row):
            values = groups[feature(name, holdout, row)]
            return (sum(values) + penalty*ybar)/(len(values) + penalty) if values else ybar

    result = []
    for row in data[holdout]:
        fixed = None
        if rules == "four" or (rules == "two" and holdout in ("P5", "P1")):
            fixed = anchor(holdout, row[1])
        result.append(evaluate(row) if fixed is None else fixed)
    return result


def error(pred, truth):
    return sqrt(average((x - row[-1]) ** 2 for x, row in zip(pred, truth, strict=True)))


def main():
    data = samples()
    output = json.loads((ROOT / "results.json").read_text())
    check(sum(row[-1] < 0 for group in data.values() for row in group), output["negative_targets"])
    assert output["negative_targets"] == 170
    assert set(output["folds"]) == set(REFERENCE)
    scores = defaultdict(dict)
    subsets = ("eligible", "all")
    families = {
        "baseline": [(subset, feature_name, penalty)
                     for subset in subsets for feature_name in FEATURES
                     for penalty in ((0,) if feature_name == "mean" else (0, 20))],
        "linear": [(subset, "linear", 0) for subset in subsets],
        "isotonic": [(subset, "isotonic", 0) for subset in subsets],
    }
    for outer in REFERENCE:
        test = data[outer]
        training = [pos for pos in REFERENCE if pos != outer]
        fold = output["folds"][outer]
        check(len(test), fold["n"])
        check(sum(anchor(outer, row[1]) is None for row in test), fold["eligible_test"])
        check(sum(anchor(pos, row[1]) is None for pos in training for row in data[pos]),
              fold["eligible_train"])
        estimates = {}
        for name in ("mean", "linear", "isotonic"):
            for prefix, rules, subset in (
                ("", "none", "all"), ("two_", "two", "all"),
                ("four_", "four", "all"), ("eligible_", "four", "eligible"),
            ):
                estimates[prefix + name] = predictions(data, training, outer, subset, name, rules=rules)
        estimates["four_chemical"] = predictions(data, training, outer, "all", "chemical")
        for family, candidates in families.items():
            def cross_validated_score(candidate):
                return average(error(predictions(data, [p for p in training if p != inner],
                                                 inner, *candidate), data[inner])
                               for inner in training)
            chosen = min(candidates, key=cross_validated_score)
            label = "_".join(map(str, chosen)) if family == "baseline" else chosen[0]
            check(label, fold["selected"][family])
            estimates["nested_" + family] = predictions(data, training, outer, *chosen)
        assert set(estimates) == set(fold["rmse"])
        for name, pred in estimates.items():
            scores[name][outer] = error(pred, test)
            check(scores[name][outer], fold["rmse"][name])
    averages = {name: average(values.values()) for name, values in scores.items()}
    assert set(averages) == set(output["mean_fold_rmse"]) == set(output["pooled_rmse"])
    for name, by_fold in scores.items():
        check(averages[name], output["mean_fold_rmse"][name])
        check(sqrt(sum(len(data[p])*v**2 for p, v in by_fold.items()) / 780),
              output["pooled_rmse"][name])

    for family in ("linear", "isotonic"):
        assert averages["four_mean"] < averages["four_" + family]
        assert averages["nested_" + family] < averages["eligible_mean"] < averages["nested_baseline"]
        assert {p for p in data if scores["nested_" + family][p] < scores["nested_baseline"][p]} == {
            "P6", "P4", "P3", "P2", "P1'",
        }
        assert {p for p in data if scores["nested_" + family][p] == scores["nested_baseline"][p]} == {
            "P5", "P1",
        }
    print("Verified 840 signed ratios, all outer scores, and nested model selection.")


if __name__ == "__main__":
    main()
