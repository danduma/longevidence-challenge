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


def feature(name, position, row):
    enzyme, residue, _, _ = row
    residue_class = RESIDUE_CLASSES[residue]
    return {
        "mean": None, "enzyme": enzyme, "residue": residue,
        "chemical": residue_class == RESIDUE_CLASSES[REFERENCE[position]],
        "chemical_group": residue_class,
    }[name]


def error(pred, truth):
    return sqrt(average((x - row[-1]) ** 2 for x, row in zip(pred, truth, strict=True)))



PRIORS = dict(zip(("P6", "P5", "P4", "P1"), (0.2, 0.8, 0.2, 0.1)))


def predict(data, training, held, model, penalty, blend, weight):
    train = [(pos, row) for pos in training
             for row in data[pos] for _ in range(weight if pos in PRIORS else 4)]
    if model == "linear":
        slope, intercept = linear_regression([r[2] for _, r in train],
                                             [r[3] for _, r in train])

        def evaluate(row):
            return slope * row[2] + intercept
    elif model == "isotonic":
        curve = increasing_curve([r for _, r in train])

        def evaluate(row):
            return curve(row[2])
    else:
        center = average(r[3] for _, r in train)
        sums = Counter()
        counts = Counter()
        for pos, row in train:
            key = feature(model, pos, row)
            sums[key] += row[3]
            counts[key] += 1

        def evaluate(row):
            key = feature(model, held, row)
            return ((sums[key] + penalty * center) / (counts[key] + penalty)
                    if counts[key] else center)
    return [estimate + blend * (PRIORS.get(held, estimate) - estimate)
            for estimate in map(evaluate, data[held])]


def main():
    data = samples()
    output = json.loads((ROOT / "results.json").read_text())
    check(sum(map(len, data.values())), output["n"])
    check(sum(row[-1] < 0 for rows in data.values() for row in rows), output["negative"])
    assert set(data) == set(output["folds"])
    scores = defaultdict(list)
    for outer in REFERENCE:
        training = tuple(p for p in REFERENCE if p != outer)
        for family in ("baseline", "linear", "isotonic"):
            models = FEATURES if family == "baseline" else (family,)
            inner_loss = {}
            for model in models:
                penalties = (0, 80) if family == "baseline" and model != "mean" else (0,)
                for penalty in penalties:
                    for weight in (4, 2, 1):
                        curves = {inner: predict(data, tuple(p for p in training if p != inner),
                                                  inner, model, penalty, 0, weight)
                                  for inner in training}
                        for blend in (0, 0.25, 0.5, 0.75):
                            validation = []
                            for inner, base in curves.items():
                                adjusted = [x + blend * (PRIORS.get(inner, x) - x) for x in base]
                                validation.append(error(adjusted, data[inner]))
                            inner_loss[model, penalty, blend, weight] = average(validation)

            for regime in ("equal", "selected", "fixed"):
                eligible = (c for c in inner_loss
                            if (regime != "equal" or c[3] == 4)
                            and (regime != "fixed" or c[2] == 0.75))
                model, penalty, blend, weight = min(eligible, key=inner_loss.__getitem__)
                name = f"{family}_{regime}"
                recorded = output["folds"][outer][name]
                check(f"{model}:{penalty // 4}:{blend}:{weight}", recorded["choice"])
                outcome = error(predict(data, training, outer, model, penalty, blend, weight),
                                data[outer])
                check(outcome, recorded["rmse"])
                scores[name].append(outcome)

    summary = {key: average(values) for key, values in scores.items()}
    assert summary.keys() == output["mean_fold_rmse"].keys()
    for key, value in summary.items():
        check(value, output["mean_fold_rmse"][key])
    for reporter in ("linear", "isotonic"):
        assert summary[reporter + "_equal"] < summary["baseline_equal"]
        assert summary[reporter + "_selected"] < summary[reporter + "_equal"]
        assert summary[reporter + "_fixed"] < summary["baseline_fixed"]
        assert summary["baseline_fixed"] < summary[reporter + "_selected"]
    print("Verified signed targets, nested selections, per-position errors and headline comparisons.")


if __name__ == "__main__":
    main()
