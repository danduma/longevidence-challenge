"""Position-grouped audit of the frozen C10 matched table (standard library)."""
import csv
import json
from bisect import bisect_left
from collections import defaultdict
from math import sqrt
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parent
POSITIONS = ("P6", "P5", "P4", "P3", "P2", "P1", "P1'")
CANON = "ENLYFQS"


def load():
    with (ROOT / "sources/biochemical_common_keys.csv").open() as source:
        rows = []
        for row in csv.DictReader(source):
            motif = row["TEVs"]
            if motif == CANON:
                continue
            index, = [i for i, (a, b) in enumerate(zip(motif, CANON)) if a != b]
            rows.append((POSITIONS[index], row["TEVp"], motif[index],
                         float(row["AUC_corrected"]), float(row["vinirel"])))
        return rows


def rmse(rows, predict):
    return sqrt(mean((predict(row) - row[-1]) ** 2 for row in rows))


def isotonic(train):
    groups = defaultdict(list)
    for _, _, _, x, y in train:
        groups[x].append(y)
    knots = sorted(groups)
    blocks = []
    for i, x in enumerate(knots):
        blocks.append((i, i, sum(groups[x]), len(groups[x])))
        while len(blocks) > 1:
            left, right = blocks[-2:]
            if left[2] / left[3] <= right[2] / right[3]:
                break
            blocks[-2:] = [(left[0], right[1], left[2]+right[2], left[3]+right[3])]
    levels = [total / count for lo, hi, total, count in blocks for _ in range(lo, hi+1)]

    def predict(row):
        x = row[3]
        i = bisect_left(knots, x)
        if i == 0 or i == len(knots):
            return levels[min(i, len(knots)-1)]
        weight = (x - knots[i-1]) / (knots[i] - knots[i-1])
        return levels[i-1] + weight * (levels[i] - levels[i-1])

    return predict


def linear(train):
    x_mean = mean(row[3] for row in train)
    y_mean = mean(row[4] for row in train)
    slope = sum((row[3] - x_mean) * (row[4] - y_mean) for row in train)
    slope /= sum((row[3] - x_mean) ** 2 for row in train)
    return lambda row: y_mean + slope * (row[3] - x_mean)


# Fixed literature-derived residue classes. P5/P1 retain the earlier convention.
PREFERRED = {"P6": "DE", "P4": "AILMV"}


def fixed_prediction(row):
    position, _, residue, _, _ = row
    if position == "P5":
        return 1.0
    if position == "P1" or (position in PREFERRED and residue not in PREFERRED[position]):
        return 0.0
    return None


def context(model):
    def predict(row):
        fixed = fixed_prediction(row)
        return model(row) if fixed is None else fixed
    return predict


def pooled_mean(train):
    level = mean(row[-1] for row in train)
    return lambda row: level


CHEMICAL_GROUPS = ("DE", "AILMV", "FWY", "KRH", "NQST", "C", "G", "P")
GROUP = {aa: i for i, group in enumerate(CHEMICAL_GROUPS) for aa in group}
REFERENCE = dict(zip(POSITIONS, CANON))


def chemical_match(row):
    position, _, residue, _, _ = row
    return GROUP[REFERENCE[position]] == GROUP[residue]


FEATURES = {
    "mean": lambda row: None,
    "enzyme": lambda row: row[1],
    "residue": lambda row: row[2],
    "chemical": chemical_match,
    "chemical_group": lambda row: GROUP[row[2]],
}


def grouped_mean(train, feature, shrinkage=0):
    overall = mean(row[-1] for row in train)
    groups = defaultdict(list)
    for row in train:
        groups[feature(row)].append(row[-1])
    levels = {key: (sum(values) + shrinkage * overall) / (len(values) + shrinkage)
              for key, values in groups.items()}
    return lambda row: levels.get(feature(row), overall)


def select_model(training, candidates):
    losses = {name: [] for name in candidates}
    for position in POSITIONS:
        validation = [row for row in training if row[0] == position]
        if not validation:
            continue
        inner_train = [row for row in training if row[0] != position]
        for name, fit in candidates.items():
            losses[name].append(rmse(validation, fit(inner_train)))
    selected = min(candidates, key=lambda name: mean(losses[name]))
    return selected, candidates[selected](training)


def contextual_fit(fit, eligible_only):
    def train(rows):
        usable = [row for row in rows if not eligible_only or fixed_prediction(row) is None]
        return context(fit(usable))
    return train


def baseline_candidates():
    for subset in ("eligible", "all"):
        for name, feature in FEATURES.items():
            for shrinkage in ((0,) if name == "mean" else (0, 20)):
                def fit(rows, feature=feature, shrinkage=shrinkage):
                    return grouped_mean(rows, feature, shrinkage)
                yield f"{subset}_{name}_{shrinkage}", contextual_fit(fit, subset == "eligible")


def run(rows):
    folds = {}
    for position in POSITIONS:
        train = [row for row in rows if row[0] != position]
        test = [row for row in rows if row[0] == position]
        plain = pooled_mean(train)
        original = {"mean": plain, "linear": linear(train), "isotonic": isotonic(train)}
        models = dict(original)
        for name, model in original.items():
            models["two_" + name] = (
                (lambda row: 1.0) if position == "P5"
                else (lambda row: 0.0) if position == "P1"
                else model
            )
            models["four_" + name] = context(model)
        models["four_chemical"] = context(grouped_mean(train, chemical_match))
        eligible = [row for row in train if fixed_prediction(row) is None]
        for name, fit in (("mean", pooled_mean), ("linear", linear), ("isotonic", isotonic)):
            models["eligible_" + name] = context(fit(eligible))
        choices = {}
        for name, candidates in (
            ("baseline", dict(baseline_candidates())),
            ("linear", {subset: contextual_fit(linear, subset == "eligible")
                        for subset in ("eligible", "all")}),
            ("isotonic", {subset: contextual_fit(isotonic, subset == "eligible")
                          for subset in ("eligible", "all")}),
        ):
            choices[name], models["nested_" + name] = select_model(train, candidates)
        folds[position] = {"n": len(test), "selected": choices, "eligible_train": len(eligible),
                           "eligible_test": sum(fixed_prediction(row) is None for row in test),
                           "rmse": {name: rmse(test, model) for name, model in models.items()}}
    return {
        "folds": folds,
        "mean_fold_rmse": {name: mean(fold["rmse"][name] for fold in folds.values()) for name in models},
        "pooled_rmse": {name: sqrt(sum(fold["n"] * fold["rmse"][name] ** 2
                                       for fold in folds.values()) / len(rows)) for name in models},
        "negative_targets": sum(row[-1] < 0 for row in rows),
    }


if __name__ == "__main__":
    (ROOT / "results.json").write_text(json.dumps(run(load()), indent=2) + "\n")
