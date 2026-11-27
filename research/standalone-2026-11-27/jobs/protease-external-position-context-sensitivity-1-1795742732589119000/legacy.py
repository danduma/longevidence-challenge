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
