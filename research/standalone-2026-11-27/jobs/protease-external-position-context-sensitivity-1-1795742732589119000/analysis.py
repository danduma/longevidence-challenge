"""Retrospective comparison of calibration errors with external position priors."""
import json
from functools import partial
from statistics import mean

from legacy import (ROOT, POSITIONS, FEATURES, load, grouped_mean, linear, isotonic,
                    select_model, rmse)

# Modeling anchors, not measured rates. See report.md for source mapping.
ANCHORS = {"P6": 0.2, "P5": 0.8, "P4": 0.2, "P1": 0.1}
STRENGTHS = (0, 0.25, 0.5, 0.75)
WEIGHTS = (4, 2, 1)


def contextual_fit(fit, strength, context_weight):
    def train(rows):
        # Repeated rows implement integer case weights for every estimator.
        model = fit([row for row in rows
                     for _ in range(context_weight if row[0] in ANCHORS else 4)])

        def predict(row):
            estimate = model(row)
            if row[0] in ANCHORS:
                return (1 - strength) * estimate + strength * ANCHORS[row[0]]
            return estimate
        return predict
    return train


def candidates(method, weights, strengths=STRENGTHS):
    if method == "baseline":
        fits = {(name, pseudo): partial(grouped_mean, feature=feature, shrinkage=4 * pseudo)
                for name, feature in FEATURES.items()
                for pseudo in ((0,) if name == "mean" else (0, 20))}
    else:
        fits = {(method, 0): {"linear": linear, "isotonic": isotonic}[method]}
    return {f"{name}:{pseudo}:{strength}:{weight}": contextual_fit(fit, strength, weight)
            for (name, pseudo), fit in fits.items()
            for strength in strengths for weight in weights}


def analyze(rows):
    folds = {}
    for position in POSITIONS:
        train = [row for row in rows if row[0] != position]
        test = [row for row in rows if row[0] == position]
        scores = {}
        for method in ("baseline", "linear", "isotonic"):
            for regime, weights, strengths in (("equal", (4,), STRENGTHS),
                                               ("selected", WEIGHTS, STRENGTHS),
                                               ("fixed", WEIGHTS, (0.75,))):
                choice, model = select_model(train, candidates(method, weights, strengths))
                scores[f"{method}_{regime}"] = {"choice": choice, "rmse": rmse(test, model)}
        folds[position] = scores
    return {
        "n": len(rows),
        "negative": sum(row[-1] < 0 for row in rows),
        "folds": folds,
        "mean_fold_rmse": {key: mean(scores[key]["rmse"] for scores in folds.values())
                           for key in scores},
    }


if __name__ == "__main__":
    (ROOT / "results.json").write_text(json.dumps(analyze(load()), indent=2) + "\n")
