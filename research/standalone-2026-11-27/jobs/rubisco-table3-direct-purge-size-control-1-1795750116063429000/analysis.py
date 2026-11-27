"""Fixed-fold Table 3 audit of direct exclusion versus transitive purging."""
import json
import random
from pathlib import Path

import numpy as np

import prior_grouped as base
from prior_component import components, identities, neighbors, purge

DRAWS = 64
MODELS = ("embedding", "trimer", "identity_knn", "mean")


def sample(eligible, count, repeat, fold):
    seed = 20261127 + 100 * repeat + fold
    return np.array(sorted(random.Random(seed).sample(eligible.tolist(), count)))


def evaluate(x, y, exact, genus, groups, gram, identity, arm, repeat=None):
    adjacency = (genus[:, None] == genus) | (identity >= .95)

    def retain(train, held):
        if arm == "exact":
            return purge(train, held, exact)
        if arm == "component":
            return purge(train, held, groups)
        return train[~adjacency[np.ix_(train, held)].any(axis=1)]

    prediction = {name: np.full(len(y), np.nan) for name in MODELS}
    learners = [
        ("embedding", base.embedding_fit, x,
         [(d, a) for d in base.DIMENSIONS for a in base.ALPHAS]),
        ("trimer", base.kernel_fit, gram, [(a,) for a in base.KERNEL_ALPHAS]),
        ("identity_knn", neighbors, identity, [(k,) for k in (1, 3, 5, 10)]),
    ]
    details = []
    for fold, (candidate, test) in enumerate(base.folds(exact, 5, "grouped", 0)):
        train = retain(candidate, test)
        if repeat is not None:
            train = sample(train, len(purge(candidate, test, groups)), repeat, fold)
        inner = [(retain(train[a], train[b]), train[b])
                 for a, b in base.folds(exact[train], 4, "grouped", 0)]
        for a, b in [(train, test), *inner]:
            assert len(a) >= 2 and not np.isin(exact[a], exact[b]).any()
            if arm != "exact":
                assert not adjacency[np.ix_(a, b)].any()
        chosen = {}
        for name, fit, features, settings in learners:
            scores = []
            for setting in settings:
                oof = np.full(len(y), np.nan)
                for a, b in inner:
                    oof[b] = fit(features, y, a, b, *setting)
                scores.append(base.score(y[train], oof[train]))
            choice = settings[np.argmax(scores)]
            chosen[name] = list(choice)
            prediction[name][test] = fit(features, y, train, test, *choice)
        prediction["mean"][test] = y[train].mean()
        details.append({"test": test.tolist(), "train": train.tolist(), "chosen": chosen,
                        "inner_train_counts": [len(a) for a, _ in inner],
                        "nearest": identity[np.ix_(test, train)].max(axis=1).tolist()})
    return {"folds": details, "metrics": base.metrics(y, prediction, details),
            "predictions": {name: pred.tolist() for name, pred in prediction.items()}}


def summarize(draws, component):
    result = {}
    for model in MODELS:
        result[model] = {}
        for metric, reference in component["metrics"][model].items():
            if reference is None:
                continue
            values = [draw["metrics"][model][metric] for draw in draws]
            result[model][metric] = {
                "quantiles_05_50_95": np.quantile(values, [.05, .5, .95]).tolist(),
                "higher_than_component": sum(value > reference + 1e-9 for value in values),
                "lower_than_component": sum(value < reference - 1e-9 for value in values),
            }
    return result


def run():
    rows, _, x = base.inputs()
    identity = identities(rows.sequence.to_numpy())
    exact = base.components(rows)
    genus = np.array([name.casefold().split()[0] for name in rows.species])
    groups = components(exact, genus, identity, "genus95")
    gram = base.trimer_kernel(rows.sequence)
    result = {}
    for target in ("S", "KC"):
        mask = rows[target].to_numpy() > 0
        ix = np.ix_(mask, mask)
        y = np.log10(rows.loc[mask, target].to_numpy())
        archived = base.table(f"tune_stepA/pred_{target}.csv")
        assert archived.species.tolist() == rows.loc[mask, "species"].tolist()
        np.testing.assert_allclose(y, archived.y_true, atol=1e-12, rtol=0)
        args = (x[mask], y, exact[mask], genus[mask], groups[mask], gram[ix], identity[ix])
        arms = {arm: evaluate(*args, arm) for arm in ("exact", "component", "direct")}
        draws = [evaluate(*args, "direct", repeat) for repeat in range(DRAWS)]
        arms["thin"] = draws
        arms["thin_summary"] = summarize(draws, arms["component"])
        result[target] = arms
    Path(__file__).with_name("results.json").write_text(
        json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n")


if __name__ == "__main__":
    run()
