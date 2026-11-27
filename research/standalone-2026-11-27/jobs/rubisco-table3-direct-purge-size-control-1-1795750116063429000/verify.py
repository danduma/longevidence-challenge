"""Recompute fold assignments, eligibility, all nested fits and rank contrasts."""
import json
import random
from pathlib import Path

import numpy as np
from sklearn.feature_extraction import DictVectorizer
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import normalize

from independent_helpers import Counter, alignment_matrix, equal, graph, rank, refit
import prior_grouped as base

ROOT = Path(__file__).resolve().parent
GRID = {
    "embedding": [(d, a) for d in (0, 8, 16, 32, 64) for a in (.1, 1, 10, 100)],
    "trimer": [(a,) for a in (.001, .01, .1, 1, 10, 100)],
    "identity_knn": [(k,) for k in (1, 3, 5, 10)],
}
MODELS = (*GRID, "mean")


def close(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-8)


def pools(candidate, held, labels, forbidden, kind):
    if kind == "direct":
        return np.array([i for i in candidate if all(not forbidden[i, j] for j in held)])
    other = set(labels[kind][held])
    return np.array([i for i in candidate if labels[kind][i] not in other])


def recompute(y, inputs, labels, forbidden, identity, record, kind, repeat=None):
    split = list(GroupKFold(5).split(y, groups=labels["exact"]))
    assert len(record["folds"]) == len(split)
    predictions = {name: np.full(len(y), np.nan) for name in MODELS}
    for fold, (candidate, test) in enumerate(split):
        saved = record["folds"][fold]
        train = pools(candidate, test, labels, forbidden, kind)
        component = pools(candidate, test, labels, forbidden, "component")
        if kind == "direct":
            assert set(component) <= set(train)
        if repeat is not None:
            sampler = random.Random(20261127 + 100 * repeat + fold)
            train = np.sort(sampler.sample(list(train), len(component)))
        assert saved["test"] == list(test)
        assert saved["train"] == list(train)
        close(saved["nearest"], identity[np.ix_(test, train)].max(axis=1))
        inner = []
        for a, b in GroupKFold(4).split(y[train], groups=labels["exact"][train]):
            fit = pools(train[a], train[b], labels, forbidden, kind)
            inner.append((fit, train[b]))
        assert saved["inner_train_counts"] == [len(a) for a, _ in inner]
        for a, b in [(train, test), *inner]:
            assert len(a) >= 2 and set(a).isdisjoint(b)
            assert set(labels["exact"][a]).isdisjoint(labels["exact"][b])
            if kind != "exact":
                assert not forbidden[np.ix_(a, b)].any()
        for model, grid in GRID.items():
            order = np.concatenate([b for _, b in inner])
            assert set(order) == set(train) and len(order) == len(train)
            scores = []
            for parameter in grid:
                oof = np.concatenate([refit(model, inputs[model], y, a, b, parameter)
                                      for a, b in inner])
                rho = rank(y[order], oof)
                scores.append(rho if rho is not None else -np.inf)
            # Spearman implementations can differ at floating-point ties.
            selected = tuple(saved["chosen"][model])
            assert scores[grid.index(selected)] >= max(scores) - 1e-13
            predictions[model][test] = refit(model, inputs[model], y, train, test, selected)
        predictions["mean"][test] = y[train].mean()
    metrics = {}
    for model, pred in predictions.items():
        close(pred, record["predictions"][model])
        fold_ranks = [rank(y[test], pred[test]) for _, test in split]
        metrics[model] = {
            "spearman": rank(y, pred),
            "rmse": np.linalg.norm(y - pred) / np.sqrt(len(y)),
            "mean_within_fold_spearman": None if model == "mean" else np.mean(fold_ranks),
        }
        for metric, value in metrics[model].items():
            expected = record["metrics"][model][metric]
            if expected is None:
                assert value is None
            else:
                close(value, expected)
    return metrics


def main():
    for line in (ROOT / "sources.sha256").read_text().splitlines():
        from hashlib import sha256
        digest, path = line.split(maxsplit=1)
        assert sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
    rows, _, embedding = base.inputs()
    sequences = rows.sequence.tolist()
    identity = alignment_matrix(sequences)
    names = [s.casefold().split() for s in rows.species]
    exact = graph(equal([" ".join(s[:2]) for s in names]) | equal(sequences))
    forbidden = equal([s[0] for s in names]) | (identity >= .95)
    component = graph(forbidden | equal(exact))
    counts = DictVectorizer().fit_transform(
        [Counter(s[i:i+3] for i in range(len(s)-2)) for s in sequences])
    normalized = normalize(counts)
    trimer = (normalized @ normalized.T).toarray()
    outcome = json.loads((ROOT / "results.json").read_text())
    prior = json.loads((ROOT / "sources/prior-component-results.json").read_text())
    for target in ("S", "KC"):
        mask = rows[target].to_numpy() > 0
        ix = np.ix_(mask, mask)
        y = np.log10(rows.loc[mask, target].to_numpy())
        archived = base.table(f"tune_stepA/pred_{target}.csv")
        assert archived.species.tolist() == rows.loc[mask, "species"].tolist()
        close(y, archived.y_true)
        args = (y, dict(embedding=embedding[mask], trimer=trimer[ix], identity_knn=identity[ix]),
                dict(exact=exact[mask], component=component[mask]), forbidden[ix], identity[ix])
        result = outcome[target]
        for arm in ("exact", "component"):
            previous = prior["targets"][target]["genus95" if arm == "component" else arm]
            for f, old in zip(result[arm]["folds"], previous["folds"], strict=True):
                assert f["test"] == old["test"] and f["train"] == old["train"]
            for model in MODELS:
                close(result[arm]["predictions"][model], previous["predictions"][model])
        scores = {arm: recompute(*args, result[arm], arm)
                  for arm in ("exact", "component", "direct")}
        assert len(result["thin"]) == 64
        draws = [recompute(*args, draw, "direct", repeat)
                 for repeat, draw in enumerate(result["thin"])]
        for model in MODELS:
            for metric, summary in result["thin_summary"][model].items():
                values = np.array([draw[model][metric] for draw in draws])
                reference = scores["component"][model][metric]
                close(np.quantile(values, [.05, .5, .95]), summary["quantiles_05_50_95"])
                assert np.sum(values > reference + 1e-9) == summary["higher_than_component"]
                assert np.sum(values < reference - 1e-9) == summary["lower_than_component"]
        print(f"{target}: verified three fixed arms and 64 matched draws")
    print("Verified pairwise eligibility, matched counts, nested fits, predictions and rank contrasts.")


if __name__ == "__main__":
    main()
