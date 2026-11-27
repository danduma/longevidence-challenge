"""Matched Table 3 train purges with nested train-only tuning."""
import json
from itertools import combinations

from Bio.Align import PairwiseAligner
import numpy as np

import prior_grouped as base

ARMS = {"exact": None, "genus": None, "genus99": .99, "genus95": .95}


def identities(sequences):
    unique, inverse = np.unique(sequences, return_inverse=True)
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                              open_gap_score=-5, extend_gap_score=-1)
    matrix = np.eye(len(unique))
    for left, right in combinations(range(len(unique)), 2):
        alignment = aligner.align(unique[left], unique[right])[0]
        matrix[left, right] = matrix[right, left] = (
            alignment.counts().identities / alignment.length
        )
    return matrix[np.ix_(inverse, inverse)]


def components(exact, genus, identity, arm):
    groups = exact.copy()
    cutoff = ARMS[arm]
    for left, right in combinations(range(len(exact)), 2):
        if arm != "exact" and (genus[left] == genus[right] or
                               (cutoff is not None and identity[left, right] >= cutoff)):
            groups[groups == groups[right]] = groups[left]
    return groups


def neighbors(identity, y, train, test, k):
    distances = identity[np.ix_(test, train)]
    k = min(k, len(train))
    kth = np.sort(distances, axis=1)[:, -k]
    selected = distances >= kth[:, None]
    return np.round((selected @ y[train]) / selected.sum(axis=1), 10)


def purge(train, test, groups):
    return train[~np.isin(groups[train], groups[test])]


def evaluate(x, y, exact, genus, groups, gram, identity, arm):
    predictions = {name: np.full(len(y), np.nan) for name in
                   ("embedding", "trimer", "identity_knn", "mean")}
    learners = [
        ("embedding", base.embedding_fit, x, [(d, a) for d in base.DIMENSIONS for a in base.ALPHAS]),
        ("trimer", base.kernel_fit, gram, [(a,) for a in base.KERNEL_ALPHAS]),
        ("identity_knn", neighbors, identity, [(k,) for k in (1, 3, 5, 10)]),
    ]
    details = []
    choices = {name: [] for name, *_ in learners}
    unavailable = False
    for candidate, test in base.folds(exact, 5, "grouped", 0):
        train = purge(candidate, test, groups)
        inner = []
        if len(set(exact[train])) >= 4:
            for a, b in base.folds(exact[train], 4, "grouped", 0):
                inner.append((purge(train[a], train[b], groups), train[b]))
        for a, b in [(train, test), *inner]:
            assert not np.isin(groups[a], groups[b]).any()
            if arm != "exact":
                assert not np.isin(genus[a], genus[b]).any()
            if ARMS[arm] is not None:
                assert not (identity[np.ix_(a, b)] >= ARMS[arm]).any()
        details.append({
            "test": test.tolist(), "train": train.tolist(), "unpurged": len(candidate),
            "inner": [[a.tolist(), b.tolist()] for a, b in inner],
            "nearest": identity[np.ix_(test, train)].max(axis=1).tolist() if len(train) else [],
            "same_genus_test_rows": int(np.isin(genus[test], genus[train]).sum()),
        })
        if not inner or min(len(a) for a, _ in inner) < 2:
            unavailable = True
            continue
        for name, fit, features, settings in learners:
            scores = []
            for setting in settings:
                fitted = np.full(len(y), np.nan)
                for a, b in inner:
                    fitted[b] = fit(features, y, a, b, *setting)
                scores.append(base.score(y[train], fitted[train]))
            choice = settings[np.argmax(scores)]
            predictions[name][test] = fit(features, y, train, test, *choice)
            choices[name].append(list(choice))
        predictions["mean"][test] = y[train].mean()
    result = {"folds": details, "status": "unavailable" if unavailable else "fit"}
    if unavailable:
        return result
    metrics = base.metrics(y, predictions, details)
    # Constant within-fold predictions have undefined Spearman.
    for name in metrics:
        for key, value in metrics[name].items():
            if value is not None and not np.isfinite(value):
                metrics[name][key] = None
    return result | {"metrics": metrics, "choices": choices,
                     "predictions": {k: v.tolist() for k, v in predictions.items()}}


def run():
    rows, _, x = base.inputs()
    identity = identities(rows.sequence.to_numpy())
    exact = base.components(rows)
    genus = np.array([name.casefold().split()[0] for name in rows.species])
    gram = base.trimer_kernel(rows.sequence)
    groupings = {arm: components(exact, genus, identity, arm) for arm in ARMS}
    result = {"components": {arm: sorted(np.unique(g, return_counts=True)[1].tolist(), reverse=True)
                             for arm, g in groupings.items()}, "targets": {}}
    for target in ("S", "KC"):
        mask = rows[target].to_numpy() > 0
        ix = np.ix_(mask, mask)
        y = np.log10(rows.loc[mask, target].to_numpy())
        archived = base.table(f"tune_stepA/pred_{target}.csv")
        assert archived.species.tolist() == rows.loc[mask, "species"].tolist()
        np.testing.assert_allclose(y, archived.y_true, atol=1e-12, rtol=0)
        result["targets"][target] = {
            arm: evaluate(x[mask], y, exact[mask], genus[mask], g[mask], gram[ix], identity[ix], arm)
            for arm, g in groupings.items()
        }
    (base.ROOT / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    run()
