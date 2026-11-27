"""Independently reconstruct the graphs, purged folds, tuning and OOF metrics."""
from collections import Counter
from functools import cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from Bio.Align import PairwiseAligner
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import rankdata
from sklearn.feature_extraction import DictVectorizer
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import normalize

import prior_grouped as replay


def alignment_matrix(sequences):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 2, -1
    aligner.open_gap_score, aligner.extend_gap_score = -5, -1

    @cache
    def measure(a, b):
        if a == b:
            return 1.
        path = aligner.align(a, b)[0].coordinates.T.tolist()
        matches, columns = 0, 0
        for (i, j), (u, v) in zip(path, path[1:]):
            columns += max(u - i, v - j)
            matches += sum(x == y for x, y in zip(a[i:u], b[j:v]))
        return matches / columns

    return np.array([[measure(*sorted((a, b))) for b in sequences] for a in sequences])


def graph(adjacency):
    return connected_components(csr_matrix(adjacency), directed=False)[1]


def equal(values):
    return np.equal.outer(values, values)


def rank(y, p):
    if np.ptp(y) == 0 or np.ptp(p) == 0:
        return None
    return float(np.corrcoef(rankdata(y), rankdata(p))[0, 1])


def refit(name, data, y, train, test, setting):
    if name == "embedding":
        return replay.embedding_fit(data, y, train, test, *setting)
    block = data[np.ix_(test, train)]
    target = y[train]
    if name == "trimer":
        normal = data[np.ix_(train, train)] + setting[0] * np.eye(len(train))
        pred = target.mean() + block @ np.linalg.solve(normal, target - target.mean())
    else:
        predictions = []
        for row in block:
            kth = sorted(row, reverse=True)[min(setting[0], len(row)) - 1]
            predictions.append(target[row >= kth].mean())
        pred = np.array(predictions)
    return pred.round(10)


def retained(candidate, held, adjacency, groups):
    rejected = {groups[i] for i in held}
    remaining = np.array([i for i in candidate if groups[i] not in rejected], dtype=int)
    assert not adjacency[np.ix_(remaining, held)].any()
    return remaining


def verify_arm(record, y, model_data, exact, groups, adjacent, names, identity, prior):
    grids = {
        "embedding": list(product((0, 8, 16, 32, 64), (.1, 1, 10, 100))),
        "trimer": [(a,) for a in (.001, .01, .1, 1, 10, 100)],
        "identity_knn": [(k,) for k in (1, 3, 5, 10)],
    }
    preds = {name: np.full(len(y), np.nan) for name in (*grids, "mean")}
    unavailable = False
    splits = list(GroupKFold(5).split(y, groups=exact))
    assert len(record["folds"]) == len(splits)
    for fold, (candidate, test) in enumerate(splits):
        saved = record["folds"][fold]
        assert prior["folds"][fold]["test"] == saved["test"] == test.tolist()
        train = retained(candidate, test, adjacent, groups)
        assert saved["train"] == train.tolist()
        assert saved["unpurged"] == len(candidate)
        nearest = identity[np.ix_(test, train)].max(axis=1) if len(train) else []
        np.testing.assert_allclose(saved["nearest"], nearest, rtol=0, atol=1e-15)
        assert saved["same_genus_test_rows"] == len([i for i in test if names[i] in names[train]])
        inner = []
        if len(set(exact[train])) >= 4:
            for a, b in GroupKFold(4).split(y[train], groups=exact[train]):
                inner.append((retained(train[a], train[b], adjacent, groups), train[b]))
        assert saved["inner"] == [[a.tolist(), b.tolist()] for a, b in inner]
        if not inner or min(len(a) for a, _ in inner) < 2:
            unavailable = True
            continue
        if record["status"] == "unavailable":
            continue
        for name, grid in grids.items():
            scores = []
            for setting in grid:
                inner_pred = np.full(len(y), np.nan)
                for a, b in inner:
                    inner_pred[b] = refit(name, model_data[name], y, a, b, setting)
                value = rank(y[train], inner_pred[train])
                scores.append(float("-inf") if value is None else value)
            # Numerical roundoff on ties can vary between the two rank implementations.
            chosen = tuple(record["choices"][name][fold])
            assert scores[grid.index(chosen)] >= max(scores) - 1e-12
            preds[name][test] = refit(name, model_data[name], y, train, test, chosen)
        preds["mean"][test] = y[train].mean()
    assert record["status"] == ("unavailable" if unavailable else "fit")
    if unavailable:
        assert "metrics" not in record
        return
    for name, pred in preds.items():
        np.testing.assert_allclose(record["predictions"][name], pred, rtol=0, atol=2e-8)
        # Rank archived rounded scores after independent refit asserts their values.
        # This preserves exact ties at the recorded rounding boundary.
        p = np.array(record["predictions"][name])
        recomputed = rank(y, pred)
        if recomputed is not None:
            np.testing.assert_allclose(record["metrics"][name]["spearman"], recomputed,
                                       rtol=0, atol=0.001)
        within = [rank(y[test], p[test]) for _, test in splits]
        values = {"spearman": rank(y, p),
                  "rmse": np.sqrt(np.mean((p - y) ** 2)),
                  "mean_within_fold_spearman": None if None in within else np.mean(within)}
        if name == "mean":
            values["mean_within_fold_spearman"] = None
        for key, actual in values.items():
            saved = record["metrics"][name][key]
            if actual is None:
                assert saved is None
            else:
                np.testing.assert_allclose(saved, actual, rtol=0, atol=1e-12)


def main():
    root = Path(__file__).resolve().parent
    for line in (root / "sources.sha256").read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        assert sha256((root / path).read_bytes()).hexdigest() == digest, path
    report = json.loads((root / "results.json").read_text())
    prior = json.loads((root / "sources/prior_grouped_results.json").read_text())
    rows, _, features = replay.inputs()
    seq = rows.sequence.tolist()
    names = np.array([r.casefold().split()[0] for r in rows.species])
    species = [" ".join(r.casefold().split()[:2]) for r in rows.species]
    identity = alignment_matrix(seq)
    exact_adj = equal(seq) | equal(species)
    genus_adj = exact_adj | equal(names)
    adjacency = {"exact": exact_adj, "genus": genus_adj,
                 "genus99": genus_adj | (identity >= .99),
                 "genus95": genus_adj | (identity >= .95)}
    trigrams = [Counter(s[i:i+3] for i in range(len(s)-2)) for s in seq]
    matrix = normalize(DictVectorizer().fit_transform(trigrams))
    kernel = (matrix @ matrix.T).toarray()
    exact = graph(exact_adj)
    for arm, adjacent in adjacency.items():
        groups = graph(adjacent)
        assert report["components"][arm] == sorted(np.bincount(groups).tolist(), reverse=True)
        for target in ("S", "KC"):
            mask = rows[target].to_numpy() > 0
            ix = np.ix_(mask, mask)
            archived = replay.table(f"tune_stepA/pred_{target}.csv")
            assert archived.species.tolist() == rows.loc[mask, "species"].tolist()
            y = np.log10(rows.loc[mask, target].to_numpy())
            np.testing.assert_allclose(y, archived.y_true, rtol=0, atol=1e-12)
            data = {"embedding": features[mask], "trimer": kernel[ix], "identity_knn": identity[ix]}
            old = prior[target]["runs"]["grouped-0"]
            record = report["targets"][target][arm]
            verify_arm(record, y, data, exact[mask], groups[mask], adjacent[ix],
                       names[mask], identity[ix], old)
            if arm == "exact":
                for name in ("embedding", "trimer", "mean"):
                    np.testing.assert_allclose(record["predictions"][name], old["predictions"][name],
                                               rtol=0, atol=2e-8)
    print("Verified archived target rows, fixed test folds, exclusion graphs, training coverage, nested fits and OOF metrics.")


if __name__ == "__main__":
    main()
