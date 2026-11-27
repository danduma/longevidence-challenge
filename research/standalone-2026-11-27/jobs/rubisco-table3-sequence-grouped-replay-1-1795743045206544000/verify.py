"""Independently refit the archived learner and check every saved prediction."""
import ast
from collections import Counter
import hashlib
import itertools
import json
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, KFold
from sklearn.preprocessing import StandardScaler

from analysis import ROOT, deposited, inputs, table


def archived_learner():
    # Use only the inspected pure fitting functions, never the archive CLI.
    tree = ast.parse(deposited("tune_flamholz_ridge_pca.py"))
    names = {"fit_predict_fold", "inner_oof_spearman", "nested_cv"}
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    namespace = globals().copy()
    exec(compile(tree, "Version3:tune_flamholz_ridge_pca.py", "exec"), namespace)
    return namespace["nested_cv"]


def verify_hashes():
    for line in (ROOT / "sources.sha256").read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest
    assert hashlib.md5((ROOT / "sources/version3.zip").read_bytes()).hexdigest() == (
        "621cd098a0b64908063c4a6cbfec5573"
    )


def graph_and_gram(rows):
    keys = [" ".join(s.lower().split()[:2]) for s in rows.species]
    sequences = list(rows.sequence)
    edges = np.equal.outer(keys, keys) | np.equal.outer(sequences, sequences)
    _, graph_groups = connected_components(edges, directed=False)
    counts = [Counter(map("".join, zip(s, s[1:], s[2:]))) for s in sequences]
    norms = [sum(n * n for n in count.values()) ** 0.5 for count in counts]
    gram = np.eye(len(rows))
    for i, j in itertools.combinations(range(len(rows)), 2):
        gram[i, j] = gram[j, i] = sum(
            count * counts[j][word] for word, count in counts[i].items()
        ) / (norms[i] * norms[j])
    return graph_groups, gram


def splitter(grouped, n, seed):
    if grouped:
        return GroupKFold(n_splits=n, shuffle=seed != 0, random_state=seed or None)
    return KFold(n_splits=n, shuffle=True, random_state=seed)


def trimer_predict(gram, y, train, test, regularization):
    center = np.mean(y[train])
    coeff = np.linalg.solve(gram[train][:, train] + regularization * np.eye(train.size),
                            y[train] - center)
    return np.round(gram[test][:, train] @ coeff + center, 10)


def verify():
    verify_hashes()
    nested = archived_learner()
    rows, original, repaired = inputs()
    group_ids, gram = graph_and_gram(rows)
    saved = json.loads((ROOT / "results.json").read_text())
    options = ([0, 8, 16, 32, 64], [0.1, 1, 10, 100], 0)
    for trait in ("S", "KC"):
        archived = table(f"tune_stepA/pred_{trait}.csv")
        keep = rows[trait].gt(0).to_numpy()
        assert archived.species.tolist() == list(rows.loc[keep, "species"])
        y = np.log10(rows.loc[keep, trait].to_numpy())
        np.testing.assert_allclose(y, archived.y_true, atol=1e-12, rtol=0)
        groups = group_ids[keep]
        local_gram = gram[keep][:, keep]
        target = saved[trait]
        assert target["n"] == len(y) and target["groups"] == np.unique(groups).size
        assert set(target["runs"]) == {f"{mode}-{seed}" for mode in ("random", "grouped") for seed in range(5)}
        original_pred, original_choices = nested(
            original[keep], y, None, splitter(False, 5, 0), splitter(False, 4, 0), *options,
        )
        archive_error = max(abs(original_pred - archived.y_pred))
        assert archive_error < 2e-5
        np.testing.assert_allclose(archive_error, target["archive_max_error"], rtol=0, atol=1e-8)
        old_choices = table(f"tune_stepA/chosen_params_{trait}.csv")
        np.testing.assert_array_equal(original_choices[["pca_dim", "alpha"]],
                                      old_choices[["pca_dim", "alpha"]])

        for key, record in target["runs"].items():
            mode, seed = key.split("-")
            grouped, seed = mode == "grouped", int(seed)
            outer = splitter(grouped, 5, seed)
            inner = splitter(grouped, 4, seed)
            split_groups = groups if grouped else None
            embedding, chosen = nested(repaired[keep], y, split_groups, outer, inner, *options)
            assert chosen[["pca_dim", "alpha"]].to_numpy().tolist() == record["choices"]["embedding"]
            means, trimers = np.full(y.size, np.nan), np.full(y.size, np.nan)
            outer_folds = list(outer.split(y, groups=split_groups))
            for index, (train, test) in enumerate(outer_folds):
                assert test.tolist() == record["folds"][index]["test"]
                inner_folds = [(train[a], train[b]) for a, b in inner.split(y[train], groups=None if split_groups is None else split_groups[train])]
                outer_exposure = sum(g in groups[train] for g in groups[test])
                inner_exposure = sum(g in groups[a] for a, b in inner_folds for g in groups[b])
                assert outer_exposure == record["folds"][index]["exposed"]
                assert inner_exposure == record["folds"][index]["inner_exposed"]
                if grouped:
                    assert outer_exposure == inner_exposure == 0
                ranked = []
                for alpha in (0.001, 0.01, 0.1, 1, 10, 100):
                    pred = np.full(y.size, np.nan)
                    for a, b in inner_folds:
                        pred[b] = trimer_predict(local_gram, y, a, b, alpha)
                    ranked.append((spearmanr(y[train], pred[train]).statistic, alpha))
                _, best = max(ranked, key=lambda item: item[0])
                assert [best] == record["choices"]["trimer"][index]
                trimers[test] = trimer_predict(local_gram, y, train, test, best)
                means[test] = np.mean(y[train])
            for name, pred in {"embedding": embedding, "mean": means, "trimer": trimers}.items():
                expected = record["metrics"][name]
                saved_pred = np.asarray(record["predictions"][name])
                np.testing.assert_allclose(pred, saved_pred, atol=1e-8, rtol=0)
                # The independent kernel solver can differ at a rounding boundary
                # for tied scores. Verify the exact published rank from the
                # validated predictions and bound the independent rank difference.
                rank = spearmanr(y, saved_pred).statistic
                np.testing.assert_allclose(rank, expected["spearman"], atol=1e-10, rtol=0)
                tolerance = 0.001 if name == "trimer" else 1e-10
                np.testing.assert_allclose(spearmanr(y, pred).statistic, rank, atol=tolerance, rtol=0,
                                           err_msg=f"{trait} {key} {name}")
                np.testing.assert_allclose(np.linalg.norm(y - pred) / np.sqrt(y.size), expected["rmse"],
                                           atol=1e-8, rtol=0)
                if name == "mean":
                    assert expected["mean_within_fold_spearman"] is None
                else:
                    within = np.mean([spearmanr(y[b], saved_pred[b]).statistic for _, b in outer_folds])
                    np.testing.assert_allclose(within, expected["mean_within_fold_spearman"], atol=1e-10, rtol=0)
            if grouped:
                assert 0 < spearmanr(y, embedding).statistic < spearmanr(y, trimers).statistic
                assert np.linalg.norm(y - trimers) < np.linalg.norm(y - embedding) < np.linalg.norm(y - means)
    print("Verified both released rows and archived fits, all 20 repaired fits, split isolation and metrics.")


if __name__ == "__main__":
    verify()
