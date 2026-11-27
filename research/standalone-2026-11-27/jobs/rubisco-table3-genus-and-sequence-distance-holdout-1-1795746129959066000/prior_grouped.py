"""Rebuild Table 3 and compare nested CV designs on the repaired features."""
from collections import Counter
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy.linalg import solve
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, KFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"
FOCAL = "Rhodopseudomonas_palustris"
DIMENSIONS = (0, 8, 16, 32, 64)
ALPHAS = (0.1, 1, 10, 100)
KERNEL_ALPHAS = (0.001, 0.01, *ALPHAS)


def deposited(suffix):
    with ZipFile(SOURCES / "version3.zip") as archive:
        member, = [name for name in archive.namelist() if name.endswith("/" + suffix)]
        return archive.read(member)


def table(suffix):
    return pd.read_csv(BytesIO(deposited(suffix)))


def inputs():
    meta = pd.read_csv(SOURCES / "71534b19c-flamholz_rubisco_wt_meta.csv")
    species = dict(zip(meta.label, meta.species))
    records = []
    for record in deposited("flamholz/flamholz_rubisco_wt.fasta").decode().split(">")[1:]:
        header, *lines = record.splitlines()
        label = header.split()[0]
        records.append((label, species[label], "".join(lines)))
    fasta = pd.DataFrame(records, columns=["label", "species", "sequence"])
    observations = table("inputs/flamholz_dataset_S1.csv")
    medians = observations.query("primary == 1 and mutant == 0").groupby("species")[["S", "KC"]].median()
    rows = medians.reset_index().merge(fasta, on="species")
    repaired_sequence = json.loads((SOURCES / "Q6N0W9.json").read_text())["sequence"]["value"]
    assert (rows.label == FOCAL).sum() == 1
    rows.loc[rows.label == FOCAL, "sequence"] = repaired_sequence
    with np.load(SOURCES / "t12.npz") as features:
        original = np.stack([features[label] for label in rows.label]).astype(float)
        repaired = original.copy()
        repaired[rows.label == FOCAL] = features["replacement"]
    return rows, original, repaired


def components(rows):
    parents = list(range(len(rows)))

    def root(index):
        while index != parents[index]:
            index = parents[index]
        return index

    seen = {}
    for index, row in enumerate(rows.itertuples(index=False)):
        for key in [("species", tuple(row.species.casefold().split()[:2])), ("sequence", row.sequence)]:
            match = seen.setdefault(key, index)
            left, right = sorted((root(index), root(match)))
            parents[right] = left
    return np.array([root(index) for index in range(len(rows))])


def trimer_kernel(sequences):
    counts = [Counter(seq[i:i + 3] for i in range(len(seq) - 2)) for seq in sequences]
    vocabulary = sorted(set().union(*counts))
    matrix = np.array([[row[word] for word in vocabulary] for row in counts], dtype=float)
    matrix /= np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix @ matrix.T


def folds(groups, k, mode, seed):
    if mode == "random":
        return list(KFold(k, shuffle=True, random_state=seed).split(groups))
    splitter = GroupKFold(k, shuffle=bool(seed), random_state=seed if seed else None)
    result = list(splitter.split(groups, groups=groups))
    for train, test in result:
        assert not set(groups[train]) & set(groups[test])
    return result


def score(y, pred):
    value = spearmanr(y, pred).statistic
    return float(value) if np.isfinite(value) else -np.inf


def embedding_fit(x, y, train, test, dimension, alpha):
    scaler = StandardScaler()
    x_train = scaler.fit_transform(x[train])
    x_test = scaler.transform(x[test])
    if dimension:
        pca = PCA(n_components=min(dimension, len(train) - 1, x.shape[1]),
                  random_state=0, svd_solver="randomized")
        x_train = pca.fit_transform(x_train)
        x_test = pca.transform(x_test)
    return Ridge(alpha=alpha, random_state=0).fit(x_train, y[train]).predict(x_test)


def kernel_fit(kernel, y, train, test, alpha):
    mean = y[train].mean()
    gram = kernel[np.ix_(train, train)] + alpha * np.eye(len(train))
    weights = solve(gram, y[train] - mean, assume_a="pos")
    return np.round(mean + kernel[np.ix_(test, train)] @ weights, 10)


def exposure(groups, train, test):
    return int(np.isin(groups[test], groups[train]).sum())


def fit_nested(x, y, groups, kernel, mode, seed):
    outer = folds(groups, 5, mode, seed)
    predictions = {name: np.full(len(y), np.nan) for name in ("embedding", "mean", "trimer")}
    choices = {"embedding": [], "trimer": []}
    fold_details = []
    for train, test in outer:
        inner = [(train[a], train[b]) for a, b in folds(groups[train], 4, mode, seed)]
        learners = [("embedding", [(d, a) for d in DIMENSIONS for a in ALPHAS], embedding_fit, x)]
        if kernel is not None:
            learners.append(("trimer", [(a,) for a in KERNEL_ALPHAS], kernel_fit, kernel))
        for name, settings, learner, features in learners:
            best_score = -np.inf
            best = None
            for setting in settings:
                inner_predictions = np.full(len(y), np.nan)
                for a, b in inner:
                    inner_predictions[b] = learner(features, y, a, b, *setting)
                value = score(y[train], inner_predictions[train])
                if value > best_score:
                    best_score, best = value, setting
            predictions[name][test] = learner(features, y, train, test, *best)
            choices[name].append(list(best))
        predictions["mean"][test] = y[train].mean()
        fold_details.append({
            "test": test.tolist(), "exposed": exposure(groups, train, test),
            "inner_exposed": sum(exposure(groups, a, b) for a, b in inner),
        })
    if kernel is None:
        del predictions["trimer"]
    assert all(np.isfinite(pred).all() for pred in predictions.values())
    return predictions, choices, fold_details


def metrics(y, predictions, details):
    return {
        name: {
            "spearman": score(y, pred),
            "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
            "mean_within_fold_spearman": None if name == "mean" else float(
                np.mean([score(y[d["test"]], pred[d["test"]]) for d in details])
            ),
        }
        for name, pred in predictions.items()
    }


def run():
    rows, original, repaired = inputs()
    groups = components(rows)
    kernel = trimer_kernel(rows.sequence)
    output = {}
    for target in ("S", "KC"):
        mask = rows[target].to_numpy() > 0
        observed = table(f"tune_stepA/pred_{target}.csv")
        y = np.log10(rows.loc[mask, target].to_numpy())
        assert observed.species.tolist() == rows.loc[mask, "species"].tolist()
        np.testing.assert_allclose(y, observed.y_true, rtol=0, atol=1e-12)
        group = groups[mask]
        replay, choices, _ = fit_nested(original[mask], y, group, None, "random", 0)
        error = max(abs(replay["embedding"] - observed.y_pred))
        assert error < 2e-5
        archived = table(f"tune_stepA/chosen_params_{target}.csv")
        assert choices["embedding"] == archived[["pca_dim", "alpha"]].to_numpy().tolist()
        record = {"n": len(y), "groups": len(set(group)), "archive_max_error": error, "runs": {}}
        for mode in ("random", "grouped"):
            for seed in range(5):
                predictions, tuned, details = fit_nested(
                    repaired[mask], y, group, kernel[np.ix_(mask, mask)], mode, seed,
                )
                record["runs"][f"{mode}-{seed}"] = {
                    "metrics": metrics(y, predictions, details),
                    "choices": tuned, "folds": details,
                    "predictions": {name: pred.tolist() for name, pred in predictions.items()},
                }
        output[target] = record
    (ROOT / "results.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    run()
