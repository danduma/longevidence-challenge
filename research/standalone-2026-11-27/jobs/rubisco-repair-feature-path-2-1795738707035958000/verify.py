"""Verify replay through independent row reconstruction and CV orchestration."""
import csv
import argparse
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from statistics import median
from zipfile import ZipFile

import numpy as np
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"
TARGETS = ("S", "KC", "vC", "KO", "vO", "KRuBP")
FOCAL = "Rhodopseudomonas_palustris"


def csv_rows(archive, filename):
    name, = (name for name in archive.namelist() if name.endswith("/" + filename))
    return list(csv.DictReader(archive.read(name).decode().splitlines()))


def archive_fasta(archive):
    path, = [path for path in archive.namelist()
             if path.endswith("/flamholz_rubisco_wt.fasta")]
    return archive.read(path)


def old_metadata():
    return list(csv.DictReader(
        (SOURCES / "71534b19c-flamholz_rubisco_wt_meta.csv").read_text().splitlines()
    ))


def dataset(archive, target, restore):
    observations = defaultdict(list)
    for row in csv_rows(archive, "inputs/flamholz_dataset_S1.csv"):
        if row["primary"] == "1" and row["mutant"] == "0":
            try:
                value = float(row[target])
            except ValueError:
                continue
            if math.isfinite(value):
                observations[row["species"]].append(value)
    metadata = old_metadata() if restore else csv_rows(archive, "flamholz_rubisco_wt_meta.csv")
    name_by_label = {row["label"]: row["species"] for row in metadata}
    mapping = defaultdict(list)
    for line in archive_fasta(archive).decode().splitlines():
        if line.startswith(">"):
            label = line[1:].split()[0]
            mapping[name_by_label[label]].append(label)
    return [(name, label, math.log10(value))
            for name, values in sorted(observations.items())
            if (value := median(values)) > 0
            for label in mapping[name]]


def predict(X, y, train, test, dim, alpha):
    scale = StandardScaler().fit(X[train])
    A, B = scale.transform(X[train]), scale.transform(X[test])
    if dim:
        pca = PCA(min(dim, A.shape[1], A.shape[0] - 1),
                  svd_solver="randomized", random_state=0)
        A = pca.fit_transform(A)
        B = pca.transform(B)
    return Ridge(alpha=alpha, random_state=0).fit(A, y[train]).predict(B)


def cross_validate(X, y, tune):
    oof = np.empty_like(y)
    choices = []
    for train, test in KFold(5, shuffle=True, random_state=0).split(X):
        dim, alpha = 16, 10.0
        if tune:
            inner_X, inner_y = X[train], y[train]
            splits = list(KFold(4, shuffle=True, random_state=0).split(inner_X))
            fitted = np.empty_like(inner_y)
            candidates = []
            for components in (0, 8, 16, 32, 64):
                for penalty in (0.1, 1, 10, 100):
                    for tr, te in splits:
                        fitted[te] = predict(inner_X, inner_y, tr, te, components, penalty)
                    candidates.append((spearmanr(fitted, inner_y).statistic, components, penalty))
            score, dim, alpha = max(candidates, key=lambda candidate: candidate[0])
            choices.append((dim, alpha, score))
        oof[test] = predict(X, y, train, test, dim, alpha)
    return oof, choices


def check(archive, vectors, branch, target, result):
    prefix = f"flamholz/{branch}/"
    archive_rows = csv_rows(archive, prefix + f"pred_{target}.csv")
    data = dataset(archive, target, restore=True)
    names, labels, y = zip(*data, strict=True)
    assert names == tuple(row["species"] for row in archive_rows)
    np.testing.assert_allclose(y, [float(row["y_true"]) for row in archive_rows],
                               atol=1e-12, rtol=0)
    assert result["literal_rows"] == len(dataset(archive, target, restore=False))
    assert result["restored_rows"] == len(y)
    y = np.asarray(y)
    X = np.stack([vectors[label] for label in labels]).astype(float)
    oof, choices = cross_validate(X, y, tune=branch == "tune_stepA")
    expected = np.array([float(row["y_pred"]) for row in archive_rows])
    max_error = np.max(np.abs(oof - expected))
    np.testing.assert_allclose(max_error, result["max_abs_archived_prediction_error"],
                               atol=1e-9, rtol=0)
    if choices:
        params = csv_rows(archive, prefix + f"chosen_params_{target}.csv")
        archive_params = [(int(row["pca_dim"]), float(row["alpha"])) for row in params]
        assert result["chosen_parameters_match"] == ([choice[:2] for choice in choices] == archive_params)
        max_score_error = max(abs(score - float(row["inner_spearman"]))
                              for (_, _, score), row in zip(choices, params, strict=True))
        np.testing.assert_allclose(max_score_error, result["max_abs_archived_inner_score_error"],
                                   atol=1e-12, rtol=0)
        archived_choice_oof = np.empty_like(y)
        splits = KFold(5, shuffle=True, random_state=0).split(X)
        for (train, test), (dim, alpha) in zip(splits, archive_params, strict=True):
            archived_choice_oof[test] = predict(X, y, train, test, dim, alpha)
        np.testing.assert_allclose(np.max(np.abs(archived_choice_oof - expected)),
                                   result["max_abs_archived_prediction_error_with_archived_choices"],
                                   atol=1e-9, rtol=0)
    match = bool(max_error < 2e-5 and result["chosen_parameters_match"] is not False)
    assert result["baseline_match"] == match
    paired = result["paired"]
    if not match:
        assert paired is None
        return

    replace = np.array(labels) == FOCAL
    X[replace] = vectors["replacement"]
    corrected, after_choices = cross_validate(X, y, tune=branch == "tune_stepA")
    if choices:
        assert paired["tuning_changes"] == [
            {"fold": fold, "before": list(before[:2]), "after": list(after[:2])}
            for fold, (before, after) in enumerate(zip(choices, after_choices, strict=True))
            if before[:2] != after[:2]
        ]
    check_pair(y, oof, corrected, replace, paired)
    if choices:
        frozen = np.empty_like(y)
        splits = KFold(5, shuffle=True, random_state=0).split(X)
        for (train, test), (dim, alpha) in zip(splits, archive_params, strict=True):
            frozen[test] = predict(X, y, train, test, dim, alpha)
        check_pair(y, oof, frozen, replace, paired["frozen_hyperparameters"])
        np.testing.assert_allclose(
            np.max(np.abs(corrected - frozen)), paired["max_abs_retuning_change"],
            atol=1e-9, rtol=0,
        )


def check_pair(y, oof, corrected, replace, paired):
    assert [row["index"] for row in paired["focal"]] == list(np.flatnonzero(replace))
    for row in paired["focal"]:
        i = row["index"]
        np.testing.assert_allclose([row["y_true"], row["before"], row["after"]],
                                   [y[i], oof[i], corrected[i]], atol=1e-9, rtol=0)
    for name, prediction in (("before", oof), ("after", corrected)):
        assert abs(spearmanr(y, prediction).statistic - paired[f"{name}_spearman"]) < 1e-12
        assert abs(np.linalg.norm(y - prediction) / np.sqrt(len(y)) - paired[f"{name}_rmse"]) < 1e-9
    difference = np.abs(corrected - oof)
    assert abs(np.max(difference) - paired["max_abs_change"]) < 1e-9
    assert abs(np.median(difference) - paired["median_abs_change"]) < 1e-9


def check_focal_embeddings(archive, layer, vectors):
    import esm
    import torch

    runtime = ROOT / "runtime"
    cache = runtime / "hub33" if layer == "t33" and (runtime / "hub33").exists() else runtime / "hub"
    torch.hub.set_dir(os.environ.get("RUBISCO_TORCH_HOME", str(cache)))
    os.environ["TORCH_FORCE_WEIGHTS_ONLY_LOAD"] = "1"
    torch.serialization.add_safe_globals([argparse.Namespace])
    torch.set_num_threads(4)
    depth = int(layer[1:])
    name = {12: "esm2_t12_35M_UR50D", 33: "esm2_t33_650M_UR50D"}[depth]
    model, alphabet = esm.pretrained.load_model_and_alphabet(name)
    model.eval()
    records = dict((lines[0], "".join(lines[1:]))
                   for block in archive_fasta(archive).decode().split(">")[1:]
                   if (lines := block.splitlines()))
    repaired = json.loads((SOURCES / "Q6N0W9.json").read_text())["sequence"]["value"]
    for label, seq in ((FOCAL, records[FOCAL]), ("replacement", repaired)):
        seq = seq.replace(" ", "").replace("-", "")
        _, _, tokens = alphabet.get_batch_converter()([(label, seq)])
        with torch.inference_mode():
            output = model(tokens, repr_layers=[depth], return_contacts=False)
            actual = output["representations"][depth][0, 1:len(seq) + 1].mean(0).numpy()
        np.testing.assert_allclose(actual, vectors[label], atol=2e-6, rtol=0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-embedding", action="store_true")
    args = parser.parse_args()
    result = json.loads((ROOT / "results.json").read_text())
    assert set(result["replays"]) == {"t12", "t33"}
    for target in ("S", "KC"):
        assert result["replays"]["t12"]["tune_stepA"][target]["baseline_match"]
        assert not result["replays"]["t33"]["tune_stepA"][target]["baseline_match"]
    archive_path = SOURCES / "version3.zip"
    deposit = json.loads((SOURCES / "zenodo.json").read_text())
    item, = (item for item in deposit["files"] if item["key"].endswith("-Version3.zip"))
    assert "md5:" + hashlib.md5(archive_path.read_bytes()).hexdigest() == item["checksum"]
    tree = json.loads((SOURCES / "historical-tree.json").read_text())
    blobs = {item["path"]: item["sha"] for item in tree["tree"] if item["type"] == "blob"}

    def verify_historical_blob(path, data):
        prefix = f"blob {len(data)}\0".encode()
        assert hashlib.sha1(prefix + data).hexdigest() == blobs[path], path

    with ZipFile(archive_path) as archive:
        original = old_metadata()
        current = csv_rows(archive, "flamholz_rubisco_wt_meta.csv")
        assert current == [{**row, "species": row["species"].rstrip()} for row in original]
        assert result["metadata_rows_with_trailing_whitespace_removed"] == sum(
            a["species"] != b["species"] for a, b in zip(original, current, strict=True)
        )
        assert archive_fasta(archive) == (SOURCES / "71534b19c-flamholz_rubisco_wt.fasta").read_bytes()
        for full in archive.namelist():
            path = full.partition("/")[2]
            if path.endswith("/flamholz_rubisco_wt_meta.csv"):
                verify_historical_blob(path, (SOURCES / "71534b19c-flamholz_rubisco_wt_meta.csv").read_bytes())
            elif path.endswith(("/flamholz_rubisco_wt.fasta", "/flamholz_dataset_S1.csv")) or (
                ("/tune_stepA/" in path or "/embedding_vs_kinetics/" in path) and path.endswith(".csv")
            ):
                verify_historical_blob(path, archive.read(full))
        for layer, run in result["replays"].items():
            with np.load(SOURCES / f"{layer}.npz") as vectors:
                width = {"t12": 480, "t33": 1280}[layer]
                labels = {line[1:].split()[0] for line in archive_fasta(archive).decode().splitlines()
                          if line.startswith(">")}
                assert set(vectors.files) == labels | {"replacement"}
                assert all(vectors[label].shape == (width,) and np.isfinite(vectors[label]).all()
                           for label in vectors.files)
                assert not np.array_equal(vectors[FOCAL], vectors["replacement"])
                if args.verify_embedding:
                    check_focal_embeddings(archive, layer, vectors)
                for branch, targets in (("embedding_vs_kinetics", TARGETS), ("tune_stepA", TARGETS[:4])):
                    assert set(targets) == set(run[branch])
                    for target in targets:
                        check(archive, vectors, branch, target, run[branch][target])
    print("Verified target rows, archived replay residuals and any matched paired reruns.")


if __name__ == "__main__":
    main()
