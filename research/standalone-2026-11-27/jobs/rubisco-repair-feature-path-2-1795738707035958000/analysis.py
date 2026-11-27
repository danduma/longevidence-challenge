"""Rebuild the archived target rows and replay the released predictors."""
import argparse
import importlib.metadata
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"
RUNTIME = ROOT / "runtime"
TARGETS = ["S", "KC", "vC", "KO", "vO", "KRuBP"]
BASE = "flamholz/"
FOCAL = "Rhodopseudomonas_palustris"
PREDICTION_TOLERANCE = 2e-5
HISTORICAL_META = SOURCES / "71534b19c-flamholz_rubisco_wt_meta.csv"


def source(suffix):
    with zipfile.ZipFile(SOURCES / "version3.zip") as archive:
        name, = (name for name in archive.namelist() if name.endswith("/" + suffix))
        return archive.read(name)


def table(suffix):
    return pd.read_csv(io.BytesIO(source(suffix)))


def release_script(name):
    RUNTIME.mkdir(exist_ok=True)
    path = RUNTIME / name
    path.write_bytes(source(name))
    return path


def fasta(text):
    for record in text.split(">")[1:]:
        header, *lines = record.splitlines()
        yield header.split()[0], "".join(lines)


def embed(layer):
    """Run the original extractor once per distinct sequence, with batch size one."""
    records = {}
    for label, sequence in fasta(source(BASE + "flamholz_rubisco_wt.fasta").decode()):
        assert records.setdefault(label, sequence) == sequence
    replacement = json.loads((SOURCES / "Q6N0W9.json").read_text())["sequence"]["value"]
    assert FOCAL in records and replacement != records[FOCAL]
    records["replacement"] = replacement
    representatives = {}
    for label, sequence in records.items():
        representatives.setdefault(sequence, label)
    unique = RUNTIME / "unique.fasta"
    unique.write_text("".join(f">{label}\n{seq}\n" for seq, label in representatives.items()))
    raw = RUNTIME / f"t{layer}-unique.npy"
    model = {12: "esm2_t12_35M_UR50D", 33: "esm2_t33_650M_UR50D"}[layer]
    extractor = (
        "import argparse, os, runpy, sys, torch;"
        "torch.set_num_threads(4);"
        "torch.hub.set_dir(os.environ['TORCH_HOME']);"
        "torch.serialization.add_safe_globals([argparse.Namespace]);"
        "sys.argv.pop(0);"
        "runpy.run_path(sys.argv[0], run_name='__main__')"
    )
    subprocess.run([
        sys.executable, "-B", "-c", extractor, str(release_script("embed_rbcL_fasta_esm2.py")),
        "--fasta", str(unique), "--out_npy", str(raw), "--model", model,
        "--repr_layer", str(layer), "--pool", "mean", "--device", "cpu",
    ], env={**os.environ, "TORCH_HOME": os.environ.get("RUBISCO_TORCH_HOME", str(RUNTIME / "hub")),
            "TORCH_FORCE_WEIGHTS_ONLY_LOAD": "1"}, check=True)
    # Object data here were generated locally by the preceding process.
    obj = np.load(raw, allow_pickle=True).item()
    vectors = dict(zip(obj["ids"], obj["emb"], strict=True))
    np.savez_compressed(SOURCES / f"t{layer}.npz", **{
        label: vectors[representatives[seq]] for label, seq in records.items()
    })


def rows(restore=True):
    raw = table("inputs/flamholz_dataset_S1.csv")
    raw = raw[(raw.primary == 1) & (raw.mutant == 0)]
    meta = pd.read_csv(HISTORICAL_META) if restore else table(BASE + "flamholz_rubisco_wt_meta.csv")
    label_to_species = dict(zip(meta.label, meta.species, strict=True))
    ids = [label for label, _ in fasta(source(BASE + "flamholz_rubisco_wt.fasta").decode())]
    labels = pd.DataFrame({"label": ids, "species": [label_to_species[label] for label in ids]})
    medians = raw.groupby("species", as_index=False)[TARGETS].median()
    return medians.merge(labels, on="species", how="inner")


def load_script(name):
    spec = importlib.util.spec_from_file_location(name[:-3], release_script(name))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def paired_metrics(y, before, after, focal):
    delta = after - before
    focal_indices = np.flatnonzero(focal)
    return {
        "focal": [{"index": int(i), "y_true": float(y[i]),
                   "before": float(before[i]), "after": float(after[i])}
                  for i in focal_indices],
        "before_spearman": float(spearmanr(y, before).statistic),
        "after_spearman": float(spearmanr(y, after).statistic),
        "before_rmse": float(np.sqrt(np.mean((y - before) ** 2))),
        "after_rmse": float(np.sqrt(np.mean((y - after) ** 2))),
        "max_abs_change": float(np.max(np.abs(delta))),
        "median_abs_change": float(np.median(np.abs(delta))),
    }


def replay(layer):
    os.environ.setdefault("MPLCONFIGDIR", str(RUNTIME / "matplotlib"))
    nested = load_script("tune_flamholz_ridge_pca.py")
    fixed = load_script("analyze_flamholz_embeddings_vs_kinetics.py")
    vectors = np.load(SOURCES / f"t{layer}.npz")
    literal = rows(False)
    restored = rows()
    outcome = {}
    for branch, targets in [
        ("embedding_vs_kinetics", TARGETS), ("tune_stepA", TARGETS[:4])
    ]:
        outcome[branch] = {}
        for target in targets:
            old = table(f"{BASE}{branch}/pred_{target}.csv")
            sub = restored[restored[target] > 0]
            y = np.log10(sub[target].to_numpy())
            assert old.species.tolist() == sub.species.tolist()
            np.testing.assert_allclose(y, old.y_true, rtol=0, atol=1e-12)
            X = np.stack([vectors[label] for label in sub.label]).astype(float)
            if branch == "tune_stepA":
                def fit(matrix):
                    return nested.nested_cv(
                        matrix, y, None,
                        KFold(5, shuffle=True, random_state=0),
                        KFold(4, shuffle=True, random_state=0),
                        [0, 8, 16, 32, 64], [0.1, 1, 10, 100], 0,
                    )
                predictions, parameters = fit(X)
                expected = table(f"{BASE}{branch}/chosen_params_{target}.csv")
                params_match = bool((parameters[["pca_dim", "alpha"]] ==
                                     expected[["pca_dim", "alpha"]]).all().all())
                score_error = float(np.max(np.abs(parameters.inner_spearman - expected.inner_spearman)))
                archived_choice_predictions = np.empty_like(y)
                folds = KFold(5, shuffle=True, random_state=0).split(X)
                for (train, test), choice in zip(folds, expected.itertuples(), strict=True):
                    archived_choice_predictions[test] = nested.fit_predict_fold(
                        X[train], y[train], X[test], choice.pca_dim, choice.alpha, 0
                    )
                choice_error = float(np.max(np.abs(archived_choice_predictions - old.y_pred)))
            else:
                def fit(matrix):
                    return fixed.cv_predict_ridge_pca(matrix, y, 5, 0, 16, 10.0), None
                predictions, _ = fit(X)
                params_match = None
                score_error = None
                choice_error = None
            baseline_error = float(np.max(np.abs(predictions - old.y_pred)))
            match = baseline_error < PREDICTION_TOLERANCE and params_match is not False
            focal = sub.label.to_numpy() == FOCAL
            paired = None
            if match:
                replaced = X.copy()
                replaced[focal] = vectors["replacement"]
                after, after_parameters = fit(replaced)
                paired = paired_metrics(y, predictions, after, focal)
                if after_parameters is not None:
                    frozen_predictions = np.empty_like(y)
                    folds = KFold(5, shuffle=True, random_state=0).split(replaced)
                    for (train, test), choice in zip(folds, expected.itertuples(), strict=True):
                        frozen_predictions[test] = nested.fit_predict_fold(
                            replaced[train], y[train], replaced[test], choice.pca_dim, choice.alpha, 0
                        )
                    paired["frozen_hyperparameters"] = paired_metrics(
                        y, predictions, frozen_predictions, focal
                    )
                    paired["max_abs_retuning_change"] = float(np.max(np.abs(after - frozen_predictions)))
                    before_choices = parameters[["pca_dim", "alpha"]].to_numpy().tolist()
                    after_choices = after_parameters[["pca_dim", "alpha"]].to_numpy().tolist()
                    paired["tuning_changes"] = [
                        {"fold": fold, "before": before, "after": changed}
                        for fold, (before, changed) in enumerate(zip(before_choices, after_choices, strict=True))
                        if before != changed
                    ]
            outcome[branch][target] = {
                "literal_rows": int((literal[target] > 0).sum()),
                "restored_rows": len(y),
                "max_abs_archived_prediction_error": baseline_error,
                "chosen_parameters_match": params_match,
                "max_abs_archived_inner_score_error": score_error,
                "max_abs_archived_prediction_error_with_archived_choices": choice_error,
                "baseline_match": match,
                "paired": paired,
            }
            print(branch, target, outcome[branch][target], flush=True)
    (RUNTIME / f"results-t{layer}.json").write_text(json.dumps(outcome, indent=2) + "\n")
    packages = ["fair-esm", "torch", "numpy", "pandas", "scipy", "scikit-learn", "matplotlib"]
    versions = {pkg: importlib.metadata.version(pkg) for pkg in packages}
    (RUNTIME / f"versions-t{layer}.json").write_text(json.dumps(versions, indent=2) + "\n")


def summarize():
    historical = pd.read_csv(HISTORICAL_META)
    deposited = table(BASE + "flamholz_rubisco_wt_meta.csv")
    layers = [path.stem.removeprefix("results-") for path in sorted(RUNTIME.glob("results-t*.json"))]
    if not layers:
        raise SystemExit("Replay an extracted checkpoint first.")
    result = {
        "metadata_rows_with_trailing_whitespace_removed": int((historical.species != deposited.species).sum()),
        "matching_rule": {"max_abs_prediction_error_below": PREDICTION_TOLERANCE,
                          "require_matching_nested_pca_dimensions_and_penalties": True},
        "cv": {"outer_folds": 5, "inner_folds": 4, "seed": 0},
        "replays": {layer: json.loads((RUNTIME / f"results-{layer}.json").read_text()) for layer in layers},
        "software": {layer: json.loads((RUNTIME / f"versions-{layer}.json").read_text()) for layer in layers},
    }
    (ROOT / "results.json").write_text(json.dumps(result, indent=2) + "\n")


def export(layer):
    """Write local feature files in the released scripts' input format."""
    ids = np.array([label for label, _ in fasta(source(BASE + "flamholz_rubisco_wt.fasta").decode())],
                   dtype=object)
    assert np.count_nonzero(ids == FOCAL) == 1
    with np.load(SOURCES / f"t{layer}.npz") as vectors:
        original = np.stack([vectors[label] for label in ids])
        corrected = original.copy()
        corrected[ids == FOCAL] = vectors["replacement"]
        RUNTIME.mkdir(exist_ok=True)
        for name, matrix in (("original", original), ("corrected", corrected)):
            np.save(RUNTIME / f"t{layer}-{name}.npy", {"ids": ids, "emb": matrix})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["embed", "replay", "summarize", "export"])
    parser.add_argument("layer", type=int, choices=[12, 33], nargs="?")
    args = parser.parse_args()
    if args.stage == "summarize":
        summarize()
    elif args.layer is None:
        parser.error("embed, replay and export require a layer")
    else:
        {"embed": embed, "replay": replay, "export": export}[args.stage](args.layer)
