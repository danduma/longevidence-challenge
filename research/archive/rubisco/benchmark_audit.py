#!/usr/bin/env python3
"""Reconcile deposited predictions and reconstruct, rather than assume, overlap."""
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median

from input_audit import canonical_species, fasta_records, positive


def cross_fold_exposure(groups, folds):
    group_folds = defaultdict(set)
    for group, fold in zip(groups, folds, strict=True):
        if group is not None:
            group_folds[group].add(fold)
    return [i for i, group in enumerate(groups) if group is not None and len(group_folds[group]) > 1]


def implied_join_count(rows, metadata, fasta_labels, target):
    """Target row count implied by the shipped species merge and supplied IDs.

    Treat positive-target selection after aggregation as the shipped script does.
    Embedding values cannot change this count. Missing embedding IDs can, which
    is why we explicitly state the assumption that all deposited FASTA IDs were used.
    """
    grouped = defaultdict(list)
    for row in rows:
        if row.get("primary") == "1" and row.get("mutant") == "0":
            try:
                value = float(row.get(target, ""))
                if math.isfinite(value):
                    grouped[row["species"]].append(value)
            except (ValueError, TypeError):
                pass
    label_species = {r["label"]: r["species"] for r in metadata}
    return sum(bool(grouped.get(label_species.get(label))) and median(grouped[label_species[label]]) > 0
               for label in fasta_labels)


def metrics(y, predictions):
    from scipy.stats import spearmanr
    from sklearn.metrics import r2_score, mean_squared_error
    return {"n": len(y), "r2": float(r2_score(y, predictions)),
            "rmse_log10": float(mean_squared_error(y, predictions) ** 0.5),
            "spearman": float(spearmanr(y, predictions).statistic)}


def main():
    import numpy as np
    from sklearn.model_selection import KFold, GroupKFold
    root = Path(__file__).resolve().parents[1]
    source = root / "sources/rubisco-active-learning/supplementary/doubling_time_package_v2"
    rows = list(csv.DictReader((source / "inputs/flamholz_dataset_S1.csv").open()))
    meta = list(csv.DictReader((source / "manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv").open()))
    records = fasta_records((source / "manuscript_assets/flamholz/flamholz_rubisco_wt.fasta").read_text())
    sequences = dict(records)
    species_sequences = {canonical_species(r["species"]): sequences[r["label"]]
                         for r in meta if r["label"] in sequences}
    out = root / "rubisco/results"
    out.mkdir(exist_ok=True)
    summaries, exposures = [], []
    for path in sorted((source / "manuscript_assets/flamholz/embedding_vs_kinetics").glob("pred_*.csv")):
        target = path.stem.removeprefix("pred_")
        predictions = list(csv.DictReader(path.open()))
        names = [canonical_species(r["species"]) for r in predictions]
        sequence_keys = [hashlib.sha256(species_sequences[s].encode()).hexdigest() if s in species_sequences else None for s in names]
        fold = np.full(len(names), -1, dtype=int)
        for index, (_, test) in enumerate(KFold(5, shuffle=True, random_state=0).split(names)):
            fold[test] = index
        exact_exposure = cross_fold_exposure(names, fold)
        sequence_exposure = cross_fold_exposure(sequence_keys, fold)
        # This is a split-integrity check, not a refitted prediction experiment.
        corrected_fold = np.full(len(names), -1, dtype=int)
        safe_groups = [s if s is not None else "missing:" + name for s, name in zip(sequence_keys, names, strict=True)]
        for index, (_, test) in enumerate(GroupKFold(5).split(names, groups=safe_groups)):
            corrected_fold[test] = index
        assert not cross_fold_exposure(sequence_keys, corrected_fold)
        y = [float(r["y_true"]) for r in predictions]
        yh = [float(r["y_pred"]) for r in predictions]
        by_species = defaultdict(list)
        for name, row in zip(names, predictions, strict=True):
            by_species[name].append(row)
        if any(len({r["y_true"] for r in group}) != 1 for group in by_species.values()):
            raise ValueError(f"Inconsistent duplicated truth values for {target}")
        # Reweight already stored predictions. This does NOT remove training leakage.
        species_y = [float(group[0]["y_true"]) for group in by_species.values()]
        species_yh = [sum(float(r["y_pred"]) for r in group) / len(group) for group in by_species.values()]
        summaries.append({"target": target, "deposited_metrics_recomputed": metrics(y, yh),
                          "species_reweighted_stored_predictions_NOT_refitted": metrics(species_y, species_yh),
                          "unique_species": len(by_species),
                          "deposited_duplicate_rows": sum(len(v) for v in by_species.values() if len(v) > 1),
                          "implied_target_count_using_all_released_fasta_ids": implied_join_count(rows, meta, [r[0] for r in records], target),
                          "default_seed_0_reconstructed_species_exposed_rows": len(exact_exposure),
                          "default_seed_0_reconstructed_sequence_exposed_rows": len(sequence_exposure),
                          "species_without_sequence_mapping": sorted({n for n, s in zip(names, sequence_keys, strict=True) if s is None}),
                          "sequence_grouped_split_exposed_rows": len(cross_fold_exposure(sequence_keys, corrected_fold))})
        for i, row in enumerate(predictions):
            exposures.append(dict(target=target, deposited_row=i, species=row["species"],
                                  reconstructed_default_fold=int(fold[i]), sequence_sha256=sequence_keys[i],
                                  same_species_in_training=i in exact_exposure,
                                  same_released_sequence_in_training=i in sequence_exposure,
                                  grouped_fold=int(corrected_fold[i])))
    report = {"repository_commit": "ade87f378aa25c24a24c19915722aa04be92936e",
              "fold_assumption": "Shipped KFold(n_splits=5, shuffle=True, random_state=0) applied to deposited prediction row order. Original split assignments are not independently recovered.",
              "join_assumption": "All released FASTA IDs used, including repeated IDs. The missing embedding file prevents establishing the actual embedding-ID list.",
              "metric_caveat": "Recomputed/reweighted stored predictions do not estimate the performance of refitted, leakage-free models.",
              "targets": summaries}
    (out / "benchmark_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    with (out / "fold_exposure.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(exposures[0]))
        writer.writeheader()
        writer.writerows(exposures)
    for summary in summaries:
        print(summary["target"], "rows", summary["deposited_metrics_recomputed"]["n"],
              "unique", summary["unique_species"], "species exposure", summary["default_seed_0_reconstructed_species_exposed_rows"],
              "sequence exposure", summary["default_seed_0_reconstructed_sequence_exposed_rows"],
              "implied released-ID rows", summary["implied_target_count_using_all_released_fasta_ids"])


if __name__ == "__main__":
    main()
