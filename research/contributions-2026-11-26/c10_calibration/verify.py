#!python3 -B
"""C10. Recompute tables with --write, or verify the distributed results offline."""
import argparse
import hashlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import isotonic_regression

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
CANON = "ENLYFQS"
POSITIONS = ["P6", "P5", "P4", "P3", "P2", "P1", "P1'"]


def load_pairs():
    for line in (ROOT / "frozen.sha256").read_text().splitlines():
        digest, file = line.split(maxsplit=1)
        assert hashlib.sha256((ROOT / file).read_bytes()).hexdigest() == digest, file
    pairs = pd.read_csv(INPUT / "biochemical_common_keys.csv", dtype={"TEVp": str})
    keys = ["TEVp", "TEVs"]
    assert not pairs.duplicated(keys).any()
    assert (pairs[["min_raw", "min_corrected"]] >= 100).all().all()
    assert np.isfinite(pairs[["AUC_corrected", "vinirel"]]).all().all()
    rel = pd.read_table(INPUT / "protease_vitro_relative.tsv", dtype={"TEVp": str})
    initial = pd.read_table(INPUT / "data_viniabs_Luzi.txt", dtype={"TEVp": str})
    data = rel.merge(initial, on=["pos", "TEVs_AA", "TEVp"], validate="1:1")
    assert len(data) == len(initial) == len(rel)
    data["TEVs"] = [CANON[:i] + aa + CANON[i+1:] for i, aa in
                    zip(data.pos.map(POSITIONS.index), data.TEVs_AA)]
    canonical = data[data.TEVs == CANON]
    blocks = ["TEVp", "pos"]
    data = data.merge(canonical[blocks + ["viniabs"]].rename(columns={"viniabs": "control"}),
                      on=blocks, validate="m:1")
    data["ratio_error"] = abs(data.viniabs / data.control - data.vinirel)
    audit = data.groupby(blocks).agg(control=("control", "first"), source_rows=("pos", "size"),
                                     max_ratio_error=("ratio_error", "max")).reset_index()
    assert canonical.vinirel.eq(1).all()
    assert len(canonical) == 42 and audit.max_ratio_error.max() < 1e-10
    labels = data.groupby(keys).vinirel.mean()
    np.testing.assert_allclose(labels.loc[pd.MultiIndex.from_frame(pairs[keys])], pairs.vinirel,
                               rtol=0, atol=1e-13)
    position = data[data.TEVs != CANON].set_index("TEVs").pos.to_dict()
    pairs["position"] = pairs.TEVs.map(position)
    assert pairs.position.isna().equals(pairs.TEVs.eq(CANON))
    return pairs, audit


def predict(train, test, design):
    x, y = train.AUC_corrected.to_numpy(), train.vinirel.to_numpy()
    xt = test.AUC_corrected.to_numpy()
    points, ids, counts = np.unique(x, return_inverse=True, return_counts=True)
    grouped_y = np.bincount(ids, weights=y) / counts
    monotone = isotonic_regression(grouped_y, weights=counts).x
    slope, intercept = np.polyfit(x, y, 1)
    predictions = {"mean": np.full(len(test), y.mean()), "linear": intercept + slope * xt,
                   "isotonic": np.interp(xt, points, monotone)}
    if design == "position":
        predictions["panel_mean"] = test.TEVp.map(train.groupby("TEVp").vinirel.mean())
    return predictions


def evaluate(pairs):
    rows = []
    for design, column in [("enzyme", "TEVp"), ("position", "position")]:
        for variant in ["primary", "include_control"]:
            selected = pairs if variant == "include_control" else pairs[pairs.TEVs != CANON]
            for fold in sorted(selected[column].dropna().unique()):
                held = selected[column].eq(fold).fillna(False)
                train, test = selected[~held], selected[held]
                assert set(train[column].dropna()).isdisjoint(test[column])
                if design == "position":
                    assert set(train.TEVs).isdisjoint(test.TEVs)
                for model, values in predict(train, test, design).items():
                    frame = test[["TEVp", "TEVs", "vinirel"]].copy()
                    rows.append(frame.assign(design=design, variant=variant, fold=fold, model=model,
                                             train_n=len(train), prediction=values))
    return pd.concat(rows, ignore_index=True)


def metrics(g):
    y, p = g.vinirel.to_numpy(), g.prediction.to_numpy()
    e = p - y
    return dict(n=len(e), rmse=np.sqrt(np.mean(e**2)), mae=np.mean(abs(e)), bias=np.mean(e))


def compute():
    pairs, audit = load_pairs()
    predictions = evaluate(pairs)
    folds, summary = [], []
    for (design, variant, model), group in predictions.groupby(["design", "variant", "model"]):
        key = dict(design=design, variant=variant, model=model)
        expected = len(pairs) if (design, variant) == ("enzyme","include_control") else (pairs.TEVs != CANON).sum()
        assert len(group) == expected and not group.duplicated(["TEVp","TEVs"]).any()
        scores = [dict(**key, fold=f, **metrics(g)) for f, g in group.groupby("fold")]
        folds.extend(scores)
        slope, intercept = np.polyfit(group.prediction, group.vinirel, 1)
        summary.append(dict(**key, **metrics(group), macro_rmse=np.mean([s["rmse"] for s in scores]),
                            calibration_intercept=intercept, calibration_slope=slope))
    return {"control_audit": audit, "predictions": predictions, "fold_errors": pd.DataFrame(folds),
            "summary": pd.DataFrame(summary)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    for name, table in compute().items():
        path = ROOT / "results" / (name + ".csv")
        buffer = io.StringIO()
        table.to_csv(buffer, index=False, float_format="%.12g")
        if args.write:
            path.write_text(buffer.getvalue())
        else:
            expected = pd.read_csv(path, dtype={"fold":str, "TEVp":str})
            rebuilt = pd.read_csv(io.StringIO(buffer.getvalue()), dtype={"fold":str, "TEVp":str})
            pd.testing.assert_frame_equal(expected, rebuilt, rtol=1e-10, atol=1e-12)
        print(f"C10 {name}: {len(table)} rows " + ("written" if args.write else "verified"))


if __name__ == "__main__":
    main()
