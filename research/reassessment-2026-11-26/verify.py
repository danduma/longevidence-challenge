#!/usr/bin/env python3
"""Reproduce the surviving checks and the omitted-baseline counterexample.

Run --download-inputs once, or supply --rubisco-inputs with the two original
CSVs. No models, fitting services or private workspace paths are required.
"""
import argparse
import csv
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
BATCH = ROOT.parent / "contributions-2026-11-26"


def checked(path, digest):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, path.name
    return path


def csv_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def c07():
    folder = BATCH / "c07_mapping"
    hashes = dict(line.split(maxsplit=1)[::-1]
                  for line in (folder / "manifest.sha256").read_text().splitlines())

    def load(name):
        return json.loads(checked(folder / name, hashes[name]).read_text())

    old = load("evidence/uniprot/Q07N61.json")
    new = load("evidence/uniprot/Q6N0W9.json")
    pdb = load("evidence/pdb_4LF1.json")
    assert "BisA53" in old["organism"]["scientificName"]
    assert "CGA009" in new["organism"]["scientificName"]
    assert pdb["entity_src_gen"][0]["gene_src_strain"] == "CGA009"
    assert pdb["entity_poly"]["pdbx_seq_one_letter_code_can"] == (
        "MGSSHHHHHHSSGLVPRGSH" + new["sequence"]["value"])
    assert old["sequence"]["value"] != new["sequence"]["value"]
    return dict(old_accession="Q07N61", old_strain="BisA53", old_length=460,
                replacement_accession="Q6N0W9", replacement_strain="CGA009",
                replacement_length=len(new["sequence"]["value"]),
                deposition_native_core_exact_match=True,
                scope="Accession/deposition consistency; original assay identity is linked in README, not re-extracted here.",
                upstream_acceptance=False, downstream_model_improvement_tested=False)


def c08(inputs, download):
    sources = json.loads((BATCH / "c08_selection/sources.json").read_text())["inputs"]
    paths = {}
    for key in ("supp", "library"):
        spec = sources[key]
        path = inputs / Path(spec["upstream_path"]).name
        if download:
            inputs.mkdir(parents=True, exist_ok=True)
            data = urlopen(spec["download_url"], timeout=60).read()
            assert hashlib.sha256(data).hexdigest() == spec["sha256"], key
            path.write_bytes(data)
        paths[key] = checked(path, spec["sha256"])
    master = {r["internal_id"]: r for r in csv_rows(paths["library"])}
    rows = csv_rows(paths["supp"])
    assert len(rows) == 469 and len({r["internal_id"] for r in rows}) == 469
    groups = {g: [r for r in rows if r["form"] == g] for g in ("Form 1", "Form 2")}
    assert [len(groups[g]) for g in groups] == [144, 178]
    results = []
    for threshold in (10, 15):
        # Counts can change only at these exact rational breakpoints. Include
        # every breakpoint and an interior point of every intervening interval.
        points = {F("1.7"), F("2.8")}
        for group in groups.values():
            for row in group:
                if row["rate_lab_corrected"]:
                    raw = F(master[row["internal_id"]]["rate_lab"])
                    if raw > 0 and F("1.7") <= F(threshold) / raw <= F("2.8"):
                        points.add(F(threshold) / raw)
        ordered = sorted(points)
        tested = ordered + [(a + b) / 2 for a, b in zip(ordered, ordered[1:])]
        gaps = []
        for factor in tested:
            high = {g: sum(bool(r["rate_lab_corrected"]) and
                          F(master[r["internal_id"]]["rate_lab"]) * factor >= threshold
                          for r in group) for g, group in groups.items()}
            unknown_i = sum(not r["rate_lab_corrected"] and not r["active<0.5s-1"]
                            for r in groups["Form 1"])
            gap = F(high["Form 2"], 178) - F(high["Form 1"] + unknown_i, 144)
            gaps.append(gap)
        results.append(dict(threshold_s_minus1=threshold,
                            factor_interval=[1.7, 2.8],
                            exact_breakpoints_and_midpoints=len(tested),
                            minimum_lower_II_minus_upper_I_pp=float(min(gaps)) * 100,
                            order_survives_all_tested_intervals=min(gaps) > 0))
    assert results[0]["order_survives_all_tested_intervals"]
    assert not results[1]["order_survives_all_tested_intervals"]
    return dict(comparisons=results,
                scope="Selected finite library; comparable latent endpoints and one common calibration multiplier assumed. No enzyme-specific uncertainty or population inference.")


def c10():
    path = checked(BATCH / "c10_calibration/inputs/biochemical_common_keys.csv",
                   "e2c93066bc21ba5e05d823152a0f8e370af7adead7d54d9178a4efe0594edaf8")
    data = pd.read_csv(path, dtype={"TEVp": str})
    data = data[data.TEVs != "ENLYFQS"].copy()
    assert len(data) == 780 and (data.vinirel < 0).sum() == 170
    folds = []
    for enzyme, test in data.groupby("TEVp", sort=True):
        train = data[data.TEVp != enzyme]
        assert set(train.TEVp).isdisjoint(test.TEVp)
        x, y = train.AUC_corrected.to_numpy(), train.vinirel.to_numpy()
        slope = np.sum((x - x.mean()) * (y - y.mean())) / np.sum((x - x.mean()) ** 2)
        linear = y.mean() + slope * (test.AUC_corrected - x.mean())
        motif_mean = test.TEVs.map(train.groupby("TEVs").vinirel.mean())
        motif_median = test.TEVs.map(train.groupby("TEVs").vinirel.median())
        assert motif_mean.notna().all() and motif_median.notna().all()
        rmse = lambda prediction: float(np.sqrt(np.mean((prediction - test.vinirel) ** 2)))
        folds.append(dict(held_enzyme=enzyme, n=len(test),
                          pooled_mean_RMSE=rmse(y.mean()),
                          AUC_linear_RMSE=rmse(linear),
                          motif_mean_RMSE=rmse(motif_mean),
                          motif_median_RMSE=rmse(motif_median)))
    averages = {key: float(np.mean([fold[key] for fold in folds]))
                for key in folds[0] if key.endswith("RMSE")}
    assert all(f["motif_mean_RMSE"] < f["AUC_linear_RMSE"] for f in folds)
    assert abs(averages["motif_mean_RMSE"] - .1212559393396278) < 1e-12
    return dict(noncanonical_pairs=780, negative_estimates_retained=170,
                folds=folds, mean_fold_RMSE=averages,
                scope="Exploratory audit comparator for known motifs in held-enzyme folds. Does not refute position-held-out analysis or test unseen motifs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rubisco-inputs", type=Path, default=ROOT / "downloaded-inputs")
    parser.add_argument("--download-inputs", action="store_true")
    parser.add_argument("--write", action="store_true", help="regenerate results.json")
    args = parser.parse_args()
    output = dict(C07=c07(), C08=c08(args.rubisco_inputs, args.download_inputs), C10=c10())
    target = ROOT / "results.json"
    if args.write:
        target.write_text(json.dumps(output, indent=2) + "\n")
    else:
        expected = json.loads(target.read_text())

        def compare(a, b):
            if isinstance(a, dict):
                assert a.keys() == b.keys()
                for key in a:
                    compare(a[key], b[key])
            elif isinstance(a, list):
                assert len(a) == len(b)
                for left, right in zip(a, b):
                    compare(left, right)
            elif isinstance(a, float):
                assert abs(a - b) < 1e-10
            else:
                assert a == b, (a, b)
        compare(output, expected)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
