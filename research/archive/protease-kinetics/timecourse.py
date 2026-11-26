#!/usr/bin/env python3
"""Locked, target-free timecourse-score comparison on selected assay panels.

Scores are empirical recorder quantities, not measured catalytic rate constants.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

AA = "ACDEFGHIKLMNPQRSTVWY"
WT = "ENLYFQS"
NEG = ["(1-13)-mCherry", "I C151A 30*31*", "no TEVp"]
PAIR = ["TEVp", "TEVs"]
COUNTS = [f"{a}{i}" for a in ["F", "N"] for i in range(1, 6)]


def hazard_slope(flipped, unflipped):
    f, n = np.asarray(flipped, dtype=float), np.asarray(unflipped, dtype=float)
    if f.shape != (5,) or n.shape != (5,) or not np.isfinite(f).all() or not np.isfinite(n).all() or (f < 0).any() or (n < 0).any():
        raise ValueError("Expected five finite nonnegative flipped/unflipped counts")
    h = -np.log((n + 0.5) / (f + n + 1))
    t = np.arange(5)
    return float(np.dot(t, h - h[0]) / np.dot(t, t))


def pin_inputs(root):
    inputs = root / "inputs"
    inputs.mkdir(exist_ok=True)
    project = root.parent
    sources = [project / "agents/protease-corrected-sol/results/controls_corrected.tsv.gz",
               project / "agents/scout-data-sol/data/protease_vitro_relative.tsv"]
    records = []
    for path in sources:
        content = path.read_bytes()
        target = inputs / path.name
        if target.exists() and target.read_bytes() != content:
            raise ValueError("Pinned source changed; explicitly version rather than silently replace")
        target.write_bytes(content)
        records.append({"source": str(path.relative_to(project)), "file": str(target.relative_to(root)),
                        "sha256": hashlib.sha256(content).hexdigest()})
    dms = project / "agents/protease-corrected-sol/results/DMS_corrected.tsv.gz"
    before = hashlib.sha256(dms.read_bytes()).hexdigest()
    data = pd.read_csv(dms, sep="\t")
    assert hashlib.sha256(dms.read_bytes()).hexdigest() == before
    parent = data[data.TEVp.eq("wt")].copy()
    parent.TEVp = "0"
    parent.to_csv(inputs / "parent_from_corrected_DMS.csv", index=False)
    records.append({"source": str(dms.relative_to(project)), "source_sha256": before,
                    "selection": "TEVp == wt, renamed 0 as in authors' control comparison",
                    "file": "inputs/parent_from_corrected_DMS.csv",
                    "sha256": hashlib.sha256((inputs / "parent_from_corrected_DMS.csv").read_bytes()).hexdigest()})
    (inputs / "manifest.json").write_text(json.dumps(records, indent=2) + "\n")


def rho(actual, predicted):
    if len(actual) < 3 or np.std(actual) == 0 or np.std(predicted) == 0:
        return float("nan")
    return float(spearmanr(actual, predicted).statistic)


def main():
    root = Path(__file__).resolve().parent
    if not (root / "inputs/manifest.json").exists():
        pin_inputs(root)
    source = root / "inputs"
    controls = pd.read_csv(source / "controls_corrected.tsv.gz", sep="\t", dtype={"TEVp": str})
    parent = pd.read_csv(source / "parent_from_corrected_DMS.csv", dtype={"TEVp": str})
    data = pd.concat([controls[~controls.TEVp.eq("0")], parent], ignore_index=True)
    assert not data.duplicated(PAIR).any()
    all_n = len(data)
    data = data[data.TEVs.str.fullmatch(f"[{AA}]{{7}}") & data["min"].ge(100)].copy()
    data["hazard"] = [hazard_slope(r[:5], r[5:]) for r in data[COUNTS].to_numpy()]
    f, n = data[[f"F{i}" for i in range(1, 6)]].to_numpy(), data[[f"N{i}" for i in range(1, 6)]].to_numpy()
    frac = (f + 0.5) / (f + n + 1)
    data["endpoint_lift"] = frac[:, -1] - frac[:, 0]
    negative = data[data.TEVp.isin(NEG)]
    complete = negative.groupby("TEVs").TEVp.nunique()
    complete = complete[complete.eq(len(NEG))].index
    background = negative[negative.TEVs.isin(complete)].groupby("TEVs")[["AUC", "hazard", "endpoint_lift"]].mean().add_prefix("negative_")
    active = data[~data.TEVp.isin(NEG)].merge(background, on="TEVs", validate="many_to_one")
    active["background_subtracted_AUC"] = active.AUC - active.negative_AUC
    active["net_hazard"] = active.hazard - active.negative_hazard
    active["net_endpoint_lift"] = active.endpoint_lift - active.negative_endpoint_lift
    rates = pd.read_csv(source / "protease_vitro_relative.tsv", sep="\t", dtype={"TEVp": str})
    positions = {name: i for i, name in enumerate(["P6", "P5", "P4", "P3", "P2", "P1", "P1'"])}
    rates["TEVs"] = [WT[:positions[p]] + a + WT[positions[p] + 1:] for p, a in zip(rates.pos, rates.TEVs_AA, strict=True)]
    canonical = rates[rates.TEVs.eq(WT)]
    assert canonical.groupby("TEVp").vinirel.nunique().eq(1).all()
    rates = rates.groupby(PAIR, as_index=False).vinirel.mean()
    matched = active.merge(rates, on=PAIR, validate="one_to_one")
    scores = ["AUC", "background_subtracted_AUC", "hazard", "net_hazard", "endpoint_lift", "net_endpoint_lift"]
    metrics = []
    for enzyme, panel in matched.groupby("TEVp", sort=True):
        high = panel.vinirel >= panel.vinirel.quantile(0.9)
        for score in scores:
            selected = panel[score] >= panel[score].quantile(0.9)
            metrics.append({"enzyme": enzyme, "score": score, "n": len(panel), "rho": rho(panel.vinirel, panel[score]),
                            "flagged": int(selected.sum()), "true_upper_decile": int(high.sum()),
                            "overlap": int((selected & high).sum()), "precision": float(high[selected].mean()), "prevalence": float(high.mean())})
    metrics = pd.DataFrame(metrics)
    differences = metrics.pivot(index="enzyme", columns="score", values="rho")
    observed = differences.net_hazard - differences.AUC
    motifs = sorted(matched.TEVs.unique())
    lookup = {s: i for i, s in enumerate(motifs)}
    panels = [(p, np.array([lookup[s] for s in p.TEVs])) for _, p in matched.groupby("TEVp", sort=True)]
    rng = np.random.default_rng(20261126)
    samples = []
    for _ in range(3000):
        weights = np.bincount(rng.integers(len(motifs), size=len(motifs)), minlength=len(motifs))
        values = []
        for panel, indices in panels:
            sample_indices = np.repeat(np.arange(len(panel)), weights[indices])
            p = panel.iloc[sample_indices]
            values.append(rho(p.vinirel, p.net_hazard) - rho(p.vinirel, p.AUC))
        samples.append(float(np.mean(values)))
    assert np.isfinite(samples).all()
    interval = np.quantile(samples, [0.025, 0.975])
    result = {"input_control_plus_parent_rows": all_n, "depth_canonical_rows": len(data),
              "complete_negative_motifs": len(complete), "matched_pairs": len(matched),
              "selected_enzyme_panels": int(matched.TEVp.nunique()),
              "matched_motifs": len(motifs), "primary_mean_rho_difference": float(observed.mean()),
              "panels_with_positive_primary_difference": int(observed.gt(0).sum()),
              "shared_motif_bootstrap_interval": interval.tolist(),
              "primary_support_criterion_met": bool(observed.mean() > 0 and observed.gt(0).sum() >= 5 and interval[0] > 0),
              "primary_vs_background_AUC_mean_rho_difference": float((differences.net_hazard - differences.background_subtracted_AUC).mean()),
              "limitations": ["Selected related enzyme panels and shared motifs do not establish population transport.",
                              "Recorder score is not measured kcat or a validated mechanistic cleavage rate.",
                              "Bootstrap is conditional panel sensitivity, not independent biological replication.",
                              "Corrected-population input reconstruction must be independently verified against primary source data."]}
    out = root / "results"
    out.mkdir(exist_ok=True)
    metrics.to_csv(out / "metrics.csv", index=False)
    matched.to_csv(out / "matched_scores.csv", index=False)
    pd.DataFrame({"primary_mean_rho_difference": samples}).to_csv(out / "bootstrap.csv", index=False)
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(metrics.to_string(index=False), flush=True)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
