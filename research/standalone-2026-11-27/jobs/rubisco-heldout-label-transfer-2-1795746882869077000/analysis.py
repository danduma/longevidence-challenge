"""Leave out whole candidate–label cells and score predictions on log raw rates."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from math import log, sqrt
from pathlib import Path
from statistics import mean

from read_assays import observations

ROOT = Path(__file__).resolve().parent
FORMATS = {
    "Form 1": "lysate", "Form 1'": "lysate", "Form 3": "lysate",
    "Form 2": "purified", "Form 2/3": "purified",
}
MODELS = ("no_label", "published_reference", "rr_ratio", "peer", "rr_peer_midpoint", "rr_inclusive")


def load():
    with (ROOT / "sources/Supplementary_Data_3.csv").open() as stream:
        metadata = {r["internal_id"]: r for r in csv.DictReader(stream)
                    if r["rate_lab_corrected"] and not r["internal_id"].startswith("ANC_")}
    formats = {c: FORMATS[r["form"]] for c, r in metadata.items()}
    pairs = defaultdict(list)
    counts = Counter()
    for accession, c, raw, normalized, a in observations():
        if c in metadata:
            assert metadata[c]["identifier"] == accession
            counts[c] += 1
            pairs[c, a].append((log(raw), log(normalized)))
    assert counts == {c: int(float(r["N"])) for c, r in metadata.items()}
    cells = {key: tuple(map(mean, zip(*values))) for key, values in pairs.items()}
    reference = {a: rates for (c, a), rates in cells.items() if c == "Rr"}
    return formats, {key: y for key, y in cells.items() if key[0] != "Rr"}, reference


def connected_fold(cells, held):
    c, a = held
    adjacency = defaultdict(list)
    for i, j in cells:
        if (i, j) != held:
            adjacency[0, i].append((1, j))
            adjacency[1, j].append((0, i))
    queue = [(0, c)]
    seen = set(queue)
    for node in queue:
        for neighbor in adjacency[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    if (1, a) in seen:
        return [(i, j) for i, j in cells if (i, j) != held and (0, i) in seen]
    return []


def cg(multiply, rhs):
    x = [0.] * len(rhs)
    r = rhs[:]
    d = r[:]
    energy = sum(v * v for v in r)
    for _ in range(6 * len(rhs) + 1):
        if energy < 1e-22:
            return x
        product = multiply(d)
        step = energy / sum(a * b for a, b in zip(d, product))
        x = [a + step * b for a, b in zip(x, d)]
        r = [a - step * b for a, b in zip(r, product)]
        updated = sum(v * v for v in r)
        d = [a + updated / energy * b for a, b in zip(r, d)]
        energy = updated
    raise ArithmeticError("Least-squares solve failed to converge")


def peer_fit(train, held, outcomes):
    # Eliminate label effects by centering within each training label.
    c, a = held
    names = [c] + sorted({i for i, _ in train} - {c})
    index = {name: i for i, name in enumerate(names)}
    blocks = defaultdict(list)
    for key in train:
        i, j = key
        blocks[j].append((index[i], outcomes[key]))
    rhs = [0.] * len(names)
    for block in blocks.values():
        center = sum(y for _, y in block) / len(block)
        for i, y in block:
            rhs[i] += y - center
    rhs[0] = 0.

    def multiply(x):
        result = [0.] * len(names)
        for block in blocks.values():
            center = sum(x[i] for i, _ in block) / len(block)
            for i, _ in block:
                result[i] += x[i] - center
        result[0] = 0.
        return result

    effects = cg(multiply, rhs)
    return mean(y - effects[i] for i, y in blocks[a])


def score(folds):
    report = {}
    for name, subset in (
        ("all", folds),
        ("lysate", [r for r in folds if r["format"] == "lysate"]),
        ("purified", [r for r in folds if r["format"] == "purified"]),
    ):
        by_candidate = defaultdict(list)
        for r in subset:
            by_candidate[r["candidate"]].append(r)
        metrics = {}
        for model in MODELS:
            def square(row):
                return (row[model] - row["target"]) ** 2

            metrics[model] = {
                "rmse": sqrt(mean(map(square, subset))),
                "mae": mean(abs(r[model] - r["target"]) for r in subset),
                "candidate_rmse": sqrt(mean(mean(map(square, rows))
                                            for rows in by_candidate.values())),
            }
        report[name] = {"cells": len(subset), "candidates": len(by_candidate),
                        "labels": len({r["label"] for r in subset}), "metrics": metrics}
    return report


def label_sensitivity(folds):
    """Omit one label from scoring while leaving all trained predictions fixed."""
    result = {}
    for group in ("all", "lysate", "purified"):
        rows = [r for r in folds if group == "all" or r["format"] == group]
        differences = {}
        for label in sorted({r["label"] for r in rows}):
            retained = [r for r in rows if r["label"] != label]
            rmse = {model: sqrt(mean((r[model] - r["target"]) ** 2 for r in retained))
                    for model in ("rr_inclusive", "rr_ratio")}
            differences[label] = rmse["rr_inclusive"] - rmse["rr_ratio"]
        result[group] = differences
    return result


def run():
    formats, cells, reference = load()
    raw = {key: rates[0] for key, rates in cells.items()}
    with_reference = {**raw, **{("Rr", a): y[0] for a, y in reference.items()}}
    ratio = {key: rate - reference[key[1]][0]
             for key, rate in raw.items() if key[1] in reference}
    output = {}
    control_cells = {key: cells[key] for key in ratio}
    for mode in ("pooled", "within_format"):
        fits = []
        for key, target in sorted(ratio.items()):
            c, a = key
            pool = control_cells if mode == "pooled" else {
                k: v for k, v in control_cells.items() if formats[k[0]] == formats[c]
            }
            train = connected_fold(pool, key)
            if not train:
                continue
            candidate_train = [k for k in train if k[0] == c]
            rr_raw, rr_norm = reference[a]
            peer = peer_fit(train, key, raw)
            reference_train = train + [("Rr", label) for label in sorted({j for _, j in train})]
            rr_inclusive = peer_fit(reference_train, key, with_reference)
            rr_prediction = mean(ratio[k] for k in candidate_train) + rr_raw
            fits.append({
                "candidate": c, "label": a, "format": formats[c],
                "target": raw[key],
                "no_label": mean(raw[k] for k in candidate_train),
                "published_reference": mean(cells[k][1] for k in candidate_train) + rr_raw - rr_norm,
                "rr_ratio": rr_prediction,
                "peer": peer,
                "rr_inclusive": rr_inclusive,
                "rr_peer_midpoint": (peer + rr_prediction) / 2,
            })
        output[mode] = {"summary": score(fits), "folds": fits}
    restricted = {(r["candidate"], r["label"]) for r in output["within_format"]["folds"]}
    output["pooled_on_within_format_support"] = score(
        [r for r in output["pooled"]["folds"] if (r["candidate"], r["label"]) in restricted]
    )
    output["label_scoring_sensitivity"] = label_sensitivity(output["pooled"]["folds"])
    return {
        "inputs": {name: hashlib.sha256((ROOT / "sources" / name).read_bytes()).hexdigest()
                   for name in ("assays.xlsx", "Supplementary_Data_3.csv")},
        "nonreference_cells": len(cells), "with_rr": len(ratio), **output,
    }


if __name__ == "__main__":
    (ROOT / "results.json").write_text(json.dumps(run(), indent=2) + "\n")
