"""Compute candidate contrasts by profiling assay effects out of least squares."""
import csv
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

from read_assays import observations

ROOT = Path(__file__).resolve().parent
LEFT, RIGHT = "RBC2_36", "RBC8_39"


def inputs():
    with (ROOT / "sources/Supplementary_Data_3.csv").open(newline="") as source:
        eligible = {r["internal_id"]: (r["identifier"], int(float(r["N"])))
                    for r in csv.DictReader(source)
                    if r["rate_lab_corrected"] and not r["internal_id"].startswith("ANC_")}
    rows = []
    counts = defaultdict(int)
    for accession, candidate, raw, normalized, label in observations():
        if candidate in eligible:
            assert eligible[candidate][0] == accession
            counts[candidate] += 1
            rows.append((candidate, label, raw, normalized))
    assert dict(counts) == {c: n for c, (_, n) in eligible.items()}
    return rows


def graph(rows):
    adjacency = defaultdict(set)
    for candidate, label, *_ in rows:
        adjacency[candidate].add(label)
        adjacency[label].add(candidate)
    return {k: sorted(v) for k, v in adjacency.items()}


def search(adjacency, omit=()):
    path = {LEFT: [LEFT]}
    queue = [LEFT]
    for node in queue:
        for neighbor in adjacency[node]:
            if neighbor not in path and neighbor not in omit:
                path[neighbor] = path[node] + [neighbor]
                queue.append(neighbor)
    return path


def fit(rows, component, column, log_scale):
    blocks = defaultdict(list)
    names = sorted({c for c, a, *_ in rows if c in component and a in component} - {RIGHT})
    index = {c: i for i, c in enumerate(names)}
    anchor = len(names)
    for row in rows:
        c, a = row[:2]
        if a in component:
            value = math.log(row[column]) if log_scale else row[column]
            blocks[a].append((index.get(c, anchor), value))
    groups = [block for block in blocks.values() if len(block) > 1]
    rhs = [0.] * (anchor + 1)
    for block in groups:
        center = sum(y for _, y in block) / len(block)
        for i, y in block:
            rhs[i] += y - center
    rhs[anchor] = 0.

    def multiply(vector):
        result = [0.] * len(vector)
        for block in groups:
            center = sum(vector[i] for i, _ in block) / len(block)
            for i, _ in block:
                result[i] += vector[i] - center
        result[anchor] = 0.
        return result

    solution = [0.] * len(rhs)
    residual = rhs[:]
    direction = rhs[:]
    energy = sum(v * v for v in residual)
    for _ in range(4 * len(rhs)):
        if energy < 1e-20:
            break
        product = multiply(direction)
        step = energy / sum(a * b for a, b in zip(direction, product))
        solution = [x + step * d for x, d in zip(solution, direction)]
        residual = [r - step * p for r, p in zip(residual, product)]
        next_energy = sum(r * r for r in residual)
        direction = [r + (next_energy / energy) * d for r, d in zip(residual, direction)]
        energy = next_energy
    else:
        raise ArithmeticError("Candidate effects failed to converge")
    return solution[index[LEFT]]


MODELS = {
    "log_raw": (2, True),
    "log_normalized": (3, True),
    "raw": (2, False),
    "normalized": (3, False),
}


def estimate(all_rows, adjacency, omit=()):
    rows = [r for r in all_rows if r[1] not in omit]
    component = search(adjacency, omit)
    baseline = mean(r[3] for r in rows if r[0] == LEFT) - mean(r[3] for r in rows if r[0] == RIGHT)
    controls = defaultdict(list)
    for c, label, raw, _ in rows:
        if c == "Rr":
            controls[label].append(raw)
    paired = None
    if controls:
        paired = (mean(math.log(raw / mean(controls[label]))
                       for c, label, raw, _ in rows if c == LEFT)
                  - mean(math.log(raw / mean(controls[label]))
                         for c, label, raw, _ in rows if c == RIGHT))
    return {"normalized_mean": baseline, "paired_reference_log_ratio": paired, **{
        name: fit(rows, component, column, log_scale) if RIGHT in component else None
        for name, (column, log_scale) in MODELS.items()}}


def run():
    extant = inputs()
    peers = [r for r in extant if r[0] != "Rr"]
    output = {"workbook_sha256": hashlib.sha256((ROOT / "sources/assays.xlsx").read_bytes()).hexdigest()}
    for scope, rows in [("peers", peers), ("with_reference", extant)]:
        adjacency = graph(rows)
        labels = sorted({r[1] for r in rows})
        paths = search(adjacency)
        single = {a: estimate(rows, adjacency, (a,)) for a in labels}
        cuts = [list(pair) for pair in itertools.combinations(labels, 2)
                if RIGHT not in search(adjacency, pair)]
        output[scope] = {
            "rows": len(rows), "candidates": len({r[0] for r in rows}), "labels": len(labels),
            "component_candidates": len(set(paths) - set(labels)),
            "component_labels": len(set(paths) & set(labels)),
            "direct_shared_labels": sorted(set(adjacency[LEFT]) & set(adjacency[RIGHT])),
            "shortest_path": paths[RIGHT],
            "full": estimate(rows, adjacency),
            "omit_one": single,
            "disconnecting_pairs": cuts,
        }
    return output


if __name__ == "__main__":
    (ROOT / "results.json").write_text(json.dumps(run(), indent=2) + "\n")
