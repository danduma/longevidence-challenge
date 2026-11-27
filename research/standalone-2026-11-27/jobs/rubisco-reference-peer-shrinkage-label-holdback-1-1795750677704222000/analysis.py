"""Nested label-block tuning for purified candidate repeat-cell prediction.

Workbook reader and endpoint solver adapted from the preceding standalone transfer audit.
"""
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


def endpoint_rows(pool, formats, reference):
    outcomes = {key: rates[0] for key, rates in pool.items()}
    outcomes.update({("Rr", label): rates[0] for label, rates in reference.items()})
    rows = []
    for candidate, label in sorted(pool):
        if formats[candidate] != "purified":
            continue
        held = candidate, label
        train = connected_fold(pool, held)
        if not train:
            continue
        repeats = [key for key in train if key[0] == candidate]
        rr_raw, rr_normalized = reference[label]
        direct = rr_raw + mean(pool[key][0] - reference[key[1]][0] for key in repeats)
        published = rr_raw - rr_normalized + mean(pool[key][1] for key in repeats)
        controls = [("Rr", j) for j in sorted({j for _, j in train})]
        inclusive = peer_fit(train + controls, held, outcomes)
        rows.append(dict(candidate=candidate, label=label, target=outcomes[held],
                         ratio=direct, inclusive=inclusive, published=published,
                         midpoint=(direct + inclusive) / 2))
    return rows


def by_candidate(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["candidate"]].append(row)
    return groups.values()


def fit_weight(rows):
    groups = list(by_candidate(rows))
    cross = 0.
    squared = 0.
    for group in groups:
        shifts = [r["inclusive"] - r["ratio"] for r in group]
        cross += mean(d * (r["target"] - r["ratio"]) for d, r in zip(shifts, group))
        squared += mean(d * d for d in shifts)
    return min(1., max(0., cross / squared)) if squared else 0.


def scores(rows):
    return {model: sqrt(mean(mean((r[model] - r["target"]) ** 2 for r in group)
                             for group in by_candidate(rows)))
            for model in ("ratio", "inclusive", "published", "midpoint", "tuned", "oracle")}


def run():
    formats, cells, reference = load()
    pool = {k: v for k, v in cells.items() if k[1] in reference}
    evaluation = endpoint_rows(pool, formats, reference)
    labels = sorted({k[1] for k in pool}, key=lambda s: hashlib.sha256(s.encode()).digest())
    fold = {label: index % 5 for index, label in enumerate(labels)}
    tuning = []
    for index in range(5):
        # Purge the entire evaluation block before making any tuning predictions.
        training_pool = {k: v for k, v in pool.items() if fold[k[1]] != index}
        training_reference = {j: v for j, v in reference.items() if j in fold and fold[j] != index}
        inner = endpoint_rows(training_pool, formats, training_reference)
        weight = fit_weight(inner)
        test = [row for row in evaluation if fold[row["label"]] == index]
        for row in test:
            row["fold"] = index
            row["tuned"] = row["ratio"] + weight * (row["inclusive"] - row["ratio"])
        tuning.append(dict(weight=weight, tuning_cells=len(inner),
                           tuning_candidates=len(list(by_candidate(inner))), evaluation_cells=len(test)))
    oracle = fit_weight(evaluation)
    for row in evaluation:
        row["oracle"] = row["ratio"] + oracle * (row["inclusive"] - row["ratio"])
    return dict(inputs={name: hashlib.sha256((ROOT / "sources" / name).read_bytes()).hexdigest()
                        for name in ("assays.xlsx", "Supplementary_Data_3.csv")},
                tuning=tuning, oracle_weight=oracle, folds=evaluation, rmse=scores(evaluation),
                support=dict(cells=len(evaluation), candidates=len(list(by_candidate(evaluation))),
                             labels=len({r["label"] for r in evaluation})))


if __name__ == "__main__":
    (ROOT / "results.json").write_text(json.dumps(run(), indent=2) + "\n")
