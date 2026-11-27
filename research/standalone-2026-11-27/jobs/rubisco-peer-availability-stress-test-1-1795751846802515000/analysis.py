"""Exact test-label peer masks. Data loading and fit adapted from the preceding holdback audit."""
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



def replay():
    formats, cells, reference = load()
    pool = {key: values for key, values in cells.items() if key[1] in reference}
    outcomes = {key: rates[0] for key, rates in pool.items()}
    outcomes.update({("Rr", label): rates[0] for label, rates in reference.items()})
    rows = []
    for held in sorted(pool):
        candidate, label = held
        if formats[candidate] != "purified":
            continue
        history = [key for key in pool if key[0] == candidate and key != held]
        if not history:
            continue
        direct = reference[label][0] + mean(
            pool[key][0] - reference[key[1]][0] for key in history)
        published = reference[label][0] - reference[label][1] + mean(
            pool[key][1] for key in history)
        peers = sorted(key for key in pool if key[1] == label and key != held)
        elsewhere = [key for key in pool if key[1] != label]

        def fit(keys):
            train = connected_fold(keys, held)
            if not train:
                return None
            assert all(key in train for key in history)
            controls = [("Rr", j) for j in sorted({j for _, j in train})]
            return peer_fit(train + controls, held, outcomes)

        assert fit(elsewhere) is None
        rows.append(dict(candidate=candidate, label=label, target=outcomes[held],
                         ratio=direct, published=published,
                         history=mean(outcomes[key] for key in history),
                         full=fit(elsewhere + peers),
                         one={key[0]: fit(elsewhere + [key]) for key in peers}))
    return rows


def summarize(rows):
    full = [row for row in rows if row["full"] is not None]
    match = [row for row in full if any(v is not None for v in row["one"].values())]

    def rmse(frame, predict, connected=False):
        groups = defaultdict(list)
        for row in frame:
            predictions = predict(row)
            if connected:
                predictions = [p for p in predictions if p is not None]
            groups[row["candidate"]].append(mean((p - row["target"]) ** 2 for p in predictions))
        return sqrt(mean(mean(losses) for losses in groups.values()))

    def blend(row, endpoint):
        return (row["ratio"] + endpoint) / 2

    base = {model: rmse(rows, lambda r: [r[model]]) for model in ("ratio", "published", "history")}
    scores = {model: rmse(match, lambda r: [r[model]]) for model in ("ratio", "published", "history")}
    scores["full_midpoint"] = rmse(match, lambda r: [blend(r, r["full"])])
    scores["one_peer_midpoint"] = rmse(
        match, lambda r: [blend(r, p) if p is not None else None for p in r["one"].values()], True)
    scores["one_peer_fit"] = rmse(match, lambda r: list(r["one"].values()), True)
    for scope, frame in [("prior_cohort", full), ("all_ratio_eligible", rows)]:
        scores[scope] = {
            "ratio": rmse(frame, lambda r: [r["ratio"]]),
            "one_peer_with_fallback": rmse(
                frame, lambda r: [blend(r, p) if p is not None else r["ratio"]
                                  for p in r["one"].values()] or [r["ratio"]]),
        }

    return dict(eligible_cells=len(rows), prior_cells=len(full), matched_cells=len(match),
                eligible_candidates=len({r["candidate"] for r in rows}),
                prior_candidates=len({r["candidate"] for r in full}),
                eligible_labels=len({r["label"] for r in rows}),
                prior_labels=len({r["label"] for r in full}),
                singleton_masks=sum(len(r["one"]) for r in full),
                singleton_disconnected=sum(p is None for r in full for p in r["one"].values()),
                singleton_disconnection_probability=mean(
                    mean(p is None for p in r["one"].values()) for r in full),
                vulnerable_cells=sum(any(p is None for p in r["one"].values()) for r in full),
                zero_disconnected=len(full), rmse=scores, all_eligible_baselines=base)


if __name__ == "__main__":
    rows = replay()
    inputs = {name: hashlib.sha256((ROOT / "sources" / name).read_bytes()).hexdigest()
              for name in ("assays.xlsx", "Supplementary_Data_3.csv")}
    (ROOT / "results.json").write_text(json.dumps(dict(inputs=inputs, summary=summarize(rows),
                                                      rows=rows), indent=2) + "\n")
    print(json.dumps(summarize(rows), indent=2))
