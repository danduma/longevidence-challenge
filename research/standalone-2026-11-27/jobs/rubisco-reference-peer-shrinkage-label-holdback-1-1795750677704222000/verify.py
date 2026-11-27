"""Independently recompute nested tuning with candidate-centered label equations.

Workbook reader and endpoint solver adapted from the preceding standalone transfer audit.
"""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean
from zipfile import ZipFile
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parent


def workbook(path):
    with ZipFile(path) as book:
        strings = ["".join(n.itertext()) for n in ET.fromstring(book.read("xl/sharedStrings.xml"))]
        rows = ET.fromstring(book.read("xl/worksheets/sheet1.xml")).find("{*}sheetData")
    for row in list(rows)[3:]:
        data = {}
        for cell in row:
            value = cell.findtext("{*}v")
            if value is not None:
                data[cell.attrib["r"].rstrip("0123456789")] = (
                    strings[int(value)] if cell.get("t") == "s" else value
                )
        if all(k in data for k in "ABCDE"):
            yield [data[k] for k in "ABCDE"]


def component(keys, held):
    links = {}

    def root(x):
        while x in links:
            x = links[x]
        return x

    for c, a in keys:
        if (c, a) != held:
            u, v = root(("candidate", c)), root(("label", a))
            if u != v:
                links[u] = v
    target = root(("candidate", held[0]))
    if target == root(("label", held[1])):
        return [(c, a) for c, a in keys
                if (c, a) != held and root(("candidate", c)) == target]
    return []


def prepare_fit(train, held):
    blocks = defaultdict(list)
    for c, a in train:
        blocks[c].append(a)
    assays = sorted({a for _, a in train} - {held[1]})
    indices = {a: i for i, a in enumerate(assays)}
    anchor = len(assays)
    indices[held[1]] = anchor
    groups = [([indices[a] for a in labels], [(c, a) for a in labels])
              for c, labels in blocks.items() if len(labels) > 1]
    # Profile out candidate intercepts. Assemble an explicit sparse label matrix.
    matrix = [defaultdict(float) for _ in assays]
    for group, _ in groups:
        weight = 1 / len(group)
        for i in group:
            if i == anchor:
                continue
            matrix[i][i] += 1
            for j in group:
                if j != anchor:
                    matrix[i][j] -= weight
    sparse = [tuple(row.items()) for row in matrix]

    def predict(values):
        b = [0.] * anchor
        for group, keys in groups:
            center = mean(values[k] for k in keys)
            for i, k in zip(group, keys):
                if i != anchor:
                    b[i] += values[k] - center

        # Jacobi-preconditioned conjugate gradient on the grounded label system.
        x = [0.] * anchor
        r = b[:]
        d = [v / matrix[i][i] for i, v in enumerate(r)]
        delta = sum(u * v for u, v in zip(d, r))
        for _ in range(8 * anchor + 1):
            if max(map(abs, r), default=0.) < 1e-11:
                break
            product = [sum(w * d[j] for j, w in row) for row in sparse]
            step = delta / sum(u * v for u, v in zip(d, product))
            x = [u + step * v for u, v in zip(x, d)]
            r = [u - step * v for u, v in zip(r, product)]
            z = [v / matrix[i][i] for i, v in enumerate(r)]
            next_delta = sum(u * v for u, v in zip(z, r))
            d = [u + next_delta / delta * v for u, v in zip(z, d)]
            delta = next_delta
        else:
            raise ArithmeticError("Independent label fit failed to converge")
        return mean(values[held[0], a] - x[indices[a]] for a in blocks[held[0]])

    return predict


def close(actual, saved):
    assert math.isclose(actual, saved, abs_tol=2e-8, rel_tol=0), (actual, saved)


def main():
    recorded = json.loads((ROOT / "results.json").read_text())
    sources = ROOT / "sources"
    for name in ("assays.xlsx", "Supplementary_Data_3.csv"):
        assert hashlib.sha256((sources / name).read_bytes()).hexdigest() == recorded["inputs"][name]
    with (sources / "Supplementary_Data_3.csv").open() as stream:
        meta = {r["internal_id"]: r for r in csv.DictReader(stream)
                if r["rate_lab_corrected"] and not r["internal_id"].startswith("ANC_")}
    measurements = defaultdict(list)
    for accession, candidate, rate, normalized, label in workbook(sources / "assays.xlsx"):
        if candidate in meta:
            assert accession == meta[candidate]["identifier"]
            measurements[candidate, label].append((math.log(float(rate)), math.log(float(normalized))))
    for name, info in meta.items():
        assert sum(len(v) for (c, _), v in measurements.items() if c == name) == float(info["N"])
    raw, norm = ({k: mean(v[i] for v in values) for k, values in measurements.items()}
                 for i in range(2))
    eligible = [k for k in raw if k[0] != "Rr" and ("Rr", k[1]) in raw]
    purified = {c for c, r in meta.items() if r["form"] in {"Form 2", "Form 2/3"}}

    def replay(keys):
        output = {}
        for held in keys:
            c, a = held
            if c not in purified:
                continue
            train = component(keys, held)
            if not train:
                continue
            controls = [("Rr", j) for j in sorted({j for _, j in train})]
            inclusive = prepare_fit(train + controls, held)(raw)
            repeats = [j for i, j in train if i == c]
            rr = ("Rr", a)
            ratio = raw[rr] + mean(raw[c, j] - raw["Rr", j] for j in repeats)
            published = raw[rr] - norm[rr] + mean(norm[c, j] for j in repeats)
            output[held] = dict(target=raw[held], ratio=ratio, inclusive=inclusive,
                                published=published, midpoint=(ratio + inclusive) / 2)
        return output

    def candidate_weights(rows):
        counts = defaultdict(int)
        for c, _ in rows:
            counts[c] += 1
        return {key: 1 / (counts[key[0]] * len(counts)) for key in rows}

    def optimum(rows):
        # Solve the same bounded quadratic by derivative bisection.
        weights = candidate_weights(rows)
        lo, hi = 0., 1.
        for _ in range(48):
            mid = (lo + hi) / 2
            gradient = 0.
            for key, row in rows.items():
                shift = row["inclusive"] - row["ratio"]
                error = row["ratio"] + mid * shift - row["target"]
                gradient += weights[key] * shift * error
            if gradient > 0:
                hi = mid
            else:
                lo = mid
        return (lo + hi) / 2

    labels = sorted({a for _, a in eligible}, key=lambda a: hashlib.sha256(a.encode()).hexdigest())
    actual = replay(eligible)
    expected = {(r["candidate"], r["label"]): r for r in recorded["folds"]}
    assert actual.keys() == expected.keys() and len(expected) == len(recorded["folds"])
    assert recorded["support"] == dict(cells=len(actual), candidates=len({c for c, _ in actual}),
                                       labels=len({a for _, a in actual}))
    assert len(recorded["tuning"]) == 5
    for fold, saved in enumerate(recorded["tuning"]):
        held_labels = set(labels[fold::5])
        training = replay([k for k in eligible if k[1] not in held_labels])
        weight = optimum(training)
        close(weight, saved["weight"])
        test = [k for k in actual if k[1] in held_labels]
        assert (len(training), len({c for c, _ in training}), len(test)) == (
            saved["tuning_cells"], saved["tuning_candidates"], saved["evaluation_cells"])
        for k in test:
            assert expected[k]["fold"] == fold
            row = actual[k]
            row["tuned"] = (1 - weight) * row["ratio"] + weight * row["inclusive"]
    oracle = optimum(actual)
    close(oracle, recorded["oracle_weight"])
    for k, row in actual.items():
        row["oracle"] = (1 - oracle) * row["ratio"] + oracle * row["inclusive"]
        for key, value in row.items():
            close(value, expected[k][key])
    weights = candidate_weights(actual)
    scores = {}
    for model in ("ratio", "inclusive", "published", "midpoint", "tuned", "oracle"):
        scores[model] = math.sqrt(sum(weights[k] * (r[model] - r["target"]) ** 2
                                      for k, r in actual.items()))
        close(scores[model], recorded["rmse"][model])
    print("Verified all endpoint and blend predictions, five tuning weights and candidate RMSEs.")
    print(json.dumps(scores, sort_keys=True))


if __name__ == "__main__":
    main()
