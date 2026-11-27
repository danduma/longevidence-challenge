"""Rebuild each fold with union-find and candidate-centered normal equations."""
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
MODELS = ["no_label", "published_reference", "rr_ratio", "peer", "rr_peer_midpoint", "rr_inclusive"]


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


def check_summary(rows, saved):
    for group in ("all", "lysate", "purified"):
        subset = [r for r in rows if group == "all" or r["format"] == group]
        candidates = {r["candidate"] for r in subset}
        assert (len(subset), len(candidates), len({r["label"] for r in subset})) == (
            saved[group]["cells"], saved[group]["candidates"], saved[group]["labels"]
        )
        for model in MODELS:
            errors = defaultdict(list)
            for r in subset:
                errors[r["candidate"]].append(abs(r["target"] - r[model]))
            flat = [v for values in errors.values() for v in values]
            metrics = saved[group]["metrics"][model]
            close(math.sqrt(sum(v*v for v in flat) / len(flat)), metrics["rmse"])
            close(sum(flat) / len(flat), metrics["mae"])
            close(math.sqrt(mean(mean(v*v for v in values) for values in errors.values())),
                  metrics["candidate_rmse"])


def check_label_sensitivity(rows, saved):
    for group in ("all", "lysate", "purified"):
        totals = defaultdict(lambda: [0, 0., 0.])
        for row in rows:
            if group != "all" and row["format"] != group:
                continue
            block = totals[row["label"]]
            block[0] += 1
            block[1] += (row["rr_inclusive"] - row["target"]) ** 2
            block[2] += (row["rr_ratio"] - row["target"]) ** 2
        assert totals.keys() == saved[group].keys()
        n, inclusive, ratio = map(sum, zip(*totals.values()))
        for label, (count, a, b) in totals.items():
            close(math.sqrt((inclusive - a) / (n - count))
                  - math.sqrt((ratio - b) / (n - count)), saved[group][label])


def main():
    saved = json.loads((ROOT / "results.json").read_text())
    sources = ROOT / "sources"
    for name in ("assays.xlsx", "Supplementary_Data_3.csv"):
        assert hashlib.sha256((sources / name).read_bytes()).hexdigest() == saved["inputs"][name]
    with (sources / "Supplementary_Data_3.csv").open() as stream:
        meta = {r["internal_id"]: r for r in csv.DictReader(stream)
                if r["rate_lab_corrected"] and not r["internal_id"].startswith("ANC_")}
    cells = defaultdict(list)
    for accession, c, raw, norm, a in workbook(sources / "assays.xlsx"):
        if c in meta:
            assert accession == meta[c]["identifier"]
            cells[c, a].append((float(raw), float(norm)))
    for c, row in meta.items():
        assert sum(len(v) for (i, _), v in cells.items() if i == c) == float(row["N"])
    raw = {k: mean(math.log(r) for r, _ in v) for k, v in cells.items()}
    norm = {k: mean(math.log(n) for _, n in v) for k, v in cells.items()}
    ratio = {k: raw[k] - raw["Rr", k[1]] for k in cells
             if k[0] != "Rr" and ("Rr", k[1]) in cells}
    fmt = {c: "purified" if r["form"] in {"Form 2", "Form 2/3"} else "lysate"
           for c, r in meta.items()}
    assert saved["nonreference_cells"] == sum(c != "Rr" for c, _ in cells)
    assert saved["with_rr"] == len(ratio)

    predictions = {}
    for mode in ("pooled", "within_format"):
        expected = {(r["candidate"], r["label"]): r for r in saved[mode]["folds"]}
        assert len(expected) == len(saved[mode]["folds"])
        actual = {}
        for held in sorted(ratio):
            c, a = held
            allowed = ratio if mode == "pooled" else [k for k in ratio if fmt[k[0]] == fmt[c]]
            train = component(allowed, held)
            if not train:
                continue
            predict = prepare_fit(train, held)
            repeats = [k for k in train if k[0] == c]
            rr = ("Rr", a)
            peer = predict(raw)
            labels = {j for _, j in train}
            inclusive_train = [k for k in raw if k[0] == "Rr" and k[1] in labels]
            inclusive_train.extend(train)
            inclusive = prepare_fit(inclusive_train, held)(raw)
            rr_prediction = mean(ratio[k] for k in repeats) + raw[rr]
            row = dict(candidate=c, label=a, format=fmt[c], target=raw[held],
                       peer=peer, rr_inclusive=inclusive, rr_peer_midpoint=(peer + rr_prediction) / 2,
                       no_label=mean(raw[k] for k in repeats),
                       published_reference=mean(norm[k] for k in repeats) + raw[rr] - norm[rr],
                       rr_ratio=rr_prediction)
            for field, value in row.items():
                if isinstance(value, str):
                    assert value == expected[held][field]
                else:
                    close(value, expected[held][field])
            actual[held] = row
        assert actual.keys() == expected.keys()
        check_summary(list(actual.values()), saved[mode]["summary"])
        predictions[mode] = actual
    check_summary([predictions["pooled"][k] for k in predictions["within_format"]],
                  saved["pooled_on_within_format_support"])
    check_label_sensitivity(predictions["pooled"].values(), saved["label_scoring_sensitivity"])
    print("Independently recomputed fold support, every prediction and all error summaries.")


if __name__ == "__main__":
    main()
