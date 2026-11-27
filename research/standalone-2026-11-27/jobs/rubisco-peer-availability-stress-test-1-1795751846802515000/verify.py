"""Independent replay with union-find and candidate-centered label equations, adapted from the holdback verifier."""
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



def rebuild():
    sources = ROOT / "sources"
    with (sources / "Supplementary_Data_3.csv").open() as stream:
        meta = {r["internal_id"]: r for r in csv.DictReader(stream)
                if r["rate_lab_corrected"] and not r["internal_id"].startswith("ANC_")}
    raw_log, norm_log = defaultdict(list), defaultdict(list)
    for accession, candidate, rate, normalized, label in workbook(sources / "assays.xlsx"):
        if candidate in meta:
            assert accession == meta[candidate]["identifier"]
            raw_log[candidate, label].append(math.log(float(rate)))
            norm_log[candidate, label].append(math.log(float(normalized)))
    for candidate, info in meta.items():
        assert sum(len(v) for (c, _), v in raw_log.items() if c == candidate) == float(info["N"])
    raw = {k: mean(v) for k, v in raw_log.items()}
    norm = {k: mean(v) for k, v in norm_log.items()}
    eligible = {k for k in raw if k[0] != "Rr" and ("Rr", k[1]) in raw}
    output = {}
    for held in sorted(eligible):
        c, a = held
        if meta[c]["form"] not in {"Form 2", "Form 2/3"}:
            continue
        repeats = sorted(j for i, j in eligible if i == c and j != a)
        if not repeats:
            continue

        def fit(keep):
            train = component(sorted(keep), held)
            if train:
                controls = [("Rr", j) for j in sorted({j for _, j in train})]
                return prepare_fit(train + controls, held)(raw)
            return None

        same = {key for key in eligible if key[1] == a} - {held}
        rest = {key for key in eligible if key[1] != a}
        assert fit(rest) is None
        output[held] = dict(
            target=raw[held],
            ratio=mean(raw[c, j] + raw["Rr", a] - raw["Rr", j] for j in repeats),
            published=mean(norm[c, j] + raw["Rr", a] - norm["Rr", a] for j in repeats),
            history=mean(raw[c, j] for j in repeats),
            full=fit(rest | same),
            one={i: fit(rest | {(i, a)}) for i, _ in sorted(same)})
    return output


def check_summary(rows, saved):
    prior = {k for k, r in rows.items() if r["full"] is not None}
    matched = {k for k, r in rows.items() if any(p is not None for p in r["one"].values())}
    assert prior == matched  # Any peer path ends at one eligible same-label peer.
    assert saved["eligible_cells"] == len(rows)
    assert saved["prior_cells"] == saved["matched_cells"] == saved["zero_disconnected"] == len(prior)
    for prefix, keys in [("eligible", rows), ("prior", prior)]:
        assert saved[prefix + "_candidates"] == len({c for c, a in keys})
        assert saved[prefix + "_labels"] == len({a for c, a in keys})
    mask_sizes = [len(rows[k]["one"]) for k in prior]
    failures = [list(rows[k]["one"].values()).count(None) for k in prior]
    assert sum(mask_sizes) == saved["singleton_masks"]
    assert sum(failures) == saved["singleton_disconnected"]
    assert sum(n > 0 for n in failures) == saved["vulnerable_cells"]
    close(mean(n / d for n, d in zip(failures, mask_sizes)),
          saved["singleton_disconnection_probability"])

    def score(keys, model):
        counts = defaultdict(int)
        for c, a in keys:
            counts[c] += 1
        total = 0.
        for key in keys:
            r = rows[key]
            points = [p for p in r["one"].values() if p is not None]
            if model == "full_midpoint":
                predictions = [(r["ratio"] + r["full"]) * 0.5]
            elif model == "one_peer_midpoint":
                predictions = [(r["ratio"] + p) * 0.5 for p in points]
            elif model == "one_peer_fit":
                predictions = points
            elif model == "one_peer_with_fallback":
                predictions = [(r["ratio"] + (r["ratio"] if p is None else p)) * 0.5
                               for p in r["one"].values()] or [r["ratio"]]
            else:
                predictions = [r[model]]
            total += sum((p - r["target"]) ** 2 for p in predictions) / (
                len(predictions) * counts[key[0]] * len(counts))
        return math.sqrt(total)

    scores = saved["rmse"]
    for model in ("ratio", "published", "history", "full_midpoint",
                  "one_peer_midpoint", "one_peer_fit"):
        close(score(matched, model), scores[model])
    for scope, keys in [("prior_cohort", prior), ("all_ratio_eligible", rows)]:
        for model in ("ratio", "one_peer_with_fallback"):
            close(score(keys, model), scores[scope][model])
    for model in ("ratio", "published", "history"):
        close(score(rows, model), saved["all_eligible_baselines"][model])
    return scores


def main():
    saved = json.loads((ROOT / "results.json").read_text())
    for name in ("assays.xlsx", "Supplementary_Data_3.csv"):
        assert hashlib.sha256((ROOT / "sources" / name).read_bytes()).hexdigest() == saved["inputs"][name]
    actual = rebuild()
    recorded = {(r["candidate"], r["label"]): r for r in saved["rows"]}
    assert actual.keys() == recorded.keys() and len(recorded) == len(saved["rows"])
    for key, row in actual.items():
        expected = recorded[key]
        for model in ("target", "ratio", "published", "history", "full"):
            if row[model] is None:
                assert expected[model] is None
            else:
                close(row[model], expected[model])
        assert row["one"].keys() == expected["one"].keys()
        for peer, prediction in row["one"].items():
            if prediction is None:
                assert expected["one"][peer] is None
            else:
                close(prediction, expected["one"][peer])
    scores = check_summary(actual, saved["summary"])
    prior = json.loads((ROOT / "sources/prior_results.json").read_text())
    prior_rows = {(r["candidate"], r["label"]): r for r in prior["folds"]}
    assert prior_rows.keys() == {k for k, r in actual.items() if r["full"] is not None}
    for key, row in prior_rows.items():
        close(actual[key]["full"], row["inclusive"])
        close(actual[key]["ratio"], row["ratio"])
    print(json.dumps(scores, indent=2))


if __name__ == "__main__":
    main()
