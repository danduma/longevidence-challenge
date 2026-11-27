"""Independently profile candidate effects and solve for assay offsets.

Uses union-find for connectivity and dense Cholesky for the reduced assay
system. The producing analysis uses breadth-first search and candidate effects.
"""
import csv
import hashlib
from itertools import combinations
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

from read_assays import observations

ROOT = Path(__file__).resolve().parent
A, B = "RBC2_36", "RBC8_39"


def close(actual, expected):
    if actual is None or expected is None:
        assert actual is expected, (actual, expected)
    else:
        assert math.isclose(actual, expected, abs_tol=2e-8, rel_tol=1e-8), (actual, expected)


def components(edges, removed):
    parent = {}

    def root(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for candidate, label in edges:
        if label in removed:
            continue
        parent.setdefault(candidate, candidate)
        parent.setdefault(label, label)
        parent[root(candidate)] = root(label)
    if A not in parent:
        return set()
    return {vertex for vertex in parent if root(vertex) == root(A)}


def cholesky(matrix):
    for i, row in enumerate(matrix):
        for j in range(i):
            previous = matrix[j]
            row[j] = (row[j] - sum(a*b for a, b in zip(row[:j], previous[:j]))) / previous[j]
        row[i] = math.sqrt(row[i] - sum(x*x for x in row[:i]))


def solve(matrix, rhs):
    x = []
    for i, row in enumerate(matrix):
        x.append((rhs[i] - sum(a*b for a, b in zip(row[:i], x))) / row[i])
    for i in reversed(range(len(x))):
        x[i] = (x[i] - sum(matrix[j][i]*x[j] for j in range(i+1, len(x)))) / matrix[i][i]
    return x


def estimates(rows, component):
    groups = defaultdict(list)
    for c, label, raw, normalized in rows:
        groups[c].append((label, raw, normalized))
    values = {"normalized_mean": mean(row[2] for row in groups[A])
              - mean(row[2] for row in groups[B])}
    values["paired_reference_log_ratio"] = None
    if "Rr" in groups:
        reference = defaultdict(list)
        for label, raw, _ in groups["Rr"]:
            reference[label].append(raw)
        values["paired_reference_log_ratio"] = (
            mean(math.log(raw / mean(reference[label])) for label, raw, _ in groups[A])
            - mean(math.log(raw / mean(reference[label])) for label, raw, _ in groups[B]))
    models = {"raw": (1, float), "log_raw": (1, math.log),
              "normalized": (2, float), "log_normalized": (2, math.log)}
    if B not in component:
        return {**values, **dict.fromkeys(models)}
    assays = sorted({label for c in groups if c in component for label, *_ in groups[c]})
    index = {label: i for i, label in enumerate(assays[:-1])}
    n = len(index)
    matrix = [[0.] * n for _ in range(n)]
    active = [group for c, group in groups.items() if c in component]
    for group in active:
        counts = Counter(index[row[0]] for row in group if row[0] in index)
        for i, weight in counts.items():
            matrix[i][i] += weight
            for j, weight_j in counts.items():
                matrix[i][j] -= weight * weight_j / len(group)
    cholesky(matrix)
    for model, (column, transform) in models.items():
        rhs = [0.] * n
        for group in active:
            y = [transform(row[column]) for row in group]
            center = mean(y)
            for (label, *_), rate in zip(group, y):
                if label in index:
                    rhs[index[label]] += rate - center
        offsets = dict(zip(assays, solve(matrix, rhs) + [0.]))
        values[model] = (mean(transform(row[column]) - offsets[row[0]] for row in groups[A])
                         - mean(transform(row[column]) - offsets[row[0]] for row in groups[B]))
    return values


def main():
    recorded = json.loads((ROOT / "results.json").read_text())
    assert recorded["workbook_sha256"] == hashlib.sha256((ROOT / "sources/assays.xlsx").read_bytes()).hexdigest()
    with (ROOT / "sources/Supplementary_Data_3.csv").open(newline="") as source:
        expected_counts = {(row["identifier"], row["internal_id"]): int(float(row["N"]))
                           for row in csv.DictReader(source)
                           if row["rate_lab_corrected"] and not row["internal_id"].startswith("ANC_")}
    data = list(observations())
    assert Counter((a, c) for a, c, *_ in data if (a, c) in expected_counts) == expected_counts
    for scope, with_reference in [("peers", False), ("with_reference", True)]:
        reported = recorded[scope]
        rows = [(c, label, raw, norm) for accession, c, raw, norm, label in data
                if (accession, c) in expected_counts and (with_reference or c != "Rr")]
        edges = {(c, label) for c, label, *_ in rows}
        labels = sorted({label for _, label in edges})
        candidates = {c for c, _ in edges}
        full = components(edges, ())
        assert reported["rows"] == len(rows)
        assert reported["candidates"] == len(candidates)
        assert reported["labels"] == len(labels)
        assert reported["component_candidates"] == len(full & candidates)
        assert reported["component_labels"] == len(full.intersection(labels))
        shared = sorted({label for c, label in edges if c == A} & {label for c, label in edges if c == B})
        assert reported["direct_shared_labels"] == shared == []
        path = reported["shortest_path"]
        assert path[0] == A and path[-1] == B and len(path) % 2 == 1
        assert all((c, label) in edges and (d, label) in edges
                   for c, label, d in zip(path[::2], path[1::2], path[2::2]))
        reached = {A}
        for distance in range(1, len(path)):
            touched = {label for c, label in edges if c in reached}
            reached = {c for c, label in edges if label in touched}
            if B in reached:
                break
        assert 2 * distance == len(path) - 1
        pairs = [list(pair) for pair in combinations(labels, 2) if B not in components(edges, pair)]
        assert pairs == reported["disconnecting_pairs"] == [["CA187", "CA188"]]
        assert set(reported["omit_one"]) == set(labels)
        for label in [None, *labels]:
            removed = () if label is None else (label,)
            component = components(edges, removed)
            assert B in component
            computed = estimates([row for row in rows if row[1] not in removed], component)
            expected = reported["full"] if label is None else reported["omit_one"][label]
            assert computed.keys() == expected.keys()
            for name, value in computed.items():
                close(value, expected[name])
            assert (computed["log_raw"] > 0) == with_reference
            assert (computed["log_normalized"] > 0) == with_reference
            if label is None:
                assert computed["normalized_mean"] > 0
    print("Verified paths, all single- and two-label connectivity checks, and all contrast estimates.")


if __name__ == "__main__":
    main()
