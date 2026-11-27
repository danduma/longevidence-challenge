"""Independently recover the result from raw-scale master rows."""
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
from statistics import median, stdev

from read_assays import observations

ROOT = Path(__file__).resolve().parent


def load(name):
    with (ROOT / "sources" / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def state(rows):
    summary = dict.fromkeys(("total", "numeric", "with_dispersion", "low_flag", "unresolved"), 0)
    for row in rows:
        summary["total"] += 1
        if row["rate_lab"]:
            summary["numeric"] += 1
            if row["std_rate_lab"]:
                assert int(float(row["N_rate_lab"])) >= 2
                summary["with_dispersion"] += 1
        elif row["active<0.5s-1"]:
            assert float(row["active<0.5s-1"]) == 1
            summary["low_flag"] += 1
        else:
            summary["unresolved"] += 1
    return summary


def verify():
    saved = json.loads((ROOT / "results.json").read_text())
    master = {r["internal_id"]: r for r in load("rubisco_library.csv") if r["from_lab"] == "1.0"}
    supplement = load("Supplementary_Data_3.csv")
    added = [r["internal_id"] for r in load("fast_formII_clusters_expansion.csv") if r["expanded"] == "1"]
    assert len(added) == len(set(added)) == 44 and set(added) <= master.keys()
    assert len(supplement) == len(master) == 469
    assert {r["internal_id"] for r in supplement} == master.keys()
    for r in supplement:
        m = master[r["internal_id"]]
        for field, other in (("identifier", "identifier"), ("form", "form"),
                             ("rate_lab_corrected", "rate_lab_corrected"),
                             ("N", "N_rate_lab"), ("stdev", "std_rate_lab_corrected"),
                             ("active<0.5s-1", "active<0.5s-1")):
            assert r[field] == m[other], (field, r["internal_id"])

    factor = median(float(r["rate_lab_corrected"]) / float(r["rate_lab"]) for r in master.values()
                    if r["rate_lab"] and float(r["rate_lab"]) > 0)
    for r in master.values():
        for field in ("rate_lab", "std_rate_lab"):
            if r[field]:
                assert math.isclose(float(r[field]) * factor, float(r[field + "_corrected"]))
    ancestors = {i for i in master if i.startswith("ANC_")}
    cohorts = {"representative": master.keys() - set(added) - ancestors,
               "expansion": set(added), "ancestors": ancestors}
    for name, ids in cohorts.items():
        assert state(master[i] for i in ids) == saved["cohorts"][name], name
    assert saved["cohorts"]["representative"]["numeric"] == 221

    extant = master.keys() - ancestors
    pools = {"representative": cohorts["representative"], "extant": extant,
             "extant_form_II": {i for i in extant if master[i]["form"] == "Form 2"}, "all": master}
    for pool, ids in pools.items():
        expected = saved["pools"][pool]
        assert state(master[i] for i in ids) == expected["counts"], pool
        numeric = {i: float(master[i]["rate_lab"]) for i in ids if master[i]["rate_lab"]}
        leader = sorted(numeric, key=numeric.get)[-1]
        assert expected["point_leader"] == leader
        for model in ("SD", "SD_over_sqrt_N"):
            widths = {}
            for i in numeric:
                row = master[i]
                if row["std_rate_lab"]:
                    scale = math.sqrt(float(row["N_rate_lab"])) if model != "SD" else 1
                    widths[i] = float(row["std_rate_lab"]) / scale
            verify_corners(numeric, widths, leader, factor, expected[model])
    verify_replay(master, pools, saved["assay_replay"])
    assert saved["pools"]["extant"]["SD"]["minimax_choice"] != saved["pools"]["extant"]["SD_over_sqrt_N"]["minimax_choice"]
    print("Verified cohort membership, calibration, leader breakpoints and adversarial choices for all pools.")


def verify_replay(master, pools, expected):
    assays = defaultdict(list)
    for identifier, candidate, raw, normalized, label in observations():
        assert identifier == master[candidate]["identifier"]
        assays[candidate].append((round(raw * 10), round(normalized * 10), label))
    for candidate, entries in assays.items():
        row = master[candidate]
        if not row["rate_lab"]:
            assert row["active<0.5s-1"] == "1.0"
            continue
        assert len(entries) == float(row["N_rate_lab"])
        values = [entry[1] / 10 for entry in entries]
        assert abs(sum(values) / len(values) - float(row["rate_lab"])) <= 0.10000001
        if row["std_rate_lab"]:
            assert abs(stdev(values) - float(row["std_rate_lab"])) <= 0.10000001
    for key, pool in pools.items():
        expected_pool = expected[key]
        ids = sorted(i for i in pool if master[i]["rate_lab"])
        counts = {i: len(assays[i]) for i in ids}
        raw = {i: sum(r for r, _, _ in assays[i]) for i in ids}
        norm = {i: sum(n for _, n, _ in assays[i]) for i in ids}
        reductions = defaultdict(Counter)
        numbers = defaultdict(Counter)
        for i in ids:
            for _, n, block in assays[i]:
                reductions[block][i] += n
                numbers[block][i] += 1
        winner = max(ids, key=lambda i: norm[i] / counts[i])
        assert winner == expected_pool["normalized_leader"]
        assert max(ids, key=lambda i: raw[i] / counts[i]) == expected_pool["raw_leader"]
        changed = {}
        for block, removed in reductions.items():
            remaining = {i: counts[i] - numbers[block][i] for i in ids}
            leading = max((i for i in ids if remaining[i]),
                          key=lambda i: (norm[i] - removed[i]) / remaining[i])
            if leading != winner:
                changed[block] = leading
        assert len(reductions) == expected_pool["assay_labels"]
        assert changed == expected_pool["changed_after_omission"], key
    a, b = "RBC2_36", "RBC8_39"
    shared = {row[2] for row in assays[a]}.intersection(row[2] for row in assays[b])
    gap = (sum(x[1] for x in assays[a]) / len(assays[a])
           - sum(x[1] for x in assays[b]) / len(assays[b])) / 10
    assert sorted(shared) == expected["extant_pair"]["shared_assays"]
    assert abs(gap - expected["extant_pair"]["normalized_gap"]) < 1e-9


def verify_corners(scores, widths, leader, factor, expected):
    def rivals_at(width):
        return [j for j in widths if j != leader and scores[j] + widths[j] * width
                >= max(0, scores[leader] - widths[leader] * width)]

    low, high = 0, 2
    assert not rivals_at(low) and rivals_at(high)
    for _ in range(45):
        middle = (low + high) / 2
        if rivals_at(middle):
            high = middle
        else:
            low = middle
    assert abs(high - expected["first_tie_width"]) < 1e-9
    assert rivals_at(high) == [expected["first_tie_rival"]]
    possible, losses = [], {}
    for candidate in widths:
        # Turn every rival down and this candidate up, then reverse the corner.
        favorable = {j: max(0, scores[j] + (w if j == candidate else -w))
                     for j, w in widths.items()}
        unfavorable = {j: max(0, scores[j] + (-w if j == candidate else w))
                       for j, w in widths.items()}
        if favorable[candidate] + 1e-10 >= max(favorable.values()):
            possible.append(candidate)
        losses[candidate] = (max(unfavorable.values()) - unfavorable[candidate]) * factor
    ranking = sorted(losses, key=losses.get)
    assert ranking[0] == expected["minimax_choice"] and losses[ranking[0]] < losses[ranking[1]]
    assert sorted(possible) == expected["possible_leaders"]
    assert abs(losses[ranking[0]] - expected["minimax_regret"]) < 1e-9
    assert abs(losses[leader] - expected["point_choice_regret"]) < 1e-9


if __name__ == "__main__":
    verify()
