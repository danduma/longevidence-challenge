"""Candidate-level sensitivity analysis using the frozen author supplement."""
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
from statistics import mean

from read_assays import observations

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"


def read(name):
    with (SOURCES / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main():
    expanded = {row["internal_id"] for row in read("fast_formII_clusters_expansion.csv")
                if row["expanded"] == "1"}
    rows = read("Supplementary_Data_3.csv")
    cohorts = {"representative": [], "expansion": [], "ancestors": []}
    for row in rows:
        key = row["internal_id"]
        cohort = "ancestors" if key.startswith("ANC_") else (
            "expansion" if key in expanded else "representative")
        cohorts[cohort].append(row)
    pools = {
        "representative": cohorts["representative"],
        "extant": cohorts["representative"] + cohorts["expansion"],
        "extant_form_II": [row for row in rows if row["form"] == "Form 2"
                          and not row["internal_id"].startswith("ANC_")],
        "all": rows,
    }
    output = {"cohorts": {name: counts(group) for name, group in cohorts.items()},
              "pools": {name: evaluate(group) for name, group in pools.items()},
              "assay_replay": replay(pools)}
    (ROOT / "results.json").write_text(json.dumps(output, indent=2) + "\n")


def replay(pools):
    data = defaultdict(list)
    for _, candidate, raw, normalized, assay in observations():
        data[candidate].append((raw, normalized, assay))
    result = {}
    for name, rows in pools.items():
        eligible = sorted(row["internal_id"] for row in rows if row["rate_lab_corrected"])
        normalized = {i: mean(x[1] for x in data[i]) for i in eligible}
        point = max(normalized, key=normalized.get)
        changed = {}
        blocks = sorted({x[2] for i in eligible for x in data[i]})
        for block in blocks:
            kept = {i: [x[1] for x in data[i] if x[2] != block] for i in eligible}
            new = max((i for i in eligible if kept[i]), key=lambda i: mean(kept[i]))
            if new != point:
                changed[block] = new
        result[name] = {
            "assay_labels": len(blocks),
            "raw_leader": max(eligible, key=lambda i: mean(x[0] for x in data[i])),
            "normalized_leader": point,
            "changed_after_omission": changed,
        }
    a, b = "RBC2_36", "RBC8_39"
    result["extant_pair"] = {
        "shared_assays": sorted({x[2] for x in data[a]} & {x[2] for x in data[b]}),
        "normalized_gap": round(mean(x[1] for x in data[a]) - mean(x[1] for x in data[b]), 9),
    }
    return result


def counts(rows):
    return {
        "total": len(rows),
        "numeric": sum(bool(row["rate_lab_corrected"]) for row in rows),
        "with_dispersion": sum(bool(row["rate_lab_corrected"] and row["stdev"]) for row in rows),
        "low_flag": sum(row["active<0.5s-1"] == "1.0" for row in rows),
        "unresolved": sum(not (row["rate_lab_corrected"] or row["active<0.5s-1"]) for row in rows),
    }


def evaluate(rows):
    measured = [row for row in rows if row["rate_lab_corrected"]]
    leader = max(measured, key=lambda row: float(row["rate_lab_corrected"]))["internal_id"]
    eligible = [row for row in measured if row["stdev"]]
    return {
        "counts": counts(rows),
        "point_leader": leader,
        "SD": stress(eligible, leader, False),
        "SD_over_sqrt_N": stress(eligible, leader, True),
    }


def stress(rows, leader, divide_by_sqrt_n):
    values = {}
    for row in rows:
        width = float(row["stdev"])
        if divide_by_sqrt_n:
            width /= math.sqrt(float(row["N"]))
        values[row["internal_id"]] = (float(row["rate_lab_corrected"]), width)
    mean, width = values[leader]
    tipping, rival = min(
        ((mean - other_mean) / (width + other_width), key)
        for key, (other_mean, other_width) in values.items()
        if key != leader and width + other_width > 0
    )
    lower = {key: max(0, m - w) for key, (m, w) in values.items()}
    upper = {key: m + w for key, (m, w) in values.items()}
    regret = {key: max(0, max(upper[j] for j in values if j != key) - lower[key])
              for key in values}
    best = min(regret, key=regret.get)
    return {
        "first_tie_width": round(tipping, 9),
        "first_tie_rival": rival,
        "possible_leaders": sorted(key for key in values if regret_possible(key, lower, upper)),
        "minimax_choice": best,
        "minimax_regret": round(regret[best], 9),
        "point_choice_regret": round(regret[leader], 9),
    }


def regret_possible(key, lower, upper):
    return upper[key] + 1e-10 >= max(lower[j] for j in lower if j != key)


if __name__ == "__main__":
    main()
