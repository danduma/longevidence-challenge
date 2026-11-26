#!/usr/bin/env python3
"""Apply a common, inspectable decision rule to three recorded scout judgments."""
import csv
import json
from pathlib import Path
from statistics import median


def main():
    root = Path(__file__).resolve().parent
    config = json.loads((root / "ranking_inputs.json").read_text())
    data = {int(r["challenge_id"]): r for r in csv.DictReader((root / "agents/scout-data-sol/results/challenge_ranking.csv").open())}
    mechanism = {int(r["challenge_id"]): r for r in csv.DictReader((root / "agents/scout-mechanisms-sol/outputs/challenge_ranking.csv").open())}
    assert set(data) == set(mechanism) == set(range(1, 13))
    rows = []
    for cid in range(1, 13):
        d, m = data[cid], mechanism[cid]
        name = d["challenge_name"]
        assert name == m["exact_title"]
        first = [float(d[k]) for k in ["data_availability", "computational_tractability", "impact", "validation_path"]]
        second = [float(m[k]) for k in config["dimensions"]]
        third = config["astra_scores_by_official_id"][str(cid)]
        scores = [median(v) for v in zip(first, second, third, strict=True)]
        assert all(1 <= x <= 5 for x in scores)
        row = dict(challenge_id=cid, challenge=name, **dict(zip(config["dimensions"], scores, strict=True)))
        for scenario, weights in config["weight_scenarios"].items():
            assert abs(sum(weights) - 1) < 1e-12
            row[scenario + "_score"] = round(sum(a * b for a, b in zip(scores, weights, strict=True)), 6)
        rows.append(row)
    scenarios = {}
    for scenario in config["weight_scenarios"]:
        ranked = sorted(rows, key=lambda r: (-r[scenario + "_score"], -r["validation"], -r["tractability"], r["challenge_id"]))
        for row in ranked:
            row[scenario + "_rank"] = 1 + sum(other[scenario + "_score"] > row[scenario + "_score"] for other in rows)
        cutoff = ranked[2][scenario + "_score"]
        scenarios[scenario] = {"top_three_cutoff_score": cutoff,
                               "at_or_above_third_place": [r["challenge_id"] for r in ranked if r[scenario + "_score"] >= cutoff],
                               "tie_spanning_third_place": sum(r[scenario + "_score"] >= cutoff for r in ranked) > 3}
    ranked = sorted(rows, key=lambda r: (r["primary_rank"], -r["validation"], -r["tractability"], r["challenge_id"]))
    with (root / "ranking.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(ranked[0]))
        writer.writeheader()
        writer.writerows(ranked)
    (root / "ranking_sensitivity.json").write_text(json.dumps(scenarios, indent=2) + "\n")
    for row in ranked:
        print(row["primary_rank"], row["challenge_id"], row["challenge"], row["primary_score"])
    print(json.dumps(scenarios, indent=2))


if __name__ == "__main__":
    main()
