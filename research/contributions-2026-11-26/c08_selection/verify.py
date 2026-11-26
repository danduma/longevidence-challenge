#!/usr/bin/env python3
"""C08 finite-library identification analysis. Default mode checks frozen results."""
import argparse
import csv
import hashlib
import io
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
PRIOR = ROOT.parents[2] / "millennium/agents/rubisco-current-sol"
FLAG = "active<0.5s-1"
RATE = "rate_lab_corrected"
STATES = ("numeric_corrected_rate", "reported_below_0.5", "unreported_or_unsuccessful")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def records(inputs):
    data = read_csv(inputs["supp"])
    master = {r["internal_id"]: r for r in read_csv(inputs["library"]) if r["from_lab"] == "1.0"}
    assert len(master) == len(data) == len({r["internal_id"] for r in data}) == 469
    assert len({r["identifier"] for r in data}) == len(data)
    ledger = read_csv(inputs["ledger"]) if "ledger" in inputs else None
    assert ledger is None or len(ledger) == len(data)
    cells = json.loads(inputs["notebook"].read_text())["cells"]
    code = "\n".join("".join(c["source"]) for c in cells if c["cell_type"] == "code")
    assert "library_df[library_df['from_lab'] == 1]" in code
    assert "library_df['rate_lab'] * 10**b" in code
    factors = []
    for i, row in enumerate(data):
        original = master[row["internal_id"]]
        for k, v in row.items():
            if k == "new_to_this_study":
                continue
            field = {"stdev": "std_rate_lab_corrected", "N": "N_rate_lab"}.get(k, k)
            assert original[field] == v, (i, k)
        assert row[FLAG] in ("", "1.0") and not (row[RATE] and row[FLAG])
        state = STATES[0] if row[RATE] else STATES[1] if row[FLAG] else STATES[2]
        if ledger is not None:
            assert all(ledger[i][k] == v for k, v in row.items())
            assert ledger[i]["measurement_state"] == state and int(ledger[i]["source_row"]) == i + 2
        if row[RATE]:
            factors.append(float(row[RATE]) / float(original["rate_lab"]))
        row["state"] = state
    assert math.isclose(min(factors), max(factors), rel_tol=1e-12)
    return data, max(factors)


def summarize(rows, factor, stratum, group):
    rates = [float(r[RATE]) for r in rows if r[RATE]]
    c = sum(r["state"] == STATES[1] for r in rows)
    n, m = len(rows), len(rates)
    u = n - m - c
    ascertained = m + c
    med = lambda values: median(values) if values else ""
    for t in (10, 15):
        h = sum(v >= t for v in rates)
        yield dict(stratum=stratum, group=group, threshold=t, eligible=n, numeric=m, censored=c,
                   unknown=u, high=h, numeric_coverage=m/n, threshold_coverage=ascertained/n,
                   numeric_share=h/m if m else "", ascertained_share=h/ascertained if ascertained else "",
                   confirmed_yield=h/n, full_share_low=h/n, full_share_high=(h+u)/n,
                   numeric_median=med(rates), ascertained_median_low=med(rates+[0]*c),
                   ascertained_median_sup=med(rates+[0.5]*c),
                   full_median_low=med(rates+[0]*(c+u)),
                   full_median_sup=med(rates+[0.5]*c+[math.inf]*u),
                   raw_flag_sensitivity_ascertained_sup=med(rates+[0.5*factor]*c),
                   raw_flag_sensitivity_full_sup=med(rates+[0.5*factor]*c+[math.inf]*u))


def comparisons(summary):
    for t in (10, 15):
        forms = [s for s in summary if s["stratum"] == "form" and s["threshold"] == t]
        for a, b in itertools.combinations(forms, 2):
            def possibilities(s):
                return {Fraction(h, s["eligible"]) for h in range(s["high"], s["high"]+s["unknown"]+1)}
            left, right = possibilities(a), possibilities(b)
            yield dict(threshold=t, form_a=a["group"], form_b=b["group"],
                       difference_low=float(min(left)-max(right)), difference_high=float(max(left)-min(right)),
                       a_greater_possible=int(max(left)>min(right)), tie_possible=int(bool(left & right)),
                       b_greater_possible=int(max(right)>min(left)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--public-inputs", type=Path, help="directory with the three named upstream files")
    p.add_argument("--write", action="store_true", help="regenerate result CSVs")
    args = p.parse_args()
    manifest = json.loads((ROOT/"sources.json").read_text())
    inputs = {}
    for name, spec in manifest["inputs"].items():
        if args.public_inputs and not spec.get("upstream_path"):
            continue
        path = (args.public_inputs / Path(spec["upstream_path"]).name if args.public_inputs
                else PRIOR / spec["local_path"])
        assert digest(path) == spec["sha256"], name
        inputs[name] = path
    data, factor = records(inputs)
    summary = list(summarize(data, factor, "all", "all"))
    for k in ("form", "prok_euk", "new_to_this_study"):
        for g in sorted({r[k] for r in data}):
            summary.extend(summarize([r for r in data if r[k] == g], factor, k, g))
    for form, new in sorted({(r["form"], r["new_to_this_study"]) for r in data}):
        subset = [r for r in data if (r["form"], r["new_to_this_study"]) == (form, new)]
        summary.extend(summarize(subset, factor, "form_by_new", f"{form} | {new}"))
    states = [dict(supplement_csv_row=i, form=r["form"], prok_euk=r["prok_euk"],
                   new_to_this_study=r["new_to_this_study"], state=r["state"])
              for i, r in enumerate(data, 2)]
    for name, table in [("states", states), ("denominators", summary), ("form_comparisons", list(comparisons(summary)))]:
        stream = io.StringIO(newline="")
        w = csv.DictWriter(stream, fieldnames=table[0], lineterminator="\n")
        w.writeheader()
        w.writerows(table)
        path = ROOT / (name + ".csv")
        output = stream.getvalue().encode()
        if args.write:
            path.write_bytes(output)
        else:
            assert path.read_bytes() == output, name
    print("C08 verified 469 source records, 32 denominator rows and 20 form contrasts.")


if __name__ == "__main__":
    main()
