"""Reproduce the three evidence audits. Python standard library, run from this directory."""
import csv
import io
import itertools
import json
import math
import statistics as st
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

DATA = Path("data")


def read_csv(name):
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def rubisco():
    rows = [r for r in read_csv("rubisco.csv") if
            (r["primary"], r["mutant"], r["temp_C"], r["isoform"]) == ("1", "0", "25", "1")]
    paired = [r for r in rows if r["S"] and r["vC"]]
    red = [r for r in rows if r["species"].startswith("Galdieria partita") and r["S"]]
    maize = [r for r in rows if r["species"].startswith("Zea mays") and r["vC"]]
    red_min = min(float(r["S"]) for r in red)
    maize_min = min(float(r["vC"]) for r in maize)
    studies = lambda rr: {r["pmid_or_doi"] for r in rr}
    # Species spellings retained, trailing whitespace removed, cultivar labels preserved.
    groups = defaultdict(list)
    for r in rows:
        groups[r["species"].strip()].append(r)
    joint, same_ref = 0, 0
    for rr in groups.values():
        a, b = ([r for r in rr if r[col]] for col in ("S", "vC"))
        joint += bool(a and b)
        same_ref += bool(studies(a) & studies(b))
    fastest_specific = max((float(r["S"]),r["species"]) for r in paired if float(r["vC"]) >= maize_min)
    # Restore within-paper joins split over two assay rows. Use per-label/paper medians.
    stratified = defaultdict(list)
    for r in rows:
        stratified[r["species"].strip(), r["short_ref"]].append(r)
    joint_estimates = []
    for (species, source), rr in stratified.items():
        specs, rates = ([float(r[c]) for r in rr if r[c]] for c in ("S", "vC"))
        if specs and rates:
            joint_estimates.append((st.median(specs), st.median(rates), species, source))
    passing_rate = [x for x in joint_estimates if x[1] >= maize_min]
    return dict(study_label_pairs=len(joint_estimates), study_label_joint_hits=sum(x[0]>=red_min for x in passing_rate),
                study_label_largest_S_above_maize_floor=max(passing_rate),
                runner_up_study=max(x for x in passing_rate if x[3] != max(passing_rate)[3]),
                rows=len(rows), same_row=len(paired), labels=len(groups), joint_any_source=joint,
                joint_same_citation=same_ref, red=[(r["short_ref"],r["S"],r["vC"]) for r in red],
                maize=[(r["short_ref"],r["vC"]) for r in maize],
                control_source_overlap=sorted(studies(red)&studies(maize)),
                largest_S_above_maize_floor=fastest_specific,
                same_row_hits=sum(float(r["S"]) >= red_min and float(r["vC"]) >= maize_min for r in paired))


def sheet_cells(archive, idx):
    xml = ET.fromstring(archive.read(f"xl/worksheets/sheet{idx}.xml"))
    return {c.get("r"): c.findtext("{*}v") for c in xml.findall(".//{*}c")}


def cryo():
    with zipfile.ZipFile(next(DATA.glob("41467*.xlsx"))) as archive:
        cool, warm = [sheet_cells(archive, i) for i in (1, 4)]
    def cols(cells, row, names):
        return [float(cells[f"{c}{row}"]) for c in names]
    cooling = []
    for row in (5, 6, 7):
        measured = cols(cool, row, "STU")
        model = cols(cool, row, "V")[0]
        cooling.append(dict(mean=st.mean(measured), sd=st.stdev(measured), model=model,
                            model_relative_error=model/st.mean(measured)-1))
    warming = []
    for row, iron in [(6, 10.7), (7, 4.6)]:
        # Matched probe readings in each experiment, not nine independent bags.
        bags = [st.mean(cols(warm, row, names)) for names in ("MPS", "NQT", "ORU")]
        warming.append(dict(mean=st.mean(bags), sd=st.stdev(bags), per_iron=st.mean(bags)/iron))
    return dict(cooling=cooling, warming=warming, raw_ratio=warming[1]["mean"]/warming[0]["mean"],
                normalized_ratio=warming[1]["per_iron"]/warming[0]["per_iron"],
                all_successes_needed=math.floor(math.log(.05)/math.log(.99))+1)


def quadruplet():
    rows = read_csv("elife-76941-fig6-data2.csv")
    codons = [c for c in rows[0] if c]
    assert all(r[""].endswith(c) for r,c in zip(rows,codons))
    edges = Counter()
    for a,b in itertools.combinations(codons,2):
        diff = [i+1 for i in range(4) if a[i] != b[i]]
        if len(diff)==1:
            edges[diff[0]] += 1
    maxima = [max(codons,key=lambda c:float(r[c])) for r in rows]
    # All 4^4 codons, 3 substitutions at each position. Divide by 2 for undirected edges.
    total_per_position=4**4*3//2
    # Independent report with fourth-base and triplet controls (2021 Fig. 4a).
    with zipfile.ZipFile(DATA / "quad2021.zip") as archive:
        with zipfile.ZipFile(io.BytesIO(archive.read("DeBenedictis2021_SourceDataFile.xlsx"))) as workbook:
            cells = sheet_cells(workbook, 14)
    ratios = []
    for row in range(5, 10):
        def row_mean(columns):
            return st.mean(float(cells[f"{c}{row}"]) for c in columns if cells.get(f"{c}{row}"))
        ratios.append(row_mean("KLMN")/row_mean("CDEF"))
    return dict(fourth_base_mismatch_over_cognate_range=[min(ratios),max(ratios)], correct_row_max=sum(a==b for a,b in zip(codons,maxima)), rows=len(rows),
                undirected_edges=[edges[i] for i in range(1,5)], possible_per_position=total_per_position,
                unique_first_three=len({c[:3] for c in codons}),
                diagonal_range=[f(float(r[c]) for r,c in zip(rows,codons)) for f in (min,max)],
                balanced_random_change=1-(16*13*12+4*12*11)/(256*255))


if __name__ == "__main__":
    print(json.dumps(dict(rubisco=rubisco(),cryo=cryo(),quadruplet=quadruplet()), indent=2))
