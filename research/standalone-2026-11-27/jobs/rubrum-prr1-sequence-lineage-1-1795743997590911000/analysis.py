"""Compare the documented parent and descendant with the released assay row."""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"


def load(name):
    return (SOURCES / name).read_text()


def translation(record):
    return "".join(re.search(r'/translation="([^"]+)"', record)[1]
                   .replace("\nFT", "\n").split())


def fasta(name):
    return "".join("".join(load(name).splitlines()[1:]).split())


def run():
    with (SOURCES / "rr.csv").open(newline="") as handle:
        row, = csv.DictReader(handle)
    rr = row["large_subunit"]
    fusion = fasta("frey.faa")
    suffix = "TS" + "R" * 10 + "LE" + "H" * 6
    sequences = {
        "K01999.1": translation(load("K01999.gb")),
        "X00286.1": translation(load("X00286.embl")),
        "P04718": "".join(load("P04718.txt").split("SQ   SEQUENCE")[1]
                          .splitlines()[1:-1]).replace(" ", ""),
        "WP_011390153.1": fasta("WP_011390153.1.fasta"),
        "Frey_R10_core": fusion.removesuffix(suffix),
    }
    comparisons = {}
    for name, seq in sequences.items():
        if len(seq) != len(rr):
            raise ValueError(f"Unexpected core length in {name}")
        comparisons[name] = {
            "identical_residues": sum(a == b for a, b in zip(seq, rr)),
            "mismatches": [[i, a, b] for i, (a, b) in
                           enumerate(zip(seq, rr), 1) if a != b],
        }
    return {"reference": row["identifier"], "length": len(rr),
            "removed_descendant_suffix": suffix, "comparisons": comparisons}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
