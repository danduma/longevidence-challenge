"""Audit the positional coverage and early slopes in PERRC supplement 2."""
import json
import re
from pathlib import Path
from statistics import linear_regression, mean
from xml.etree import ElementTree as ET
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent


def read_cells():
    with ZipFile(ROOT / "sources/sb5c00154_si_002.xlsx") as book:
        strings = [
            "".join(node.itertext())
            for node in ET.fromstring(book.read("xl/sharedStrings.xml"))
        ]
        sheet = ET.fromstring(book.read("xl/worksheets/sheet1.xml"))
    return {
        cell.get("r"): strings[int(value.text)] if cell.get("t") == "s" else value.text
        for cell in sheet.findall(".//{*}c")
        if (value := cell.find("{*}v")) is not None
    }


def main():
    cells = read_cells()
    blocks = {}
    for address, label in cells.items():
        match = re.fullmatch(r"B(\d+)", address)
        if not match or not re.fullmatch(r"\w+ on [A-Z]{7}", label):
            continue
        top = int(match[1])
        header = top + 2
        first = header + 1
        stop = first
        while cells.get(f"B{stop}", "").isdigit():
            stop += 1
        times = [int(cells[f"B{row}"]) for row in range(first, stop)]
        slopes = {}
        for col in "CDEFGHIJKLMNOPQRST":
            concentration = cells[f"{col}{header}"]
            values = [float(cells[f"{col}{row}"]) for row in range(first, stop)]
            slopes.setdefault(concentration, []).append(linear_regression(times, values).slope)
        blocks[label] = times, slopes

    # The canonical label is identified in the article, and index 0 is P6.
    canonical = "ENLYFQS"
    positions = sorted({
        6 - index
        for label in blocks
        for index, (base, other) in enumerate(zip(canonical, label.split()[-1]))
        if base != other
    })
    ratios = {}
    for enzyme in sorted({label.split()[0] for label in blocks}):
        native = blocks[f"{enzyme} on {canonical}"][1]
        alternate = blocks[f"{enzyme} on ENLYFES"][1]
        ratios[enzyme] = {
            concentration: round(mean(alternate[concentration]) / mean(native[concentration]), 8)
            for concentration in ("50", "25")
        }
    results = {
        "kinetic_pair_count": len(blocks),
        "noncanonical_pair_count": sum(not label.endswith(canonical) for label in blocks),
        "substitution_positions": positions,
        "maximum_nonempty_position_folds": len(positions),
        "trainable_held_position_folds": len(positions) if len(positions) > 1 else 0,
        "trace_time_seconds": sorted({tuple(t) for t, _ in blocks.values()}),
        "trace_columns": sum(len(v) for _, slopes in blocks.values() for v in slopes.values()),
        "alternate_over_canonical_rfu_slope": ratios,
        "archived_slope_count": sum(
            f"{col}{row}" in cells for row in (20, 42, 64, 86) for col in "CDEFGHIJKLMNOPQRST"
        ),
    }
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
