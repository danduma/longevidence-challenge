"""Independently verify PERRC coverage and RFU slope ratios from the source workbook."""
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path
from xml.dom import minidom
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"


def get_text(node):
    return "".join(n.firstChild.data for n in node.getElementsByTagName("t") if n.firstChild)


def node_text(node):
    if node.nodeType == node.TEXT_NODE:
        return node.data
    return "".join(node_text(child) for child in node.childNodes)


def rows_from_workbook(path):
    with ZipFile(path) as book:
        strings = [
            get_text(node)
            for node in minidom.parseString(book.read("xl/sharedStrings.xml")).getElementsByTagName("si")
        ]
        sheet = minidom.parseString(book.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in sheet.getElementsByTagName("row"):
        parts = {}
        for cell in row.getElementsByTagName("c"):
            values = cell.getElementsByTagName("v")
            if not values or not values[0].firstChild:
                continue
            value = values[0].firstChild.data
            col = cell.getAttribute("r").rstrip("0123456789")
            parts[col] = strings[int(value)] if cell.getAttribute("t") == "s" else value
        rows.append(parts)
    return rows


def covariance_slope(trace):
    count = len(trace)
    sx = sum(x for x, _ in trace)
    sy = sum(y for _, y in trace)
    sxx = sum(x * x for x, _ in trace)
    sxy = sum(x * y for x, y in trace)
    return (count * sxy - sx * sy) / (count * sxx - sx * sx)


def main():
    manifest = json.loads((SOURCES / "figshare-29104880.json").read_text())["files"][0]
    path = SOURCES / manifest["name"]
    assert hashlib.md5(path.read_bytes()).hexdigest() == manifest["computed_md5"]
    rows = rows_from_workbook(path)
    published = minidom.parse(str(SOURCES / "perrc.xml"))
    table = published.getElementsByTagName("tbody")[0]
    pairs = set()
    enzyme = None
    for row in table.getElementsByTagName("tr"):
        values = [node_text(cell).strip() for cell in row.getElementsByTagName("td")[:3]]
        if len(row.getElementsByTagName("td")) == 7:
            enzyme = values[0]
            motif = values[2]
        else:
            motif = values[0]
        pairs.add((enzyme, motif))

    parsed = {}
    time_grids = set()
    column_count = 0
    saved_count = 0
    for index, header in enumerate(rows):
        if header.get("B") != "Time (sec)":
            continue
        enzyme, motif = rows[index - 2]["B"].split(" on ")
        columns = [col for col in header if col != "B"]
        assert set(Counter(header[col] for col in columns).values()) == {2}
        trace = []
        for row in rows[index + 1:]:
            if not row.get("B", "").isdigit():
                break
            trace.append(row)
        time_grids.add(tuple(int(row["B"]) for row in trace))
        column_count += len(columns)
        archived = next(row for row in rows[index + 1:] if row.get("B") == "Slope")
        saved_count += sum(col in archived for col in columns)
        slopes = {}
        for col in columns:
            value = covariance_slope([(Fraction(row["B"]), Fraction(row[col])) for row in trace])
            slopes.setdefault(header[col], []).append(value)
        assert (enzyme, motif) not in parsed
        parsed[enzyme, motif] = {dose: sum(values) / len(values) for dose, values in slopes.items()}

    assert set(parsed) == pairs
    canonical = "ENLYFQS"
    assert all(len(motif) == len(canonical) for _, motif in pairs)
    noncanonical = {pair for pair in pairs if pair[1] != canonical}
    positions = {
        distance
        for _, motif in noncanonical
        for distance, (a, b) in enumerate(zip(reversed(canonical), reversed(motif)), start=0)
        if a != b
    }
    # P1 is one residue before the terminal P1-prime residue.
    assert positions == {1}
    folds = sum(bool(positions - {position}) for position in positions)
    assert folds == 0

    result = json.loads((ROOT / "results.json").read_text())
    assert result["kinetic_pair_count"] == len(pairs)
    assert result["noncanonical_pair_count"] == len(noncanonical)
    assert result["substitution_positions"] == sorted(positions)
    assert result["maximum_nonempty_position_folds"] == len(positions)
    assert result["trainable_held_position_folds"] == folds
    assert result["trace_time_seconds"] == [list(grid) for grid in sorted(time_grids)]
    assert result["trace_columns"] == column_count
    assert result["archived_slope_count"] == saved_count
    reported = result["alternate_over_canonical_rfu_slope"]
    assert set(reported) == {enzyme for enzyme, _ in noncanonical}
    for enzyme, motif in noncanonical:
        assert set(reported[enzyme]) == {"25", "50"}
        for dose in reported[enzyme]:
            ratio = parsed[enzyme, motif][dose] / parsed[enzyme, canonical][dose]
            assert math.isclose(reported[enzyme][dose], float(ratio), rel_tol=0, abs_tol=5.1e-9)
    print("Verified workbook–Table 1 coverage, slope counts and raw-trace ratios.")


if __name__ == "__main__":
    main()
