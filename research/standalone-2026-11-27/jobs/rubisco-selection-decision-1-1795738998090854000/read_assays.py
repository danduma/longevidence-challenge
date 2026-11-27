"""Read published S02 assay observations using only the standard library."""
from pathlib import Path
import xml.etree.ElementTree as ET
from zipfile import ZipFile

SOURCE = Path(__file__).resolve().parent / "sources" / "pnas.2501433122.sd02.xlsx"
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def observations():
    with ZipFile(SOURCE) as book:
        strings = ["".join(item.itertext())
                   for item in ET.fromstring(book.read("xl/sharedStrings.xml"))]
        sheet = ET.fromstring(book.read("xl/worksheets/sheet1.xml"))
    for row in sheet.findall(".//m:row", NS):
        if int(row.get("r")) <= 3:
            continue
        cells = {}
        for cell in row.findall("m:c", NS):
            value = cell.findtext("m:v", namespaces=NS)
            if value is not None:
                cells[cell.get("r")[0]] = strings[int(value)] if cell.get("t") == "s" else value
        if set("ABCDE") <= cells.keys():
            yield cells["A"], cells["B"], float(cells["C"]), float(cells["D"]), cells["E"]
