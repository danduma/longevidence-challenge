"""Independently translate the archived CDSs and check the recorded comparison."""
import csv
from itertools import product
import json
from pathlib import Path
import re

root = Path(__file__).resolve().parent
sources = root / "sources"
result = json.loads((root / "results.json").read_text())
codons = dict(zip(map("".join, product("TCAG", repeat=3)),
                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def translate(name, marker):
    record = (sources / name).read_text()
    coordinates = re.search(r"\bCDS\s+(\d+)\.\.(\d+)", record)
    start, end = map(int, coordinates.groups())
    dna = "".join(re.findall(r"[acgt]+", record.split(marker, 1)[1]))
    cds = dna[start - 1:end].upper()
    protein = "".join(codons[cds[i:i + 3]] for i in range(0, len(cds), 3))
    assert protein.startswith("M") and protein.count("*") == 1 and protein.endswith("*")
    return protein[:-1], cds


def fasta(name):
    return "".join("".join((sources / name).read_text().splitlines()[1:]).split())


with (sources / "rr.csv").open() as handle:
    row, = csv.DictReader(handle)
rr = row["large_subunit"]
assert row["internal_id"] == "Rr"
assert row["identifier"] == result["reference"] == "WP_011390153.1"
assert len(rr) == result["length"] == 466
old, old_cds = translate("K01999.gb", "\n        1 ")
new, new_cds = translate("X00286.embl", "0 other;")
assert old == rr == fasta("WP_011390153.1.fasta")
assert new[:90] == rr[:90] and new[91:] == rr[91:] and (new[90], rr[90]) == ("H", "D")
assert [i for i, (a, b) in enumerate(zip(old_cds, new_cds), 1) if a != b] == [271]
fusion = fasta("frey.faa")
assert fusion[:466] == rr
assert fusion[466:] == result["removed_descendant_suffix"] == "TSRRRRRRRRRRLEHHHHHH"
uniprot = (sources / "P04718.txt").read_text().splitlines()
assert "".join(line.strip().replace(" ", "") for line in uniprot if line.startswith("     ") and line.strip().replace(" ", "").isalpha()) == new
assert "plasmid pRR2119" in (sources / "K01999.gb").read_text()
assert "replaced by X00286.1" in (sources / "K01999.gb").read_text()
assert "GenBank" in (sources / "atcc.html").read_text() and "K01999" in (sources / "atcc.html").read_text()
for name, seq in {"K01999.1": old, "X00286.1": new, "P04718": new,
                  "WP_011390153.1": rr, "Frey_R10_core": fusion[:466]}.items():
    differences = [[i + 1, seq[i], residue] for i, residue in enumerate(rr) if seq[i] != residue]
    assert result["comparisons"][name] == {
        "identical_residues": len(rr) - len(differences), "mismatches": differences}
print("Verified both CDS translations, the residue-91 discrepancy, and the released Rr and descendant matches.")
