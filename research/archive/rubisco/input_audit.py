#!/usr/bin/env python3
"""Audit the pinned 2026 Rubisco release without silently dropping alias records.

Compatible metadata is not proof of paired enzyme assays. The conservative
pair table is for prioritization only, never a declaration of challenge success.
"""
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
from statistics import median


def canonical_species(value):
    return " ".join(value.split())


def fasta_records(text):
    records = []
    seen = {}
    for block in text.split(">")[1:]:
        lines = block.splitlines()
        label, sequence = lines[0].strip(), "".join(lines[1:]).strip().upper()
        if label in seen and seen[label] != sequence:
            raise ValueError(f"Conflicting sequences for FASTA header {label}")
        seen[label] = sequence
        records.append((label, sequence))
    return records


def positive(value):
    try:
        x = float(value)
        return x if 0 < x < float("inf") else None
    except (ValueError, TypeError):
        return None


def paired_measurements(rows):
    groups = defaultdict(list)
    direct = []
    for index, row in enumerate(rows):
        if row.get("primary") != "1" or row.get("mutant") != "0":
            continue
        speed, specificity = positive(row.get("vC")), positive(row.get("S"))
        identity = {k: row.get(k, "") for k in ["species", "temp_C", "pH", "isoform", "variant", "heterologous_expression", "pmid_or_doi"]}
        identity["species"] = canonical_species(identity["species"])
        if not identity["pmid_or_doi"] or not identity["temp_C"] or not identity["pH"]:
            if speed and specificity:
                direct.append(dict(identity, vC=speed, S=specificity, row_indices=[index], provenance="single_source_row"))
            continue
        groups[tuple(identity.items())].append((index, row))
    output = direct
    for key, observations in groups.items():
        speeds = [positive(r.get("vC")) for _, r in observations]
        specificities = [positive(r.get("S")) for _, r in observations]
        speeds, specificities = [x for x in speeds if x], [x for x in specificities if x]
        if speeds and specificities:
            output.append(dict(dict(key), vC=median(speeds), S=median(specificities),
                               row_indices=[i for i, _ in observations],
                               provenance="single_source_row" if len(observations) == 1 else "metadata_compatible_rows_not_verified_paired_assay"))
    return output


def main():
    root = Path(__file__).resolve().parents[1]
    source = root / "sources/rubisco-active-learning/supplementary/doubling_time_package_v2"
    rows = list(csv.DictReader((source / "inputs/flamholz_dataset_S1.csv").open()))
    meta = list(csv.DictReader((source / "manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv").open()))
    records = fasta_records((source / "manuscript_assets/flamholz/flamholz_rubisco_wt.fasta").read_text())
    wt = [r for r in rows if r["primary"] == "1" and r["mutant"] == "0"]
    labels = Counter(label for label, _ in records)
    original_species = {r["species"] for r in wt}
    canonical_names = {canonical_species(s) for s in original_species}
    meta_species = {r["species"] for r in meta}
    alias_groups = defaultdict(list)
    for s in sorted(original_species):
        alias_groups[canonical_species(s)].append(s)
    pairs = paired_measurements(rows)
    exact_sequences = defaultdict(list)
    unique_records = dict(records)
    for label, sequence in unique_records.items():
        exact_sequences[sequence].append(label)
    target_overlap = {}
    for target in ["vC", "S", "KC", "KO", "vO", "KRuBP"]:
        exact = {r["species"] for r in wt if positive(r[target])}
        cleaned = {canonical_species(r["species"]) for r in wt if positive(r[target])}
        target_overlap[target] = {
            "exact_species_with_labels": len(exact & meta_species),
            "canonical_species_with_labels": len(cleaned & meta_species),
            "source_exact_species": len(exact), "source_canonical_species": len(cleaned)}
    out = root / "rubisco/results"
    out.mkdir(exist_ok=True)
    summary = {
        "source_rows": len(rows), "primary_wild_type_rows": len(wt),
        "original_species_keys": len(original_species), "whitespace_normalized_species_keys": len(canonical_names),
        "metadata_rows": len(meta), "metadata_unique_species": len(meta_species),
        "fasta_records": len(records), "fasta_unique_headers": len(unique_records),
        "fasta_duplicate_headers": {k: v for k, v in labels.items() if v > 1},
        "unique_sequences": len(exact_sequences),
        "distinct_species_with_identical_sequence_groups": [v for v in exact_sequences.values() if len(v) > 1],
        "whitespace_aliases": {k: v for k, v in alias_groups.items() if len(v) > 1},
        "target_coverage": target_overlap,
        "compatible_speed_specificity_pairs": len(pairs),
        "compatible_pairs_at_25C": sum(positive(r["temp_C"]) == 25 for r in pairs),
        "caveats": ["Sequence identity is an attribute of this downloaded release; assay-sample sequence identity is unverified.",
                    "RbcL-only equality does not imply equal holoenzymes: small subunits and assay context can differ.",
                    "Compatible metadata is not proof of the paired control assays required by the challenge."]}
    (out / "input_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "compatible_pairs.json").write_text(json.dumps(pairs, indent=2) + "\n")
    manifest = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [source / "inputs/flamholz_dataset_S1.csv", source / "manuscript_assets/flamholz/flamholz_rubisco_wt_meta.csv", source / "manuscript_assets/flamholz/flamholz_rubisco_wt.fasta"]}
    (out / "input_hashes.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ["whitespace_aliases", "distinct_species_with_identical_sequence_groups"]}, indent=2))


if __name__ == "__main__":
    main()
