import unittest

from input_audit import canonical_species, fasta_records, paired_measurements


class InputAuditTests(unittest.TestCase):
    def test_whitespace_aliases_join_but_cultivars_remain_distinct(self):
        self.assertEqual(canonical_species(" Zea  mays \n"), "Zea mays")
        self.assertNotEqual(canonical_species("Zea mays"), canonical_species("Zea mays cv. Carella"))

    def test_duplicate_fasta_headers_are_retained_and_conflicts_rejected(self):
        records = fasta_records(">A\nAA\nAA\n>A\nAAAA\n>B\nCCCC\n")
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0][1], "AAAA")
        with self.assertRaises(ValueError):
            fasta_records(">A\nAAAA\n>A\nAAAC\n")

    def test_pairs_may_join_split_rows_only_within_same_source_and_conditions(self):
        base = dict(species="Zea mays", primary="1", mutant="0", temp_C="25", pH="8", isoform="1", variant="", heterologous_expression="", pmid_or_doi="paper-a")
        rows = [dict(base, vC="4", S=""), dict(base, vC="", S="80"),
                dict(base, pmid_or_doi="paper-b", vC="", S="90"),
                dict(base, temp_C="30", vC="6", S=""),
                dict(base, mutant="1", vC="99", S="999")]
        pairs = paired_measurements(rows)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0]["vC"], 4)
        self.assertEqual(pairs[0]["S"], 80)

    def test_unreported_condition_does_not_authorize_joining_different_assays(self):
        rows = [dict(species="A", primary="1", mutant="0", temp_C="25", pH="", vC="4", S="", pmid_or_doi="paper"),
                dict(species="A", primary="1", mutant="0", temp_C="25", pH="", vC="", S="80", pmid_or_doi="paper")]
        self.assertEqual(paired_measurements(rows), [])


if __name__ == "__main__":
    unittest.main()
