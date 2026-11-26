import unittest

from benchmark_audit import cross_fold_exposure, implied_join_count


class BenchmarkAuditTests(unittest.TestCase):
    def test_identical_group_in_another_fold_is_exposed(self):
        self.assertEqual(cross_fold_exposure(["a", "b", "a", "c"], [0, 0, 1, 1]), [0, 2])
        self.assertEqual(cross_fold_exposure(["a", "a", "b", "c"], [0, 0, 1, 1]), [])

    def test_missing_identity_is_not_counted_as_a_shared_group(self):
        self.assertEqual(cross_fold_exposure([None, None, "a"], [0, 1, 0]), [])

    def test_alias_metadata_duplicates_reweight_a_species_but_do_not_add_measurements(self):
        rows = [dict(species="A", primary="1", mutant="0", vC="2"),
                dict(species="A ", primary="1", mutant="0", vC="4"),
                dict(species="B", primary="1", mutant="0", vC=""),
                dict(species="C", primary="1", mutant="1", vC="9")]
        metadata = [dict(label="A", species="A"), dict(label="A", species="A")]
        labels = ["A", "A"]
        self.assertEqual(implied_join_count(rows, metadata, labels, "vC"), 2)
        self.assertEqual(implied_join_count(rows, metadata[:1], labels[:1], "vC"), 1)


if __name__ == "__main__":
    unittest.main()
