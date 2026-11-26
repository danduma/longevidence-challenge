import unittest

import numpy as np
from timecourse import hazard_slope


class TimecourseTests(unittest.TestCase):
    def test_constant_fraction_has_zero_slope_after_initial_subtraction(self):
        self.assertAlmostEqual(hazard_slope(np.full(5, 20), np.full(5, 80)), 0.0)

    def test_linear_cumulative_hazard_recovers_sampling_index_slope(self):
        t = np.arange(5)
        total = 10**8
        unflipped = np.round(total * np.exp(-0.2 * t))
        flipped = total - unflipped
        self.assertAlmostEqual(hazard_slope(flipped, unflipped), 0.2, places=6)

    def test_invalid_count_or_curve_length_is_rejected(self):
        with self.assertRaises(ValueError):
            hazard_slope([-1, 1, 1, 1, 1], [100] * 5)
        with self.assertRaises(ValueError):
            hazard_slope([1, 2], [100, 100])


if __name__ == "__main__":
    unittest.main()
