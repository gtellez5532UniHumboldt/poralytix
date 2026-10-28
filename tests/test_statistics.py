"""
Pruebas unitarias de métodos bioestadísticos: Wilson CI, Cohen's h, Cochran-Armitage y Tablas 2x2.
"""
import unittest
import numpy as np
from cytoworkstation.core.statistics import (
    wilson_binomial_ci,
    cohen_h_proportions,
    contingency_2x2_p_and_or,
    cochran_armitage_trend_test,
    classify_affectation_status
)


class TestStatistics(unittest.TestCase):
    def test_wilson_ci(self):
        lo, hi = wilson_binomial_ci(50, 100)
        self.assertAlmostEqual((lo + hi) / 2.0, 50.0, delta=1.0)
        self.assertTrue(0.0 <= lo < 50.0 < hi <= 100.0)

        # Caso extremo 0%
        lo0, hi0 = wilson_binomial_ci(0, 100)
        self.assertEqual(lo0, 0.0)
        self.assertGreater(hi0, 0.0)

    def test_cohen_h(self):
        # Mismas proporciones -> h = 0
        h_zero = cohen_h_proportions(50.0, 50.0)
        self.assertAlmostEqual(h_zero, 0.0, places=5)

        # Gran diferencia (10% vs 90%) -> h grande (> 0.8)
        h_large = cohen_h_proportions(90.0, 10.0)
        self.assertGreater(h_large, 0.8)

    def test_contingency_2x2(self):
        # Afectados altos en EXP vs bajos en CTRL
        chi2, p_val, or_val, or_lo, or_hi, rr_val = contingency_2x2_p_and_or(
            aff_exp=80, unaff_exp=20,
            aff_ref=5, unaff_ref=95
        )
        self.assertLess(p_val, 0.001)
        self.assertGreater(or_val, 10.0)
        self.assertGreater(rr_val, 5.0)

    def test_cochran_armitage(self):
        doses = [0.5, 1.0, 2.0, 4.0]
        aff_counts = [5, 15, 60, 95]
        tot_counts = [100, 100, 100, 100]
        z_stat, p_val, slope = cochran_armitage_trend_test(doses, aff_counts, tot_counts)
        self.assertGreater(z_stat, 0.0)
        self.assertLess(p_val, 0.001)
        self.assertGreater(slope, 0.0)

    def test_classify_status(self):
        self.assertIn("Severa", classify_affectation_status(0.0001, 0.95, 85.0))
        self.assertIn("No Afectado", classify_affectation_status(0.5, 0.05, 2.0))


if __name__ == "__main__":
    unittest.main()
