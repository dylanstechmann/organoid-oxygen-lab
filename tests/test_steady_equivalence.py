import tempfile
import unittest
from pathlib import Path

import numpy as np
from oxygenlab.cli import equivalence_demo
from oxygenlab.model import Parameters, zero_order_analytic
from oxygenlab.steady_equivalence import equivalent_parameters, steady_equivalence


class EquivalenceTests(unittest.TestCase):
    def test_mm_profiles_match_while_uptake_and_diffusion_scale_change(self):
        report = steady_equivalence(Parameters(shells=40), 3)
        self.assertLess(report["max_profile_difference_mol_m3"], 1e-9)
        self.assertAlmostEqual(report["uptake_ratio"], 3)
        outputs = report["rate_outputs"]
        self.assertAlmostEqual(outputs["scaled"]["diffusion_time_scale_s"] * 3, outputs["baseline"]["diffusion_time_scale_s"])
        self.assertLess(outputs["scaled"]["relative_mass_balance_error"], 1e-7)
        self.assertFalse(report["model_fitted"])

    def test_independent_zero_order_formula_retains_same_profile(self):
        p = Parameters(radius_um=100, kinetics="zero_order", shells=40)
        scaled = equivalent_parameters(p, 2.5)
        radii = np.linspace(0, p.radius_um, 12)
        np.testing.assert_allclose(zero_order_analytic(radii, p), zero_order_analytic(radii, scaled), atol=1e-14)
        self.assertLess(steady_equivalence(p, 2.5)["max_profile_difference_mol_m3"], 1e-9)

    def test_fixed_surface_and_zero_uptake_keep_explicit_degenerate_cases(self):
        report = steady_equivalence(Parameters(transfer_m_s=None, vmax_mol_m3_s=0), 2)
        self.assertIsNone(report["dimensionless_groups"]["scaled"]["biot_number"])
        self.assertIsNone(report["uptake_ratio"])
        self.assertEqual(report["rate_outputs"]["scaled"]["total_uptake_mol_s"], 0)

    def test_invalid_or_overflowing_scale_is_rejected(self):
        for scale in (True, 0, -1, float("nan"), float("inf"), "2"):
            with self.subTest(scale=scale), self.assertRaises(ValueError):
                equivalent_parameters(Parameters(), scale)
        with self.assertRaises(ValueError):
            equivalent_parameters(Parameters(vmax_mol_m3_s=1e100), 1e300)

    def test_report_refuses_existing_output_and_keeps_units_and_nonclaim(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report"
            report = equivalence_demo(path, Parameters(shells=20))
            self.assertFalse(report["biological_validation_performed"])
            self.assertIn("mol_m3", (path / "profiles.csv").read_text(encoding="utf-8"))
            self.assertIn("no measurements fitted", (path / "REPORT.md").read_text(encoding="utf-8"))
            with self.assertRaises(FileExistsError):
                equivalence_demo(path)


if __name__ == "__main__":
    unittest.main()
