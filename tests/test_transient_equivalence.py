import json
import tempfile
import unittest
from pathlib import Path

from oxygenlab.cli import transient_equivalence_demo
from oxygenlab.model import Parameters
from oxygenlab.transient_equivalence import transient_rate_equivalence


class TransientEquivalenceTests(unittest.TestCase):
    def test_joint_rate_scaling_preserves_path_and_scales_model_time(self):
        parameters = Parameters(radius_um=150, shells=40)
        report = transient_rate_equivalence(
            parameters, rate_scale=3, total_time_s=180, time_steps=61,
            initial_oxygen_mol_m3=0.04,
        )
        self.assertLess(report["max_dimensionless_time_grid_difference"], 1e-14)
        self.assertLess(report["max_transient_profile_difference_mol_m3"], 1e-9)
        self.assertLess(report["max_steady_profile_difference_mol_m3"], 1e-9)
        self.assertAlmostEqual(report["time_to_95pct_core_change_s"]["baseline_over_scaled"], 3)
        self.assertAlmostEqual(report["final_rate_outputs"]["uptake_scaled_over_baseline"], 3)
        self.assertAlmostEqual(report["final_rate_outputs"]["influx_scaled_over_baseline"], 3)
        self.assertLess(report["max_relative_transient_mass_balance_error"]["baseline"], 1e-8)
        self.assertLess(report["max_relative_transient_mass_balance_error"]["scaled"], 1e-8)
        self.assertFalse(report["model_fitted"])
        self.assertFalse(report["biological_validation_performed"])

    def test_fixed_surface_boundary_obeys_the_same_transient_time_scaling(self):
        report = transient_rate_equivalence(
            Parameters(radius_um=120, shells=24, transfer_m_s=None),
            rate_scale=2.5, total_time_s=100, time_steps=31,
        )
        self.assertIsNone(report["scaled_parameters"]["transfer_m_s"])
        self.assertLess(report["max_transient_profile_difference_mol_m3"], 1e-9)
        self.assertAlmostEqual(report["time_to_95pct_core_change_s"]["baseline_over_scaled"], 2.5)

    def test_invalid_rate_scale_is_rejected(self):
        for scale in (True, 0, -1, float("nan"), float("inf"), "2"):
            with self.subTest(scale=scale), self.assertRaises(ValueError):
                transient_rate_equivalence(Parameters(), rate_scale=scale)

    def test_report_writes_reproducible_json_csv_and_caveats(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "transient-equivalence"
            report = transient_equivalence_demo(
                output, Parameters(radius_um=120, shells=24), rate_scale=2,
                total_time_s=80, time_steps=21,
            )
            self.assertEqual(json.loads((output / "transient_equivalence.json").read_text(encoding="utf-8")), report)
            self.assertIn("dimensionless_time_tau", (output / "matched_trajectories.csv").read_text(encoding="utf-8"))
            self.assertIn("No biological measurements", (output / "REPORT.md").read_text(encoding="utf-8"))
            with self.assertRaises(FileExistsError):
                transient_equivalence_demo(output)


if __name__ == "__main__":
    unittest.main()
