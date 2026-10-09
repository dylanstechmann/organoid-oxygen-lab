"""Independent geometric references; no synthetic source promoted to measured biology."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from oxygenlab.cli import main
from oxygenlab.spatial_coverage import radial_coverage


class SpatialCoverageTests(unittest.TestCase):
    def test_volume_regions_use_spherical_geometry_not_linear_radius(self):
        result = radial_coverage([20., 50., 80.], 100.)
        self.assertAlmostEqual(result["inward_of_first_sample_sphere_volume_fraction"], .008)
        self.assertAlmostEqual(result["outward_of_last_sample_sphere_volume_fraction"], .488)
        self.assertAlmostEqual(result["largest_internal_radial_gap_fraction"], .3)
        self.assertFalse(result["center_sample_present"])
        self.assertFalse(result["surface_sample_present"])

    def test_boundary_points_do_not_eliminate_internal_gaps(self):
        result = radial_coverage([0., 50., 100.], 100.)
        self.assertTrue(result["center_sample_present"])
        self.assertTrue(result["surface_sample_present"])
        self.assertEqual(result["outward_of_last_sample_sphere_volume_fraction"], 0)
        self.assertEqual(result["largest_internal_radial_gap_fraction"], .5)

    def test_single_point_has_no_internal_gap_statistic(self):
        result = radial_coverage([50.], 100.)
        self.assertIsNone(result["largest_internal_radial_gap_fraction"])
        self.assertAlmostEqual(result["inward_of_first_sample_sphere_volume_fraction"] + result["outward_of_last_sample_sphere_volume_fraction"], 1)

    def test_outside_nonfinite_boolean_and_duplicate_radii_fail(self):
        for radii, radius in [([0.,101.],100.),([0.,float('nan')],100.),([0.,True],100.),([0.,50.,50.],100.),([0.],0.),([0.],True),([],100.)]:
            with self.assertRaises(ValueError):
                radial_coverage(radii, radius)

    def test_cli_demo_is_explicit_and_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "new"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["coverage-demo","--out",str(out)]),0)
            report = json.loads((out / "coverage.json").read_text())
            self.assertFalse(report["biological_data_used"])
            self.assertFalse(report["oxygen_or_viability_inferred"])
            before = (out / "coverage.json").read_bytes()
            with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(SystemExit) as raised:
                    main(["coverage-demo","--out",str(out)])
                self.assertEqual(raised.exception.code,2)
            self.assertEqual((out / "coverage.json").read_bytes(),before)


if __name__ == "__main__":
    unittest.main()
