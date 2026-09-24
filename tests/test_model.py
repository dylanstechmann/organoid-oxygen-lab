from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import solve_bvp

from oxygenlab.cli import save
from oxygenlab.model import Parameters, solve, zero_order_analytic


class OxygenTests(unittest.TestCase):
    def test_no_consumption_returns_bath_concentration(self):
        for transfer in [None, 1e-5]:
            result = solve(Parameters(vmax_mol_m3_s=0, transfer_m_s=transfer))
            np.testing.assert_array_equal(result.concentration_mol_m3, np.full(160, .2))
            self.assertEqual(result.summary["total_uptake_mol_s"], 0)
            self.assertEqual(result.summary["surface_influx_mol_s"], 0)

    def test_constant_uptake_matches_analytic_and_refines_quadratically(self):
        for transfer in [None, 2e-5]:
            errors = []
            for shells in [20, 40, 80]:
                p = Parameters(radius_um=100, shells=shells, transfer_m_s=transfer, kinetics="zero_order")
                result = solve(p)
                errors.append(np.max(np.abs(result.concentration_mol_m3 - zero_order_analytic(result.radius_um, p))))
                expected_uptake = 4 * np.pi / 3 * (p.radius_um * 1e-6) ** 3 * p.vmax_mol_m3_s
                self.assertAlmostEqual(result.summary["total_uptake_mol_s"] / expected_uptake, 1, places=12)
            self.assertGreater(errors[0] / errors[1], 3.9)
            self.assertGreater(errors[1] / errors[2], 3.9)
            self.assertLess(errors[-1], 1e-6)

    def test_nonlinear_profile_agrees_with_independent_boundary_value_solver(self):
        p = Parameters(radius_um=300, shells=240)
        result = solve(p)
        radius = p.radius_um * 1e-6
        da = p.vmax_mol_m3_s * radius ** 2 / (p.diffusivity_m2_s * p.bulk_oxygen_mol_m3)
        km = p.km_mol_m3 / p.bulk_oxygen_mol_m3
        bi = p.transfer_m_s * radius / p.diffusivity_m2_s
        def rhs(x, y): return np.vstack((y[1], da * y[0] / (km + y[0])))
        def boundary(a, b): return np.array([a[1], b[0] + b[1] / bi - 1])
        x = np.linspace(0, 1, 100)
        initial = np.vstack((np.interp(x, result.radius_um / p.radius_um, result.concentration_mol_m3 / p.bulk_oxygen_mol_m3), np.zeros(len(x))))
        reference = solve_bvp(rhs, boundary, x, initial, S=np.array([[0., 0.], [0., -2.]]), tol=1e-8, max_nodes=10000)
        self.assertTrue(reference.success, reference.message)
        expected = reference.sol(result.radius_um / p.radius_um)[0] * p.bulk_oxygen_mol_m3
        np.testing.assert_allclose(result.concentration_mol_m3, expected, atol=2e-6, rtol=1e-4)

    def test_positive_monotone_profiles_and_conserved_flux(self):
        for radius in [50, 300, 800]:
            result = solve(Parameters(radius_um=radius))
            c = result.concentration_mol_m3
            self.assertTrue(np.all(c >= 0)); self.assertTrue(np.all(c <= .2))
            self.assertTrue(np.all(np.diff(c) >= -1e-12))
            self.assertLess(result.summary["relative_mass_balance_error"], 1e-7)
            self.assertAlmostEqual(result.shell_volume_fraction.sum(), 1, places=12)

    def test_boundary_resistance_lowers_oxygen(self):
        fixed = solve(Parameters(transfer_m_s=None))
        limited = solve(Parameters(transfer_m_s=1e-5))
        self.assertTrue(np.all(limited.concentration_mol_m3 < fixed.concentration_mol_m3))
        self.assertLess(limited.summary["surface_oxygen_mol_m3"], .2)

    def test_threshold_fraction_is_volume_weighted(self):
        base = Parameters(radius_um=300)
        initial = solve(base)
        threshold = float(initial.concentration_mol_m3[80])
        result = solve(replace(base, threshold_mol_m3=threshold))
        self.assertAlmostEqual(result.summary["fraction_volume_below_threshold"], .5 ** 3, places=12)
        self.assertEqual(solve(replace(base, threshold_mol_m3=0)).summary["fraction_volume_below_threshold"], 0)

    def test_invalid_parameters_and_infeasible_constant_uptake(self):
        for field, value in [("radius_um", 0), ("diffusivity_m2_s", -1), ("vmax_mol_m3_s", float("nan")),
                             ("transfer_m_s", 0), ("km_mol_m3", 0), ("shells", True), ("shells", 4),
                             ("kinetics", "unknown"), ("radius_um", "300")]:
            with self.assertRaises(ValueError): solve(replace(Parameters(), **{field: value}))
        with self.assertRaisesRegex(ValueError, "negative oxygen"):
            solve(Parameters(radius_um=1000, kinetics="zero_order"))
        with self.assertRaises(RuntimeError): solve(Parameters(), max_iterations=1)

    def test_configuration_and_outputs_preserve_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); config = root / "p.json"
            config.write_text('{"radius_um":100,"transfer_m_s":null}')
            p = Parameters.from_json(config)
            result = solve(p)
            report = save(result, root / "result", input_sha256="a" * 64)
            self.assertEqual(report["input_sha256"], "a" * 64)
            self.assertEqual(len(report["parameters_sha256"]), 64)
            self.assertEqual(json.loads((root / "result/summary.json").read_text())["parameters"]["radius_um"], 100)
            with self.assertRaises(FileExistsError): save(result, root / "result")
            config.write_text('{"radius_microns_typo":100}')
            with self.assertRaises(ValueError): Parameters.from_json(config)


if __name__ == "__main__":
    unittest.main()
