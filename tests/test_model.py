from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import solve_bvp

from oxygenlab.cli import save, transfer_sweep_rows, vmax_sweep_rows
from oxygenlab.model import (
    Parameters,
    critical_radius,
    solve,
    zero_order_analytic,
    zero_order_critical_radius,
)


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

    def test_higher_vmax_does_not_raise_core_oxygen(self):
        rows = vmax_sweep_rows(Parameters(radius_um=300, shells=40))
        self.assertGreater(rows[0]["minimum_oxygen_mol_m3"], rows[-1]["minimum_oxygen_mol_m3"])
        self.assertTrue(all(row["core_did_not_rise"] for row in rows))
        self.assertTrue(all(row["relative_mass_balance_error"] < 1e-6 for row in rows))

    def test_higher_surface_transfer_raises_core_toward_fixed_surface_limit(self):
        rows = transfer_sweep_rows(Parameters(radius_um=300, shells=80))
        self.assertEqual(rows[-1]["boundary"], "fixed_surface")
        self.assertTrue(all(row["core_did_not_fall"] for row in rows))
        self.assertGreater(rows[-1]["minimum_oxygen_mol_m3"], rows[0]["minimum_oxygen_mol_m3"])
        self.assertTrue(all(row["relative_mass_balance_error"] < 1e-6 for row in rows))
        for values in [(1e-5, float("nan")), (0,), (1e-5, 1e-6), (None,)]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                transfer_sweep_rows(transfer_values_m_s=values)

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

    def test_transient_ramp_up_converges_to_steady(self):
        from oxygenlab.model import solve_transient
        p = Parameters(radius_um=150, shells=40)
        steady = solve(p)
        # 150 um sphere has tau_diff = (1.5e-4)^2 / 2e-9 = 11.25 s
        # In 120 s (10+ time constants), it reaches steady state
        sol = solve_transient(p, total_time_s=120.0, time_steps=40, initial_oxygen_mol_m3=0.0)
        self.assertEqual(len(sol.time_s), 40)
        self.assertEqual(sol.concentration_history_mol_m3.shape, (40, 40))
        self.assertEqual(sol.core_oxygen_history_mol_m3[0], 0.0)
        # Monotonically increasing core oxygen during ramp up from anoxia
        self.assertTrue(np.all(np.diff(sol.core_oxygen_history_mol_m3) >= -1e-12))
        # Long-time final profile converges to steady-state solution
        final_profile = sol.concentration_history_mol_m3[-1]
        np.testing.assert_allclose(final_profile, steady.concentration_mol_m3, atol=1e-4)
        self.assertIsNotNone(sol.summary["time_to_half_steady_core_s"])
        self.assertIsNotNone(sol.summary["time_to_95pct_steady_core_s"])
        self.assertLess(sol.summary["time_to_half_steady_core_s"], sol.summary["time_to_95pct_steady_core_s"])
        self.assertLess(sol.summary["max_scaled_step_residual"], 1e-8)
        self.assertLess(sol.summary["max_relative_transient_mass_balance_error"], 1e-7)

    def test_transient_rejects_unconverged_controls_and_unmodeled_initial_state(self):
        from oxygenlab.model import solve_transient
        p = Parameters(radius_um=150, shells=40)
        for kwargs in [
            {"max_iterations": 0},
            {"max_iterations": 1},
            {"initial_oxygen_mol_m3": p.bulk_oxygen_mol_m3 + 0.01},
            {"tolerance": float("nan")},
        ]:
            with self.subTest(kwargs=kwargs), self.assertRaises((ValueError, RuntimeError)):
                solve_transient(p, total_time_s=120, time_steps=2, **kwargs)

    def test_transient_timestep_refinement_is_reported_and_stable(self):
        from oxygenlab.model import solve_transient
        p = Parameters(radius_um=150, shells=40)
        coarse = solve_transient(p, total_time_s=120, time_steps=21)
        fine = solve_transient(p, total_time_s=120, time_steps=41)
        self.assertAlmostEqual(coarse.summary["final_core_oxygen_mol_m3"],
                               fine.summary["final_core_oxygen_mol_m3"], delta=2e-4)
        self.assertLess(coarse.summary["max_relative_transient_mass_balance_error"], 1e-6)
        self.assertLess(fine.summary["max_relative_transient_mass_balance_error"], 1e-6)

    def test_transient_parameter_validation(self):
        from oxygenlab.model import solve_transient
        p = Parameters()
        with self.assertRaises(ValueError):
            solve_transient(p, total_time_s=0)
        with self.assertRaises(ValueError):
            solve_transient(p, total_time_s=float("nan"))
        with self.assertRaises(ValueError):
            solve_transient(p, time_steps=1)
        with self.assertRaises(ValueError):
            solve_transient(p, initial_oxygen_mol_m3=-0.01)


if __name__ == "__main__":
    unittest.main()


class CriticalRadiusTests(unittest.TestCase):
    """The threshold-limited radius is checked against an independent closed form."""

    def test_closed_form_matches_bisection_for_both_surface_conditions(self):
        for label, parameters in (
            ("fixed surface", Parameters(kinetics="zero_order", transfer_m_s=None, shells=2000)),
            ("finite transfer", Parameters(kinetics="zero_order", shells=2000)),
        ):
            with self.subTest(surface=label):
                analytic = zero_order_critical_radius(parameters)
                numeric = critical_radius(parameters, tolerance_um=0.005, shells=2000)
                self.assertEqual(numeric["status"], "bracketed")
                # The bisection brackets the closed-form root within its tolerance.
                low, high = numeric["bracket_um"]
                self.assertLessEqual(low, analytic)
                self.assertLessEqual(analytic, high + 0.01)
                self.assertAlmostEqual(numeric["critical_radius_um"], analytic, delta=0.02)

    def test_closed_form_satisfies_the_balance_it_solves(self):
        parameters = Parameters(kinetics="zero_order", shells=400)
        radius_um = zero_order_critical_radius(parameters)
        at_radius = replace(parameters, radius_um=radius_um)
        # c(0) must equal the threshold at the critical radius, by construction.
        self.assertAlmostEqual(float(zero_order_analytic(0.0, at_radius)),
                               parameters.threshold_mol_m3, places=12)
        # A slightly larger sphere must fall below it.
        larger = replace(parameters, radius_um=radius_um * 1.01)
        self.assertLess(float(zero_order_analytic(0.0, larger)), parameters.threshold_mol_m3)

    def test_michaelis_menten_radius_is_at_least_the_zero_order_radius(self):
        result = critical_radius(Parameters(shells=400), tolerance_um=0.05, shells=400)
        self.assertEqual(result["status"], "bracketed")
        self.assertGreaterEqual(result["critical_radius_um"], result["zero_order_analytic_radius_um"])
        self.assertAlmostEqual(result["sampled_core_at_critical_radius_mol_m3"],
                               Parameters().threshold_mol_m3, delta=5e-4)

    def test_larger_uptake_shrinks_the_radius(self):
        radii = [critical_radius(Parameters(vmax_mol_m3_s=vmax, shells=200),
                                 tolerance_um=0.1, shells=200)["critical_radius_um"]
                 for vmax in (0.01, 0.02, 0.04)]
        self.assertTrue(radii[0] > radii[1] > radii[2], radii)

    def test_degenerate_cases_are_named_rather_than_guessed(self):
        self.assertEqual(critical_radius(Parameters(vmax_mol_m3_s=0.0))["status"],
                         "unbounded_without_uptake")
        self.assertIsNone(critical_radius(Parameters(vmax_mol_m3_s=0.0))["critical_radius_um"])
        at_threshold = critical_radius(Parameters(threshold_mol_m3=0.2))
        self.assertEqual(at_threshold["status"], "no_radius_qualifies_bath_at_or_below_threshold")
        self.assertIsNone(at_threshold["critical_radius_um"])
        self.assertIsNone(zero_order_critical_radius(Parameters(vmax_mol_m3_s=0.0)))
        self.assertIsNone(zero_order_critical_radius(Parameters(threshold_mol_m3=0.2)))

    def test_search_bound_is_reported_rather_than_extrapolated(self):
        result = critical_radius(Parameters(vmax_mol_m3_s=1e-6, shells=64), max_radius_um=50.0)
        self.assertEqual(result["status"], "exceeds_search_bound")
        self.assertIsNone(result["critical_radius_um"])
        self.assertIn("raise max_radius_um", " ".join(result["interpretation"]))

    def test_invalid_search_settings_are_rejected(self):
        for kwargs in ({"tolerance_um": 0}, {"tolerance_um": -1}, {"tolerance_um": float("nan")},
                       {"max_radius_um": 0}, {"max_radius_um": -5}, {"tolerance_um": True}):
            with self.subTest(**kwargs):
                with self.assertRaises(ValueError):
                    critical_radius(Parameters(), **kwargs)

    def test_infeasible_probe_radii_are_recorded_with_reasons(self):
        result = critical_radius(Parameters(kinetics="zero_order", shells=200),
                                 tolerance_um=0.05, shells=200)
        self.assertEqual(result["status"], "bracketed")
        self.assertTrue(result["infeasible_searched_radii"])
        for entry in result["infeasible_searched_radii"]:
            self.assertGreater(entry["radius_um"], result["critical_radius_um"])
            self.assertTrue(entry["reason"])

    def test_fine_mesh_zero_order_solves_instead_of_failing_its_residual_gate(self):
        # An exact linear solve carries backward error that grows with refinement.
        # These feasible configurations were rejected before the acceptance floor.
        for shells in (160, 2000, 5000):
            for radius_um in (50, 100, 150):
                with self.subTest(shells=shells, radius_um=radius_um):
                    parameters = Parameters(kinetics="zero_order", shells=shells, radius_um=radius_um)
                    summary = solve(parameters).summary
                    self.assertAlmostEqual(summary["minimum_sampled_oxygen_mol_m3"],
                                           float(zero_order_analytic(0.0, parameters)), places=9)
                    self.assertLessEqual(summary["scaled_residual"],
                                         summary["residual_acceptance_tolerance"])
        # The nonlinear branch keeps the stricter tolerance.
        nonlinear = solve(Parameters(shells=2000)).summary
        self.assertEqual(nonlinear["residual_acceptance_tolerance"], 1e-8)

    def test_negative_zero_order_oxygen_still_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "negative oxygen"):
            solve(Parameters(kinetics="zero_order", radius_um=600, shells=200))
