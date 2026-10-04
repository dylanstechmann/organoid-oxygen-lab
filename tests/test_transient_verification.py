"""Independent PDE reference and meaningful transient-exposure summaries."""

import unittest

import numpy as np

from oxygenlab.model import Parameters, solve_transient


class TransientVerificationTests(unittest.TestCase):
    def test_early_diffusion_refines_against_independent_spherical_series(self):
        # With fixed surface and no uptake, w=x*(1-c/c_bulk) solves the
        # one-dimensional heat equation with w=0 at x=0,1 and w(x,0)=x.
        # Its sine coefficients are 2*(-1)^(n+1)/(n*pi), independent of the FV code.
        p = Parameters(radius_um=150, shells=160, transfer_m_s=None, vmax_mol_m3_s=0)
        diffusion_time = (p.radius_um * 1e-6) ** 2 / p.diffusivity_m2_s
        tau = 0.1  # Early dynamics, well before convergence to steady oxygen.
        modes = np.arange(1, 101)
        coefficients = 2 * (-1.0) ** (modes + 1) / (modes * np.pi)
        coefficients *= np.exp(-modes ** 2 * np.pi ** 2 * tau)
        errors = []
        for time_points in [33, 65, 129]:
            solution = solve_transient(p, total_time_s=tau * diffusion_time, time_steps=time_points)
            x = solution.radius_um / p.radius_um
            deficit = coefficients @ (np.sin(modes[:, None] * np.pi * x[None, :]) / x[None, :])
            reference = p.bulk_oxygen_mol_m3 * (1 - deficit)
            errors.append(float(np.max(np.abs(solution.concentration_history_mol_m3[-1] - reference))))
            self.assertLess(solution.summary["max_relative_transient_mass_balance_error"], 1e-8)
        # Backward Euler is first order in time, unlike the spatial O(dx^2) check.
        self.assertGreater(errors[0] / errors[1], 1.9)
        self.assertGreater(errors[1] / errors[2], 1.9)
        self.assertLess(errors[-1], 4e-4)

    def test_core_change_times_describe_both_rise_and_relaxation_from_above(self):
        p = Parameters(radius_um=150, shells=40)
        for initial in [0.05, p.bulk_oxygen_mol_m3]:
            solution = solve_transient(p, total_time_s=120, time_steps=40,
                                       initial_oxygen_mol_m3=initial)
            summary = solution.summary
            times = [summary[f"time_to_{fraction}pct_initial_to_steady_core_change_s"]
                     for fraction in [50, 95]]
            with self.subTest(initial=initial):
                self.assertGreater(times[0], 0)
                self.assertLess(times[0], times[1])
                self.assertLess(times[1], 120)
                for fraction, time in zip([50, 95], times):
                    target = initial + fraction / 100 * (summary["steady_core_oxygen_mol_m3"] - initial)
                    self.assertAlmostEqual(summary["core_change_targets_mol_m3"][f"{fraction}pct"], target)
                    self.assertAlmostEqual(np.interp(time, solution.time_s,
                                                     solution.core_oxygen_history_mol_m3), target)

    def test_core_already_at_steady_value_has_zero_change_time(self):
        p = Parameters(shells=20, vmax_mol_m3_s=0)
        summary = solve_transient(p, total_time_s=1, time_steps=3,
                                   initial_oxygen_mol_m3=p.bulk_oxygen_mol_m3).summary
        for fraction in [50, 95]:
            self.assertEqual(summary[f"time_to_{fraction}pct_initial_to_steady_core_change_s"], 0)

    def test_core_change_targets_outside_horizon_are_not_reported_as_reached(self):
        p = Parameters(radius_um=150, shells=40)
        summary = solve_transient(p, total_time_s=0.1, time_steps=4,
                                   initial_oxygen_mol_m3=p.bulk_oxygen_mol_m3).summary
        for fraction in [50, 95]:
            self.assertIsNone(summary[f"time_to_{fraction}pct_initial_to_steady_core_change_s"])

    def test_duplicate_config_keys_cannot_hide_parameter_edits(self):
        with self.assertRaisesRegex(ValueError, "duplicate configuration key"):
            Parameters.from_json_bytes(b'{"radius_um":100,"radius_um":200}')
        self.assertEqual(Parameters.from_json_bytes(b'\xef\xbb\xbf{"radius_um":100}').radius_um, 100)


if __name__ == "__main__":
    unittest.main()
