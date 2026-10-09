"""Constructed transient symmetry for jointly scaled oxygen transport rates."""

from dataclasses import asdict

import numpy as np

from oxygenlab.model import Parameters, solve_transient
from oxygenlab.steady_equivalence import equivalent_parameters


def transient_rate_equivalence(
    parameters=None,
    rate_scale=2.0,
    total_time_s=1200.0,
    time_steps=120,
    initial_oxygen_mol_m3=0.0,
):
    """Compare two constructed trajectories at the same dimensionless times.

    Scaling D, Vmax and finite surface transfer by ``rate_scale`` maps the
    original trajectory at t onto the scaled trajectory at t/rate_scale.
    This is a model symmetry, not an estimate of biological timing.
    """
    baseline_parameters = Parameters() if parameters is None else parameters
    scaled_parameters = equivalent_parameters(baseline_parameters, rate_scale)
    baseline = solve_transient(
        baseline_parameters,
        total_time_s=total_time_s,
        time_steps=time_steps,
        initial_oxygen_mol_m3=initial_oxygen_mol_m3,
    )
    scaled = solve_transient(
        scaled_parameters,
        total_time_s=total_time_s / rate_scale,
        time_steps=time_steps,
        initial_oxygen_mol_m3=initial_oxygen_mol_m3,
    )

    baseline_tau = baseline.time_s / baseline.summary["diffusion_time_scale_s"]
    scaled_tau = scaled.time_s / scaled.summary["diffusion_time_scale_s"]
    max_time_grid_difference = float(np.max(np.abs(baseline_tau - scaled_tau)))
    max_profile_difference = float(
        np.max(np.abs(baseline.concentration_history_mol_m3 - scaled.concentration_history_mol_m3))
    )
    max_steady_profile_difference = float(
        np.max(np.abs(baseline.steady_solution.concentration_mol_m3
                      - scaled.steady_solution.concentration_mol_m3))
    )

    baseline_95 = baseline.summary["time_to_95pct_initial_to_steady_core_change_s"]
    scaled_95 = scaled.summary["time_to_95pct_initial_to_steady_core_change_s"]
    t95_ratio = baseline_95 / scaled_95 if baseline_95 is not None and scaled_95 else None
    baseline_uptake = baseline.summary["final_total_uptake_mol_s"]
    scaled_uptake = scaled.summary["final_total_uptake_mol_s"]
    baseline_influx = baseline.summary["final_surface_influx_mol_s"]
    scaled_influx = scaled.summary["final_surface_influx_mol_s"]

    trajectory = [
        {
            "dimensionless_time_tau": float(tau),
            "baseline_time_s": float(t_base),
            "scaled_time_s": float(t_scaled),
            "baseline_core_oxygen_mol_m3": float(c_base),
            "scaled_core_oxygen_mol_m3": float(c_scaled),
            "baseline_volume_mean_oxygen_mol_m3": float(m_base),
            "scaled_volume_mean_oxygen_mol_m3": float(m_scaled),
        }
        for tau, t_base, t_scaled, c_base, c_scaled, m_base, m_scaled in zip(
            baseline_tau,
            baseline.time_s,
            scaled.time_s,
            baseline.core_oxygen_history_mol_m3,
            scaled.core_oxygen_history_mol_m3,
            baseline.volume_mean_history_mol_m3,
            scaled.volume_mean_history_mol_m3,
        )
    ]
    profiles = [
        {
            "radius_um": float(radius),
            "baseline_oxygen_mol_m3": float(c_base),
            "scaled_oxygen_mol_m3": float(c_scaled),
        }
        for radius, c_base, c_scaled in zip(
            baseline.radius_um,
            baseline.concentration_history_mol_m3[-1],
            scaled.concentration_history_mol_m3[-1],
        )
    ]
    return {
        "schema_version": 1,
        "origin": "constructed_transient_rate_equivalence",
        "rate_scale": float(rate_scale),
        "baseline_parameters": asdict(baseline_parameters),
        "scaled_parameters": asdict(scaled_parameters),
        "total_time_s": float(total_time_s),
        "time_steps": int(time_steps),
        "initial_oxygen_mol_m3": float(initial_oxygen_mol_m3),
        "time_relation": "scaled time = baseline time / rate_scale at matched dimensionless time",
        "max_dimensionless_time_grid_difference": max_time_grid_difference,
        "max_transient_profile_difference_mol_m3": max_profile_difference,
        "max_steady_profile_difference_mol_m3": max_steady_profile_difference,
        "time_to_95pct_core_change_s": {
            "baseline": baseline_95,
            "scaled": scaled_95,
            "baseline_over_scaled": t95_ratio,
            "expected_baseline_over_scaled": float(rate_scale) if t95_ratio is not None else None,
        },
        "final_rate_outputs": {
            "total_uptake_mol_s": {"baseline": float(baseline_uptake), "scaled": float(scaled_uptake)},
            "surface_influx_mol_s": {"baseline": float(baseline_influx), "scaled": float(scaled_influx)},
            "uptake_scaled_over_baseline": float(scaled_uptake / baseline_uptake) if baseline_uptake else None,
            "influx_scaled_over_baseline": float(scaled_influx / baseline_influx) if baseline_influx else None,
        },
        "max_relative_transient_mass_balance_error": {
            "baseline": float(baseline.summary["max_relative_transient_mass_balance_error"]),
            "scaled": float(scaled.summary["max_relative_transient_mass_balance_error"]),
        },
        "trajectory": trajectory,
        "final_matched_time_profiles": profiles,
        "model_fitted": False,
        "biological_validation_performed": False,
        "limits": [
            "For the same geometry, bath, uptake law, Km, initial profile and dimensionless time grid, jointly scaling diffusivity, maximal uptake and finite surface transfer (when present) preserves the backward-Euler transient concentration path when physical time is divided by that factor; a fixed-surface boundary remains fixed.",
            "The characteristic diffusion scale R²/D and model-derived core-change times scale inversely with the rate factor. These are mathematical model times, not measured equilibration or biological response times.",
            "Absolute uptake and surface influx scale with the rate factor at matched states; mass-balance checks are reported separately for both integrations.",
            "No measurements were fitted. The constructed symmetry establishes no organoid-specific rates, oxygen response, viability, culture recommendation or biological validation.",
        ],
    }
