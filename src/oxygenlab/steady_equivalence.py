"""An exact steady-model symmetry: concentration alone cannot fix absolute rates."""

import math
from dataclasses import asdict, replace

import numpy as np

from oxygenlab.model import Parameters, solve


def equivalent_parameters(parameters, rate_scale):
    """Scale D, maximal uptake and finite film transfer together; keep R/bath/Km."""
    parameters.validate()
    if isinstance(rate_scale, bool) or not isinstance(rate_scale, (int, float)) or not math.isfinite(rate_scale) or rate_scale <= 0:
        raise ValueError("rate_scale must be positive and finite")
    return replace(parameters,
                   diffusivity_m2_s=parameters.diffusivity_m2_s * rate_scale,
                   vmax_mol_m3_s=parameters.vmax_mol_m3_s * rate_scale,
                   transfer_m_s=None if parameters.transfer_m_s is None else parameters.transfer_m_s * rate_scale).validate()


def steady_equivalence(parameters=None, rate_scale=2.0):
    """Compare constructed settings, not measurements; no parameter calibration."""
    parameters = Parameters() if parameters is None else parameters
    changed = equivalent_parameters(parameters, rate_scale)
    baseline, scaled = solve(parameters), solve(changed)
    uptake = baseline.summary["total_uptake_mol_s"]
    report = {
        "schema_version": 1, "origin": "constructed_model_equivalence",
        "rate_scale": rate_scale, "baseline_parameters": asdict(parameters), "scaled_parameters": asdict(changed),
        "max_profile_difference_mol_m3": float(np.max(np.abs(scaled.concentration_mol_m3 - baseline.concentration_mol_m3))),
        "profiles": [{"radius_um": float(r), "baseline_oxygen_mol_m3": float(a), "scaled_oxygen_mol_m3": float(b)}
                     for r, a, b in zip(baseline.radius_um, baseline.concentration_mol_m3, scaled.concentration_mol_m3)],
        "dimensionless_groups": {label: {"damkohler_number": s.summary["damkohler_number"],
                                         "biot_number": s.summary["biot_number"],
                                         "km_over_bulk": p.km_mol_m3 / p.bulk_oxygen_mol_m3}
                                 for label, s, p in (("baseline", baseline, parameters), ("scaled", scaled, changed))},
        "rate_outputs": {label: {"total_uptake_mol_s": s.summary["total_uptake_mol_s"],
                                  "surface_influx_mol_s": s.summary["surface_influx_mol_s"],
                                  "relative_mass_balance_error": s.summary["relative_mass_balance_error"],
                                  "diffusion_time_scale_s": (p.radius_um * 1e-6) ** 2 / p.diffusivity_m2_s}
                         for label, s, p in (("baseline", baseline, parameters), ("scaled", scaled, changed))},
        "uptake_ratio": None if uptake == 0 else scaled.summary["total_uptake_mol_s"] / uptake,
        "model_fitted": False, "biological_validation_performed": False,
        "limits": ["Within this homogeneous steady sphere, scaling D, Vmax and finite k together preserves Da, Bi and Km/c_bulk, hence the concentration profile at the same physical radii.",
                   "Absolute uptake and influx scale; R^2/D changes. That diffusion scale is not a measured equilibration time.",
                   "More precise steady concentrations alone do not remove this symmetry. Independent rate, diffusivity or temporal information could constrain it only under valid model/measurement assumptions.",
                   "No unique fitted parameters, experimental identifiability, viability or recommended culture conditions are established."]}
    return report
