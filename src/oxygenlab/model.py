"""Spherical finite-volume diffusion with oxygen-dependent uptake.

Concentrations are mol/m^3 (= mM); lengths are converted to SI once.
The surrounding bath concentration is fixed. No flow field or cell death is
modeled. A finite surface mass-transfer coefficient represents a boundary
resistance, not a simulation of a pump or its perfusion rate.
"""

from dataclasses import asdict, dataclass, replace
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded


@dataclass(frozen=True)
class Parameters:
    radius_um: float = 300.0
    diffusivity_m2_s: float = 2e-9
    bulk_oxygen_mol_m3: float = 0.2
    vmax_mol_m3_s: float = 0.02
    km_mol_m3: float = 0.01
    transfer_m_s: float | None = 2e-5
    threshold_mol_m3: float = 0.02
    shells: int = 160
    kinetics: str = "michaelis_menten"

    def validate(self):
        for name in ["radius_um", "diffusivity_m2_s", "bulk_oxygen_mol_m3", "km_mol_m3"]:
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a positive finite number")
        for name in ["vmax_mol_m3_s", "threshold_mol_m3"]:
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be a nonnegative finite number")
        if self.transfer_m_s is not None and (isinstance(self.transfer_m_s, bool)
                or not isinstance(self.transfer_m_s, (int, float))
                or not np.isfinite(self.transfer_m_s) or self.transfer_m_s <= 0):
            raise ValueError("transfer_m_s must be positive and finite, or null for a fixed surface concentration")
        if isinstance(self.shells, bool) or not isinstance(self.shells, int) or not 8 <= self.shells <= 5000:
            raise ValueError("shells must be an integer between 8 and 5000")
        if not isinstance(self.kinetics, str) or self.kinetics not in {"michaelis_menten", "zero_order"}:
            raise ValueError("kinetics must be michaelis_menten or zero_order")
        return self

    @classmethod
    def from_json(cls, path):
        return cls.from_json_bytes(Path(path).read_bytes())

    @classmethod
    def from_json_bytes(cls, data):
        """Parse one UTF-8 snapshot that can also be hashed by a caller."""
        def unique_object(pairs):
            result = {}
            for key, value in pairs:
                if key in result:
                    raise ValueError(f"duplicate configuration key: {key}")
                result[key] = value
            return result

        values = json.loads(data.decode("utf-8-sig"), object_pairs_hook=unique_object)
        if not isinstance(values, dict) or set(values) - set(cls.__dataclass_fields__):
            raise ValueError("configuration must be an object containing only recognized parameter names")
        return cls(**values).validate()


@dataclass
class Solution:
    radius_um: np.ndarray
    concentration_mol_m3: np.ndarray
    shell_volume_fraction: np.ndarray
    summary: dict


@dataclass
class TransientSolution:
    time_s: np.ndarray
    radius_um: np.ndarray
    concentration_history_mol_m3: np.ndarray
    core_oxygen_history_mol_m3: np.ndarray
    surface_oxygen_history_mol_m3: np.ndarray
    volume_mean_history_mol_m3: np.ndarray
    shell_volume_fraction: np.ndarray
    steady_solution: Solution
    summary: dict


def _system(parameters):
    p = parameters
    radius = p.radius_um * 1e-6
    edges = np.linspace(0, 1, p.shells + 1)
    x = (edges[:-1] + edges[1:]) / 2
    dx = 1 / p.shells
    volumes = np.diff(edges ** 3) / 3
    conductance = edges[1:-1] ** 2 / dx
    boundary_resistance = 0 if p.transfer_m_s is None else p.diffusivity_m2_s / (p.transfer_m_s * radius)
    boundary_g = 1 / (dx / 2 + boundary_resistance)
    diagonal = np.zeros(p.shells)
    diagonal[:-1] += conductance
    diagonal[1:] += conductance
    diagonal[-1] += boundary_g
    banded = np.zeros((3, p.shells))
    banded[0, 1:] = -conductance
    banded[1] = diagonal
    banded[2, :-1] = -conductance
    boundary = np.zeros(p.shells); boundary[-1] = boundary_g
    damkohler = p.vmax_mol_m3_s * radius ** 2 / (p.diffusivity_m2_s * p.bulk_oxygen_mol_m3)
    if not np.isfinite(damkohler) or not np.isfinite(banded).all() or boundary_g <= 0:
        raise ValueError("parameters exceed the representable dimensionless range")
    return x, volumes, conductance, banded, boundary, boundary_g, damkohler


def solve(parameters=Parameters(), *, tolerance=1e-8, max_iterations=300):
    p = parameters.validate()
    if not np.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError("tolerance must lie in (0, 1)")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    x, volumes, conductance, banded, boundary, boundary_g, da = _system(p)
    km = p.km_mol_m3 / p.bulk_oxygen_mol_m3

    def residual(u):
        # Flux differences avoid subtracting large diagonal terms in flat profiles.
        balance = np.zeros_like(u)
        flux = conductance * (u[1:] - u[:-1])
        balance[:-1] -= flux; balance[1:] += flux
        balance[-1] += boundary_g * (u[-1] - 1)
        rate = np.ones_like(u) if p.kinetics == "zero_order" else u / (km + u)
        return balance + da * volumes * rate

    def norm(u):
        return float(np.max(np.abs(residual(u)) / volumes) / max(1, da))

    u = np.ones(p.shells)
    iterations = 0
    # An exact linear solve still carries backward error of order eps*||A||, which grows
    # with mesh refinement. Judging it by the Newton iteration's absolute tolerance
    # rejected feasible fine-mesh configurations, so the linear branch accepts a
    # round-off floor. The nonlinear branch keeps the stricter tolerance it enforces above.
    accept_tolerance = tolerance
    if p.vmax_mol_m3_s == 0:
        pass
    elif p.kinetics == "zero_order":
        accept_tolerance = max(
            tolerance,
            float(np.finfo(float).eps * np.max(np.abs(banded)) / np.min(volumes)) / max(1, da),
        )
        u = solve_banded((1, 1), banded, boundary - da * volumes)
        iterations = 1
        if np.any(u < 0):
            raise ValueError("zero-order uptake predicts negative oxygen; this assumption is infeasible for these parameters")
    else:
        # q(c) <= Vmax*c/Km: the first-order uptake solution is a positive
        # lower bound. Newton steps from below avoid exhausting positivity at
        # an almost anoxic center while the outer shells are still changing.
        lower_system = banded.copy()
        lower_system[1] += da * volumes / km
        u = solve_banded((1, 1), lower_system, boundary)
        for iterations in range(1, max_iterations + 1):
            old_norm = norm(u)
            if old_norm <= tolerance:
                break
            jacobian = banded.copy()
            jacobian[1] += da * volumes * km / (km + u) ** 2
            step = solve_banded((1, 1), jacobian, -residual(u))
            negative = step < 0
            alpha = min(1.0, 0.99 * float(np.min(-u[negative] / step[negative]))) if negative.any() else 1.0
            accepted = False
            for _ in range(60):
                candidate = u + alpha * step
                if np.all(candidate >= 0) and norm(candidate) < old_norm:
                    u = candidate; accepted = True; break
                alpha *= 0.5
            if not accepted:
                raise RuntimeError("nonlinear line search failed; refine parameters or mesh")
        if norm(u) > tolerance:
            raise RuntimeError("nonlinear solver did not meet the requested residual tolerance")
    if norm(u) > accept_tolerance or not np.isfinite(u).all() or np.any(u > 1 + 1e-8):
        raise RuntimeError("solution failed residual or physical-bound checks")
    radius = p.radius_um * 1e-6
    rate = np.ones_like(u) if p.kinetics == "zero_order" else u / (km + u)
    uptake = float(4 * np.pi * radius ** 3 * p.vmax_mol_m3_s * np.sum(volumes * rate))
    influx = float(4 * np.pi * radius * p.diffusivity_m2_s * p.bulk_oxygen_mol_m3 * boundary_g * (1 - u[-1]))
    surface = p.bulk_oxygen_mol_m3 if p.transfer_m_s is None else p.bulk_oxygen_mol_m3 - influx / (4 * np.pi * radius ** 2 * p.transfer_m_s)
    concentrations = u * p.bulk_oxygen_mol_m3
    summary = {"parameters": asdict(p), "iterations": iterations, "scaled_residual": norm(u),
        "residual_acceptance_tolerance": float(accept_tolerance),
        "damkohler_number": float(da), "biot_number": None if p.transfer_m_s is None else float(p.transfer_m_s * radius / p.diffusivity_m2_s),
        "minimum_sampled_oxygen_mol_m3": float(concentrations.min()),
        "surface_oxygen_mol_m3": float(surface), "volume_mean_oxygen_mol_m3": float(np.sum(3 * volumes * concentrations)),
        "fraction_volume_below_threshold": float(np.sum(3 * volumes[concentrations < p.threshold_mol_m3])),
        "total_uptake_mol_s": uptake, "surface_influx_mol_s": influx,
        "relative_mass_balance_error": float(abs(influx - uptake) / max(abs(influx), abs(uptake), 1e-30)),
        "interpretation": ["Illustrative model, not calibrated to a cell line or organoid experiment.",
            "The minimum is sampled at the innermost shell center, not exactly at r=0.",
            "Below-threshold volume uses piecewise-constant shell values; the threshold is user supplied, not a viability cutoff.",
            "Spherical, homogeneous, steady state with a fixed bath; no vascularization, growth, heterogeneous uptake or cell death.",
            "Surface mass transfer is a boundary resistance and cannot be converted to pump flow without another model."]}
    return Solution(x * p.radius_um, concentrations, 3 * volumes, summary)


def zero_order_analytic(radius_um, parameters):
    """Independent closed-form reference at requested radial positions."""
    p = parameters.validate()
    r = np.asarray(radius_um, dtype=float) * 1e-6
    radius = p.radius_um * 1e-6
    if not np.isfinite(r).all() or np.any(r < 0) or np.any(r > radius):
        raise ValueError("reference radii must lie within the sphere")
    surface_drop = 0 if p.transfer_m_s is None else p.vmax_mol_m3_s * radius / (3 * p.transfer_m_s)
    return p.bulk_oxygen_mol_m3 - surface_drop - p.vmax_mol_m3_s * (radius ** 2 - r ** 2) / (6 * p.diffusivity_m2_s)


def zero_order_critical_radius(parameters):
    """Largest radius whose exact zero-order core concentration reaches the threshold.

    Setting c(0) = threshold in the closed-form zero-order solution gives a
    quadratic in R. With a finite surface transfer coefficient the surface drop
    contributes the linear term; a fixed surface concentration drops it:

        (vmax / 6D) R^2 + (vmax / 3k) R - (c_bulk - threshold) = 0

    Returns ``None`` when the balance has no positive root, which happens when
    uptake is zero (no radius is limited) or when the bath concentration is
    already at or below the threshold (no radius qualifies).
    """
    p = parameters.validate()
    headroom = p.bulk_oxygen_mol_m3 - p.threshold_mol_m3
    if p.vmax_mol_m3_s == 0 or headroom <= 0:
        return None
    quadratic = p.vmax_mol_m3_s / (6 * p.diffusivity_m2_s)
    linear = 0.0 if p.transfer_m_s is None else p.vmax_mol_m3_s / (3 * p.transfer_m_s)
    radius_m = (math.sqrt(linear ** 2 + 4 * quadratic * headroom) - linear) / (2 * quadratic)
    return float(radius_m * 1e6)


def critical_radius(parameters=Parameters(), *, tolerance_um=0.01, max_radius_um=20000.0,
                    shells=None):
    """Find the largest sphere whose sampled core oxygen stays at the threshold.

    The numerical search bisects on the solver's minimum sampled concentration,
    which is monotone decreasing in radius for fixed uptake. That minimum sits at
    the innermost shell center rather than exactly at r=0, so the reported radius
    is the largest radius whose *sampled* core reaches the threshold; refining
    ``shells`` tightens the gap. The threshold is a user-selected reporting level,
    not a viability, hypoxia or potency cutoff.
    """
    p = parameters.validate()
    if isinstance(tolerance_um, bool) or not isinstance(tolerance_um, (int, float)) \
            or not math.isfinite(tolerance_um) or tolerance_um <= 0:
        raise ValueError("tolerance_um must be a positive finite number")
    if isinstance(max_radius_um, bool) or not isinstance(max_radius_um, (int, float)) \
            or not math.isfinite(max_radius_um) or max_radius_um <= 0:
        raise ValueError("max_radius_um must be a positive finite number")
    search_shells = p.shells if shells is None else shells
    analytic = zero_order_critical_radius(replace(p, kinetics="zero_order"))

    infeasible_radii = []

    def sampled_core(radius_um):
        """Sampled core oxygen, treating an infeasible zero-order radius as below threshold.

        Zero-order uptake predicts negative oxygen beyond a finite radius and the
        solver fails closed there. Such a radius cannot hold the core at any
        positive threshold, so the search reads it as below threshold and records it.
        """
        candidate = replace(p, radius_um=radius_um, shells=search_shells).validate()
        try:
            return float(solve(candidate).summary["minimum_sampled_oxygen_mol_m3"])
        except (ValueError, RuntimeError) as exc:
            # The solver fails closed when a configuration predicts negative oxygen or
            # misses its residual/physical bounds. Record which radius and why rather
            # than swallowing it.
            infeasible_radii.append({"radius_um": float(radius_um),
                                     "reason": f"{type(exc).__name__}: {exc}"})
            return float("-inf")

    result = {
        "schema_version": 1,
        "analysis": "largest radius whose sampled core oxygen reaches the reporting threshold",
        "parameters": asdict(p),
        "threshold_mol_m3": p.threshold_mol_m3,
        "search_shells": search_shells,
        "tolerance_um": tolerance_um,
        "zero_order_analytic_radius_um": analytic,
        "critical_radius_um": None,
        "status": None,
        "bracket_um": None,
        "bisection_iterations": 0,
        "sampled_core_at_critical_radius_mol_m3": None,
        "infeasible_searched_radii": infeasible_radii,
        "interpretation": [],
    }
    if p.bulk_oxygen_mol_m3 <= p.threshold_mol_m3:
        result["status"] = "no_radius_qualifies_bath_at_or_below_threshold"
        result["interpretation"].append(
            "The bath concentration is already at or below the reporting threshold, so no sphere qualifies.")
        return result
    if p.vmax_mol_m3_s == 0:
        result["status"] = "unbounded_without_uptake"
        result["interpretation"].append(
            "Uptake is zero, so the sphere never falls below the bath concentration and no radius is limiting.")
        return result
    if sampled_core(max_radius_um) >= p.threshold_mol_m3:
        result["status"] = "exceeds_search_bound"
        result["bracket_um"] = [max_radius_um, None]
        result["interpretation"].append(
            f"Even a {max_radius_um:g} um sphere keeps its sampled core at the threshold under these "
            "illustrative parameters; raise max_radius_um to search further.")
        return result

    low = min(p.radius_um, max_radius_um) / 2
    while low > tolerance_um and sampled_core(low) < p.threshold_mol_m3:
        low /= 2
    if sampled_core(low) < p.threshold_mol_m3:
        result["status"] = "below_resolution"
        result["bracket_um"] = [0.0, low]
        result["interpretation"].append(
            "Even the smallest searched sphere falls below the threshold at this uptake; the critical radius "
            "is smaller than the requested tolerance.")
        return result
    high = low * 2
    while high <= max_radius_um and sampled_core(high) >= p.threshold_mol_m3:
        low, high = high, high * 2
    high = min(high, max_radius_um)

    iterations = 0
    while high - low > tolerance_um:
        middle = 0.5 * (low + high)
        if sampled_core(middle) >= p.threshold_mol_m3:
            low = middle
        else:
            high = middle
        iterations += 1
    result.update({
        "critical_radius_um": float(low),
        "status": "bracketed",
        "bracket_um": [float(low), float(high)],
        "bisection_iterations": iterations,
        "sampled_core_at_critical_radius_mol_m3": sampled_core(low),
    })
    result["interpretation"].extend([
        "Bisection on the solver's minimum sampled concentration, which decreases monotonically with radius "
        "at fixed uptake.",
        "The sampled minimum lies at the innermost shell center, not exactly at r=0, so this is the largest "
        "radius whose sampled core reaches the threshold; refine shells to tighten that gap.",
        "Michaelis-Menten uptake is at most vmax, so its critical radius is at least the zero-order value; "
        "the closed-form zero-order radius is reported beside it as an independent check.",
        "The threshold is a user-selected reporting level. It is not a hypoxia, death, viability or potency "
        "cutoff, and these parameters are illustrative rather than measured.",
    ])
    if infeasible_radii:
        result["interpretation"].append(
            "Some searched radii were infeasible for this configuration: the solver fails closed when it "
            "would predict negative oxygen or miss its residual bounds, which zero-order uptake does beyond "
            "a finite radius. Those radii were read as below threshold and are listed with their reasons.")
    return result


def solve_transient(
    parameters=Parameters(),
    total_time_s=1200.0,
    time_steps=120,
    initial_oxygen_mol_m3=0.0,
    *,
    tolerance=1e-8,
    max_iterations=50,
):
    """Solve the time-dependent PDE in (r, t) from initial conditions to steady state.

    PDE: dc/dt = (D/r^2) d/dr(r^2 dc/dr) - R(c)
    Boundary conditions:
      r=0: dc/dr = 0
      r=R: -D dc/dr = k_transfer (c - c_bulk) [or c=c_bulk if transfer_m_s is None]
    Initial condition:
      c(r, 0) = initial_oxygen_mol_m3 (default 0.0, representing anoxia ramp-up)
    """
    p = parameters.validate()
    if isinstance(total_time_s, bool) or not isinstance(total_time_s, (int, float)) or not np.isfinite(total_time_s) or total_time_s <= 0:
        raise ValueError("total_time_s must be a positive finite number")
    if isinstance(time_steps, bool) or not isinstance(time_steps, int) or time_steps < 2 or time_steps > 50000:
        raise ValueError("time_steps must be an integer between 2 and 50000")
    if isinstance(initial_oxygen_mol_m3, bool) or not isinstance(initial_oxygen_mol_m3, (int, float)) or not np.isfinite(initial_oxygen_mol_m3) or initial_oxygen_mol_m3 < 0:
        raise ValueError("initial_oxygen_mol_m3 must be a nonnegative finite number")
    if initial_oxygen_mol_m3 > p.bulk_oxygen_mol_m3:
        raise ValueError("initial oxygen cannot exceed the fixed bath concentration in this model")
    if isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)) or not np.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError("tolerance must be positive, finite, and less than 1")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")

    steady = solve(p)
    x, volumes, conductance, banded, boundary, boundary_g, da = _system(p)
    km = p.km_mol_m3 / p.bulk_oxygen_mol_m3
    radius = p.radius_um * 1e-6
    t_scale = radius ** 2 / p.diffusivity_m2_s

    time_s = np.linspace(0.0, float(total_time_s), time_steps)
    tau = time_s / t_scale

    u_history = np.zeros((time_steps, p.shells))
    u_init = float(initial_oxygen_mol_m3) / p.bulk_oxygen_mol_m3
    u_history[0] = u_init

    u_cur = np.full(p.shells, u_init)
    step_residuals = []
    max_step_iterations = 0
    max_mass_balance_error = 0.0
    for step in range(1, time_steps):
        dtau = float(tau[step] - tau[step - 1])
        if dtau <= 0:
            continue
        u_old = np.copy(u_cur)
        converged = False
        for it in range(1, max_iterations + 1):
            flux = conductance * (u_cur[1:] - u_cur[:-1])
            diff = np.zeros_like(u_cur)
            diff[:-1] -= flux
            diff[1:] += flux
            diff[-1] += boundary_g * (u_cur[-1] - 1.0)

            rate = np.ones_like(u_cur) if p.kinetics == "zero_order" else u_cur / (km + u_cur)
            f = volumes * (u_cur - u_old) / dtau + diff + da * volumes * rate
            scale = max(1.0, da, 1.0 / dtau)
            residual_norm = float(np.max(np.abs(f) / volumes) / scale)

            jac = banded.copy()
            jac[1] += volumes / dtau
            if p.kinetics != "zero_order":
                jac[1] += da * volumes * km / (km + u_cur) ** 2

            step_delta = solve_banded((1, 1), jac, -f)
            if not np.isfinite(step_delta).all():
                raise RuntimeError(f"transient solver produced a non-finite step at time index {step}")
            alpha = 1.0
            if np.any(step_delta < 0):
                alpha = min(alpha, 0.99 * float(np.min(-u_cur[step_delta < 0] / step_delta[step_delta < 0])))
            accepted = False
            for _ in range(50):
                candidate = u_cur + alpha * step_delta
                if np.all(candidate >= 0) and np.all(candidate <= 1) and np.isfinite(candidate).all():
                    candidate_flux = conductance * (candidate[1:] - candidate[:-1])
                    candidate_diff = np.zeros_like(candidate)
                    candidate_diff[:-1] -= candidate_flux
                    candidate_diff[1:] += candidate_flux
                    candidate_diff[-1] += boundary_g * (candidate[-1] - 1.0)
                    candidate_rate = np.ones_like(candidate) if p.kinetics == "zero_order" else candidate / (km + candidate)
                    candidate_f = volumes * (candidate - u_old) / dtau + candidate_diff + da * volumes * candidate_rate
                    candidate_norm = float(np.max(np.abs(candidate_f) / volumes) / scale)
                    if candidate_norm < residual_norm or candidate_norm <= tolerance:
                        accepted = True
                        break
                alpha *= 0.5
            if not accepted:
                raise RuntimeError(f"transient solver line search failed at time index {step}")
            u_cur = candidate
            if candidate_norm <= tolerance and np.max(np.abs(alpha * step_delta)) <= max(tolerance, np.sqrt(tolerance)):
                residual_norm = candidate_norm
                converged = True
                max_step_iterations = max(max_step_iterations, it)
                break
        if not converged:
            raise RuntimeError(f"transient solver did not converge at time index {step} after {max_iterations} iterations")
        step_residuals.append(residual_norm)
        old_mass = float(np.sum(volumes * u_old))
        new_mass = float(np.sum(volumes * u_cur))
        surface_supply = dtau * boundary_g * (1.0 - u_cur[-1])
        total_uptake = dtau * da * float(np.sum(volumes * (np.ones_like(u_cur) if p.kinetics == "zero_order" else u_cur / (km + u_cur))))
        balance_error = abs((new_mass - old_mass) - (surface_supply - total_uptake))
        balance_scale = max(abs(new_mass - old_mass), abs(surface_supply), abs(total_uptake), 1e-30)
        max_mass_balance_error = max(max_mass_balance_error, balance_error / balance_scale)
        u_history[step] = u_cur

    concentrations = u_history * p.bulk_oxygen_mol_m3
    core_history = concentrations[:, 0]

    influx_history = np.zeros(time_steps)
    surface_history = np.zeros(time_steps)
    uptake_history = np.zeros(time_steps)
    for step in range(time_steps):
        u_step = u_history[step]
        influx_step = float(4 * np.pi * radius * p.diffusivity_m2_s * p.bulk_oxygen_mol_m3 * boundary_g * (1.0 - u_step[-1]))
        influx_history[step] = influx_step
        if p.transfer_m_s is None:
            surface_history[step] = p.bulk_oxygen_mol_m3
        else:
            surface_history[step] = p.bulk_oxygen_mol_m3 - influx_step / (4 * np.pi * radius ** 2 * p.transfer_m_s)
        rate_step = np.ones_like(u_step) if p.kinetics == "zero_order" else u_step / (km + u_step)
        uptake_history[step] = float(4 * np.pi * radius ** 3 * p.vmax_mol_m3_s * np.sum(volumes * rate_step))

    volume_mean_history = np.sum(3 * volumes[None, :] * concentrations, axis=1)

    steady_core = float(steady.summary["minimum_sampled_oxygen_mol_m3"])
    t_half = None
    t_95 = None
    if steady_core > 1e-9:
        initial_core = float(core_history[0])

        def first_threshold_crossing(target):
            if math.isclose(initial_core, target, rel_tol=1e-12, abs_tol=1e-12):
                return 0.0
            crossed = core_history >= target if initial_core < target else core_history <= target
            indices = np.flatnonzero(crossed)
            return float(time_s[indices[0]]) if indices.size else None

        t_half = first_threshold_crossing(0.5 * steady_core)
        t_95 = first_threshold_crossing(0.95 * steady_core)

    initial_core = float(core_history[0])
    core_change = steady_core - initial_core

    def time_to_fraction_of_core_change(fraction):
        # Unlike fractions of absolute steady oxygen, these targets also describe
        # relaxation from above steady state. Interpolation resolves only the
        # sampled time grid, not the PDE's discretization error.
        if core_change == 0:
            return 0.0
        target = initial_core + fraction * core_change
        crossed = core_history >= target if core_change > 0 else core_history <= target
        indices = np.flatnonzero(crossed)
        if not indices.size:
            return None
        index = int(indices[0])
        if index == 0:
            return 0.0
        left, right = float(core_history[index - 1]), float(core_history[index])
        weight = (target - left) / (right - left)
        return float(time_s[index - 1] + weight * (time_s[index] - time_s[index - 1]))

    summary = {
        "parameters": asdict(p),
        "total_time_s": float(total_time_s),
        "time_steps": int(time_steps),
        "initial_oxygen_mol_m3": float(initial_oxygen_mol_m3),
        "diffusion_time_scale_s": float(t_scale),
        "initial_core_oxygen_mol_m3": float(core_history[0]),
        "final_core_oxygen_mol_m3": float(core_history[-1]),
        "steady_core_oxygen_mol_m3": steady_core,
        "final_vs_steady_core_error_mol_m3": float(abs(core_history[-1] - steady_core)),
        "max_step_iterations": int(max_step_iterations),
        "max_scaled_step_residual": float(max(step_residuals, default=0.0)),
        "max_relative_transient_mass_balance_error": float(max_mass_balance_error),
        "time_to_half_steady_core_s": t_half,
        "time_to_95pct_steady_core_s": t_95,
        "time_to_50pct_initial_to_steady_core_change_s": time_to_fraction_of_core_change(0.5),
        "time_to_95pct_initial_to_steady_core_change_s": time_to_fraction_of_core_change(0.95),
        "core_change_targets_mol_m3": {
            "50pct": initial_core + 0.5 * core_change,
            "95pct": initial_core + 0.95 * core_change,
        },
        "final_surface_oxygen_mol_m3": float(surface_history[-1]),
        "final_volume_mean_oxygen_mol_m3": float(volume_mean_history[-1]),
        "final_fraction_below_threshold": float(np.sum(3 * volumes[concentrations[-1] < p.threshold_mol_m3])),
        "final_total_uptake_mol_s": float(uptake_history[-1]),
        "final_surface_influx_mol_s": float(influx_history[-1]),
        "interpretation": [
            "Transient PDE mode solving del(c)/del(t) = D nabla^2(c) - R(c) from initial concentration.",
            "Illustrative numerical model; parameters need experimental calibration.",
            "Core oxygen is tracked at innermost shell center, surface at r=R with Robin/Dirichlet boundary.",
            "Absolute steady-core threshold times are first sampled crossings of 0.5 or 0.95 times the steady concentration; they may be unreachable from above.",
            "Initial-to-steady core-change times use targets c_initial + fraction*(c_steady-c_initial), interpolated between the first bracketing samples in the direction of that change; null means not reached during the simulated interval.",
            "An initial core already equal to steady core has zero core-change time; this does not establish that the whole initial radial profile is steady.",
            "Accepted steps meet the scaled backward-Euler residual; the reported conservation error compares accumulation with surface influx minus uptake. Refine both time and radial grids before interpreting exposure times.",
            "Spherical, homogeneous organoid model with fixed bath oxygen; no vascularization or cell death."
        ]
    }

    return TransientSolution(
        time_s=time_s,
        radius_um=x * p.radius_um,
        concentration_history_mol_m3=concentrations,
        core_oxygen_history_mol_m3=core_history,
        surface_oxygen_history_mol_m3=surface_history,
        volume_mean_history_mol_m3=volume_mean_history,
        shell_volume_fraction=3 * volumes,
        steady_solution=steady,
        summary=summary,
    )
