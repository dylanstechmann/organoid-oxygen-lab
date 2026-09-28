"""Spherical finite-volume diffusion with oxygen-dependent uptake.

Concentrations are mol/m^3 (= mM); lengths are converted to SI once.
The surrounding bath concentration is fixed. No flow field or cell death is
modeled. A finite surface mass-transfer coefficient represents a boundary
resistance, not a simulation of a pump or its perfusion rate.
"""

from dataclasses import asdict, dataclass
import json
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
        values = json.loads(Path(path).read_text(encoding="utf-8"))
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
    if p.vmax_mol_m3_s == 0:
        pass
    elif p.kinetics == "zero_order":
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
    if norm(u) > tolerance or not np.isfinite(u).all() or np.any(u > 1 + 1e-8):
        raise RuntimeError("solution failed residual or physical-bound checks")
    radius = p.radius_um * 1e-6
    rate = np.ones_like(u) if p.kinetics == "zero_order" else u / (km + u)
    uptake = float(4 * np.pi * radius ** 3 * p.vmax_mol_m3_s * np.sum(volumes * rate))
    influx = float(4 * np.pi * radius * p.diffusivity_m2_s * p.bulk_oxygen_mol_m3 * boundary_g * (1 - u[-1]))
    surface = p.bulk_oxygen_mol_m3 if p.transfer_m_s is None else p.bulk_oxygen_mol_m3 - influx / (4 * np.pi * radius ** 2 * p.transfer_m_s)
    concentrations = u * p.bulk_oxygen_mol_m3
    summary = {"parameters": asdict(p), "iterations": iterations, "scaled_residual": norm(u),
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

    steady = solve(p)
    x, volumes, conductance, banded, boundary, boundary_g, da = _system(p)
    km = p.km_mol_m3 / p.bulk_oxygen_mol_m3
    radius = p.radius_um * 1e-6
    t_scale = radius ** 2 / p.diffusivity_m2_s

    time_s = np.linspace(0.0, float(total_time_s), time_steps)
    tau = time_s / t_scale

    u_history = np.zeros((time_steps, p.shells))
    u_init = min(1.0, float(initial_oxygen_mol_m3) / p.bulk_oxygen_mol_m3)
    u_history[0] = u_init

    u_cur = np.full(p.shells, u_init)
    for step in range(1, time_steps):
        dtau = float(tau[step] - tau[step - 1])
        if dtau <= 0:
            continue
        u_old = np.copy(u_cur)
        for it in range(max_iterations):
            flux = conductance * (u_cur[1:] - u_cur[:-1])
            diff = np.zeros_like(u_cur)
            diff[:-1] -= flux
            diff[1:] += flux
            diff[-1] += boundary_g * (u_cur[-1] - 1.0)

            rate = np.ones_like(u_cur) if p.kinetics == "zero_order" else u_cur / (km + u_cur)
            f = volumes * (u_cur - u_old) / dtau + diff + da * volumes * rate

            jac = banded.copy()
            jac[1] += volumes / dtau
            if p.kinetics != "zero_order":
                jac[1] += da * volumes * km / (km + u_cur) ** 2

            step_delta = solve_banded((1, 1), jac, -f)
            u_cur += step_delta
            u_cur = np.maximum(u_cur, 0.0)
            if np.max(np.abs(step_delta)) < tolerance:
                break
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
        half_mask = core_history >= 0.5 * steady_core
        if np.any(half_mask):
            t_half = float(time_s[np.argmax(half_mask)])
        pct95_mask = core_history >= 0.95 * steady_core
        if np.any(pct95_mask):
            t_95 = float(time_s[np.argmax(pct95_mask)])

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
        "time_to_half_steady_core_s": t_half,
        "time_to_95pct_steady_core_s": t_95,
        "final_surface_oxygen_mol_m3": float(surface_history[-1]),
        "final_volume_mean_oxygen_mol_m3": float(volume_mean_history[-1]),
        "final_fraction_below_threshold": float(np.sum(3 * volumes[concentrations[-1] < p.threshold_mol_m3])),
        "final_total_uptake_mol_s": float(uptake_history[-1]),
        "final_surface_influx_mol_s": float(influx_history[-1]),
        "interpretation": [
            "Transient PDE mode solving del(c)/del(t) = D nabla^2(c) - R(c) from initial concentration.",
            "Illustrative numerical model; parameters need experimental calibration.",
            "Core oxygen is tracked at innermost shell center, surface at r=R with Robin/Dirichlet boundary.",
            "Time to 50% and 95% of steady-state core oxygen indicate penetration time from anoxia.",
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

