"""Configuration-driven numerical experiments and reproducible reports."""

import argparse
from contextlib import contextmanager
import csv
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import shutil
from tempfile import TemporaryDirectory

import numpy as np
import scipy

from oxygenlab import __version__
from oxygenlab.model import Parameters, solve, solve_transient, zero_order_analytic


@contextmanager
def _report_directory(output):
    """Prepare a complete report before publishing it to a new directory."""
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as temporary:
        staged = Path(temporary)
        yield staged
        # Reserve the destination exclusively. A directory rename can replace
        # another writer's empty directory on POSIX, so publish its prepared
        # entries only after mkdir has established that this destination is ours.
        output.mkdir()
        try:
            for entry in staged.iterdir():
                entry.rename(output / entry.name)
        except BaseException:
            shutil.rmtree(output)
            raise


def save(solution, output, *, input_sha256=None):
    with _report_directory(output) as staged:
        return _write_solution(solution, staged, input_sha256=input_sha256)


def _write_solution(solution, output, *, input_sha256=None):
    report = {**solution.summary, "schema_version": 1,
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "oxygenlab": __version__}}
    canonical = json.dumps(report["parameters"], sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    report["parameters_sha256"] = hashlib.sha256(canonical).hexdigest()
    report["input_sha256"] = input_sha256
    (output / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (output / "profile.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["radius_um", "oxygen_mol_m3", "shell_volume_fraction"])
        writer.writerows(zip(solution.radius_um, solution.concentration_mol_m3, solution.shell_volume_fraction))
    lines = ["# Spherical oxygen model", "", "Illustrative numerical model; parameters need experimental calibration.", "",
        f"Radius: {report['parameters']['radius_um']:g} µm. Concentrations are mol/m³ (= mM).", "",
        "| Quantity | Result |", "|---|---:|",
        f"| Minimum sampled oxygen (mol/m³) | {report['minimum_sampled_oxygen_mol_m3']:.6g} |",
        f"| Surface oxygen (mol/m³) | {report['surface_oxygen_mol_m3']:.6g} |",
        f"| Volume below user threshold | {report['fraction_volume_below_threshold']:.4f} |",
        f"| Integrated uptake (mol/s) | {report['total_uptake_mol_s']:.6g} |",
        f"| Relative mass-balance error | {report['relative_mass_balance_error']:.3g} |",
        f"| Scaled equation residual | {report['scaled_residual']:.3g} |", "",
        "## Interpretation", "", *[f"- {item}" for item in report["interpretation"]], ""]
    (output / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    return report


def demo(output, *, plot=False):
    with _report_directory(output) as staged:
        return _write_demo(staged, plot=plot)


def _write_demo(output, *, plot=False):
    parameters = Parameters()
    (output / "parameters.json").write_text(json.dumps(asdict(parameters), indent=2) + "\n", encoding="utf-8")
    save(solve(parameters), output / "single_sphere")
    rows, profiles = [], {}
    for condition, transfer in [("fixed_surface", None), ("finite_transfer", parameters.transfer_m_s)]:
        for radius in [50, 100, 200, 300, 400, 600, 800]:
            result = solve(replace(parameters, radius_um=radius, transfer_m_s=transfer))
            s = result.summary
            rows.append({"boundary": condition, "radius_um": radius,
                         "minimum_oxygen_mol_m3": s["minimum_sampled_oxygen_mol_m3"],
                         "surface_oxygen_mol_m3": s["surface_oxygen_mol_m3"],
                         "fraction_below_threshold": s["fraction_volume_below_threshold"],
                         "mass_balance_error": s["relative_mass_balance_error"]})
            if condition == "finite_transfer" and radius in [100, 300, 800]:
                profiles[radius] = result
    with (output / "radius_sweep.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    convergence = []
    for shells in [20, 40, 80, 160]:
        p = replace(parameters, radius_um=100, shells=shells, kinetics="zero_order")
        solution = solve(p)
        error = float(np.max(np.abs(solution.concentration_mol_m3 - zero_order_analytic(solution.radius_um, p))))
        convergence.append({"shells": shells, "max_abs_error_mol_m3": error})
    validation = {"kind": "numerical verification, not biological validation", "analytic_zero_order": convergence,
                  "max_sweep_mass_balance_error": max(row["mass_balance_error"] for row in rows),
                  "parameter_source": "illustrative values selected for software demonstration; not fitted or copied as a calibrated organoid parameter set"}
    (output / "validation.json").write_text(json.dumps(validation, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = ["# Oxygen model numerical demonstration", "", "All parameters are illustrative. No biological measurements are included.", "",
             "## Boundary resistance and size", "", "The bath is fixed at 0.2 mol/m³. The user-selected reporting threshold is 0.02 mol/m³.",
             "Finite-transfer scenarios use k = 2 × 10⁻⁵ m/s; fixed-surface scenarios have no exterior resistance.", "",
             "| Boundary | Radius (µm) | Minimum sampled oxygen (mol/m³) | Volume below threshold |",
             "|---|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['boundary']} | {row['radius_um']} | {row['minimum_oxygen_mol_m3']:.5g} | {row['fraction_below_threshold']:.3f} |")
    lines += ["", "## Analytic reference and mesh refinement", "",
              "The constant-uptake reference has a closed-form solution. Halving shell width should reduce",
              "the maximum concentration error by about a factor of four in this smooth case.", "",
              "| Shells | Maximum absolute error (mol/m³) |", "|---|---:|"]
    for item in convergence:
        lines.append(f"| {item['shells']} | {item['max_abs_error_mol_m3']:.5g} |")
    lines += ["", f"Maximum relative mass-balance error across the size scenarios: {validation['max_sweep_mass_balance_error']:.3g}.", "",
              "Unit tests additionally compare the nonlinear profile against an independent SciPy boundary-value solver.",
              "Equation residuals and numerical convergence do not establish accuracy for a particular organoid.", ""]
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
        for radius, solution in profiles.items():
            axes[0].plot(solution.radius_um / radius, solution.concentration_mol_m3, label=f"R = {radius} µm")
        axes[0].axhline(parameters.threshold_mol_m3, color="#666666", linestyle="--", linewidth=1, label="User reporting threshold")
        axes[0].set(xlabel="Relative radius r/R", ylabel="Oxygen (mol/m³ = mM)", title="Finite surface-transfer resistance")
        axes[0].legend(fontsize=8); axes[0].grid(alpha=0.2)
        for condition, label in [("fixed_surface", "Fixed surface concentration"), ("finite_transfer", "Finite surface transfer")]:
            subset = [row for row in rows if row["boundary"] == condition]
            axes[1].plot([row["radius_um"] for row in subset], [row["fraction_below_threshold"] for row in subset], "o-", label=label)
        axes[1].set(xlabel="Sphere radius (µm)", ylabel="Volume fraction below threshold", ylim=(-0.03, 1.03), title="Illustrative size scenarios")
        axes[1].legend(fontsize=8); axes[1].grid(alpha=0.2)
        fig.suptitle("Spherical oxygen transport • illustrative parameters, no biological calibration", fontsize=12)
        try:
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            fig.savefig(output / "oxygen_demo.png", dpi=160)
        finally:
            plt.close(fig)
        lines += ["![Oxygen profiles and size scenarios](oxygen_demo.png)", ""]
    (output / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    return validation


def vmax_sweep_rows(parameters=None, vmax_values=(0.0, 0.005, 0.01, 0.02, 0.04)):
    """Illustrative uptake sweep at fixed geometry. Not a fitted rate."""
    parameters = Parameters() if parameters is None else parameters
    rows = []
    previous = None
    for vmax in vmax_values:
        result = solve(replace(parameters, vmax_mol_m3_s=float(vmax)))
        summary = result.summary
        core = float(summary["minimum_sampled_oxygen_mol_m3"])
        rows.append({
            "radius_um": float(parameters.radius_um),
            "vmax_mol_m3_s": float(vmax),
            "minimum_oxygen_mol_m3": core,
            "fraction_below_threshold": float(summary["fraction_volume_below_threshold"]),
            "relative_mass_balance_error": float(summary["relative_mass_balance_error"]),
            "core_did_not_rise": previous is None or core <= previous + 1e-9,
        })
        previous = core
    return rows


def sweep_vmax(output, parameters=None, vmax_values=(0.0, 0.005, 0.01, 0.02, 0.04)):
    rows = vmax_sweep_rows(parameters, vmax_values)
    with _report_directory(output) as staged:
        with (staged / "vmax_sweep.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        lines = [
            "# Illustrative uptake sweep",
            "",
            "Question: at fixed radius, does the sampled core oxygen fall as `vmax_mol_m3_s` rises,",
            "and does the relative mass-balance error stay small?",
            "",
            "These `vmax` values are not measured uptake. A per-cell rate still needs a cell density",
            "before it can be entered as mol/(m³·s).",
            "",
            "| vmax (mol/m³/s) | Minimum sampled oxygen (mol/m³) | Volume below threshold | Mass-balance error |",
            "|---:|---:|---:|---:|",
        ]
        for row in rows:
            lines.append(
                f"| {row['vmax_mol_m3_s']:g} | {row['minimum_oxygen_mol_m3']:.6g} | "
                f"{row['fraction_below_threshold']:.4f} | {row['relative_mass_balance_error']:.3g} |"
            )
        monotone = all(row["core_did_not_rise"] for row in rows)
        worst = max(row["relative_mass_balance_error"] for row in rows)
        lines += [
            "",
            f"Core oxygen was nonincreasing across this grid: {monotone}.",
            f"Worst relative mass-balance error: {worst:.3g}.",
            "A monotone illustrative grid is not a calibration and not a hypoxia threshold.",
            "",
        ]
        (staged / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    return {"rows": rows, "core_nonincreasing": monotone, "worst_relative_mass_balance_error": worst}


def transfer_sweep_rows(parameters=None, transfer_values_m_s=(2e-6, 5e-6, 1e-5, 2e-5, 5e-5)):
    """Vary only the illustrative surface mass-transfer coefficient (m/s)."""
    parameters = (Parameters() if parameters is None else parameters).validate()
    values = tuple(transfer_values_m_s)
    if not values:
        raise ValueError("transfer_values_m_s must contain at least one value")
    for value in values:
        if value is None:
            raise ValueError("transfer_values_m_s must contain finite coefficients in m/s")
        replace(parameters, transfer_m_s=value).validate()
    if any(left >= right for left, right in zip(values, values[1:])):
        raise ValueError("transfer_values_m_s must be strictly increasing")

    rows = []
    previous = None
    for transfer in (*values, None):
        summary = solve(replace(parameters, transfer_m_s=transfer)).summary
        core = float(summary["minimum_sampled_oxygen_mol_m3"])
        rows.append({
            "boundary": "fixed_surface" if transfer is None else "finite_transfer",
            "radius_um": float(parameters.radius_um),
            "transfer_m_s": transfer,
            "minimum_oxygen_mol_m3": core,
            "surface_oxygen_mol_m3": float(summary["surface_oxygen_mol_m3"]),
            "fraction_below_threshold": float(summary["fraction_volume_below_threshold"]),
            "relative_mass_balance_error": float(summary["relative_mass_balance_error"]),
            "core_did_not_fall": previous is None or core >= previous - 1e-9,
        })
        previous = core
    return rows


def sweep_transfer(output, parameters=None, transfer_values_m_s=(2e-6, 5e-6, 1e-5, 2e-5, 5e-5), *, plot=False):
    """Report one boundary-resistance sensitivity question at fixed uptake and size."""
    parameters = Parameters() if parameters is None else parameters
    values = tuple(transfer_values_m_s)
    rows = transfer_sweep_rows(parameters, values)
    finite = rows[:-1]
    fixed = rows[-1]
    monotone = all(row["core_did_not_fall"] for row in rows)
    worst = max(row["relative_mass_balance_error"] for row in rows)
    with _report_directory(output) as staged:
        settings = {
            "base_parameters": asdict(parameters),
            "sweep_parameter": "transfer_m_s",
            "sweep_unit": "m/s",
            "transfer_values_m_s": list(values),
            "fixed_surface_reference": "transfer_m_s = null; no exterior resistance",
            "provenance": "illustrative values; no measured transfer coefficient or pump-flow conversion",
        }
        (staged / "settings.json").write_text(json.dumps(settings, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        with (staged / "transfer_sweep.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        lines = [
            "# Illustrative surface-transfer sweep", "",
            "Question: at fixed radius and uptake, does sampled core oxygen rise as",
            "the surface mass-transfer coefficient `transfer_m_s` rises?", "",
            "These coefficients are illustrative m/s values, not measured for an",
            "organoid. The fixed-surface row removes exterior resistance and is",
            "a limiting reference, not a finite transfer coefficient or pump flow.", "",
            "| Boundary | transfer (m/s) | Minimum sampled oxygen (mol/m³) | Surface oxygen (mol/m³) | Mass-balance error |",
            "|---|---:|---:|---:|---:|",
        ]
        for row in rows:
            transfer = "none" if row["transfer_m_s"] is None else f"{row['transfer_m_s']:.2g}"
            lines.append(
                f"| {row['boundary']} | {transfer} | {row['minimum_oxygen_mol_m3']:.6g} | "
                f"{row['surface_oxygen_mol_m3']:.6g} | {row['relative_mass_balance_error']:.3g} |"
            )
        lines += [
            "", f"Core oxygen was nondecreasing across the grid and fixed limit: {monotone}.",
            f"Worst relative mass-balance error: {worst:.3g}.",
            "The user-selected threshold is a descriptive concentration cut, not a viability limit.",
            "The sphere is homogeneous and steady; this grid is not a calibration.", "",
        ]
        if plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.semilogx([row["transfer_m_s"] for row in finite],
                        [row["minimum_oxygen_mol_m3"] for row in finite], "o-", label="Sampled core")
            ax.semilogx([row["transfer_m_s"] for row in finite],
                        [row["surface_oxygen_mol_m3"] for row in finite], "s-", label="Surface")
            ax.axhline(fixed["minimum_oxygen_mol_m3"], color="#555555", linestyle="--",
                       label="Fixed-surface core limit")
            ax.set(xlabel="Surface mass-transfer coefficient k (m/s)",
                   ylabel="Oxygen concentration (mol/m³)",
                   title="Illustrative boundary-resistance sensitivity")
            ax.grid(alpha=0.2)
            ax.legend()
            fig.tight_layout()
            try:
                fig.savefig(staged / "transfer_sweep.png", dpi=160)
            finally:
                plt.close(fig)
            lines += ["![Illustrative surface-transfer sensitivity](transfer_sweep.png)", ""]
        (staged / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    return {"rows": rows, "core_nondecreasing": monotone, "worst_relative_mass_balance_error": worst}


def save_transient(transient_solution, output, *, input_sha256=None, plot=False):
    with _report_directory(output) as staged:
        return _write_transient(transient_solution, staged, input_sha256=input_sha256, plot=plot)


def _write_transient(sol, output, *, input_sha256=None, plot=False):
    report = {
        **sol.summary,
        "schema_version": 1,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "oxygenlab": __version__,
        },
    }
    canonical = json.dumps(report["parameters"], sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    report["parameters_sha256"] = hashlib.sha256(canonical).hexdigest()
    report["input_sha256"] = input_sha256
    (output / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    with (output / "transient_profile.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["time_s", "core_oxygen_mol_m3", "surface_oxygen_mol_m3", "volume_mean_oxygen_mol_m3"])
        for t, c, s, vm in zip(sol.time_s, sol.core_oxygen_history_mol_m3, sol.surface_oxygen_history_mol_m3, sol.volume_mean_history_mol_m3):
            writer.writerow([f"{t:.4f}", f"{c:.6g}", f"{s:.6g}", f"{vm:.6g}"])

    n_times = len(sol.time_s)
    snapshot_indices = sorted(list(set([0, n_times // 4, n_times // 2, 3 * n_times // 4, n_times - 1])))
    snap_headers = ["radius_um"] + [f"t_{sol.time_s[idx]:.1f}s" for idx in snapshot_indices] + ["steady_state"]
    with (output / "radial_snapshots.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(snap_headers)
        for r_idx, r in enumerate(sol.radius_um):
            row = [f"{r:.2f}"]
            for t_idx in snapshot_indices:
                row.append(f"{sol.concentration_history_mol_m3[t_idx, r_idx]:.6g}")
            row.append(f"{sol.steady_solution.concentration_mol_m3[r_idx]:.6g}")
            writer.writerow(row)

    p = report["parameters"]
    t_half_str = f"{report['time_to_half_steady_core_s']:.2f} s" if report["time_to_half_steady_core_s"] is not None else "N/A"
    t_95_str = f"{report['time_to_95pct_steady_core_s']:.2f} s" if report["time_to_95pct_steady_core_s"] is not None else "N/A"
    change_times = [report[f"time_to_{fraction}pct_initial_to_steady_core_change_s"] for fraction in [50, 95]]
    change_strings = [f"{value:.2f} s" if value is not None else "not reached" for value in change_times]
    lines = [
        "# Spherical Transient Oxygen Diffusion Model", "",
        "Illustrative time-dependent numerical simulation from the recorded uniform initial oxygen.", "",
        f"Radius: {p['radius_um']:g} µm. Total simulated time: {report['total_time_s']:g} s ({report['time_steps']} saved time points).",
        f"Initial oxygen: {report['initial_oxygen_mol_m3']:g} mol/m³. Bulk oxygen: {p['bulk_oxygen_mol_m3']:g} mol/m³.", "",
        "| Quantity | Value |", "|---|---:|",
        f"| Diffusion time scale R²/D (s) | {report['diffusion_time_scale_s']:.2f} |",
        f"| Initial core oxygen (mol/m³) | {report['initial_core_oxygen_mol_m3']:.6g} |",
        f"| Final core oxygen at {report['total_time_s']}s (mol/m³) | {report['final_core_oxygen_mol_m3']:.6g} |",
        f"| Steady-state core oxygen (mol/m³) | {report['steady_core_oxygen_mol_m3']:.6g} |",
        f"| Difference vs steady state (mol/m³) | {report['final_vs_steady_core_error_mol_m3']:.4e} |",
        f"| Time to 50% steady-state core | {t_half_str} |",
        f"| Time to 95% steady-state core | {t_95_str} |",
        f"| Time to 50% of initial-to-steady core change | {change_strings[0]} |",
        f"| Time to 95% of initial-to-steady core change | {change_strings[1]} |",
        f"| Maximum scaled time-step residual | {report['max_scaled_step_residual']:.3g} |",
        f"| Maximum relative dynamic mass-balance error | {report['max_relative_transient_mass_balance_error']:.3g} |",
        f"| Final surface oxygen (mol/m³) | {report['final_surface_oxygen_mol_m3']:.6g} |",
        f"| Final volume-mean oxygen (mol/m³) | {report['final_volume_mean_oxygen_mol_m3']:.6g} |",
        f"| Final below-threshold fraction | {report['final_fraction_below_threshold']:.4f} |", "",
        "## Interpretation", "", *[f"- {item}" for item in report["interpretation"]], "",
    ]

    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        try:
            ax1.plot(sol.time_s, sol.core_oxygen_history_mol_m3, label="Sampled Core", color="#e11d48", lw=2)
            ax1.plot(sol.time_s, sol.volume_mean_history_mol_m3, label="Volume Mean", color="#2563eb", lw=1.5, ls="--")
            ax1.plot(sol.time_s, sol.surface_oxygen_history_mol_m3, label="Surface", color="#059669", lw=1.5)
            ax1.axhline(report["steady_core_oxygen_mol_m3"], color="#64748b", ls=":", label="Steady Core")
            ax1.set(xlabel="Time (s)", ylabel="Oxygen (mol/m³)", title="Illustrative Oxygen Transient")
            ax1.grid(alpha=0.3)
            ax1.legend()

            for idx in snapshot_indices:
                ax2.plot(sol.radius_um, sol.concentration_history_mol_m3[idx], label=f"t = {sol.time_s[idx]:.0f} s")
            ax2.plot(sol.radius_um, sol.steady_solution.concentration_mol_m3, "k--", label="Steady State", lw=1.8)
            ax2.set(xlabel="Radius (µm)", ylabel="Oxygen (mol/m³)", title="Radial Profiles over Time")
            ax2.grid(alpha=0.3)
            ax2.legend()

            fig.tight_layout()
            fig.savefig(output / "transient_ramp.png", dpi=160)
        finally:
            plt.close(fig)
        lines += ["![Illustrative transient oxygen ramp](transient_ramp.png)", ""]

    (output / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description="Conservative spherical oxygen transport model")
    commands = parser.add_subparsers(dest="command", required=True)
    d = commands.add_parser("demo", help="illustrative size sweep and analytic convergence checks")
    d.add_argument("--out", required=True); d.add_argument("--plot", action="store_true")
    s = commands.add_parser("solve", help="solve one JSON parameter configuration")
    s.add_argument("config"); s.add_argument("--out", required=True)
    sweep = commands.add_parser("sweep-vmax", help="illustrative uptake sweep at the default radius")
    sweep.add_argument("--out", required=True)
    transfer = commands.add_parser("sweep-transfer", help="illustrative surface-transfer sweep (m/s)")
    transfer.add_argument("--out", required=True)
    transfer.add_argument("--plot", action="store_true")
    trans = commands.add_parser("transient", help="transient PDE simulation from uniform initial oxygen")
    trans.add_argument("config", nargs="?", default=None, help="optional path to JSON parameter configuration")
    trans.add_argument("--out", required=True, help="output directory for transient report")
    trans.add_argument("--total-time-s", type=float, default=1200.0, help="total simulation time in seconds (default: 1200)")
    trans.add_argument("--time-steps", type=int, default=120, help="number of saved time points including t=0 (default: 120)")
    trans.add_argument("--initial-oxygen", type=float, default=0.0, help="initial uniform oxygen in mol/m3 (default: 0.0, anoxic)")
    trans.add_argument("--plot", action="store_true", help="generate transient ramp plot")
    args = parser.parse_args(argv)
    try:
        if args.command == "demo":
            report = demo(args.out, plot=args.plot)
        elif args.command == "sweep-vmax":
            report = sweep_vmax(args.out)
        elif args.command == "sweep-transfer":
            report = sweep_transfer(args.out, plot=args.plot)
        elif args.command == "transient":
            if args.config:
                config_bytes = Path(args.config).read_bytes()
                p = Parameters.from_json_bytes(config_bytes)
                sha = hashlib.sha256(config_bytes).hexdigest()
            else:
                p = Parameters()
                sha = None
            sol = solve_transient(p, total_time_s=args.total_time_s, time_steps=args.time_steps,
                                  initial_oxygen_mol_m3=args.initial_oxygen)
            report = save_transient(sol, args.out, input_sha256=sha, plot=args.plot)
        else:
            config_bytes = Path(args.config).read_bytes()
            report = save(solve(Parameters.from_json_bytes(config_bytes)), args.out,
                          input_sha256=hashlib.sha256(config_bytes).hexdigest())
        print(json.dumps(report, indent=2, allow_nan=False))
    except (ValueError, OSError, RuntimeError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

