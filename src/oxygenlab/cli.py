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
from oxygenlab.model import Parameters, solve, zero_order_analytic


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


def main(argv=None):
    parser = argparse.ArgumentParser(description="Conservative spherical oxygen transport model")
    commands = parser.add_subparsers(dest="command", required=True)
    d = commands.add_parser("demo", help="illustrative size sweep and analytic convergence checks")
    d.add_argument("--out", required=True); d.add_argument("--plot", action="store_true")
    s = commands.add_parser("solve", help="solve one JSON parameter configuration")
    s.add_argument("config"); s.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "demo":
            report = demo(args.out, plot=args.plot)
        else:
            config = Path(args.config)
            report = save(solve(Parameters.from_json(config)), args.out, input_sha256=hashlib.sha256(config.read_bytes()).hexdigest())
        print(json.dumps(report, indent=2, allow_nan=False))
    except (ValueError, OSError, RuntimeError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
