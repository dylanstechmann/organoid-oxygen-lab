import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from oxygenlab.cli import (
    constructed_inverse_demo,
    demo,
    local_sensitivity,
    main,
    save,
    save_transient,
    sweep_km,
    sweep_transfer,
    write_critical_radius,
)
from oxygenlab.model import Parameters, solve, solve_transient


class ReportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.output = self.root / "report"
        self.solution = solve(Parameters(radius_um=100))

    def assert_no_report_or_staging(self):
        self.assertEqual(list(self.root.iterdir()), [])

    def test_failed_serialization_leaves_no_output(self):
        self.solution.summary["surface_oxygen_mol_m3"] = float("nan")
        with self.assertRaises(ValueError):
            save(self.solution, self.output)
        self.assert_no_report_or_staging()

    def test_failed_last_write_rolls_back_and_allows_retry(self):
        write_text = Path.write_text

        def fail_report(path, *args, **kwargs):
            if path.name == "REPORT.md":
                self.assertFalse(self.output.exists())
                self.assertTrue((path.parent / "summary.json").is_file())
                self.assertTrue((path.parent / "profile.csv").is_file())
                raise OSError("injected report write failure")
            return write_text(path, *args, **kwargs)

        with patch.object(Path, "write_text", fail_report):
            with self.assertRaisesRegex(OSError, "injected"):
                save(self.solution, self.output)
        self.assert_no_report_or_staging()
        report = save(self.solution, self.output)
        self.assertEqual(json.loads((self.output / "summary.json").read_text(encoding="utf-8")), report)
        self.assertIn("µm", (self.output / "REPORT.md").read_text(encoding="utf-8"))

    def test_failed_publication_removes_already_published_entries(self):
        rename = Path.rename
        calls = []

        def fail_second_move(path, target):
            calls.append(path.name)
            if len(calls) == 2:
                self.assertEqual(len(list(self.output.iterdir())), 1)
                raise OSError("injected publication failure")
            return rename(path, target)

        with patch.object(Path, "rename", fail_second_move):
            with self.assertRaisesRegex(OSError, "injected"):
                save(self.solution, self.output)
        self.assertEqual(len(calls), 2)
        self.assert_no_report_or_staging()

    def test_existing_output_files_and_directories_are_preserved(self):
        for directory in [False, True]:
            with self.subTest(directory=directory):
                output = self.root / str(directory)
                if directory:
                    output.mkdir()
                    sentinel = output / "user.txt"
                else:
                    sentinel = output
                sentinel.write_bytes(b"existing user output")
                with self.assertRaises(FileExistsError):
                    save(self.solution, output)
                with self.assertRaises(FileExistsError):
                    demo(output)
                self.assertEqual(sentinel.read_bytes(), b"existing user output")
        self.assertEqual({entry.name for entry in self.root.iterdir()}, {"False", "True"})

    def test_destination_created_during_generation_is_not_reused(self):
        write_text = Path.write_text

        def create_destination(path, *args, **kwargs):
            result = write_text(path, *args, **kwargs)
            if path.name == "REPORT.md":
                self.output.mkdir()
            return result

        with patch.object(Path, "write_text", create_destination):
            with self.assertRaises(FileExistsError):
                save(self.solution, self.output)
        self.assertEqual(list(self.root.iterdir()), [self.output])
        self.assertEqual(list(self.output.iterdir()), [])

    def test_demo_solver_failure_removes_nested_report(self):
        def fail_after_single_sphere(parameters):
            if parameters.radius_um == 50:
                self.assertFalse(self.output.exists())
                raise RuntimeError("injected solver failure")
            return solve(parameters)

        with patch("oxygenlab.cli.solve", side_effect=fail_after_single_sphere):
            with self.assertRaisesRegex(RuntimeError, "injected"):
                demo(self.output)
        self.assert_no_report_or_staging()

    def test_demo_plot_failure_removes_report_and_closes_figure(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("plot extra is not installed")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure

        before = plt.get_fignums()
        with patch.object(Figure, "savefig", side_effect=OSError("injected plot failure")):
            with self.assertRaisesRegex(OSError, "injected"):
                demo(self.output, plot=True)
        self.assertEqual(plt.get_fignums(), before)
        self.assert_no_report_or_staging()

    def test_cli_write_failure_returns_error_without_success_output(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(Path, "write_text", side_effect=OSError("injected write failure")):
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                with self.assertRaises(SystemExit) as caught:
                    main(["demo", "--out", str(self.output)])
        self.assertEqual(caught.exception.code, 2)
        self.assertEqual(stdout.getvalue(), "")
        self.assertIn("injected write failure", stderr.getvalue())
        self.assert_no_report_or_staging()

    def run_cli(self, *arguments):
        # Fail on implicit text encodings, even on machines whose normal locale
        # is UTF-8. Disable UTF-8 coercion to exercise ASCII/legacy locales too.
        environment = {**os.environ, "LC_ALL": "C", "PYTHONUTF8": "0", "PYTHONCOERCECLOCALE": "0"}
        return subprocess.run(
            [sys.executable, "-X", "warn_default_encoding", "-W", "error::EncodingWarning",
             "-m", "oxygenlab.cli", *map(str, arguments)],
            env=environment, capture_output=True, encoding="utf-8", timeout=60,
        )

    def test_cli_demo_and_solve_generate_complete_utf8_reports(self):
        completed = self.run_cli("demo", "--out", self.output)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        validation = json.loads(completed.stdout)
        self.assertEqual(json.loads((self.output / "validation.json").read_text(encoding="utf-8")), validation)
        self.assertIn("10⁻⁵", (self.output / "REPORT.md").read_text(encoding="utf-8"))
        self.assertTrue((self.output / "single_sphere/profile.csv").is_file())

        config = self.output / "parameters.json"
        solved = self.root / "solved"
        completed = self.run_cli("solve", config, "--out", solved)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads(completed.stdout)
        self.assertEqual(report["input_sha256"], hashlib.sha256(config.read_bytes()).hexdigest())
        self.assertEqual({path.name for path in solved.iterdir()}, {"summary.json", "profile.csv", "REPORT.md"})
        self.assertEqual(json.loads((solved / "summary.json").read_text(encoding="utf-8")), report)
        self.assertEqual({path.name for path in self.root.iterdir()}, {"report", "solved"})

    def test_cli_existing_output_and_invalid_config_fail_cleanly(self):
        self.output.mkdir()
        sentinel = self.output / "user.txt"
        sentinel.write_bytes(b"keep me")
        completed = self.run_cli("demo", "--out", self.output)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("already exists", completed.stderr)
        self.assertEqual(completed.stdout, "")
        self.assertEqual(sentinel.read_bytes(), b"keep me")

        config = self.root / "invalid.json"
        config.write_text('{"radius_um": -1}', encoding="utf-8")
        completed = self.run_cli("solve", config, "--out", self.root / "invalid-report")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("radius_um", completed.stderr)
        self.assertEqual(completed.stdout, "")
        self.assertEqual({entry.name for entry in self.root.iterdir()}, {"report", "invalid.json"})

    def test_km_sweep_moves_the_critical_radius_the_documented_way(self):
        result = sweep_km(self.output)
        rows = result["rows"]
        self.assertTrue(all(row["status"] == "bracketed" for row in rows))
        radii = [row["critical_radius_um"] for row in rows]
        self.assertEqual(radii, sorted(radii))
        self.assertGreater(radii[-1], radii[0])
        # Michaelis-Menten uptake never exceeds vmax, so the zero-order radius is a lower bound.
        for row in rows:
            self.assertGreaterEqual(row["critical_radius_um"] + 0.05, row["zero_order_analytic_radius_um"])
        self.assertTrue(result["nondecreasing_in_km"])
        report = (self.output / "REPORT.md").read_text(encoding="utf-8")
        self.assertIn("not measured", report)
        self.assertTrue((self.output / "km_critical_radius.csv").is_file())
        with self.assertRaises(Exception):
            sweep_km(self.output)

    def test_transfer_sweep_report_records_units_and_preserves_existing_output(self):
        result = sweep_transfer(self.output)
        self.assertTrue(result["core_nondecreasing"])
        settings = json.loads((self.output / "settings.json").read_text(encoding="utf-8"))
        self.assertEqual(settings["sweep_parameter"], "transfer_m_s")
        self.assertEqual(settings["sweep_unit"], "m/s")
        self.assertIn("not measured", (self.output / "REPORT.md").read_text(encoding="utf-8"))
        self.assertTrue((self.output / "transfer_sweep.csv").is_file())
        with self.assertRaises(FileExistsError):
            sweep_transfer(self.output)

    def test_local_sensitivity_reports_parameter_tradeoffs_without_measurement_claims(self):
        report = local_sensitivity(self.output, Parameters(shells=40))
        self.assertEqual(report["analysis"], "local finite-difference log-parameter sensitivity")
        self.assertEqual(report["radial_measurement_count"], 40)
        self.assertEqual({item["parameter"] for item in report["parameter_sensitivities"]},
                         {"radius_um", "vmax_mol_m3_s", "km_mol_m3", "transfer_m_s"})
        self.assertTrue(all(np.isfinite(item["local_profile_sensitivity_norm_mol_m3"])
                            for item in report["parameter_sensitivities"]))
        self.assertIn("not a biological identifiability conclusion",
                      report["interpretation"])
        self.assertTrue((self.output / "sensitivity.json").is_file())
        self.assertIn("Illustrative numerical sensitivity",
                      (self.output / "REPORT.md").read_text(encoding="utf-8"))

    def test_constructed_inverse_demo_is_seeded_and_explicitly_synthetic(self):
        first = constructed_inverse_demo(self.output, seed=9)
        second_output = self.root / "inverse-repeat"
        second = constructed_inverse_demo(second_output, seed=9)
        self.assertEqual(first, second)
        self.assertEqual(first["data_status"], "constructed_synthetic_profile_only")
        self.assertEqual(first["n_constructed_observations"], 18)
        self.assertTrue(first["optimizer_success"])
        self.assertTrue(np.isfinite(first["weighted_residual_sum_squares"]))
        self.assertTrue((self.output / "inverse_demo.json").is_file())
        self.assertEqual(len((self.output / "constructed_profile.csv").read_text(encoding="utf-8").splitlines()), 19)
        self.assertIn("No biological measurements were fit",
                      (self.output / "REPORT.md").read_text(encoding="utf-8"))

    def test_config_hash_and_parameters_identify_one_snapshot_for_both_cli_modes(self):
        read_bytes = Path.read_bytes
        for command in ["solve", "transient"]:
            with self.subTest(command=command):
                config = self.root / f"{command}.json"
                config.write_bytes(b'{"radius_um":100,"shells":20}')
                snapshot = read_bytes(config)
                reads = []

                def replace_after_read(path):
                    data = read_bytes(path)
                    if path == config:
                        reads.append(path)
                        path.write_bytes(b'{"radius_um":200,"shells":20}')
                    return data

                output = self.root / command
                args = [command, str(config), "--out", str(output)]
                if command == "transient":
                    args += ["--total-time-s", "60", "--time-steps", "20"]
                with patch.object(Path, "read_bytes", replace_after_read), contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(main(args), 0)
                report = json.loads((output / "summary.json").read_text(encoding="utf-8"))
                self.assertEqual(reads, [config])
                self.assertEqual(report["input_sha256"], hashlib.sha256(snapshot).hexdigest())
                self.assertEqual(report["parameters"]["radius_um"], 100)

    def test_transient_plot_failures_do_not_publish_success_or_leak_figures(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("plot extra is not installed")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure

        solution = solve_transient(Parameters(radius_um=100, shells=20), total_time_s=60, time_steps=20)
        before = plt.get_fignums()
        for failing_method in ["tight_layout", "savefig"]:
            with self.subTest(method=failing_method):
                with patch.object(Figure, failing_method, side_effect=OSError("injected transient plot failure")):
                    with self.assertRaisesRegex(OSError, "injected transient plot failure"):
                        save_transient(solution, self.output, plot=True)
                self.assertEqual(plt.get_fignums(), before)
                self.assert_no_report_or_staging()
        save_transient(solution, self.output, plot=True)
        self.assertTrue((self.output / "transient_ramp.png").is_file())
        self.assertIn("initial-to-steady core change", (self.output / "REPORT.md").read_text(encoding="utf-8"))
        self.assertEqual(plt.get_fignums(), before)

    def test_cli_transient_report_and_preserves_existing_output(self):
        transient_out = self.root / "transient-out"
        completed = self.run_cli("transient", "--out", transient_out, "--total-time-s", "60", "--time-steps", "20")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue((transient_out / "summary.json").is_file())
        self.assertTrue((transient_out / "transient_profile.csv").is_file())
        self.assertTrue((transient_out / "radial_snapshots.csv").is_file())
        self.assertTrue((transient_out / "REPORT.md").is_file())

        summary = json.loads((transient_out / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["total_time_s"], 60.0)
        self.assertEqual(summary["time_steps"], 20)
        self.assertEqual(summary["initial_core_oxygen_mol_m3"], 0.0)

        # Output cannot be overwritten
        completed_again = self.run_cli("transient", "--out", transient_out)
        self.assertEqual(completed_again.returncode, 2)
        self.assertIn("already exists", completed_again.stderr)


if __name__ == "__main__":
    unittest.main()


class CriticalRadiusReportTests(unittest.TestCase):
    def test_report_publishes_both_estimates_and_refuses_reuse(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "critical"
            report = write_critical_radius(output, Parameters(shells=200), tolerance_um=0.05, shells=200)
            saved = json.loads((output / "critical_radius.json").read_text(encoding="utf-8"))
            self.assertEqual(saved["status"], "bracketed")
            self.assertEqual(saved["critical_radius_um"], report["critical_radius_um"])
            self.assertIsNotNone(saved["zero_order_analytic_radius_um"])
            self.assertIsNone(saved["input_sha256"])
            text = (output / "REPORT.md").read_text(encoding="utf-8")
            self.assertIn("Threshold-limited sphere size", text)
            self.assertIn("Closed-form zero-order radius", text)
            self.assertIn("not a hypoxia, death, viability or potency", text)
            with self.assertRaises(FileExistsError):
                write_critical_radius(output, Parameters(shells=200))

    def test_degenerate_report_states_the_reason_without_a_radius(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "no-uptake"
            report = write_critical_radius(output, Parameters(vmax_mol_m3_s=0.0, shells=64))
            self.assertEqual(report["status"], "unbounded_without_uptake")
            text = (output / "REPORT.md").read_text(encoding="utf-8")
            self.assertIn("unavailable", text)
            self.assertIn("no radius is limiting", text)
