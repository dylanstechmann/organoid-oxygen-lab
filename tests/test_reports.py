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

from oxygenlab.cli import demo, main, save
from oxygenlab.model import Parameters, solve


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


if __name__ == "__main__":
    unittest.main()
