"""Constructed schema fixtures test eligibility; they establish no biology."""

import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from oxygenlab.cli import main
from oxygenlab.measurement_intake import qualify_profiles


class MeasuredProfileTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(); self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.csv, self.manifest = self.root / "test.csv", self.root / "manifest.json"
        self.csv.write_bytes(b"specimen_id,radius_um,oxygen_mol_m3\na,0,0.1\na,50,0.15\na,100,0.2\nb,0,0.1\nb,50,0.15\nb,100,0.2\n")
        (self.root / "native.txt").write_bytes(b"constructed-test-export")
        common = {key: "constructed-test-only" for key in ("biological_identity_record", "sensor_calibration_record", "radius_measurement_record", "bath_measurement_record", "steady_state_record", "cell_type", "culture_conditions")}
        common.update(source_file="native", geometry="sphere", steady_state=True, complete_profile=True, temperature_degC=25.0, specimen_radius_um=100.0, bath_oxygen_mol_m3=0.2)
        self.document = {"schema_version": 1, "origin": "direct_calibrated_sensor", "source_url": "https://example.invalid/test-only", "license": "test fixture", "attribution": "constructed schema test", "acquisition_record": "test-only",
                         "units": {"radius": "um", "oxygen": "mol/m3", "temperature": "degC"}, "csv_sha256": hashlib.sha256(self.csv.read_bytes()).hexdigest(),
                         "source_files": {"native": {"path": "native.txt", "sha256": hashlib.sha256(b"constructed-test-export").hexdigest()}},
                         "specimens": {key: dict(common, original_specimen_id=key, biological_unit_id=key, source_trace_id=key, role=role) for key, role in (("a", "calibration"), ("b", "validation"))}}

    def qualify(self, document=None):
        self.manifest.write_text(json.dumps(document or self.document), encoding="utf-8")
        return qualify_profiles(self.csv, self.manifest)

    def test_qualified_counts_remain_intake_and_never_biological_validation(self):
        report = self.qualify()
        self.assertEqual(report["n_biological_units"], 2)
        self.assertEqual(report["n_specimens"], 2)
        self.assertFalse(report["model_fitted"])
        self.assertFalse(report["biological_validation_performed"])
        self.assertEqual(report["profiles"][0]["sampled_oxygen_range_mol_m3"], [0.1, 0.2])

    def test_complete_declaration_does_not_hide_a_missing_center_sample(self):
        raw = self.csv.read_bytes().replace(b"a,0,0.1", b"a,20,0.1")
        self.csv.write_bytes(raw)
        self.document["csv_sha256"] = hashlib.sha256(raw).hexdigest()
        report = self.qualify()
        profile = next(p for p in report["profiles"] if p["specimen_id"] == "a")
        self.assertFalse(profile["spatial_coverage"]["center_sample_present"])
        self.assertAlmostEqual(profile["spatial_coverage"]["inward_of_first_sample_sphere_volume_fraction"], .008)
        self.assertFalse(report["biological_validation_performed"])

    def test_uncertainty_is_optional_source_bound_and_not_propagated(self):
        report = self.qualify()
        self.assertIsNone(report["profiles"][0]["measurement_uncertainty"]["quantities"]["oxygen"]["value"])
        self.document["specimens"]["a"]["measurement_uncertainty"] = {"oxygen": {
            "kind": "absolute_bound", "value": .003, "unit": "mol/m3", "source_file": "native",
            "source_record": "fixture only", "scope": "Declared sensor bound in constructed fixture"}}
        report = self.qualify()
        uncertainty = report["profiles"][0]["measurement_uncertainty"]
        self.assertEqual(uncertainty["quantities"]["oxygen"]["value"], .003)
        self.assertFalse(uncertainty["propagated_to_model"])
        self.assertFalse(report["model_fitted"])

    def test_rejects_averaged_synthetic_digitized_or_model_inferred_oxygen(self):
        for origin in ("synthetic", "figure_digitized", "averaged", "model_inferred"):
            document = copy.deepcopy(self.document); document["origin"] = origin
            with self.assertRaisesRegex(ValueError, "direct_calibrated_sensor"):
                self.qualify(document)

    def test_rejects_same_biological_unit_in_training_and_validation(self):
        document = copy.deepcopy(self.document); document["specimens"]["b"]["biological_unit_id"] = "a"
        with self.assertRaisesRegex(ValueError, "both calibration and validation"):
            self.qualify(document)

    def test_rejects_original_specimen_split_into_multiple_ids(self):
        document = copy.deepcopy(self.document); document["specimens"]["b"]["original_specimen_id"] = "a"
        with self.assertRaisesRegex(ValueError, "original specimen"):
            self.qualify(document)

    def test_rejects_source_trace_aliases_as_replicates(self):
        document = copy.deepcopy(self.document)
        document["source_files"]["alias"] = copy.deepcopy(document["source_files"]["native"])
        document["specimens"]["b"].update(source_file="alias", source_trace_id="a")
        with self.assertRaisesRegex(ValueError, "source trace"):
            self.qualify(document)

    def test_rejects_no_heldout_units_incomplete_or_unsuitable_geometry(self):
        document = copy.deepcopy(self.document); document["specimens"]["b"]["role"] = "calibration"
        with self.assertRaisesRegex(ValueError, "held-out"):
            self.qualify(document)
        for key, value in (("geometry", "hemisphere"), ("steady_state", False), ("complete_profile", False)):
            document = copy.deepcopy(self.document); document["specimens"]["b"][key] = value
            with self.assertRaisesRegex(ValueError, "steady spherical"):
                self.qualify(document)

    def test_rejects_pressure_or_gas_percent_units_and_missing_calibration(self):
        document = copy.deepcopy(self.document); document["units"]["oxygen"] = "mmHg"
        with self.assertRaisesRegex(ValueError, "canonical units"):
            self.qualify(document)
        document = copy.deepcopy(self.document); document["specimens"]["b"]["sensor_calibration_record"] = ""
        with self.assertRaisesRegex(ValueError, "sensor_calibration_record"):
            self.qualify(document)

    def test_rejects_hash_changes_missing_specimens_and_duplicate_keys(self):
        document = copy.deepcopy(self.document); del document["specimens"]["b"]
        with self.assertRaisesRegex(ValueError, "exactly all"):
            self.qualify(document)
        document = copy.deepcopy(self.document); document["csv_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "CSV hash mismatch"):
            self.qualify(document)
        (self.root / "native.txt").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "source_file hash mismatch"):
            self.qualify()
        self.manifest.write_text('{"schema_version": 1, "schema_version": 1}')
        with self.assertRaisesRegex(ValueError, "duplicate"):
            qualify_profiles(self.csv, self.manifest)

    def test_rejects_nonfinite_context_and_radii_outside_measured_specimen(self):
        for key, value in (("specimen_radius_um", 99), ("specimen_radius_um", True), ("temperature_degC", float("nan")), ("bath_oxygen_mol_m3", 0)):
            document = copy.deepcopy(self.document); document["specimens"]["b"][key] = value
            with self.assertRaises(ValueError):
                self.qualify(document)

    def test_rejects_duplicate_positions_insufficient_points_and_nonfinite_oxygen(self):
        original = self.csv.read_bytes()
        for replacement in (b"b,50,0.15", b"b,100,nan", b"b,100,-0.1", b"b,100,"):
            self.csv.write_bytes(original.replace(b"b,100,0.2", replacement))
            document = copy.deepcopy(self.document); document["csv_sha256"] = hashlib.sha256(self.csv.read_bytes()).hexdigest()
            with self.assertRaises(ValueError):
                self.qualify(document)
        self.csv.write_bytes(original.replace(b"b,100,0.2\n", b""))
        document = copy.deepcopy(self.document); document["csv_sha256"] = hashlib.sha256(self.csv.read_bytes()).hexdigest()
        with self.assertRaisesRegex(ValueError, ">=3"):
            self.qualify(document)

    def test_cli_publishes_qualification_only_and_refuses_existing_output(self):
        self.qualify()
        output = self.root / "report"
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(["qualify-profiles", str(self.csv), str(self.manifest), "--out", str(output)]), 0)
        report = json.loads((output / "qualification.json").read_text())
        self.assertFalse(report["biological_validation_performed"])
        self.assertFalse((output / "parameters.json").exists())
        with self.assertRaises(SystemExit):
            main(["qualify-profiles", str(self.csv), str(self.manifest), "--out", str(output)])

    def test_failed_intake_leaves_no_output(self):
        document = copy.deepcopy(self.document); document["origin"] = "synthetic"
        self.manifest.write_text(json.dumps(document))
        output = self.root / "report"
        with self.assertRaises(SystemExit) as error:
            main(["qualify-profiles", str(self.csv), str(self.manifest), "--out", str(output)])
        self.assertEqual(error.exception.code, 2)
        self.assertFalse(output.exists())
