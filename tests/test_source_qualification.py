"""Constructed catalog fixtures test decisions only; they add no source evidence."""

import contextlib
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from oxygenlab.cli import main, source_qualification
from oxygenlab.source_qualification import qualify_source_catalog


class SourceQualificationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.frozen_catalog = Path(__file__).parents[1] / "docs/public-source-qualification.json"
        self.frozen_bytes = self.frozen_catalog.read_bytes()

    def write_catalog(self, document):
        path = self.root / "catalog.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        return path

    def complete_constructed_document(self):
        return {
            "schema_version": 1,
            "kind": "public-source-qualification-not-a-measurement-dataset",
            "audit_date": "constructed-test-only",
            "n_eligible_raw_datasets_acquired": 1,
            "decision": "Constructed schema fixture; not a source decision.",
            "method": "Constructed schema fixture only.",
            "sources": [{
                "id": "constructed-test-only",
                "url": "https://example.invalid/constructed-test-only",
                "doi": "constructed-test-only",
                "license_observation": "Constructed reuse-terms field only.",
                "individual_observations": [{"observation_id": "constructed-observation-id-only"}],
                "calibration_record": "Constructed calibration field only.",
                "independent_unit_ids": {"calibration": ["unit-a"], "validation": ["unit-b"]},
                "eligibility": "eligible_as_acquired",
                "reasons": [],
                "raw_measurement_sha256": None,
            }],
        }

    def test_frozen_catalog_completes_with_two_exclusions_and_zero_eligible_cohorts(self):
        document = json.loads(self.frozen_bytes)
        report = qualify_source_catalog(self.frozen_catalog)
        self.assertEqual(report["status"], "completed_zero_eligible_raw_cohorts")
        self.assertEqual(report["n_candidates"], 2)
        self.assertEqual(report["n_eligible_raw_cohorts"], 0)
        self.assertEqual(report["catalog_sha256"], hashlib.sha256(self.frozen_bytes).hexdigest())
        self.assertFalse(report["external_fetch_performed"])
        self.assertFalse(report["model_fitted"])
        self.assertFalse(report["biological_validation_performed"])
        for row, archived in zip(report["sources"], document["sources"]):
            self.assertEqual(row["source_id"], archived["id"])
            self.assertEqual(row["source_url"], archived["url"])
            self.assertEqual(row["reuse_terms"], archived["license_observation"])
            self.assertIsNone(row["individual_observations"])
            self.assertIsNone(row["individual_observation_count"])
            self.assertIsNone(row["calibration_record"])
            self.assertIsNone(row["independent_unit_ids"])
            self.assertFalse(row["eligible_raw_cohort"])
            self.assertEqual(row["eligibility"], "excluded")
            for reason in archived["reasons"]:
                self.assertIn(reason, row["exclusion_reasons"])
            for required in ("individual_observations", "calibration_record", "independent_unit_ids"):
                self.assertTrue(any(required in reason for reason in row["exclusion_reasons"]))
        self.assertEqual(self.frozen_catalog.read_bytes(), self.frozen_bytes)

    def test_complete_constructed_fields_can_satisfy_eligibility_rule(self):
        report = qualify_source_catalog(self.write_catalog(self.complete_constructed_document()))
        self.assertEqual(report["n_eligible_raw_cohorts"], 1)
        self.assertEqual(report["status"], "completed_eligible_raw_cohorts_present")
        row = report["sources"][0]
        self.assertTrue(row["eligible_raw_cohort"])
        self.assertEqual(row["individual_observation_count"], 1)
        self.assertEqual(row["exclusion_reasons"], [])

    def test_each_missing_required_field_retains_an_explicit_exclusion(self):
        required_fields = (
            "id", "url", "license_observation", "individual_observations",
            "calibration_record", "independent_unit_ids",
        )
        for field in required_fields:
            with self.subTest(field=field):
                document = self.complete_constructed_document()
                del document["sources"][0][field]
                document["n_eligible_raw_datasets_acquired"] = 0
                report = qualify_source_catalog(self.write_catalog(document))
                row = report["sources"][0]
                self.assertFalse(row["eligible_raw_cohort"])
                self.assertEqual(row["eligibility"], "excluded")
                self.assertTrue(any(field in reason for reason in row["exclusion_reasons"]))
                self.assertEqual(report["n_eligible_raw_cohorts"], 0)
                self.assertEqual(report["status"], "completed_zero_eligible_raw_cohorts")

    def test_complete_fields_do_not_clear_an_archived_source_exclusion(self):
        document = self.complete_constructed_document()
        document["sources"][0]["eligibility"] = "ineligible_as_acquired"
        document["sources"][0]["reasons"] = ["Constructed geometry exclusion only."]
        document["n_eligible_raw_datasets_acquired"] = 0
        row = qualify_source_catalog(self.write_catalog(document))["sources"][0]
        self.assertFalse(row["eligible_raw_cohort"])
        self.assertIn("Constructed geometry exclusion only.", row["exclusion_reasons"])

    def test_independent_units_require_disjoint_calibration_and_validation_ids(self):
        for unit_ids in (
            {"calibration": ["same-unit"], "validation": ["same-unit"]},
            {"calibration": ["unit-a"], "validation": []},
        ):
            with self.subTest(unit_ids=unit_ids):
                document = self.complete_constructed_document()
                document["sources"][0]["independent_unit_ids"] = unit_ids
                document["n_eligible_raw_datasets_acquired"] = 0
                row = qualify_source_catalog(self.write_catalog(document))["sources"][0]
                self.assertFalse(row["eligible_raw_cohort"])
                self.assertIsNone(row["independent_unit_ids"])
                self.assertTrue(any("independent" in reason.lower() for reason in row["exclusion_reasons"]))

    def test_cli_outputs_are_deterministic_and_existing_output_is_preserved(self):
        first, second = self.root / "first", self.root / "second"
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            self.assertEqual(main(["qualify-sources", str(self.frozen_catalog), "--out", str(first)]), 0)
        source_qualification(self.frozen_catalog, second)
        names = {"qualification.json", "qualification.csv", "REPORT.md"}
        self.assertEqual({path.name for path in first.iterdir()}, names)
        self.assertEqual({path.name for path in second.iterdir()}, names)
        for name in names:
            self.assertEqual((first / name).read_bytes(), (second / name).read_bytes())
        self.assertEqual(json.loads(stdout.getvalue()), json.loads((first / "qualification.json").read_text(encoding="utf-8")))

        with (first / "qualification.csv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["eligibility"] == "excluded" for row in rows))
        self.assertTrue(all(row["individual_observations"] == "" for row in rows))
        self.assertTrue(all(row["calibration_record"] == "" for row in rows))
        self.assertTrue(all(row["independent_unit_ids"] == "" for row in rows))
        self.assertTrue(all(json.loads(row["exclusion_reasons"]) for row in rows))
        markdown = (first / "REPORT.md").read_text(encoding="utf-8")
        self.assertIn("Eligible raw cohorts in this archived search: **0**", markdown)
        self.assertIn("completed qualification outcome", markdown)
        self.assertIn("Individual observations", markdown)
        self.assertIn("https://pmc.ncbi.nlm.nih.gov/articles/PMC5332570/", markdown)
        self.assertIn("https://pmc.ncbi.nlm.nih.gov/articles/PMC6545061/", markdown)

        sentinel = first / "sentinel.txt"
        sentinel.write_text("preserve me", encoding="utf-8")
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
            main(["qualify-sources", str(self.frozen_catalog), "--out", str(first)])
        self.assertEqual(caught.exception.code, 2)
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "preserve me")
        self.assertEqual(self.frozen_catalog.read_bytes(), self.frozen_bytes)


if __name__ == "__main__":
    unittest.main()
