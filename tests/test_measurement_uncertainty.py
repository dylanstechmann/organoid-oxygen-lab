"""Declaration validation, not experimental uncertainty estimation."""

import unittest

from oxygenlab.measurement_uncertainty import uncertainty_declarations


class UncertaintyTests(unittest.TestCase):
    def setUp(self):
        self.sources = {"native": "a" * 64}
        self.entry = {"kind": "standard_deviation", "value": .002, "unit": "mol/m3",
                      "source_file": "native", "source_record": "constructed fixture",
                      "scope": "Between repeated sensor readings; constructed test only"}

    def audit(self, entry):
        return uncertainty_declarations({"measurement_uncertainty": {"oxygen": entry}}, self.sources)

    def test_absence_and_explicit_zero_stay_distinct(self):
        missing = uncertainty_declarations({}, self.sources)["quantities"]["oxygen"]
        zero = self.audit({**self.entry, "value": 0})["quantities"]["oxygen"]
        self.assertIsNone(missing["value"])
        self.assertEqual(zero["value"], 0)
        self.assertNotEqual(missing["status"], zero["status"])

    def test_sd_retains_source_scope_without_ci_or_propagation(self):
        report = self.audit(self.entry)
        self.assertEqual(report["quantities"]["oxygen"]["source_sha256"], "a" * 64)
        self.assertEqual(report["quantities"]["oxygen"]["kind"], "standard_deviation")
        self.assertFalse(report["confidence_interval_inferred"])
        self.assertFalse(report["propagated_to_model"])

    def test_expanded_uncertainty_requires_factor_and_never_infers_confidence(self):
        entry = {**self.entry, "kind": "expanded_uncertainty"}
        for factor in (None, 0, True, float("nan")):
            with self.assertRaises(ValueError):
                self.audit({**entry, "coverage_factor": factor})
        self.assertFalse(self.audit({**entry, "coverage_factor": 2})["confidence_interval_inferred"])
        with self.assertRaises(ValueError):
            self.audit({**self.entry, "coverage_factor": 2})

    def test_invalid_units_sources_scope_and_values_fail(self):
        for patch in ({"value": True}, {"value": -.1}, {"value": float("inf")},
                      {"unit": "%"}, {"kind": "95_percent_CI"}, {"source_file": "missing"},
                      {"scope": ""}, {"extra": 1}):
            with self.subTest(patch=patch), self.assertRaises(ValueError):
                self.audit({**self.entry, **patch})
        with self.assertRaises(ValueError):
            uncertainty_declarations({"measurement_uncertainty": {"unknown": self.entry}}, self.sources)


if __name__ == "__main__":
    unittest.main()
