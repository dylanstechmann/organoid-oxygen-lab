# Source-declared measurement uncertainty

Profile intake accepts an optional `measurement_uncertainty` object on each
specimen. Its four optional quantities are `specimen_radius`, `radial_position`
(both `um`), `oxygen` and `bath_oxygen` (both `mol/m3`). Omitted quantities are
reported as `not_reported_in_manifest` with a null value, rather than zero.
Existing schema-1 manifests remain valid.

Each entry needs `kind`, `value`, `unit`, `source_file`, `source_record` and
`scope`. Kinds are `standard_deviation`, `standard_uncertainty`,
`expanded_uncertainty` or `absolute_bound`. Values must be finite and
nonnegative; explicit zero remains a source declaration. Expanded uncertainty
also requires a positive `coverage_factor`. Other kinds cannot have that field.
Source files must already be hash-verified by intake.

For example, this is a **constructed schema fragment**, not measured precision:

```json
{"measurement_uncertainty":{"oxygen":{
  "kind":"absolute_bound", "value":0.003, "unit":"mol/m3",
  "source_file":"native", "source_record":"example calibration record",
  "scope":"Declared sensor bound; point applicability must be reviewed"
}}}
```

The report retains the source-file hash and scope. It does not authenticate the
record, infer a probability distribution or confidence interval, treat SD as SE,
assume independent errors, or propagate uncertainty into fitted parameters.
A specimen-level annotation does not silently supply an error bar for every
point. Sensor precision also cannot establish spatial completeness, biological
replication, steady state or homogeneous uptake. Those remain separate checks.

All 71 unit tests passed in Docker `dev`; the intake tests passed again after
an equivalent ordered-radius iteration cleanup. Ruff passed for the changed
modules and new tests. A fresh illustrative `demo --plot` retained analytic
mesh-refinement errors of 1.04e-5 to 1.63e-7 mol/m3 from 20 to 160 shells and
maximum scenario mass-balance error of 1.41e-9. The plot was inspected. These
are numerical checks, without biological calibration or uncertainty propagation.
