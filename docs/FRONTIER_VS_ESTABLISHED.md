# Frontier versus established

**Status date: 2026-11-27**

This note separates reproducible software qualification and numerical checks
from biological claims that remain unestablished. The results below are from an
AI-assisted, agent-run software analysis; they are not an independent
experimental review.

## Established: reproducible software and numerical checks

- `oxygenlab qualify-sources docs/public-source-qualification.json --out <fresh-dir>`
  deterministically produces JSON, CSV and Markdown qualification reports from
  the byte-hashed frozen catalog. Repeated fresh output directories were
  byte-identical. Missing observations, calibration records and independent-unit
  IDs remain explicit, and archived exclusion reasons are retained.
- The archived search contains **two candidates and zero eligible raw cohorts**.
  `completed_zero_eligible_raw_cohorts` is a completed bounded-search outcome,
  not evidence that no suitable public dataset exists elsewhere.
- `python -m unittest discover -s tests -v` completed **88 tests successfully**,
  with **two optional plot tests skipped** because Matplotlib was unavailable.
- The independent constant-uptake analytic regression passed its smooth-case
  mesh-refinement check. This checks numerical convergence of the implemented
  equations; it is not biological validation.
- Source-linked measurement intake checks identifier whitespace and reports
  calibration/validation counts separately for biological units, specimens and
  radial points. Constructed fixtures establish this software behavior only.

See the [generated qualification report](../examples/source-qualification/REPORT.md),
[numerical methods](METHODS.md) and [model card](MODEL_CARD.md) for scope and
limitations.

## Frontier and not established

- **Public dataset availability beyond the archived search:** unknown. The two
  archived exclusions do not establish that no suitable public cohort exists.
- **Measured oxygen uptake and calibration:** no eligible raw cohort was acquired,
  no oxygen uptake or boundary parameter was fitted, and no sensor or incubation
  condition was calibrated. Incubator gas percentage is not converted to local
  dissolved oxygen with a universal factor.
- **Biological validation:** not performed. Software intake, numerical residuals
  and constructed fixtures do not validate organoid physiology, geometry,
  survival, maturation or independent experimental prediction.
- **Hypoxia, death or potency interpretation:** not established. A user-selected
  reporting threshold and the volume below it are descriptive numerical summaries,
  not biological response cutoffs.
- **Paper reproduction:** not claimed. Leedale et al. and McMurtrey provide
  research context; the checks here do not reproduce their experimental systems,
  results or validation.

## Update when evidence changes

Update this note, the model card and the generated qualification report whenever
new evidence changes a boundary. A future measured-data claim must identify the
exact source bytes and reuse terms, individual observations, calibration and
unit conversion records, and independent calibration/held-out biological-unit
IDs. Record the fresh command, test count, skips and output hashes with the
change. Move an item only when the corresponding evidence and independent checks
exist; figure digitization, additional model simulations or a larger unit-test
count alone do not establish biological calibration or validation.
