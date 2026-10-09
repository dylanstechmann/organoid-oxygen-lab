# Source-profile sampling diagnostics — 2026-10-09

The AI assistant added spatial-coverage diagnostics to the existing source-linked
oxygen-profile intake. Reports now retain center/surface sampling, radial span,
largest internal gap and spherical geometric regions outside that span. A source
declaration of a complete profile no longer hides these sampling limitations.

The [constructed CLI example](../examples/spatial-coverage/REPORT.md) contains
no oxygen measurements. Samples at 20, 50 and 80 um in a declared 100 um sphere
leave geometric regions of 0.8% inward and 48.8% outward of the sample span.
Adding boundary points makes those regions zero but does not fill the gaps
between point measurements. No adequate-sampling threshold, hypoxia threshold,
viability result or fitted parameter is introduced.

All 66 unit tests passed in Docker `dev`, including five direct coverage tests
and an intake integration test. A fresh `demo --plot` run retained quadratic
analytic mesh refinement: maximum errors fell from approximately 1.04e-5 to
1.63e-7 mol/m3 over 20 to 160 shells. Maximum scenario mass-balance error was
approximately 1.41e-9. The illustrative plot was inspected. These are numerical
checks with illustrative settings; no biological calibration was performed.

Reproduce the new contrast into a fresh directory:

```sh
oxygenlab coverage-demo --out artifacts/new-spatial-coverage
```

The reduced sphere still has homogeneous uptake, no growth and no vascular flow.
The ResearchDesk dental/organoid tracks cite this repo as an independent methods
tool; they do not turn its illustrative settings into dental culture parameters.
