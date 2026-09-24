# Spherical oxygen model

Illustrative numerical model; parameters need experimental calibration.

Radius: 300 µm. Concentrations are mol/m³ (= mM).

| Quantity | Result |
|---|---:|
| Minimum sampled oxygen (mol/m³) | 0.00789644 |
| Surface oxygen (mol/m³) | 0.115714 |
| Volume below user threshold | 0.0610 |
| Integrated uptake (mol/s) | 1.90651e-12 |
| Relative mass-balance error | 6.73e-13 |
| Scaled equation residual | 5.59e-11 |

## Interpretation

- Illustrative model, not calibrated to a cell line or organoid experiment.
- The minimum is sampled at the innermost shell center, not exactly at r=0.
- Below-threshold volume uses piecewise-constant shell values; the threshold is user supplied, not a viability cutoff.
- Spherical, homogeneous, steady state with a fixed bath; no vascularization, growth, heterogeneous uptake or cell death.
- Surface mass transfer is a boundary resistance and cannot be converted to pump flow without another model.
