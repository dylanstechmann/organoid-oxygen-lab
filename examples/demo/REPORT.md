# Oxygen model numerical demonstration

All parameters are illustrative. No biological measurements are included.

## Boundary resistance and size

The bath is fixed at 0.2 mol/m³. The user-selected reporting threshold is 0.02 mol/m³.
Finite-transfer scenarios use k = 2 × 10⁻⁵ m/s; fixed-surface scenarios have no exterior resistance.

| Boundary | Radius (µm) | Minimum sampled oxygen (mol/m³) | Volume below threshold |
|---|---:|---:|---:|
| fixed_surface | 50 | 0.19603 | 0.000 |
| fixed_surface | 100 | 0.18417 | 0.000 |
| fixed_surface | 200 | 0.13738 | 0.000 |
| fixed_surface | 300 | 0.064072 | 0.000 |
| fixed_surface | 400 | 0.0042046 | 0.050 |
| fixed_surface | 600 | 3.2743e-06 | 0.291 |
| fixed_surface | 800 | 5.8401e-09 | 0.443 |
| finite_transfer | 50 | 0.18025 | 0.000 |
| finite_transfer | 100 | 0.15293 | 0.000 |
| finite_transfer | 200 | 0.07848 | 0.000 |
| finite_transfer | 300 | 0.0078964 | 0.061 |
| finite_transfer | 400 | 0.00016048 | 0.291 |
| finite_transfer | 600 | 2.3952e-07 | 0.536 |
| finite_transfer | 800 | 4.7156e-10 | 0.656 |

## Analytic reference and mesh refinement

The constant-uptake reference has a closed-form solution. Halving shell width should reduce
the maximum concentration error by about a factor of four in this smooth case.

| Shells | Maximum absolute error (mol/m³) |
|---|---:|
| 20 | 1.0417e-05 |
| 40 | 2.6042e-06 |
| 80 | 6.5104e-07 |
| 160 | 1.6276e-07 |

Maximum relative mass-balance error across the size scenarios: 1.41e-09.

Unit tests additionally compare the nonlinear profile against an independent SciPy boundary-value solver.
Equation residuals and numerical convergence do not establish accuracy for a particular organoid.

![Oxygen profiles and size scenarios](oxygen_demo.png)
