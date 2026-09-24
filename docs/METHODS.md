# Equations and numerical verification

Let x=r/R, u=c/c_bulk, κ=K_m/c_bulk and
Da=V_max R²/(D c_bulk). For Michaelis–Menten uptake,

$$\frac{1}{x^2}(x^2u')'=\mathrm{Da}\frac{u}{\kappa+u}.$$

The boundary conditions are u′(0)=0 and either u(1)=1 or
u′(1)=Bi(1−u(1)), where Bi=kR/D. There is no consumption outside the sphere;
the finite-transfer boundary replaces that exterior region with a resistance.

## Finite volumes

Uniform shell faces x_i divide [0,1]. Unknowns are at shell midpoints.
The dimensionless shell volumes without 4π are
v_i=(x_(i+1)³−x_i³)/3. Interface conductance is x_face²/Δx.
The outer half-cell and film resistances are added in series:

$$g_b=\left(\frac{\Delta x}{2}+\frac{1}{\mathrm{Bi}}\right)^{-1}.$$

For a fixed surface the second resistance vanishes. Integrated equations are

$$A u-b+\mathrm{Da}\,v\odot\frac{u}{\kappa+u}=0.$$

A is tridiagonal with positive diagonal and negative neighbor conductances.
The central face has zero area, so no division by zero at r=0 is needed.
The Newton Jacobian adds Da v_i κ/(κ+u_i)² to the diagonal. A banded solve
computes the Newton step. Backtracking preserves nonnegativity and reduces
the residual. Because u/(κ+u) ≤ u/κ for nonnegative u, the linear first-order
uptake solution supplies a positive lower-bound initialization.

Convergence requires max_i |residual_i|/v_i divided by max(1,Da) to be below
the requested tolerance. The default is 1e−8. A failed solve raises an error
instead of exporting a plausible-looking profile. The reported residual is a
numerical diagnostic, not an error bound for biology or mesh discretization.

## Conservation and summaries

The integrated uptake is 4πR³ V_max sum_i v_i u_i/(κ+u_i).
The surface influx is 4π D c_bulk R g_b(1−u_last). Their relative discrepancy
is reported separately from the equation residual. For finite transfer,
c_surface is reconstructed from the influx and film law.

The volume-weighted mean uses weights 3v_i, which sum to one. The volume below
a threshold sums weights for shells whose midpoint concentration is below it.
That piecewise-constant estimate changes discretely as shells cross the
threshold. Refine the mesh before interpreting a near-threshold volume fraction.
The minimum sampled concentration is at R/(2N), not exactly at the center.

## Independent references

For uniform constant uptake q, the exact concentration is

$$c(r)=c_{bulk}-\frac{qR}{3k}-\frac{q(R^2-r^2)}{6D}.$$

Omit the film term for a fixed surface. The numerical example uses a radius
where this expression remains positive and checks its error on progressively
finer grids. Negative constant-uptake solutions are rejected rather than clipped.

Tests also solve the nonlinear boundary-value problem through
[SciPy's separate `solve_bvp` implementation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_bvp.html).
That collocation solver treats the singular term 2u′/x explicitly. Profile
comparison does not reuse the finite-volume equations.

Relevant primary context: [Leedale et al. (2021)](https://doi.org/10.1371/journal.pone.0244070)
and [McMurtrey (2016)](https://doi.org/10.1089/ten.tec.2015.0375).
All equations and boundary assumptions implemented here are specified above.
