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

## Transient integration and exposure timing

The same conservative shell fluxes are used in backward Euler:

`v * (u_new - u_old) / delta_tau + A*u_new - b + Da*v*rate(u_new) = 0`.

Here `tau = t*D/R^2`. Accepted nonlinear steps preserve nonnegative oxygen and
meet a residual scaled by `max(1, Da, 1/delta_tau)`. The conservation diagnostic
compares the change in integrated shell oxygen with the step-end surface influx
minus uptake, integrated over the step. It is distinct from the distance to the
steady profile. Backward Euler is first order in time; saved points include t=0.

An independent no-uptake, fixed-surface reference follows by setting
`w = x*(1-c/c_bulk)` for an initially anoxic sphere. Then `w_tau = w_xx`,
`w(0,tau)=w(1,tau)=0` and `w(x,0)=x`. Separation of variables gives

`c/c_bulk = 1 - sum_n [2*(-1)^(n+1)/(n*pi*x)] * sin(n*pi*x) * exp(-n^2*pi^2*tau)`.

The regression test evaluates this series at positive shell midpoints at
`tau=0.1`, before steady state, and checks first-order temporal refinement on
a fixed fine radial mesh. It does not reuse finite-volume coefficients.

Absolute fractions of steady-core concentration may be unreachable when the
initial core is above steady. Separate core-change metrics use
`target = c_initial + fraction*(c_steady-c_initial)` and interpolate the first
bracketing time samples in the direction of this change. A target outside the
simulated horizon is null. A zero core difference gives zero core-change time;
this does not imply that the entire initial profile was steady. These are
numerical exposure summaries, not cell-survival or biological-response times.

## Local sensitivity and constructed inverse illustration

`oxygenlab sensitivity` perturbs each positive parameter by a relative step in
log space, solves the steady profile at both sides, and forms a centered
finite-difference sensitivity matrix over radial locations. Column cosines and
singular values reveal local tradeoffs at the chosen baseline. This profile
geometry is not a proof of identifiability from measured sensors, and the
diagnostic does not include measurement error or parameter priors.

`oxygenlab inverse-demo` creates an exact profile from the illustrative default
model, adds seeded synthetic Gaussian noise, then fits only Vmax and surface
transfer while holding other settings fixed. The observations and fitted values
are written with the report. Because generation and fitting use the same model
family, this verifies software plumbing only; it cannot calibrate organoid
oxygen uptake or a real device coefficient.
