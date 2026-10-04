# Model and parameter card

## Purpose

Transparent numerical exploration of oxygen supply and uptake in a homogeneous
sphere. The application context is organoid and tissue-engineering research.
The code does not evaluate aging, rejuvenation, lifespan, transplantation
success, potency or patient outcomes.

## Parameter provenance

Every default is an **illustrative software parameter**, selected to generate
well-resolved examples with and without appreciable oxygen depletion. The set
is not fitted to measurements, claimed as a named cell-line parameter set, or
transferred from a paper as an experimentally validated recipe. The literature
supports diffusion–uptake modeling as a research approach; it does not validate
these numbers or this reduced geometry.

The example files are outputs of `oxygenlab demo`, not experimental data.
Size and surface conditions are deliberately varied model scenarios, not an
uncertainty distribution or recommended culture settings. No statistical
confidence interval or cell-survival percentage is reported.

## Validation status

Numerical verification includes conservation, analytic convergence, positivity,
monotonicity and comparison to an independent nonlinear boundary-value solver.
Transient verification additionally compares early zero-uptake diffusion with
an independent spherical series and checks first-order temporal refinement.
Core-change timing is a numerical summary, not an assay response time.
It establishes that the specified equations are being solved consistently in
the tested cases. No biological measurements have been used for calibration
or out-of-sample validation.

## Assumptions requiring review

1. Spherical shape, homogeneous effective diffusivity and uniform maximal
   uptake density.
2. Either steady state or transient diffusion from recorded uniform initial
   oxygen, with continuously maintained bath concentration and constant boundary
   resistance.
3. Michaelis–Menten uptake with no changes in metabolic state, proliferation,
   death, density, differentiation or stress feedback.
4. No vascularization, perfused channels, necrotic core, convection inside
   tissue, other nutrients or coupled signaling.
5. A user-selected threshold used only to summarize concentration profiles.

## Collaborator handoff

Record cell type, geometry and size distribution, measurement units, culture
conditions, oxygen sensor calibration, density and uptake estimation method.
Estimate parameters from one dataset, state identifiability limitations, and
validate on independent sizes or conditions. Measure the bath near the tissue
or model its transport explicitly. A measured pump flow cannot by itself
identify the surface mass-transfer coefficient.

Keep private measurements outside Git. Record a source hash and data license
for any future calibration dataset. Generated outputs retain the complete
parameter object, its canonical SHA-256, software versions and solver checks.
