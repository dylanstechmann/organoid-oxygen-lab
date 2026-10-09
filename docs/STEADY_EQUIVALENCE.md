# Steady oxygen profiles and absolute rates

The steady model depends on `Da = Vmax R² / (D c_bulk)`, `Bi = kR/D` for a
finite surface film, and `Km/c_bulk`. At fixed radius, bath and Km, changing
`D -> aD`, `Vmax -> aVmax` and `k -> ak` preserves these groups for positive a.
The fixed-surface limit remains fixed-surface. Both uptake laws implemented here
share this symmetry. Concentrations at the same physical radii remain unchanged.

Absolute integrated uptake and influx scale by a; the characteristic diffusion
scale `R²/D` scales by 1/a. That scale is not a measured equilibration time and
does not supply a viability or culture-condition recommendation.

The [constructed example](../examples/steady-equivalence/REPORT.md) doubles
those three parameters. Its numerical profiles match exactly in this run;
uptake doubles, while `R²/D` changes from 45 to 22.5 seconds. These values are
illustrative, not measured organoid properties or biological parameter recovery.

This ambiguity applies when those rates are jointly unknown. Independently
fixing D removes this particular equivalent family; it does not guarantee all
remaining parameters are identifiable. The existing local-sensitivity command
holds D fixed, so its finite-difference matrix addresses a different question.
More precise steady concentrations or more points alone cannot remove the
joint-rate symmetry. Independent rate, diffusivity or temporal information
could constrain it only with suitable measurement and model assumptions.

No parameters were fitted. All 76 unit tests passed in Docker `dev`, including
an independent zero-order analytic check of the equivalence. A fresh `demo
--plot` retained analytic mesh refinement and maximum scenario mass-balance
error 1.41e-9. The new illustration and numerical demo were visually checked.
