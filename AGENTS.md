# Agent instructions — organoid-oxygen-lab

Work only in this repository. This is a steady homogeneous-sphere oxygen
solver with Michaelis–Menten uptake. Parameters in the demo are illustrative.

## Do not

- Fit `vmax_mol_m3_s` to a story and call it measured uptake.
- Convert incubator %O2 into dissolved mol/m³ with a universal factor.
- Treat `threshold_mol_m3` as hypoxia, death, or potency.
- Treat `transfer_m_s` as a pump flow from `open-perfusion-rig`.
- Add vascularization, growth, or CFD unless you also add a new fixture that shows the old sphere tests still pass.
- Claim the Leedale or McMurtrey papers were reproduced.

## First commands

```bash
python -m pip install -e ".[plots]"
python -m unittest discover -s tests -v
oxygenlab demo --out artifacts/agent-demo --plot
```

Use a fresh `--out` directory every time.

## Improve, in this order

1. Read `docs/METHODS.md` and `docs/MODEL_CARD.md`. If a comment contradicts the mass balance, fix the comment or the code, and add a regression test.
2. `oxygenlab sweep-vmax` already asks whether core oxygen falls as illustrative `vmax` rises. Do not add another vmax grid. A different single-parameter question (for example `km` or surface `transfer_m_s`) is allowed if the report states the question and keeps units.
3. If you add a parameter, give its unit in the CLI/JSON schema and reject non-finite values the way existing settings are rejected.
4. Do not “calibrate” against a dataset that is not in this repo.

## Done when

Unit tests and the analytic mesh-refinement check still pass, and the new figure is labeled illustrative.
