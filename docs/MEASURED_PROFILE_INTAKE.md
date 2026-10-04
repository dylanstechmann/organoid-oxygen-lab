# Public-source qualification and measured oxygen profiles

The bounded public-source audit on 2026-10-04 acquired **zero eligible raw
complete-sphere oxygen-profile cohorts**. No uptake/boundary parameter was
fitted, and no biological validation was performed.
[The catalog](public-source-qualification.json) records the exclusions.

[Murphy et al. (2017)](https://doi.org/10.1098/rsif.2016.0851) measured MSC
spheroids with an OX-10 electrode at 25 degrees C, stepping 10 micrometers
toward the center. The profile figure reports n=5 and medium oxygen of
267 micromolar (0.267 mol/m3). The published means/figures do not supply
individual radial tables and independent biological identities. The linked
Figshare collection lists three supplementary figure PDFs under CC BY 4.0;
[inventory and metadata hashes](public-supplement-inventory.json) are retained.
No figures were digitized into supposed raw data.

[Sheth and Gratzl (2019)](https://doi.org/10.1098/rspa.2018.0647) used a radial
electrode array under tumor hemi-spheroids. Its downloadable manuscript archive
is CC BY 4.0 and contains two Word documents and four presentations. Inspection
found no embedded spreadsheet/data files in those four presentations. The
archive SHA-256 is `e616d985f5e71d693d57ca1ccf9a54cae3d5b43af04d4273d52cc661a84f3d5c`.
The archive does not provide a qualified raw, calibrated, independently
identified complete-sphere cohort. A hemi-sphere might be modeled through a
reflection argument under suitable substrate/symmetry conditions, but that
requires an explicit geometry/boundary justification; it is not established by
this solver's existing sphere tests.

[Page snapshots and hashes](public-source-access-20261004.json) are retained
locally under `artifacts/public-source-audit-20261004/`; they are source-access
evidence, not measurements. Recheck pages in a fresh directory:

```bash
.venv/bin/python scripts/probe_public_sources.py --out artifacts/new-source-audit
```

For reproducible supplement inspection, download the author archive from the
[Figshare file](https://ndownloader.figshare.com/files/15200057) into a local
artifact directory, then run:

```bash
.venv/bin/python scripts/inspect_public_supplement.py /path/to/rspa20180647_si_001.zip --out artifacts/new-supplement-inventory.json
```

The inspector verifies the recorded hash before listing members and embedded
data. It executes no author code and does not turn figures into measurements.
Raw author files remain local artifacts; original code and curated audit text
remain MIT. The author archive remains CC BY 4.0 with the attribution above.

## Measurement cohort intake

The new command checks a source-linked cohort before future fitting:

```bash
oxygenlab qualify-profiles /path/to/profiles.csv /path/to/manifest.json --out artifacts/new-profile-intake
```

It writes `qualification.json` and `REPORT.md` into a fresh directory. It does
**no fitting, model scoring or biological validation**. Failure leaves no
successful report. Numeric eligibility is distinct from proof that a tissue
is homogeneous, spherical or at steady state.

The CSV header must be exactly:

```csv
specimen_id,radius_um,oxygen_mol_m3
```

Record radial position from the center in micrometers, strictly increasing
within each specimen, with at least three unique positions inside its measured
radius. Oxygen is nonnegative dissolved mol/m3. Keep original sensor outputs
and the source-supported conversion/calibration record. Gas oxygen percentage,
pressure, normalized fluorescence, hypoxia reporter expression and
model-inferred oxygen cannot silently become dissolved concentration. Negative
sensor estimates need an explicit censoring/uncertainty analysis before use;
this narrow input format rejects them rather than clipping them.

## Manifest contract (schema version 1)

| Field | Required value or meaning |
|---|---|
| `schema_version` | Integer `1` |
| `origin` | `direct_calibrated_sensor`; synthetic, averaged, digitized and model-inferred profiles are rejected |
| `source_url`, `license`, `attribution`, `acquisition_record` | Nonempty source, applicable use permission, credit and measurement/conversion record |
| `units` | Exactly `{"radius":"um","oxygen":"mol/m3","temperature":"degC"}` |
| `csv_sha256` | Hash of the exact supplied CSV bytes |
| `source_files` | Object keyed by file ID, each with relative native `path` and exact `sha256` |
| `specimens` | Object keyed by exactly the CSV's specimen IDs |

Each specimen needs `original_specimen_id`, `biological_unit_id`,
`biological_identity_record`, `source_trace_id`, `source_file`,
`sensor_calibration_record`, `radius_measurement_record`,
`bath_measurement_record`, `steady_state_record`, `cell_type`,
`culture_conditions`, `geometry: "sphere"`, `steady_state: true`,
`complete_profile: true`, finite `temperature_degC`, positive
`specimen_radius_um`, positive **measured local** `bath_oxygen_mol_m3`, and
`role` set to `calibration` or `validation`.

Define biological unit at the level intended for generalization (for example,
independent donor or independent culture preparation) using source identities.
At least one calibration and one held-out validation unit are required. One
biological unit cannot occupy both roles; one original specimen or source trace
cannot be renamed into multiple specimens. Alias source filenames are checked
using bytes hashes. Multiple genuinely separate traces may share one native
file if their original trace identities are distinct and supported. Source
paths must stay inside the manifest directory. Do not split radial positions
from one specimen into training and validation.

The full specimen records, source hashes and manifest hash are included in the
report. The gate checks declared provenance and identity consistency; it cannot
authenticate acquisition or catch deliberately falsified biological IDs.
Constructed test fixtures establish this software behavior only.

Once actual qualified measurements are available, independently constrain
transport/bath parameters, assess uptake and boundary identifiability, fit only
on the calibration units, freeze parameters, and score the held-out units with
measurement uncertainty and deviations from the sphere assumptions visible.
A pump flow is not the surface mass-transfer coefficient. This intake command
does not change the illustrative parameters or the model's validation status.
