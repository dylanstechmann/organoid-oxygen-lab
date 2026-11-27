# Archived public-source qualification

Audit date: 2026-10-04. Frozen catalog SHA-256: `e421226f34591c9d8aafcbcd64b9f59daadb9e8d14be36f2e35e638628afe515`.
Candidates reviewed: 2. Eligible raw cohorts in this archived search: **0**.
Outcome: `completed_zero_eligible_raw_cohorts`. Zero eligible data is a completed qualification outcome, not a report-generation failure.

No external data were fetched, no observations or reviewer decisions were added, and no model was fitted.

| Source identity | Source link | Reuse terms | Individual observations | Calibration | Independent-unit IDs | Eligibility |
|---|---|---|---|---|---|---|
| murphy2017-msc; DOI: 10.1098/rsif.2016.0851 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5332570/ | Figshare supplementary figures are CC BY 4.0; no individual raw profile dataset was acquired. | **missing from archived catalog** | **missing from archived catalog** | **missing from archived catalog** | excluded |
| sheth2019-hemisphere; DOI: 10.1098/rspa.2018.0647 | https://pmc.ncbi.nlm.nih.gov/articles/PMC6545061/ | Inspected Figshare manuscript archive is CC BY 4.0; no qualified raw numeric measurement cohort acquired. | **missing from archived catalog** | **missing from archived catalog** | **missing from archived catalog** | excluded |

## Exclusions and archived observations

### murphy2017-msc

Archived measurement description: Human bone-marrow MSC spheroids; Unisense OX-10 microelectrode at 25 degC, readings every 10 um toward center, n=5 in reported oxygen-profile figure. Measured medium oxygen reported as 267 uM (0.267 mol/m3).

Exclusion reasons:
- Article presents profile means/plots; individual unaveraged radial tables and biological-unit identities were not acquired.
- Supplementary materials listed by the primary article are three figure PDFs, not identified raw radial tables.
- Figure digitization would not recover raw replicates or independent held-out identities.
- Required `calibration_record` for individual observations is missing from the archived catalog.
- Required `individual_observations` records are missing or invalid in the archived catalog; aggregate descriptions and figures do not fill this field.
- Required `independent_unit_ids` must provide nonempty, disjoint calibration and validation lists.

### sheth2019-hemisphere

Archived measurement description: MCF7-R tumor hemi-spheroids measured using five gold micro-disc oxygen electrodes.

Archived file observation: Archive contains two DOCX and four PPTX files; no embedded spreadsheet/data members found in the four PPTX files.

Exclusion reasons:
- A hemi-spheroid on an electrode substrate needs an explicit geometry and boundary applicability argument before using the complete-sphere solver.
- An example of voltammetry processing is not a raw calibrated complete-sphere profile cohort.
- Direct radial per-specimen oxygen tables with independent calibration/validation biological identities were not acquired.
- Required `calibration_record` for individual observations is missing from the archived catalog.
- Required `individual_observations` records are missing or invalid in the archived catalog; aggregate descriptions and figures do not fill this field.
- Required `independent_unit_ids` must provide nonempty, disjoint calibration and validation lists.

## Archived decision and limits

No raw dataset acquired satisfies this repository measurement contract. No public-data fit or physical/biological validation was performed.

Inspect primary article/repository methods, data availability and publicly listed files; refuse figure means, synthetic/model outputs, unidentified biological replicates, restricted files and incompatible units/geometry. This is a bounded source audit, not proof that no suitable public dataset exists.

- This is a deterministic presentation of the archived search and its declared decisions, not a new source review.
- Zero eligible raw cohorts in this archived search does not prove that no suitable public dataset exists.
- Reuse-terms text is retained as an observation, not authenticated as legal permission for a raw dataset.
- Field completeness does not replace `qualify-profiles` raw-byte, calibration, geometry and independent-unit intake checks.
- No observations were digitized or constructed, no parameters were fitted, and no biological validation was performed.
