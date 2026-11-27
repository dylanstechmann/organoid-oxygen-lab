"""Deterministic review of the frozen public-source qualification catalog.

This module reorganizes archived source decisions and missing-field evidence. It
does not acquire observations, perform a new review, or qualify raw profiles for
model fitting.
"""

import hashlib
import json
from pathlib import Path


_TEXT_REQUIREMENTS = (
    ("id", "Required source identity field `id` is missing or is not a nonempty string."),
    ("url", "Required source link field `url` is missing or is not a nonempty string."),
    ("license_observation", "Required reuse-terms field `license_observation` is missing or is not a nonempty string."),
    ("calibration_record", "Required `calibration_record` for individual observations is missing from the archived catalog."),
)


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate source-catalog key: {key}")
        result[key] = value
    return result


def _reject_nonfinite(value):
    raise ValueError(f"non-finite JSON value: {value}")


def _text(value):
    return value if isinstance(value, str) and value.strip() else None


def _individual_observations(source):
    observations = source.get("individual_observations")
    if (not isinstance(observations, list) or not observations
            or any(not isinstance(item, dict) or not item for item in observations)):
        return None, ("Required `individual_observations` records are missing or invalid in the archived catalog; "
                      "aggregate descriptions and figures do not fill this field.")
    return observations, None


def _independent_unit_ids(source):
    unit_ids = source.get("independent_unit_ids")
    roles = ("calibration", "validation")
    if not isinstance(unit_ids, dict) or set(unit_ids) != set(roles):
        return None, ("Required `independent_unit_ids` must provide nonempty, disjoint calibration and validation lists.")
    normalized = {}
    for role in roles:
        values = unit_ids[role]
        if (not isinstance(values, list) or not values
                or any(not isinstance(value, str) or not value.strip() or value != value.strip() for value in values)
                or len(set(values)) != len(values)):
            return None, ("Required `independent_unit_ids` must contain unique, nonempty, unpadded IDs in both "
                          "calibration and validation lists.")
        normalized[role] = list(values)
    overlap = set(normalized["calibration"]) & set(normalized["validation"])
    if overlap:
        return None, ("Independent calibration and validation biological-unit IDs overlap; one biological unit cannot "
                      "occupy both roles.")
    return normalized, None


def _archived_reasons(source):
    reasons = source.get("reasons", [])
    if not isinstance(reasons, list) or any(not isinstance(reason, str) or not reason.strip() for reason in reasons):
        return [], ["Archived source exclusion `reasons` must be a list of nonempty strings."]
    return list(reasons), []


def _source_row(source):
    if not isinstance(source, dict):
        raise ValueError("each source-catalog candidate must be an object")
    exclusions, reason_error = _archived_reasons(source)
    exclusions += reason_error
    for key, reason in _TEXT_REQUIREMENTS:
        if _text(source.get(key)) is None:
            exclusions.append(reason)

    observations, observation_reason = _individual_observations(source)
    if observation_reason:
        exclusions.append(observation_reason)
    unit_ids, unit_reason = _independent_unit_ids(source)
    if unit_reason:
        exclusions.append(unit_reason)

    catalog_eligibility = source.get("eligibility")
    if catalog_eligibility not in {"eligible_as_acquired", "ineligible_as_acquired"}:
        exclusions.append("Archived `eligibility` must be `eligible_as_acquired` or `ineligible_as_acquired`.")
    elif catalog_eligibility == "ineligible_as_acquired" and not exclusions:
        exclusions.append("Archived catalog marks the source ineligible but supplies no exclusion reason.")
    # Recorded exclusions cannot be cleared merely by adding metadata fields.
    eligible = catalog_eligibility == "eligible_as_acquired" and not exclusions
    return {
        "source_id": _text(source.get("id")),
        "source_url": _text(source.get("url")),
        "doi": _text(source.get("doi")),
        "reuse_terms": _text(source.get("license_observation")),
        "individual_observations": observations,
        "individual_observation_count": None if observations is None else len(observations),
        "calibration_record": _text(source.get("calibration_record")),
        "independent_unit_ids": unit_ids,
        "raw_measurement_sha256": source.get("raw_measurement_sha256"),
        "catalog_eligibility": catalog_eligibility,
        "eligible_raw_cohort": eligible,
        "eligibility": "eligible_as_acquired" if eligible else "excluded",
        "exclusion_reasons": exclusions,
        "archived_record": source,
    }


def qualify_source_catalog(catalog_path):
    """Apply explicit field-presence rules without changing archived decisions."""
    catalog_bytes = Path(catalog_path).read_bytes()
    document = json.loads(
        catalog_bytes.decode("utf-8-sig"),
        object_pairs_hook=_object,
        parse_constant=_reject_nonfinite,
    )
    if not isinstance(document, dict) or type(document.get("schema_version")) is not int or document["schema_version"] != 1:
        raise ValueError("source catalog schema_version must be 1")
    if document.get("kind") != "public-source-qualification-not-a-measurement-dataset":
        raise ValueError("unexpected source catalog kind")
    for key in ("audit_date", "decision", "method"):
        if _text(document.get(key)) is None:
            raise ValueError(f"source catalog needs nonempty {key}")
    sources = document.get("sources")
    if not isinstance(sources, list):
        raise ValueError("source catalog sources must be a list")
    declared_count = document.get("n_eligible_raw_datasets_acquired")
    if type(declared_count) is not int or declared_count < 0:
        raise ValueError("source catalog needs a nonnegative integer n_eligible_raw_datasets_acquired")

    rows = [_source_row(source) for source in sources]
    source_ids = [row["source_id"] for row in rows if row["source_id"] is not None]
    if len(set(source_ids)) != len(source_ids):
        raise ValueError("source catalog candidate IDs must be unique")
    eligible_count = sum(row["eligible_raw_cohort"] for row in rows)
    if eligible_count != declared_count:
        raise ValueError("derived eligible raw-cohort count disagrees with the archived catalog")
    return {
        "schema_version": 1,
        "kind": "archived-public-source-qualification-review",
        "qualification_scope": "frozen_archived_catalog_only",
        "audit_date": document["audit_date"],
        "catalog_sha256": hashlib.sha256(catalog_bytes).hexdigest(),
        "status": "completed_zero_eligible_raw_cohorts" if eligible_count == 0 else "completed_eligible_raw_cohorts_present",
        "n_candidates": len(rows),
        "n_eligible_raw_cohorts": eligible_count,
        "sources": rows,
        "archived_decision": document["decision"],
        "archived_method": document["method"],
        "external_fetch_performed": False,
        "model_fitted": False,
        "biological_validation_performed": False,
        "limits": [
            "This is a deterministic presentation of the archived search and its declared decisions, not a new source review.",
            "Zero eligible raw cohorts in this archived search does not prove that no suitable public dataset exists.",
            "Reuse-terms text is retained as an observation, not authenticated as legal permission for a raw dataset.",
            "Field completeness does not replace `qualify-profiles` raw-byte, calibration, geometry and independent-unit intake checks.",
            "No observations were digitized or constructed, no parameters were fitted, and no biological validation was performed.",
        ],
    }


def qualification_csv_rows(report):
    """Return stable CSV-ready rows; missing evidence stays empty, never zero."""
    rows = []
    for source in report["sources"]:
        rows.append({
            "source_id": source["source_id"] or "",
            "doi": source["doi"] or "",
            "source_url": source["source_url"] or "",
            "reuse_terms": source["reuse_terms"] or "",
            "individual_observations": "" if source["individual_observations"] is None else json.dumps(
                source["individual_observations"], ensure_ascii=False, separators=(",", ":")),
            "calibration_record": source["calibration_record"] or "",
            "independent_unit_ids": "" if source["independent_unit_ids"] is None else json.dumps(
                source["independent_unit_ids"], ensure_ascii=False, separators=(",", ":")),
            "eligibility": source["eligibility"],
            "exclusion_reasons": json.dumps(source["exclusion_reasons"], ensure_ascii=False, separators=(",", ":")),
        })
    return rows


def _markdown(value):
    if value is None:
        return "**missing from archived catalog**"
    return str(value).replace("|", "\\|").replace("\n", " ")


def qualification_markdown(report):
    """Render a reviewable table plus every archived and missing-field exclusion."""
    lines = [
        "# Archived public-source qualification",
        "",
        f"Audit date: {report['audit_date']}. Frozen catalog SHA-256: `{report['catalog_sha256']}`.",
        f"Candidates reviewed: {report['n_candidates']}. Eligible raw cohorts in this archived search: **{report['n_eligible_raw_cohorts']}**.",
        f"Outcome: `{report['status']}`. Zero eligible data is a completed qualification outcome, not a report-generation failure.",
        "",
        "No external data were fetched, no observations or reviewer decisions were added, and no model was fitted.",
        "",
        "| Source identity | Source link | Reuse terms | Individual observations | Calibration | Independent-unit IDs | Eligibility |",
        "|---|---|---|---|---|---|---|",
    ]
    for source in report["sources"]:
        identity = source["source_id"]
        if source["doi"]:
            identity = f"{identity}; DOI: {source['doi']}" if identity else f"DOI: {source['doi']}"
        observations = source["individual_observations"]
        observation_text = (None if observations is None
                            else f"{len(observations)} individual record(s) retained in JSON")
        unit_ids = source["independent_unit_ids"]
        unit_text = (None if unit_ids is None else
                     f"calibration: {', '.join(unit_ids['calibration'])}; validation: {', '.join(unit_ids['validation'])}")
        lines.append(
            f"| {_markdown(identity)} | {_markdown(source['source_url'])} | {_markdown(source['reuse_terms'])} | "
            f"{_markdown(observation_text)} | {_markdown(source['calibration_record'])} | {_markdown(unit_text)} | "
            f"{source['eligibility']} |"
        )
    lines += ["", "## Exclusions and archived observations", ""]
    for source in report["sources"]:
        lines.append(f"### {source['source_id'] or 'Missing source ID'}")
        lines.append("")
        description = source["archived_record"].get("measurement_description")
        if description:
            lines += [f"Archived measurement description: {description}", ""]
        archive_observation = source["archived_record"].get("archive_observation")
        if archive_observation:
            lines += [f"Archived file observation: {archive_observation}", ""]
        lines.append("Exclusion reasons:")
        lines.extend(f"- {reason}" for reason in source["exclusion_reasons"])
        lines.append("")
    lines += [
        "## Archived decision and limits",
        "",
        report["archived_decision"],
        "",
        report["archived_method"],
        "",
        *[f"- {limit}" for limit in report["limits"]],
        "",
    ]
    return "\n".join(lines)
