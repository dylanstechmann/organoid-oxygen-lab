"""Qualify source-linked steady spherical oxygen profiles; never fit parameters."""

import csv
import hashlib
import io
import json
import math
from pathlib import Path


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate measurement-manifest key: {key}")
        result[key] = value
    return result


def _text(record, key):
    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"measurement manifest needs nonempty {key}")
    return value


def _number(record, key):
    value = record.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"measurement manifest needs finite {key}")
    return value


def qualify_profiles(csv_path, manifest_path):
    """Check measurement eligibility and source/group leakage without model fitting.

    Requires held-out biological units. An intake pass checks declared provenance
    and numeric compatibility; it cannot authenticate acquisition or validate the
    homogeneous-sphere assumptions. Do not report its counts as biological validation.
    """
    input_bytes = Path(csv_path).read_bytes()
    path = Path(manifest_path)
    manifest_bytes = path.read_bytes()
    document = json.loads(manifest_bytes.decode("utf-8-sig"), object_pairs_hook=_object)
    if not isinstance(document, dict) or type(document.get("schema_version")) is not int or document["schema_version"] != 1:
        raise ValueError("measurement manifest schema_version must be 1")
    if document.get("origin") != "direct_calibrated_sensor":
        raise ValueError("direct_calibrated_sensor required; synthetic, figure-digitized, averaged and model-inferred oxygen are ineligible")
    for key in ("source_url", "license", "attribution", "acquisition_record"):
        _text(document, key)
    if document.get("units") != {"radius": "um", "oxygen": "mol/m3", "temperature": "degC"}:
        raise ValueError("explicit canonical units required; no gas-percent or pressure conversion is assumed")
    digest = hashlib.sha256(input_bytes).hexdigest()
    if document.get("csv_sha256") != digest:
        raise ValueError("oxygen CSV hash mismatch")
    files = document.get("source_files")
    if not isinstance(files, dict) or not files:
        raise ValueError("measurement manifest needs source_files")
    verified = {}
    source_root = path.parent.resolve()
    for file_id, entry in files.items():
        if not isinstance(file_id, str) or not file_id.strip() or not isinstance(entry, dict):
            raise ValueError("invalid source_file entry")
        relative = Path(_text(entry, "path"))
        source_path = (source_root / relative).resolve()
        if relative.is_absolute() or not source_path.is_relative_to(source_root):
            raise ValueError("source_file must be relative and inside manifest directory")
        source_digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
        if entry.get("sha256") != source_digest:
            raise ValueError(f"source_file hash mismatch: {file_id}")
        verified[file_id] = source_digest
    reader = csv.DictReader(io.StringIO(input_bytes.decode("utf-8-sig")))
    expected = ["specimen_id", "radius_um", "oxygen_mol_m3"]
    if reader.fieldnames != expected:
        raise ValueError("CSV header must be specimen_id,radius_um,oxygen_mol_m3")
    profiles = {}
    for line, row in enumerate(reader, 2):
        if None in row or any(value is None for value in row.values()) or not row["specimen_id"].strip():
            raise ValueError(f"row {line}: incomplete oxygen profile")
        try:
            radius, concentration = float(row["radius_um"]), float(row["oxygen_mol_m3"])
        except ValueError as exc:
            raise ValueError(f"row {line}: invalid radius/concentration") from exc
        if not math.isfinite(radius) or not math.isfinite(concentration) or radius < 0 or concentration < 0:
            raise ValueError(f"row {line}: radius and concentration must be finite and nonnegative")
        profiles.setdefault(row["specimen_id"].strip(), []).append((radius, concentration))
    specimens = document.get("specimens")
    if not profiles or not isinstance(specimens, dict) or set(specimens) != set(profiles):
        raise ValueError("manifest must cover exactly all CSV specimens")
    original_ids, source_traces, biological_roles = set(), set(), {}
    summary = []
    for specimen_id, points in profiles.items():
        record = specimens[specimen_id]
        if not isinstance(record, dict):
            raise ValueError(f"invalid specimen record: {specimen_id}")
        for key in ("original_specimen_id", "biological_unit_id", "source_trace_id", "biological_identity_record", "sensor_calibration_record", "radius_measurement_record", "bath_measurement_record", "steady_state_record", "cell_type", "culture_conditions"):
            _text(record, key)
        if record.get("geometry") != "sphere" or record.get("steady_state") is not True or record.get("complete_profile") is not True:
            raise ValueError("complete steady spherical profile declarations required")
        role = record.get("role")
        if role not in {"calibration", "validation"}:
            raise ValueError("specimen role must be calibration or validation")
        biological_id = record["biological_unit_id"]
        if biological_id in biological_roles and biological_roles[biological_id] != role:
            raise ValueError("one biological unit cannot appear in both calibration and validation")
        biological_roles[biological_id] = role
        original_id = record["original_specimen_id"]
        if original_id in original_ids:
            raise ValueError("one original specimen cannot be split into multiple specimen_ids")
        original_ids.add(original_id)
        file_id = _text(record, "source_file")
        if file_id not in verified:
            raise ValueError("unknown specimen source_file")
        identity = (verified[file_id], record["source_trace_id"])
        if identity in source_traces:
            raise ValueError("one source trace cannot be split into replicate specimens")
        source_traces.add(identity)
        _number(record, "temperature_degC")
        specimen_radius = _number(record, "specimen_radius_um")
        bath = _number(record, "bath_oxygen_mol_m3")
        radii = [point[0] for point in points]
        if specimen_radius <= 0 or bath <= 0 or len(points) < 3 or max(radii) > specimen_radius or any(b <= a for a, b in zip(radii, radii[1:])):
            raise ValueError("profile needs >=3 strictly increasing radii inside positive measured specimen radius and positive bath oxygen")
        summary.append({"specimen_id": specimen_id, "biological_unit_id": biological_id, "role": role,
                        "n_points": len(points), "specimen_radius_um": specimen_radius,
                        "bath_oxygen_mol_m3": bath, "sampled_radius_um": [radii[0], radii[-1]],
                        "sampled_oxygen_range_mol_m3": [min(c for _, c in points), max(c for _, c in points)]})
    if set(biological_roles.values()) != {"calibration", "validation"}:
        raise ValueError("independent calibration and held-out validation biological units required")
    return {"schema_version": 1, "status": "qualified_declared_profiles_for_future_comparison",
            "csv_sha256": digest, "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
            "source_files_sha256": verified, "source_url": document["source_url"],
            "license": document["license"], "attribution": document["attribution"],
            "acquisition_record": document["acquisition_record"], "specimen_records": specimens,
            "n_specimens": len(summary), "n_biological_units": len(biological_roles), "profiles": summary,
            "model_fitted": False, "biological_validation_performed": False,
            "limits": ["Intake checks declared source identity and units; it cannot authenticate measurements or prove spherical homogeneity.",
                       "No oxygen model has been fitted or scored against these profiles.",
                       "Parameter identifiability and measurement uncertainty need separate analysis."]}
