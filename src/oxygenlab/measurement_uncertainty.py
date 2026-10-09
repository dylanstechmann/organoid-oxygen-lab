"""Retain source-declared uncertainty without inventing error bars or propagating them."""

import math

QUANTITIES = {
    "specimen_radius": "um",
    "radial_position": "um",
    "oxygen": "mol/m3",
    "bath_oxygen": "mol/m3",
}
KINDS = {"standard_deviation", "standard_uncertainty", "expanded_uncertainty", "absolute_bound"}


def uncertainty_declarations(record, verified_sources):
    """Optional per-specimen annotations; omitted quantities remain unreported.

    Values describe the scope stated in the source record, not every CSV point.
    Source hashes bind files, but do not authenticate their uncertainty claims.
    """
    declarations = record.get("measurement_uncertainty", {})
    if not isinstance(declarations, dict) or set(declarations) - set(QUANTITIES):
        raise ValueError("measurement_uncertainty must contain only known quantities")
    result = {}
    for quantity, unit in QUANTITIES.items():
        entry = declarations.get(quantity)
        if quantity not in declarations:
            result[quantity] = {"status": "not_reported_in_manifest", "value": None, "unit": unit}
            continue
        fields = {"kind", "value", "unit", "source_file", "source_record", "scope"}
        if not isinstance(entry, dict) or set(entry) - (fields | {"coverage_factor"}) or not fields <= set(entry):
            raise ValueError(f"{quantity}: uncertainty declaration needs explicit fields")
        if entry["kind"] not in KINDS or entry["unit"] != unit:
            raise ValueError(f"{quantity}: unknown uncertainty kind or noncanonical unit")
        value = entry["value"]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError(f"{quantity}: uncertainty value must be finite and nonnegative")
        for key in ("source_file", "source_record", "scope"):
            if not isinstance(entry[key], str) or not entry[key].strip():
                raise ValueError(f"{quantity}: uncertainty needs {key}")
        if entry["source_file"] not in verified_sources:
            raise ValueError(f"{quantity}: unknown uncertainty source_file")
        if entry["kind"] == "expanded_uncertainty":
            factor = entry.get("coverage_factor")
            if isinstance(factor, bool) or not isinstance(factor, (int, float)) or not math.isfinite(factor) or factor <= 0:
                raise ValueError(f"{quantity}: expanded uncertainty needs positive coverage_factor")
        elif "coverage_factor" in entry:
            raise ValueError(f"{quantity}: coverage_factor is only valid for expanded uncertainty")
        result[quantity] = {"status": "source_declared_not_authenticated", **entry,
                            "source_sha256": verified_sources[entry["source_file"]]}
    return {"quantities": result, "propagated_to_model": False,
            "confidence_interval_inferred": False,
            "limits": "SD, standard uncertainty, expanded uncertainty and bounds retain their declared meanings. No independence, distribution, pointwise applicability or adequate precision is inferred."}
