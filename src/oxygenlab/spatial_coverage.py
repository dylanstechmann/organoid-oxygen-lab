"""Describe radial sampling geometry without inventing oxygen between measurements."""

import math
from itertools import pairwise


def radial_coverage(radii_um, specimen_radius_um):
    """Geometric diagnostics, not measured volume coverage or an assay adequacy cutoff."""
    if isinstance(specimen_radius_um, bool) or not isinstance(specimen_radius_um, (float, int)) or not math.isfinite(specimen_radius_um) or specimen_radius_um <= 0:
        raise ValueError("specimen_radius_um must be positive and finite")
    radii = list(radii_um)
    if not radii or any(isinstance(r, bool) or not isinstance(r, (float, int)) or not math.isfinite(r) or r < 0 or r > specimen_radius_um for r in radii):
        raise ValueError("Sample radii must be finite and inside the measured sphere")
    if any(b <= a for a, b in pairwise(radii)):
        raise ValueError("Sample radii must be strictly increasing")
    fractions = [r / specimen_radius_um for r in radii]
    return {
        "radial_span_fraction": [fractions[0], fractions[-1]],
        "center_sample_present": radii[0] == 0,
        "surface_sample_present": radii[-1] == specimen_radius_um,
        "inward_of_first_sample_sphere_volume_fraction": fractions[0] ** 3,
        "outward_of_last_sample_sphere_volume_fraction": 1 - fractions[-1] ** 3,
        "largest_internal_radial_gap_fraction": max((b-a for a, b in pairwise(fractions)), default=None),
        "interpretation": "Volume fractions are geometric regions outside the radial sample span in a declared sphere, not directly measured tissue volumes. Point samples do not measure every location between them; no oxygen, viability or minimum sampling adequacy is inferred."
    }
