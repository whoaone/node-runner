"""v5.0.0 item 8: built-in MMPDS-traceable material library.

This module ships a curated starter set of aerospace metallics and a
couple of carbon-fiber composites with values cross-referenced against
public engineering handbooks (MIL-HDBK-5J / matweb / vendor datasheets)
that track MMPDS-16 within published rounding. Each entry carries a
``source`` field that names both the public reference and the MMPDS
section so a value-by-value reconciliation pass remains possible once a
machine-readable MMPDS source is available.

Units (all entries):

    Stress:   psi
    Density:  lbm/in^3
    Thermal:  1/degF
    Temp:     degF

These are US customary, matching the user's aerospace Nastran workflow
(Femap is unitless; we do not convert). The picker dialog surfaces these
units on every load to force the user to confirm their model is in
matching units before populating MAT1/MAT8 fields.
"""

from __future__ import annotations

from typing import Any

MATERIAL_LIBRARY_VERSION = 1

SOURCE_UNITS = {
    "stress":             "psi",
    "density":            "lbm/in^3",
    "thermal_expansion":  "1/degF",
    "temperature":        "degF",
}

SOURCE_UNITS_TEXT = (
    "psi (stress) / lbm/in^3 (density) / 1/degF (thermal) / degF (Tref)"
)

# Curated list. Each entry uses the following keys depending on type:
#
# MAT1:  name, type='MAT1', source, basis, E, G, nu, rho, alpha, Tref
# MAT8:  name, type='MAT8', source, basis, E1, E2, nu12, G12, G1z, G2z,
#                            rho, a1, a2, Tref
#
# Values are deliberately conservative single-point references. Real
# design analyses must consult the live MMPDS edition - this library is
# a productivity convenience, not an authoritative datasheet.
MATERIALS_US: list[dict[str, Any]] = [
    # ---------- Aluminums ----------
    {"name": "Al 2024-T3 sheet (<=0.249 in)",
     "type": "MAT1", "category": "Aluminum",
     "source": "MMPDS-16 Section 3.2.3 (cross-ref MIL-HDBK-5J Table 3.2.3.0(b1))",
     "basis": "S-basis, RT",
     "E": 10.5e6, "G": 4.00e6, "nu": 0.33,
     "rho": 0.100, "alpha": 12.6e-6, "Tref": 70.0},
    {"name": "Al 2024-T351 plate (0.250-2.000 in)",
     "type": "MAT1", "category": "Aluminum",
     "source": "MMPDS-16 Section 3.2.3",
     "basis": "S-basis, RT",
     "E": 10.5e6, "G": 4.00e6, "nu": 0.33,
     "rho": 0.100, "alpha": 12.6e-6, "Tref": 70.0},
    {"name": "Al 6061-T6 sheet/plate",
     "type": "MAT1", "category": "Aluminum",
     "source": "MMPDS-16 Section 3.6.2 (cross-ref MIL-HDBK-5J Table 3.6.2.0(b))",
     "basis": "S-basis, RT",
     "E": 9.9e6, "G": 3.80e6, "nu": 0.33,
     "rho": 0.098, "alpha": 13.0e-6, "Tref": 70.0},
    {"name": "Al 7050-T7451 plate (1.001-2.000 in)",
     "type": "MAT1", "category": "Aluminum",
     "source": "MMPDS-16 Section 3.7.5",
     "basis": "S-basis, RT",
     "E": 10.3e6, "G": 3.90e6, "nu": 0.33,
     "rho": 0.102, "alpha": 12.4e-6, "Tref": 70.0},
    {"name": "Al 7075-T6 sheet (<=0.249 in)",
     "type": "MAT1", "category": "Aluminum",
     "source": "MMPDS-16 Section 3.7.4 (cross-ref MIL-HDBK-5J Table 3.7.4.0(b1))",
     "basis": "S-basis, RT",
     "E": 10.3e6, "G": 3.90e6, "nu": 0.33,
     "rho": 0.101, "alpha": 12.9e-6, "Tref": 70.0},
    {"name": "Al 7075-T73 forging",
     "type": "MAT1", "category": "Aluminum",
     "source": "MMPDS-16 Section 3.7.4",
     "basis": "S-basis, RT",
     "E": 10.3e6, "G": 3.90e6, "nu": 0.33,
     "rho": 0.101, "alpha": 12.9e-6, "Tref": 70.0},
    {"name": "Al 2219-T87 sheet/plate",
     "type": "MAT1", "category": "Aluminum",
     "source": "MMPDS-16 Section 3.2.7",
     "basis": "S-basis, RT",
     "E": 10.6e6, "G": 4.00e6, "nu": 0.33,
     "rho": 0.103, "alpha": 12.4e-6, "Tref": 70.0},
    {"name": "Al-Li 2050-T84 plate",
     "type": "MAT1", "category": "Aluminum-Lithium",
     "source": "Vendor datasheet (Constellium); cross-ref MMPDS-16 Section 3.2",
     "basis": "Typical, RT",
     "E": 11.4e6, "G": 4.30e6, "nu": 0.33,
     "rho": 0.0955, "alpha": 12.4e-6, "Tref": 70.0},
    # ---------- Titaniums ----------
    {"name": "Ti-6Al-4V annealed sheet (Grade 5)",
     "type": "MAT1", "category": "Titanium",
     "source": "MMPDS-16 Section 5.4.1",
     "basis": "S-basis, RT",
     "E": 16.0e6, "G": 6.20e6, "nu": 0.31,
     "rho": 0.160, "alpha": 4.8e-6, "Tref": 70.0},
    {"name": "Ti-6Al-4V STA bar/forging",
     "type": "MAT1", "category": "Titanium",
     "source": "MMPDS-16 Section 5.4.1",
     "basis": "S-basis, RT",
     "E": 16.4e6, "G": 6.20e6, "nu": 0.31,
     "rho": 0.160, "alpha": 4.8e-6, "Tref": 70.0},
    {"name": "CP Ti Grade 2 sheet",
     "type": "MAT1", "category": "Titanium",
     "source": "MMPDS-16 Section 5.2.1",
     "basis": "S-basis, RT",
     "E": 14.9e6, "G": 6.50e6, "nu": 0.31,
     "rho": 0.163, "alpha": 5.0e-6, "Tref": 70.0},
    # ---------- Steels ----------
    {"name": "4130 normalized sheet",
     "type": "MAT1", "category": "Steel",
     "source": "MMPDS-16 Section 2.3.1 (cross-ref MIL-HDBK-5J Table 2.3.1.0(b))",
     "basis": "S-basis, RT",
     "E": 29.0e6, "G": 11.0e6, "nu": 0.32,
     "rho": 0.283, "alpha": 7.1e-6, "Tref": 70.0},
    {"name": "4340 HT 200 ksi bar",
     "type": "MAT1", "category": "Steel",
     "source": "MMPDS-16 Section 2.3.1",
     "basis": "S-basis, RT",
     "E": 29.0e6, "G": 11.0e6, "nu": 0.32,
     "rho": 0.283, "alpha": 7.1e-6, "Tref": 70.0},
    {"name": "17-4 PH H1025 bar/plate",
     "type": "MAT1", "category": "Steel (PH)",
     "source": "MMPDS-16 Section 2.6.9",
     "basis": "S-basis, RT",
     "E": 28.5e6, "G": 11.2e6, "nu": 0.27,
     "rho": 0.282, "alpha": 6.0e-6, "Tref": 70.0},
    {"name": "15-5 PH H1025 bar/plate",
     "type": "MAT1", "category": "Steel (PH)",
     "source": "MMPDS-16 Section 2.6.7",
     "basis": "S-basis, RT",
     "E": 28.5e6, "G": 11.2e6, "nu": 0.27,
     "rho": 0.283, "alpha": 6.0e-6, "Tref": 70.0},
    {"name": "A286 STA bar (AMS 5732)",
     "type": "MAT1", "category": "Steel (austenitic)",
     "source": "MMPDS-16 Section 6.2.1",
     "basis": "S-basis, RT",
     "E": 29.1e6, "G": 11.2e6, "nu": 0.31,
     "rho": 0.286, "alpha": 9.2e-6, "Tref": 70.0},
    {"name": "300M forging (HT 280 ksi)",
     "type": "MAT1", "category": "Steel (low-alloy)",
     "source": "MMPDS-16 Section 2.3.1.13",
     "basis": "S-basis, RT",
     "E": 29.0e6, "G": 11.2e6, "nu": 0.32,
     "rho": 0.283, "alpha": 6.5e-6, "Tref": 70.0},
    # ---------- Nickels / superalloys ----------
    {"name": "Inconel 718 STA bar/forging",
     "type": "MAT1", "category": "Nickel superalloy",
     "source": "MMPDS-16 Section 6.3.5",
     "basis": "S-basis, RT",
     "E": 29.6e6, "G": 11.5e6, "nu": 0.29,
     "rho": 0.297, "alpha": 7.1e-6, "Tref": 70.0},
    {"name": "Inconel 625 annealed sheet",
     "type": "MAT1", "category": "Nickel superalloy",
     "source": "MMPDS-16 Section 6.3.3",
     "basis": "S-basis, RT",
     "E": 30.0e6, "G": 11.5e6, "nu": 0.29,
     "rho": 0.305, "alpha": 7.1e-6, "Tref": 70.0},
    {"name": "Hastelloy X solution-annealed sheet",
     "type": "MAT1", "category": "Nickel superalloy",
     "source": "MMPDS-16 Section 6.3.2",
     "basis": "S-basis, RT",
     "E": 28.6e6, "G": 11.0e6, "nu": 0.32,
     "rho": 0.297, "alpha": 7.7e-6, "Tref": 70.0},
    {"name": "Waspaloy STA bar",
     "type": "MAT1", "category": "Nickel superalloy",
     "source": "MMPDS-16 Section 6.3.7",
     "basis": "S-basis, RT",
     "E": 30.4e6, "G": 11.7e6, "nu": 0.30,
     "rho": 0.296, "alpha": 7.2e-6, "Tref": 70.0},
    # ---------- Magnesium / Beryllium ----------
    {"name": "Mg AZ31B-H24 sheet",
     "type": "MAT1", "category": "Magnesium",
     "source": "MMPDS-16 Section 4.2.1",
     "basis": "S-basis, RT",
     "E": 6.5e6, "G": 2.4e6, "nu": 0.35,
     "rho": 0.0639, "alpha": 14.5e-6, "Tref": 70.0},
    # ---------- Carbon-fiber composites (MAT8 orthotropic ply) ----------
    {"name": "IM7/8552 UD ply (Hexcel)",
     "type": "MAT8", "category": "Composite (carbon UD)",
     "source": "NIAR NCAMP NPS-81226 (B-basis, RT/dry); cross-ref Hexcel datasheet",
     "basis": "B-basis, RT/dry",
     "E1": 23.8e6, "E2": 1.20e6, "nu12": 0.32,
     "G12": 0.70e6, "G1z": 0.70e6, "G2z": 0.50e6,
     "rho": 0.0578,
     "a1": -0.10e-6, "a2": 16.0e-6, "Tref": 70.0},
    {"name": "T800S/3900-2 UD ply (Toray)",
     "type": "MAT8", "category": "Composite (carbon UD)",
     "source": "Toray T800S/3900-2 datasheet (typical, RT/dry)",
     "basis": "Typical, RT/dry",
     "E1": 22.0e6, "E2": 1.20e6, "nu12": 0.32,
     "G12": 0.65e6, "G1z": 0.65e6, "G2z": 0.50e6,
     "rho": 0.0563,
     "a1": -0.10e-6, "a2": 16.5e-6, "Tref": 70.0},
]


def materials_by_category() -> dict[str, list[dict[str, Any]]]:
    """Group MATERIALS_US by category for an alphabetized picker view."""
    out: dict[str, list[dict[str, Any]]] = {}
    for m in MATERIALS_US:
        out.setdefault(m.get("category", "Other"), []).append(m)
    return out
