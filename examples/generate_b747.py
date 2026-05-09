"""Generate a Boeing 747-class BDF for Node Runner stress testing.

Builds primary + secondary structure: fuselage (skin / frames / stringers /
longerons), wings, vertical and horizontal stabilizers, floor structure
(panels, lateral cross-beams, longitudinal beams, stanchions), engine
attachment beams. No landing gear.

Streaming writer: at 1M+ nodes, holding everything in pyNastran's BDF
object would consume 1-2 GB of Python memory and take several minutes.
We keep the BDF object only for the small materials-and-properties
section, then write GRID / CQUAD4 / CBAR cards directly via pyNastran's
`print_card_8` helper.

Materials (MAT1):
  1: Aluminum 2024-T3 (skin)
  2: Aluminum 7075-T6 (structural members)
  3: Titanium Ti-6Al-4V (heavy-load fittings - declared, unused)

Properties:
  PSHELL 1: Fuselage skin (1.6 mm)        PBEAM 10: Fuselage frame
  PSHELL 2: Wing upper skin (4.0 mm)      PBEAM 11: Fuselage stringer
  PSHELL 3: Wing lower skin (3.5 mm)      PBEAM 12: Wing main spar
  PSHELL 4: Vertical stab skin (2.0 mm)   PBEAM 13: Wing rib
  PSHELL 5: Horizontal stab skin (1.8 mm) PBEAM 14: Wing stringer
  PSHELL 6: Floor panel (1.0 mm)          PBEAM 15: Empennage spar
                                          PBEAM 16: Floor cross-beam
                                          PBEAM 17: Longeron
                                          PBEAM 18: Engine pylon

Usage:
  python examples/generate_b747.py [output.bdf]
"""

from __future__ import annotations

import math
import sys
import time
from io import StringIO
from pathlib import Path

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.field_writer_8 import print_card_8


# ---------- 747-400 dimensions (approximate, meters) ----------
FUSELAGE_LEN = 68.0
FUSELAGE_RADIUS = 3.25
NOSE_CONE_LEN = 8.0
TAIL_CONE_LEN = 8.0

WING_ROOT_X = 24.0
WING_ROOT_CHORD = 12.0
WING_TIP_CHORD = 4.0
WING_HALF_SPAN = 28.75
WING_SWEEP_DEG = 37.0
WING_DIHEDRAL_DEG = 7.0
WING_ROOT_Z = -1.5

VS_ROOT_X = 56.0
VS_ROOT_CHORD = 8.0
VS_TIP_CHORD = 2.5
VS_HEIGHT = 9.0
VS_SWEEP_DEG = 35.0
VS_BASE_Z = FUSELAGE_RADIUS

HS_ROOT_X = 60.0
HS_ROOT_CHORD = 6.0
HS_TIP_CHORD = 2.0
HS_HALF_SPAN = 11.0
HS_SWEEP_DEG = 30.0
HS_Z = VS_BASE_Z + VS_HEIGHT * 0.95

# ---------- Discretization (sized for ~1M nodes) ----------
FUSELAGE_FRAMES = 1100
FUSELAGE_STRINGERS = 250
WING_RIBS = 280
WING_CHORD_NODES = 380
VS_RIBS = 200
VS_CHORD_NODES = 240
HS_RIBS = 180
HS_CHORD_NODES = 220
FLOOR_X_STATIONS = 700
FLOOR_Y_STATIONS = 100

# Beam thinning - skins are full-density, structural beams are sampled
# so the file stays in the 100-200 MiB range and parses in reasonable time.
FUSELAGE_FRAME_BEAM_EVERY_N = 4
FUSELAGE_STRINGER_BEAM_EVERY_N = 5
WING_STRINGER_EVERY_N = 1
FLOOR_LONG_BEAM_EVERY_N = 4
FLOOR_LATERAL_BEAM_EVERY_N = 2
STANCHION_EVERY_N = 4

# ---------- ID base offsets (room for 1.5M+ ids per region) ----------
FUSELAGE_NID_BASE = 1
LEFT_WING_NID_BASE = 2_000_000
RIGHT_WING_NID_BASE = 4_000_000
VS_NID_BASE = 6_000_000
LEFT_HS_NID_BASE = 7_000_000
RIGHT_HS_NID_BASE = 8_000_000
ENGINE_NID_BASE = 9_000_000
FLOOR_NID_BASE = 9_500_000

FUSELAGE_EID_BASE = 1
LEFT_WING_EID_BASE = 2_000_000
RIGHT_WING_EID_BASE = 4_000_000
VS_EID_BASE = 6_000_000
LEFT_HS_EID_BASE = 7_000_000
RIGHT_HS_EID_BASE = 8_000_000
ENGINE_EID_BASE = 9_000_000
FLOOR_EID_BASE = 9_500_000


class Counter:
    """Cumulative tally of grids and elements for end-of-run reporting."""
    __slots__ = ("nodes", "elements")

    def __init__(self) -> None:
        self.nodes = 0
        self.elements = 0


def w_grid(fh, counter, nid, x, y, z):
    fh.write(print_card_8(['GRID', nid, 0, float(x), float(y), float(z)]))
    counter.nodes += 1


def w_cquad4(fh, counter, eid, pid, n1, n2, n3, n4):
    fh.write(print_card_8(['CQUAD4', eid, pid, n1, n2, n3, n4]))
    counter.elements += 1


def w_cbar(fh, counter, eid, pid, n1, n2, x1, x2, x3):
    fh.write(print_card_8(['CBAR', eid, pid, n1, n2,
                           float(x1), float(x2), float(x3)]))
    counter.elements += 1


# ---------- Geometry helpers ----------

def fuselage_radius_at(x: float) -> float:
    if x < NOSE_CONE_LEN:
        ratio = max(x / NOSE_CONE_LEN, 0.0)
        return 0.30 + (FUSELAGE_RADIUS - 0.30) * math.sqrt(ratio)
    if x > FUSELAGE_LEN - TAIL_CONE_LEN:
        ratio = (x - (FUSELAGE_LEN - TAIL_CONE_LEN)) / TAIL_CONE_LEN
        return FUSELAGE_RADIUS - (FUSELAGE_RADIUS - 0.50) * ratio
    return FUSELAGE_RADIUS


def naca_thickness(xi: float, max_t: float = 0.12) -> float:
    if xi <= 0.0 or xi >= 1.0:
        return 0.0
    return 5.0 * max_t * (
        0.2969 * math.sqrt(xi)
        - 0.1260 * xi
        - 0.3516 * xi ** 2
        + 0.2843 * xi ** 3
        - 0.1036 * xi ** 4
    )


# ---------- Materials & properties via pyNastran (small section) ----------

def setup_materials_and_properties() -> str:
    """Return BDF text for materials and properties.

    Uses pyNastran's BDF object to get properly formatted multi-line PBEAM
    cards (which are tedious to format by hand), then strips out just the
    bulk-data section so we can splice it into our streaming output.
    """
    model = BDF(debug=False)
    model.add_mat1(1, E=72.4e9, G=None, nu=0.33, rho=2780.0,
                   comment='Aluminum 2024-T3 (skin)')
    model.add_mat1(2, E=71.7e9, G=None, nu=0.33, rho=2810.0,
                   comment='Aluminum 7075-T6 (structural members)')
    model.add_mat1(3, E=113.8e9, G=None, nu=0.342, rho=4430.0,
                   comment='Titanium Ti-6Al-4V (heavy-load fittings)')

    model.add_pshell(1, mid1=1, t=0.0016, comment='Fuselage skin (1.6 mm Al 2024-T3)')
    model.add_pshell(2, mid1=1, t=0.0040, comment='Wing upper skin (4.0 mm Al 2024-T3)')
    model.add_pshell(3, mid1=1, t=0.0035, comment='Wing lower skin (3.5 mm Al 2024-T3)')
    model.add_pshell(4, mid1=1, t=0.0020, comment='Vertical stab skin (2.0 mm)')
    model.add_pshell(5, mid1=1, t=0.0018, comment='Horizontal stab skin (1.8 mm)')
    model.add_pshell(6, mid1=1, t=0.0010, comment='Floor panel (1.0 mm)')

    def pbeam(pid, A, I1, I2, J, comment):
        model.add_pbeam(
            pid, 2, xxb=[0.0], so=['C'],
            area=[A], i1=[I1], i2=[I2], i12=[0.0], j=[J],
            comment=comment,
        )

    pbeam(10, 1.5e-4, 4.0e-8, 1.5e-8, 1.0e-9,  'Fuselage frame (Z-section)')
    pbeam(11, 5.0e-5, 5.0e-9, 5.0e-9, 1.0e-10, 'Fuselage stringer (hat-section)')
    pbeam(12, 8.0e-4, 5.0e-7, 5.0e-7, 1.0e-7,  'Wing main spar (I-section)')
    pbeam(13, 4.0e-4, 1.0e-7, 1.0e-7, 5.0e-8,  'Wing rib (truss)')
    pbeam(14, 1.0e-4, 1.0e-8, 1.0e-8, 5.0e-9,  'Wing stringer')
    pbeam(15, 2.0e-4, 5.0e-8, 5.0e-8, 2.0e-8,  'Empennage spar')
    pbeam(16, 3.0e-4, 1.0e-7, 1.0e-7, 5.0e-8,  'Floor cross-beam')
    pbeam(17, 1.5e-4, 3.0e-8, 3.0e-8, 1.0e-8,  'Longeron')
    pbeam(18, 1.0e-3, 1.0e-6, 1.0e-6, 5.0e-7,  'Engine pylon')

    buf = StringIO()
    model.write_bdf(buf, size=8, is_double=False, close=False)
    raw = buf.getvalue()

    # Extract the bulk-data section. pyNastran may emit either a full deck
    # (with BEGIN BULK / ENDDATA) or a punch deck (cards only) depending on
    # what's present. Handle both: if BEGIN BULK is in the output, use it
    # as the start sentinel; otherwise treat every non-pyNastran-header
    # line as bulk.
    has_begin_bulk = any(
        line.lstrip().upper().startswith('BEGIN BULK') for line in raw.splitlines()
    )
    out_lines = []
    in_bulk = not has_begin_bulk
    for line in raw.splitlines():
        upper = line.lstrip().upper()
        if upper.startswith('BEGIN BULK'):
            in_bulk = True
            continue
        if upper.startswith('ENDDATA'):
            break
        if not in_bulk:
            continue
        # Skip pyNastran's auto-generated banner comments (they reference
        # punch=, encoding=, etc. and are not useful inside our deck).
        stripped = line.strip()
        if stripped.startswith('$pyNastran'):
            continue
        out_lines.append(line)
    text = '\n'.join(out_lines).rstrip()
    if not text.strip():
        raise RuntimeError(
            "setup_materials_and_properties produced no bulk content; "
            "check pyNastran's write_bdf output format."
        )
    return text + '\n'


# ---------- Streaming structural builders ----------

def write_fuselage(fh, counter, eid_start: int) -> int:
    nid_base = FUSELAGE_NID_BASE
    cols = FUSELAGE_STRINGERS

    def nid(fi, sj):
        return nid_base + fi * cols + sj

    # Grids
    for fi in range(FUSELAGE_FRAMES):
        x = fi * FUSELAGE_LEN / (FUSELAGE_FRAMES - 1)
        r = fuselage_radius_at(x)
        for sj in range(cols):
            theta = 2.0 * math.pi * sj / cols
            y = r * math.cos(theta)
            z = r * math.sin(theta)
            w_grid(fh, counter, nid(fi, sj), x, y, z)

    eid = eid_start

    # Skin (full density)
    for fi in range(FUSELAGE_FRAMES - 1):
        for sj in range(cols):
            sj1 = (sj + 1) % cols
            w_cquad4(
                fh, counter, eid, 1,
                nid(fi, sj), nid(fi + 1, sj),
                nid(fi + 1, sj1), nid(fi, sj1),
            )
            eid += 1

    # Frames as CBARs (every Nth frame station, all stringers around the ring)
    for fi in range(0, FUSELAGE_FRAMES, FUSELAGE_FRAME_BEAM_EVERY_N):
        for sj in range(cols):
            sj1 = (sj + 1) % cols
            w_cbar(
                fh, counter, eid, 10,
                nid(fi, sj), nid(fi, sj1),
                1.0, 0.0, 0.0,
            )
            eid += 1

    # Stringers as CBARs (every Nth stringer column)
    for sj in range(0, cols, FUSELAGE_STRINGER_BEAM_EVERY_N):
        theta = 2.0 * math.pi * sj / cols
        cy = math.cos(theta)
        cz = math.sin(theta)
        for fi in range(FUSELAGE_FRAMES - 1):
            w_cbar(
                fh, counter, eid, 11,
                nid(fi, sj), nid(fi + 1, sj),
                0.0, cy, cz,
            )
            eid += 1

    # Longerons - 4 distinguished stringer rows with PBEAM 17
    for sj in (0, cols // 4, cols // 2, 3 * cols // 4):
        theta = 2.0 * math.pi * sj / cols
        cy = math.cos(theta)
        cz = math.sin(theta)
        for fi in range(FUSELAGE_FRAMES - 1):
            w_cbar(
                fh, counter, eid, 17,
                nid(fi, sj), nid(fi + 1, sj),
                0.0, cy, cz,
            )
            eid += 1

    return eid


def write_wing(
    fh, counter, side: str, nid_base: int, eid_start: int,
) -> int:
    sgn = +1.0 if side == 'right' else -1.0
    sweep = math.radians(WING_SWEEP_DEG)
    dihedral = math.radians(WING_DIHEDRAL_DEG)
    y_root = sgn * FUSELAGE_RADIUS
    cols = WING_CHORD_NODES

    def nid(si, cj, surf):
        return nid_base + si * (cols * 2) + surf * cols + cj

    for si in range(WING_RIBS):
        eta = si / (WING_RIBS - 1)
        chord = (1 - eta) * WING_ROOT_CHORD + eta * WING_TIP_CHORD
        x_le = WING_ROOT_X + eta * WING_HALF_SPAN * math.tan(sweep)
        y_here = y_root + sgn * eta * WING_HALF_SPAN
        z_le = WING_ROOT_Z + eta * WING_HALF_SPAN * math.tan(dihedral)
        for cj in range(cols):
            xi = cj / (cols - 1)
            x = x_le + xi * chord
            half_t = naca_thickness(xi) * chord
            for surf, dz in ((0, +half_t), (1, -half_t)):
                w_grid(fh, counter, nid(si, cj, surf), x, y_here, z_le + dz)

    eid = eid_start

    # Upper / lower skin
    for si in range(WING_RIBS - 1):
        for cj in range(cols - 1):
            w_cquad4(
                fh, counter, eid, 2,
                nid(si, cj, 0), nid(si + 1, cj, 0),
                nid(si + 1, cj + 1, 0), nid(si, cj + 1, 0),
            )
            eid += 1
    for si in range(WING_RIBS - 1):
        for cj in range(cols - 1):
            w_cquad4(
                fh, counter, eid, 3,
                nid(si, cj, 1), nid(si, cj + 1, 1),
                nid(si + 1, cj + 1, 1), nid(si + 1, cj, 1),
            )
            eid += 1

    # 3 spars (front / mid / rear)
    for cj in (int(0.15 * (cols - 1)), int(0.50 * (cols - 1)), int(0.85 * (cols - 1))):
        for si in range(WING_RIBS - 1):
            w_cbar(fh, counter, eid, 12,
                   nid(si, cj, 0), nid(si + 1, cj, 0),
                   0.0, 0.0, 1.0); eid += 1
            w_cbar(fh, counter, eid, 12,
                   nid(si, cj, 1), nid(si + 1, cj, 1),
                   0.0, 0.0, 1.0); eid += 1
        for si in range(WING_RIBS):
            w_cbar(fh, counter, eid, 12,
                   nid(si, cj, 0), nid(si, cj, 1),
                   1.0, 0.0, 0.0); eid += 1

    # Ribs (chordwise truss segments)
    rib_chord_nodes = [
        int(0.05 * (cols - 1)), int(0.25 * (cols - 1)), int(0.5 * (cols - 1)),
        int(0.75 * (cols - 1)), int(0.95 * (cols - 1)),
    ]
    for si in range(WING_RIBS):
        for k in range(len(rib_chord_nodes) - 1):
            cj0 = rib_chord_nodes[k]
            cj1 = rib_chord_nodes[k + 1]
            w_cbar(fh, counter, eid, 13,
                   nid(si, cj0, 0), nid(si, cj1, 0),
                   0.0, 1.0, 0.0); eid += 1
            w_cbar(fh, counter, eid, 13,
                   nid(si, cj0, 1), nid(si, cj1, 1),
                   0.0, 1.0, 0.0); eid += 1

    # Wing stringers
    stringer_chord = [int(p * (cols - 1)) for p in (0.20, 0.35, 0.65, 0.80)]
    for cj in stringer_chord:
        for surf in (0, 1):
            for si in range(0, WING_RIBS - 1, WING_STRINGER_EVERY_N):
                w_cbar(fh, counter, eid, 14,
                       nid(si, cj, surf), nid(si + 1, cj, surf),
                       0.0, 0.0, 1.0); eid += 1

    return eid


def write_vertical_stab(fh, counter, eid_start: int) -> int:
    nid_base = VS_NID_BASE
    sweep = math.radians(VS_SWEEP_DEG)
    cols = VS_CHORD_NODES

    def nid(si, cj, surf):
        return nid_base + si * (cols * 2) + surf * cols + cj

    for si in range(VS_RIBS):
        eta = si / (VS_RIBS - 1)
        chord = (1 - eta) * VS_ROOT_CHORD + eta * VS_TIP_CHORD
        x_le = VS_ROOT_X + eta * VS_HEIGHT * math.tan(sweep)
        z_here = VS_BASE_Z + eta * VS_HEIGHT
        for cj in range(cols):
            xi = cj / (cols - 1)
            x = x_le + xi * chord
            half_t = naca_thickness(xi) * chord
            for surf, dy in ((0, +half_t), (1, -half_t)):
                w_grid(fh, counter, nid(si, cj, surf), x, dy, z_here)

    eid = eid_start
    for surf in (0, 1):
        for si in range(VS_RIBS - 1):
            for cj in range(cols - 1):
                if surf == 0:
                    quad = (nid(si, cj, 0), nid(si + 1, cj, 0),
                            nid(si + 1, cj + 1, 0), nid(si, cj + 1, 0))
                else:
                    quad = (nid(si, cj, 1), nid(si, cj + 1, 1),
                            nid(si + 1, cj + 1, 1), nid(si + 1, cj, 1))
                w_cquad4(fh, counter, eid, 4, *quad); eid += 1

    for cj in (int(0.20 * (cols - 1)), int(0.70 * (cols - 1))):
        for si in range(VS_RIBS - 1):
            w_cbar(fh, counter, eid, 15,
                   nid(si, cj, 0), nid(si + 1, cj, 0),
                   0.0, 0.0, 1.0); eid += 1
            w_cbar(fh, counter, eid, 15,
                   nid(si, cj, 1), nid(si + 1, cj, 1),
                   0.0, 0.0, 1.0); eid += 1

    return eid


def write_horizontal_stab(
    fh, counter, side: str, nid_base: int, eid_start: int,
) -> int:
    sgn = +1.0 if side == 'right' else -1.0
    sweep = math.radians(HS_SWEEP_DEG)
    cols = HS_CHORD_NODES

    def nid(si, cj, surf):
        return nid_base + si * (cols * 2) + surf * cols + cj

    for si in range(HS_RIBS):
        eta = si / (HS_RIBS - 1)
        chord = (1 - eta) * HS_ROOT_CHORD + eta * HS_TIP_CHORD
        x_le = HS_ROOT_X + eta * HS_HALF_SPAN * math.tan(sweep)
        y_here = sgn * eta * HS_HALF_SPAN
        for cj in range(cols):
            xi = cj / (cols - 1)
            x = x_le + xi * chord
            half_t = naca_thickness(xi) * chord
            for surf, dz in ((0, +half_t), (1, -half_t)):
                w_grid(fh, counter, nid(si, cj, surf), x, y_here, HS_Z + dz)

    eid = eid_start
    for surf in (0, 1):
        for si in range(HS_RIBS - 1):
            for cj in range(cols - 1):
                if surf == 0:
                    quad = (nid(si, cj, 0), nid(si + 1, cj, 0),
                            nid(si + 1, cj + 1, 0), nid(si, cj + 1, 0))
                else:
                    quad = (nid(si, cj, 1), nid(si, cj + 1, 1),
                            nid(si + 1, cj + 1, 1), nid(si + 1, cj, 1))
                w_cquad4(fh, counter, eid, 5, *quad); eid += 1

    for cj in (int(0.25 * (cols - 1)), int(0.75 * (cols - 1))):
        for si in range(HS_RIBS - 1):
            w_cbar(fh, counter, eid, 15,
                   nid(si, cj, 0), nid(si + 1, cj, 0),
                   0.0, 0.0, 1.0); eid += 1

    return eid


def write_floor(fh, counter, eid_start: int) -> int:
    nid_base = FLOOR_NID_BASE
    floor_z = -1.0

    def nid(xi, yi):
        return nid_base + xi * FLOOR_Y_STATIONS + yi

    x_min = NOSE_CONE_LEN + 1.0
    x_max = FUSELAGE_LEN - TAIL_CONE_LEN - 1.0
    y_half = FUSELAGE_RADIUS * 0.85

    floor_x_coords = []
    for xi in range(FLOOR_X_STATIONS):
        x = x_min + xi * (x_max - x_min) / (FLOOR_X_STATIONS - 1)
        floor_x_coords.append(x)
        for yi in range(FLOOR_Y_STATIONS):
            y = -y_half + yi * (2 * y_half) / (FLOOR_Y_STATIONS - 1)
            w_grid(fh, counter, nid(xi, yi), x, y, floor_z)

    eid = eid_start

    # Floor panels
    for xi in range(FLOOR_X_STATIONS - 1):
        for yi in range(FLOOR_Y_STATIONS - 1):
            w_cquad4(fh, counter, eid, 6,
                     nid(xi, yi), nid(xi + 1, yi),
                     nid(xi + 1, yi + 1), nid(xi, yi + 1)); eid += 1

    # Lateral cross-beams
    for xi in range(0, FLOOR_X_STATIONS, FLOOR_LATERAL_BEAM_EVERY_N):
        for yi in range(FLOOR_Y_STATIONS - 1):
            w_cbar(fh, counter, eid, 16,
                   nid(xi, yi), nid(xi, yi + 1),
                   1.0, 0.0, 0.0); eid += 1

    # Longitudinal beams (every Nth y-station to keep count manageable)
    for yi in range(0, FLOOR_Y_STATIONS, FLOOR_LONG_BEAM_EVERY_N):
        for xi in range(FLOOR_X_STATIONS - 1):
            w_cbar(fh, counter, eid, 16,
                   nid(xi, yi), nid(xi + 1, yi),
                   0.0, 1.0, 0.0); eid += 1

    # Stanchions: vertical struts from floor outboard edges to nearest
    # fuselage frame node. Mapped from floor edge (y, z) to closest stringer.
    cols = FUSELAGE_STRINGERS

    def fuselage_nid(frame_i, stringer_j):
        return FUSELAGE_NID_BASE + frame_i * cols + stringer_j

    fuselage_x_stations = [
        i * FUSELAGE_LEN / (FUSELAGE_FRAMES - 1) for i in range(FUSELAGE_FRAMES)
    ]

    def nearest_frame(x):
        return min(
            range(len(fuselage_x_stations)),
            key=lambda i: abs(fuselage_x_stations[i] - x),
        )

    def nearest_stringer(target_rad):
        target = target_rad % (2 * math.pi)
        best_i = 0
        best_diff = float('inf')
        for s in range(cols):
            ang = 2.0 * math.pi * s / cols
            diff = abs(((ang - target + math.pi) % (2 * math.pi)) - math.pi)
            if diff < best_diff:
                best_diff = diff
                best_i = s
        return best_i

    right_stringer = nearest_stringer(math.atan2(floor_z, +y_half))
    left_stringer = nearest_stringer(math.atan2(floor_z, -y_half))

    for xi in range(0, FLOOR_X_STATIONS, STANCHION_EVERY_N):
        fx = floor_x_coords[xi]
        frame_i = nearest_frame(fx)
        # Right edge stanchion
        w_cbar(fh, counter, eid, 16,
               nid(xi, FLOOR_Y_STATIONS - 1),
               fuselage_nid(frame_i, right_stringer),
               1.0, 0.0, 0.0); eid += 1
        # Left edge stanchion
        w_cbar(fh, counter, eid, 16,
               nid(xi, 0),
               fuselage_nid(frame_i, left_stringer),
               1.0, 0.0, 0.0); eid += 1

    return eid


def write_engines(fh, counter, eid_start: int) -> int:
    eid = eid_start
    specs = [
        ('right', 0.30, -1.5, -2.5),
        ('right', 0.65, -1.0, -2.5),
        ('left',  0.30, -1.5, -2.5),
        ('left',  0.65, -1.0, -2.5),
    ]
    for k, (side, eta, dx, dz) in enumerate(specs):
        sgn = +1.0 if side == 'right' else -1.0
        chord = (1 - eta) * WING_ROOT_CHORD + eta * WING_TIP_CHORD
        x_le = WING_ROOT_X + eta * WING_HALF_SPAN * math.tan(math.radians(WING_SWEEP_DEG))
        y_here = sgn * (FUSELAGE_RADIUS + eta * WING_HALF_SPAN)
        z_here = WING_ROOT_Z + eta * WING_HALF_SPAN * math.tan(math.radians(WING_DIHEDRAL_DEG))
        x_eng = x_le + 0.30 * chord + dx
        z_eng = z_here + dz

        anchor_nid = ENGINE_NID_BASE + k * 10 + 1
        engine_nid = ENGINE_NID_BASE + k * 10 + 2
        w_grid(fh, counter, anchor_nid, x_le + 0.30 * chord, y_here, z_here)
        w_grid(fh, counter, engine_nid, x_eng, y_here, z_eng)
        w_cbar(fh, counter, eid, 18,
               anchor_nid, engine_nid,
               1.0, 0.0, 0.0); eid += 1

    return eid


# ---------- Main ----------

def main(out_path: Path) -> None:
    t0 = time.time()
    counter = Counter()

    print(f"Building materials and properties...")
    setup_text = setup_materials_and_properties()

    print(f"Streaming to {out_path}...")
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write("$ Boeing 747-class structural model (Node Runner)\n")
        fh.write("$ Generated by examples/generate_b747.py\n")
        fh.write("$ Discretization: "
                 f"frames={FUSELAGE_FRAMES} stringers={FUSELAGE_STRINGERS} "
                 f"wing={WING_RIBS}x{WING_CHORD_NODES} "
                 f"vstab={VS_RIBS}x{VS_CHORD_NODES} "
                 f"hstab={HS_RIBS}x{HS_CHORD_NODES} "
                 f"floor={FLOOR_X_STATIONS}x{FLOOR_Y_STATIONS}\n")
        fh.write("BEGIN BULK\n")

        fh.write("$ ---- Materials and properties ----\n")
        fh.write(setup_text)

        fh.write("$ ---- Fuselage ----\n")
        write_fuselage(fh, counter, FUSELAGE_EID_BASE)
        print(f"  fuselage done ({counter.nodes:,} nodes)")

        fh.write("$ ---- Right wing ----\n")
        write_wing(fh, counter, 'right', RIGHT_WING_NID_BASE, RIGHT_WING_EID_BASE)
        print(f"  right wing done ({counter.nodes:,} nodes)")

        fh.write("$ ---- Left wing ----\n")
        write_wing(fh, counter, 'left', LEFT_WING_NID_BASE, LEFT_WING_EID_BASE)
        print(f"  left wing done ({counter.nodes:,} nodes)")

        fh.write("$ ---- Vertical stabilizer ----\n")
        write_vertical_stab(fh, counter, VS_EID_BASE)
        print(f"  vertical stab done ({counter.nodes:,} nodes)")

        fh.write("$ ---- Right horizontal stab ----\n")
        write_horizontal_stab(fh, counter, 'right', RIGHT_HS_NID_BASE, RIGHT_HS_EID_BASE)
        fh.write("$ ---- Left horizontal stab ----\n")
        write_horizontal_stab(fh, counter, 'left',  LEFT_HS_NID_BASE,  LEFT_HS_EID_BASE)
        print(f"  horizontal stabs done ({counter.nodes:,} nodes)")

        fh.write("$ ---- Floor (panels, beams, stanchions) ----\n")
        write_floor(fh, counter, FLOOR_EID_BASE)
        print(f"  floor done ({counter.nodes:,} nodes)")

        fh.write("$ ---- Engine pylons ----\n")
        write_engines(fh, counter, ENGINE_EID_BASE)

        fh.write("ENDDATA\n")

    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"Wrote {counter.nodes:,} nodes, {counter.elements:,} elements, "
          f"{size_mb:.1f} MiB in {time.time() - t0:.1f}s")


if __name__ == '__main__':
    out = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else Path(__file__).resolve().parent / 'b747.bdf'
    )
    main(out)
