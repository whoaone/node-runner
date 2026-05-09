"""Generate a Boeing 747-class BDF for Node Runner testing.

Builds primary + secondary structure: fuselage skin / frames / stringers,
wing skin / spars / ribs / stringers, vertical and horizontal stabilizers,
floor structure, and engine attachment beams with point masses. No
landing gear.

Materials: Aluminum 2024-T3 (skin), Aluminum 7075-T6 (structural
members), Titanium Ti-6Al-4V (heavy-load fittings, defined but unused
here so the import surfaces a third material).

Properties:
  PSHELL 1: Fuselage skin           t = 1.6 mm
  PSHELL 2: Wing upper skin         t = 4.0 mm
  PSHELL 3: Wing lower skin         t = 3.5 mm
  PSHELL 4: Vertical stab skin      t = 2.0 mm
  PSHELL 5: Horizontal stab skin    t = 1.8 mm
  PSHELL 6: Floor panel             t = 1.0 mm

  PBEAM 10: Fuselage frame
  PBEAM 11: Fuselage stringer
  PBEAM 12: Wing main spar
  PBEAM 13: Wing rib
  PBEAM 14: Wing stringer
  PBEAM 15: Tail spar
  PBEAM 16: Floor beam
  PBEAM 17: Longeron
  PBEAM 18: Engine pylon

Usage:
  python examples/generate_b747.py [output.bdf]
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path

from pyNastran.bdf.bdf import BDF


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
WING_ROOT_Z = -1.5  # below fuselage centerline

VS_ROOT_X = 56.0
VS_ROOT_CHORD = 8.0
VS_TIP_CHORD = 2.5
VS_HEIGHT = 9.0
VS_SWEEP_DEG = 35.0
VS_BASE_Z = FUSELAGE_RADIUS  # vert stab grows from top of fuselage

HS_ROOT_X = 60.0
HS_ROOT_CHORD = 6.0
HS_TIP_CHORD = 2.0
HS_HALF_SPAN = 11.0
HS_SWEEP_DEG = 30.0
HS_Z = VS_BASE_Z + VS_HEIGHT * 0.95  # near top of vertical stab

# ---------- Discretization ----------
# Tuned to push the total node count above 50,000 so Node Runner's
# adaptive LOD threshold (50k) triggers when this model is opened.
FUSELAGE_FRAMES = 220
FUSELAGE_STRINGERS = 80
WING_RIBS = 60
WING_CHORD_NODES = 100
VS_RIBS = 50
VS_CHORD_NODES = 60
HS_RIBS = 40
HS_CHORD_NODES = 50
FLOOR_X_STATIONS = 140
FLOOR_Y_STATIONS = 24

# Stanchion density: one pair (left + right) every Nth floor x-station.
STANCHION_EVERY_N = 4
# Lateral floor cross-beams every Nth x-station (denser than original).
FLOOR_LATERAL_BEAM_EVERY_N = 2

# ---------- ID base offsets (keep ranges separate for readability) ----------
FUSELAGE_NID_BASE = 1
LEFT_WING_NID_BASE = 100_000
RIGHT_WING_NID_BASE = 150_000
VS_NID_BASE = 200_000
LEFT_HS_NID_BASE = 220_000
RIGHT_HS_NID_BASE = 240_000
ENGINE_NID_BASE = 260_000
FLOOR_NID_BASE = 300_000

FUSELAGE_EID_BASE = 1
LEFT_WING_EID_BASE = 100_000
RIGHT_WING_EID_BASE = 150_000
VS_EID_BASE = 200_000
LEFT_HS_EID_BASE = 220_000
RIGHT_HS_EID_BASE = 240_000
ENGINE_EID_BASE = 260_000
FLOOR_EID_BASE = 300_000


# ---------- Helpers ----------

def _setup_materials_properties(model: BDF) -> None:
    # Materials
    model.add_mat1(
        1, E=72.4e9, G=None, nu=0.33, rho=2780.0,
        comment='Aluminum 2024-T3 (skin)',
    )
    model.add_mat1(
        2, E=71.7e9, G=None, nu=0.33, rho=2810.0,
        comment='Aluminum 7075-T6 (structural members)',
    )
    model.add_mat1(
        3, E=113.8e9, G=None, nu=0.342, rho=4430.0,
        comment='Titanium Ti-6Al-4V (heavy-load fittings)',
    )

    # Skins
    model.add_pshell(1, mid1=1, t=0.0016, comment='Fuselage skin (1.6 mm Al 2024-T3)')
    model.add_pshell(2, mid1=1, t=0.0040, comment='Wing upper skin (4.0 mm Al 2024-T3)')
    model.add_pshell(3, mid1=1, t=0.0035, comment='Wing lower skin (3.5 mm Al 2024-T3)')
    model.add_pshell(4, mid1=1, t=0.0020, comment='Vertical stab skin (2.0 mm)')
    model.add_pshell(5, mid1=1, t=0.0018, comment='Horizontal stab skin (1.8 mm)')
    model.add_pshell(6, mid1=1, t=0.0010, comment='Floor panel (1.0 mm)')

    # Beams (Al 7075-T6) - approximate aerospace section properties
    def pbeam(pid, A, I1, I2, J, comment):
        model.add_pbeam(
            pid, 2,
            xxb=[0.0], so=['C'],
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


def _fuselage_radius_at(x: float) -> float:
    """External radius at axial station x along fuselage."""
    if x < NOSE_CONE_LEN:
        ratio = max(x / NOSE_CONE_LEN, 0.0)
        return 0.30 + (FUSELAGE_RADIUS - 0.30) * math.sqrt(ratio)
    if x > FUSELAGE_LEN - TAIL_CONE_LEN:
        ratio = (x - (FUSELAGE_LEN - TAIL_CONE_LEN)) / TAIL_CONE_LEN
        return FUSELAGE_RADIUS - (FUSELAGE_RADIUS - 0.50) * ratio
    return FUSELAGE_RADIUS


def _naca_thickness(xi: float, max_t: float = 0.12) -> float:
    """Half-thickness as fraction of chord (symmetric NACA 4-digit).

    For chord fraction xi in [0, 1], returns y_t/c with peak ~0.06 at
    xi ~= 0.30 for max_t=0.12.
    """
    if xi <= 0.0 or xi >= 1.0:
        return 0.0
    return 5.0 * max_t * (
        0.2969 * math.sqrt(xi)
        - 0.1260 * xi
        - 0.3516 * xi ** 2
        + 0.2843 * xi ** 3
        - 0.1036 * xi ** 4
    )


# ---------- Builders ----------

def build_fuselage(model: BDF, eid_start: int) -> int:
    nid_base = FUSELAGE_NID_BASE
    cols = FUSELAGE_STRINGERS

    def nid(fi, sj):
        return nid_base + fi * cols + sj

    # Grids
    for fi in range(FUSELAGE_FRAMES):
        x = fi * FUSELAGE_LEN / (FUSELAGE_FRAMES - 1)
        r = _fuselage_radius_at(x)
        for sj in range(cols):
            theta = 2.0 * math.pi * sj / cols
            y = r * math.cos(theta)
            z = r * math.sin(theta)
            model.add_grid(nid(fi, sj), [x, y, z])

    eid = eid_start

    # Skin (CQUAD4) wrapped around the cylinder
    for fi in range(FUSELAGE_FRAMES - 1):
        for sj in range(cols):
            sj1 = (sj + 1) % cols
            model.add_cquad4(
                eid, 1,
                [nid(fi, sj), nid(fi + 1, sj), nid(fi + 1, sj1), nid(fi, sj1)],
            )
            eid += 1

    # Frames (circumferential CBARs) at every other frame station to keep
    # element count reasonable; every station for the heavily loaded mid.
    for fi in range(FUSELAGE_FRAMES):
        x_here = fi * FUSELAGE_LEN / (FUSELAGE_FRAMES - 1)
        # In the constant-radius mid section, we put a frame every station.
        # In nose/tail cones we skip every other station.
        if (NOSE_CONE_LEN <= x_here <= (FUSELAGE_LEN - TAIL_CONE_LEN)) or fi % 2 == 0:
            for sj in range(cols):
                sj1 = (sj + 1) % cols
                model.add_cbar(
                    eid, 10,
                    [nid(fi, sj), nid(fi, sj1)],
                    x=[1.0, 0.0, 0.0], g0=None,
                )
                eid += 1

    # Stringers (longitudinal CBARs) - every third station only to avoid
    # ballooning the element count; visually still reads as a stringer grid.
    for sj in range(0, cols, 1):  # every stringer
        for fi in range(FUSELAGE_FRAMES - 1):
            theta = 2.0 * math.pi * sj / cols
            x_ori = [0.0, math.cos(theta), math.sin(theta)]
            model.add_cbar(
                eid, 11,
                [nid(fi, sj), nid(fi + 1, sj)],
                x=x_ori, g0=None,
            )
            eid += 1

    # Longerons: 4 specially marked stringers using PBEAM 17 as overlay
    longeron_indices = [0, cols // 4, cols // 2, 3 * cols // 4]
    for sj in longeron_indices:
        for fi in range(FUSELAGE_FRAMES - 1):
            theta = 2.0 * math.pi * sj / cols
            x_ori = [0.0, math.cos(theta), math.sin(theta)]
            model.add_cbar(
                eid, 17,
                [nid(fi, sj), nid(fi + 1, sj)],
                x=x_ori, g0=None,
            )
            eid += 1

    return eid


def build_wing(
    model: BDF, side: str, nid_base: int, eid_start: int,
) -> int:
    """Build one wing (side='left' or 'right')."""
    sgn = +1.0 if side == 'right' else -1.0
    sweep = math.radians(WING_SWEEP_DEG)
    dihedral = math.radians(WING_DIHEDRAL_DEG)
    y_root = sgn * FUSELAGE_RADIUS

    cols = WING_CHORD_NODES  # one row per surface (upper / lower) per spanwise station

    def nid(span_i, chord_j, surf):
        # surf 0 = upper, 1 = lower
        return nid_base + span_i * (cols * 2) + surf * cols + chord_j

    # Grids
    for si in range(WING_RIBS):
        eta = si / (WING_RIBS - 1)
        chord = (1 - eta) * WING_ROOT_CHORD + eta * WING_TIP_CHORD
        x_le = WING_ROOT_X + eta * WING_HALF_SPAN * math.tan(sweep)
        y_here = y_root + sgn * eta * WING_HALF_SPAN
        z_le = WING_ROOT_Z + eta * WING_HALF_SPAN * math.tan(dihedral)
        for cj in range(cols):
            xi = cj / (cols - 1)
            x = x_le + xi * chord
            half_t = _naca_thickness(xi) * chord
            for surf, dz in ((0, +half_t), (1, -half_t)):
                model.add_grid(nid(si, cj, surf), [x, y_here, z_le + dz])

    eid = eid_start
    upper_pid = 2
    lower_pid = 3

    # Upper skin
    for si in range(WING_RIBS - 1):
        for cj in range(cols - 1):
            n1 = nid(si,     cj,     0)
            n2 = nid(si + 1, cj,     0)
            n3 = nid(si + 1, cj + 1, 0)
            n4 = nid(si,     cj + 1, 0)
            model.add_cquad4(eid, upper_pid, [n1, n2, n3, n4])
            eid += 1
    # Lower skin
    for si in range(WING_RIBS - 1):
        for cj in range(cols - 1):
            n1 = nid(si,     cj,     1)
            n2 = nid(si,     cj + 1, 1)
            n3 = nid(si + 1, cj + 1, 1)
            n4 = nid(si + 1, cj,     1)
            model.add_cquad4(eid, lower_pid, [n1, n2, n3, n4])
            eid += 1

    # Spars (front, mid, rear) - chord positions ~15%, 50%, 85%
    for cj in (int(0.15 * (cols - 1)), int(0.50 * (cols - 1)), int(0.85 * (cols - 1))):
        for si in range(WING_RIBS - 1):
            # Connect upper and lower along span: bar from upper to upper next station
            # plus a vertical web from upper to lower at each station
            model.add_cbar(
                eid, 12,
                [nid(si, cj, 0), nid(si + 1, cj, 0)],
                x=[0.0, 0.0, 1.0], g0=None,
            )
            eid += 1
            model.add_cbar(
                eid, 12,
                [nid(si, cj, 1), nid(si + 1, cj, 1)],
                x=[0.0, 0.0, 1.0], g0=None,
            )
            eid += 1
        # Vertical webs at each rib station for this spar
        for si in range(WING_RIBS):
            model.add_cbar(
                eid, 12,
                [nid(si, cj, 0), nid(si, cj, 1)],
                x=[1.0, 0.0, 0.0], g0=None,
            )
            eid += 1

    # Ribs - chordwise CBARs at each station between upper/lower spar nodes
    rib_chord_nodes = [
        int(0.05 * (cols - 1)), int(0.25 * (cols - 1)), int(0.5 * (cols - 1)),
        int(0.75 * (cols - 1)), int(0.95 * (cols - 1)),
    ]
    for si in range(WING_RIBS):
        for k in range(len(rib_chord_nodes) - 1):
            cj0 = rib_chord_nodes[k]
            cj1 = rib_chord_nodes[k + 1]
            model.add_cbar(
                eid, 13,
                [nid(si, cj0, 0), nid(si, cj1, 0)],
                x=[0.0, 1.0, 0.0], g0=None,
            )
            eid += 1
            model.add_cbar(
                eid, 13,
                [nid(si, cj0, 1), nid(si, cj1, 1)],
                x=[0.0, 1.0, 0.0], g0=None,
            )
            eid += 1

    # Wing stringers along span on upper and lower surfaces
    stringer_chord_positions = [int(p * (cols - 1)) for p in (0.20, 0.35, 0.65, 0.80)]
    for cj in stringer_chord_positions:
        for surf in (0, 1):
            for si in range(WING_RIBS - 1):
                model.add_cbar(
                    eid, 14,
                    [nid(si, cj, surf), nid(si + 1, cj, surf)],
                    x=[0.0, 0.0, 1.0], g0=None,
                )
                eid += 1

    return eid


def build_vertical_stab(model: BDF, eid_start: int) -> int:
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
            half_t = _naca_thickness(xi) * chord
            for surf, dy in ((0, +half_t), (1, -half_t)):
                model.add_grid(nid(si, cj, surf), [x, dy, z_here])

    eid = eid_start
    for surf in (0, 1):
        for si in range(VS_RIBS - 1):
            for cj in range(cols - 1):
                if surf == 0:
                    quad = [nid(si, cj, 0), nid(si + 1, cj, 0),
                            nid(si + 1, cj + 1, 0), nid(si, cj + 1, 0)]
                else:
                    quad = [nid(si, cj, 1), nid(si, cj + 1, 1),
                            nid(si + 1, cj + 1, 1), nid(si + 1, cj, 1)]
                model.add_cquad4(eid, 4, quad)
                eid += 1

    # 2 spars
    for cj in (int(0.20 * (cols - 1)), int(0.70 * (cols - 1))):
        for si in range(VS_RIBS - 1):
            model.add_cbar(eid, 15, [nid(si, cj, 0), nid(si + 1, cj, 0)],
                           x=[0.0, 0.0, 1.0], g0=None); eid += 1
            model.add_cbar(eid, 15, [nid(si, cj, 1), nid(si + 1, cj, 1)],
                           x=[0.0, 0.0, 1.0], g0=None); eid += 1

    return eid


def build_horizontal_stab(
    model: BDF, side: str, nid_base: int, eid_start: int,
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
            half_t = _naca_thickness(xi) * chord
            for surf, dz in ((0, +half_t), (1, -half_t)):
                model.add_grid(nid(si, cj, surf), [x, y_here, HS_Z + dz])

    eid = eid_start
    for surf in (0, 1):
        for si in range(HS_RIBS - 1):
            for cj in range(cols - 1):
                if surf == 0:
                    quad = [nid(si, cj, 0), nid(si + 1, cj, 0),
                            nid(si + 1, cj + 1, 0), nid(si, cj + 1, 0)]
                else:
                    quad = [nid(si, cj, 1), nid(si, cj + 1, 1),
                            nid(si + 1, cj + 1, 1), nid(si + 1, cj, 1)]
                model.add_cquad4(eid, 5, quad)
                eid += 1

    # 2 spars
    for cj in (int(0.25 * (cols - 1)), int(0.75 * (cols - 1))):
        for si in range(HS_RIBS - 1):
            model.add_cbar(eid, 15, [nid(si, cj, 0), nid(si + 1, cj, 0)],
                           x=[0.0, 0.0, 1.0], g0=None); eid += 1

    return eid


def build_floor(model: BDF, eid_start: int) -> int:
    """Floor structure: panels, lateral and longitudinal beams, plus
    stanchions tying the floor edges to the fuselage frames."""
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
            model.add_grid(nid(xi, yi), [x, y, floor_z])

    eid = eid_start

    # Floor panels (CQUAD4)
    for xi in range(FLOOR_X_STATIONS - 1):
        for yi in range(FLOOR_Y_STATIONS - 1):
            model.add_cquad4(
                eid, 6,
                [nid(xi, yi), nid(xi + 1, yi),
                 nid(xi + 1, yi + 1), nid(xi, yi + 1)],
            )
            eid += 1

    # Lateral cross-beams: span side to side at every Nth x-station.
    for xi in range(0, FLOOR_X_STATIONS, FLOOR_LATERAL_BEAM_EVERY_N):
        for yi in range(FLOOR_Y_STATIONS - 1):
            model.add_cbar(
                eid, 16,
                [nid(xi, yi), nid(xi, yi + 1)],
                x=[1.0, 0.0, 0.0], g0=None,
            )
            eid += 1

    # Longitudinal beams: run end-to-end along each y-station.
    for yi in range(FLOOR_Y_STATIONS):
        for xi in range(FLOOR_X_STATIONS - 1):
            model.add_cbar(
                eid, 16,
                [nid(xi, yi), nid(xi + 1, yi)],
                x=[0.0, 1.0, 0.0], g0=None,
            )
            eid += 1

    # Stanchions: vertical struts from the floor's outboard edges up to
    # the nearest fuselage frame node on each side. We compute the
    # fuselage stringer index whose (y, z) is closest to the floor edge.
    cols = FUSELAGE_STRINGERS

    def fuselage_nid(frame_i, stringer_j):
        return FUSELAGE_NID_BASE + frame_i * cols + stringer_j

    # Pre-build x stations of the fuselage so we can find the nearest frame.
    fuselage_x_stations = [
        i * FUSELAGE_LEN / (FUSELAGE_FRAMES - 1) for i in range(FUSELAGE_FRAMES)
    ]

    def nearest_frame(x):
        # Linear search is fine - FUSELAGE_FRAMES is in the low hundreds.
        return min(
            range(len(fuselage_x_stations)),
            key=lambda i: abs(fuselage_x_stations[i] - x),
        )

    # Stringer angles for the floor-edge attachment.
    right_target = math.atan2(floor_z, +y_half)   # ~ -20 deg
    left_target = math.atan2(floor_z, -y_half)    # ~ -160 deg (third quadrant)

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

    right_stringer = nearest_stringer(right_target)
    left_stringer = nearest_stringer(left_target)

    for xi in range(0, FLOOR_X_STATIONS, STANCHION_EVERY_N):
        fx = floor_x_coords[xi]
        frame_i = nearest_frame(fx)

        # Right edge stanchion
        floor_node = nid(xi, FLOOR_Y_STATIONS - 1)
        fuselage_node = fuselage_nid(frame_i, right_stringer)
        model.add_cbar(
            eid, 16,  # reuse "floor cross-beam" property; visually consistent
            [floor_node, fuselage_node],
            x=[1.0, 0.0, 0.0], g0=None,
        )
        eid += 1

        # Left edge stanchion
        floor_node = nid(xi, 0)
        fuselage_node = fuselage_nid(frame_i, left_stringer)
        model.add_cbar(
            eid, 16,
            [floor_node, fuselage_node],
            x=[1.0, 0.0, 0.0], g0=None,
        )
        eid += 1

    return eid


def build_engines(model: BDF, eid_start: int) -> int:
    """Four engines as pylon CBARs anchored to wing rib nodes plus point
    nodes representing engine cores (no CONM2 to keep the import clean)."""
    eid = eid_start
    engine_specs = [
        # (side, span_eta, x_offset, z_offset)
        ('right', 0.30, -1.5, -2.5),
        ('right', 0.65, -1.0, -2.5),
        ('left',  0.30, -1.5, -2.5),
        ('left',  0.65, -1.0, -2.5),
    ]
    for k, (side, eta, dx, dz) in enumerate(engine_specs):
        sgn = +1.0 if side == 'right' else -1.0
        chord = (1 - eta) * WING_ROOT_CHORD + eta * WING_TIP_CHORD
        x_le = WING_ROOT_X + eta * WING_HALF_SPAN * math.tan(math.radians(WING_SWEEP_DEG))
        y_here = sgn * (FUSELAGE_RADIUS + eta * WING_HALF_SPAN)
        z_here = WING_ROOT_Z + eta * WING_HALF_SPAN * math.tan(math.radians(WING_DIHEDRAL_DEG))
        x_eng = x_le + 0.30 * chord + dx
        z_eng = z_here + dz

        anchor_nid = ENGINE_NID_BASE + k * 10 + 1
        engine_nid = ENGINE_NID_BASE + k * 10 + 2
        model.add_grid(anchor_nid, [x_le + 0.30 * chord, y_here, z_here])
        model.add_grid(engine_nid, [x_eng, y_here, z_eng])

        # Pylon as a single CBAR
        model.add_cbar(
            eid, 18,
            [anchor_nid, engine_nid],
            x=[1.0, 0.0, 0.0], g0=None,
        )
        eid += 1

    return eid


# ---------- Main ----------

def main(out_path: Path) -> None:
    t0 = time.time()
    model = BDF(debug=False)

    print("Setting up materials and properties...")
    _setup_materials_properties(model)

    print("Building fuselage...")
    eid = build_fuselage(model, FUSELAGE_EID_BASE)

    print("Building right wing...")
    eid = build_wing(model, 'right', RIGHT_WING_NID_BASE, RIGHT_WING_EID_BASE)
    print("Building left wing...")
    eid = build_wing(model, 'left', LEFT_WING_NID_BASE, LEFT_WING_EID_BASE)

    print("Building vertical stabilizer...")
    eid = build_vertical_stab(model, VS_EID_BASE)

    print("Building horizontal stabilizers...")
    eid = build_horizontal_stab(model, 'right', RIGHT_HS_NID_BASE, RIGHT_HS_EID_BASE)
    eid = build_horizontal_stab(model, 'left',  LEFT_HS_NID_BASE,  LEFT_HS_EID_BASE)

    print("Building floor...")
    eid = build_floor(model, FLOOR_EID_BASE)

    print("Building engine pylons...")
    eid = build_engines(model, ENGINE_EID_BASE)

    print(f"Built in {time.time() - t0:.2f}s. "
          f"Nodes={len(model.nodes):,}  Elements={len(model.elements):,}  "
          f"Properties={len(model.properties)}  Materials={len(model.materials)}")

    print(f"Writing {out_path} ...")
    t1 = time.time()
    model.write_bdf(str(out_path), size=8, is_double=False)
    print(f"Wrote in {time.time() - t1:.2f}s "
          f"({out_path.stat().st_size / (1024 * 1024):.1f} MiB).")


if __name__ == '__main__':
    out = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else Path(__file__).resolve().parent / 'b747.bdf'
    )
    main(out)
