"""Geometry fixes: mixed solid cell types, cylindrical / spherical CP, global coords in quality and mass, hex Jacobian."""
from __future__ import annotations

import numpy as np
import pytest
import pyvista as pv

from node_runner.model import NastranModelGenerator
from node_runner.scene_build import build_element_arrays_vectorized, build_node_coords_vectorized


def _read(tmp_path, bulk):
    path = tmp_path / "m.bdf"
    path.write_text("SOL 101\nCEND\nBEGIN BULK\n" + bulk + "ENDDATA\n")
    model, _ = NastranModelGenerator._read_bdf_robust(str(path))
    NastranModelGenerator._finalize_for_viewer(model)
    return model


MIXED_SOLIDS = """\
GRID    1               0.      0.      0.
GRID    2               1.      0.      0.
GRID    3               1.      1.      0.
GRID    4               0.      1.      0.
GRID    5               0.      0.      1.
GRID    6               1.      0.      1.
GRID    7               1.      1.      1.
GRID    8               0.      1.      1.
GRID    9               0.      0.      2.
CHEXA   1       1       1       2       3       4       5       6
        7       8
CTETRA  2       1       5       6       8       9
PSOLID  1       1
MAT1    1       1.+7            .3
"""


def test_mixed_solid_bucket_keeps_each_cell_type(tmp_path):
    model = _read(tmp_path, MIXED_SOLIDS)
    ids = np.array(sorted(model.nodes))
    arr = build_element_arrays_vectorized(model, ids, {n: i for i, n in enumerate(ids)})
    types = dict(zip(arr['eid'].tolist(), arr['cell_types'].tolist()))
    assert types == {1: pv.CellType.HEXAHEDRON, 2: pv.CellType.TETRA}


CURVED_CP = """\
CORD2C  2       0       0.      0.      0.      1.      0.      0.
        0.      0.      1.
CORD2S  3       0       0.      0.      0.      0.      0.      1.
        1.      0.      0.
GRID    1       2       5.      90.     20.
GRID    2       3       2.      90.     90.
GRID    3       2       5.      0.      20.
"""


def test_cylindrical_and_spherical_cp_land_where_pynastran_puts_them(tmp_path):
    model = _read(tmp_path, CURVED_CP)
    ids = [1, 2, 3]
    expected = np.array([model.nodes[n].get_position() for n in ids])
    assert np.allclose(build_node_coords_vectorized(model, ids), expected)
    assert np.allclose(expected[0], [20., -5., 0.])     # CORD2C 2: z along basic X, theta 0 along basic Z


OFFSET_CP = """\
CORD2R  5       0       100.    0.      0.      100.    0.      1.
        101.    0.      0.
GRID    1               0.      0.      0.
GRID    2               1.      0.      0.
GRID    3       5       -99.    1.      0.
GRID    4       5       -100.   1.      0.
CQUAD4  1       1       1       2       3       4
PSHELL  1       1       .1
MAT1    1       1.+7            .3      .1
CONM2   9       3               10.
"""


def test_quality_and_mass_use_basic_positions(tmp_path):
    gen = NastranModelGenerator()
    gen.model = _read(tmp_path, OFFSET_CP)    # a unit square whose nodes 3 and 4 sit in an offset CP
    assert gen.calculate_element_quality()[1]['aspect'] == pytest.approx(1.0)
    cg = gen.calculate_mass_summary()['cg']
    assert cg[0] == pytest.approx((10 * 1.0 + 0.01 * 0.5) / 10.01, rel=1e-6)


CUBE = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                 [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)


def test_hex_jacobian_flags_inverted_but_not_reversed():
    jac = NastranModelGenerator._calculate_hex_jacobian
    assert jac(CUBE) == pytest.approx(1.0)
    assert jac(CUBE[[4, 5, 6, 7, 0, 1, 2, 3]]) == pytest.approx(1.0)   # top face first: every corner flips alike
    inverted = CUBE.copy()
    inverted[6] = [0.2, 0.2, -0.5]                                  # corner 7 pushed through the bottom face
    assert jac(inverted) < 0
