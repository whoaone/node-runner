"""Shared pytest fixtures for the Node Runner test suite.

Run from the repo root:

    venv\\Scripts\\python.exe -m pytest tests -v

A small synthetic model (`tiny_gen`) is reused across most tests; UI
tests get a single shared QApplication via `qapp` (auto-supplied by
pytest-qt).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

# Make `node_runner` importable when pytest runs from the repo root.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pytest

from node_runner.model import NastranModelGenerator


@pytest.fixture
def tiny_gen():
    """A 4-node planar quad + a 2-node beam + 1 mat + 2 props."""
    gen = NastranModelGenerator()
    m = gen.model
    for nid, xyz in zip(
        (1, 2, 3, 4, 5, 6),
        ([0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [3, 0, 0]),
    ):
        m.add_grid(nid, list(xyz))
    m.add_mat1(1, E=2e11, G=None, nu=0.3, rho=7800.0)
    m.add_pshell(1, mid1=1, t=0.001)
    m.add_pbeam(2, 1, xxb=[0.0], so=['C'], area=[1e-4],
                i1=[1e-8], i2=[1e-8], i12=[0.0], j=[1e-9])
    m.add_cquad4(101, 1, [1, 2, 3, 4])
    m.add_cbar(201, 2, [5, 6], x=[0.0, 1.0, 0.0], g0=None)
    return gen


@pytest.fixture
def two_quads_gen():
    """Two adjacent CQUAD4 sharing edge 2-3, useful for refine/smooth/insert tests."""
    gen = NastranModelGenerator()
    m = gen.model
    for nid, xyz in zip(
        (1, 2, 3, 4, 5, 6),
        ([0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0]),
    ):
        m.add_grid(nid, list(xyz))
    m.add_mat1(1, E=2e11, G=None, nu=0.3, rho=7800.0)
    m.add_pshell(1, mid1=1, t=0.001)
    m.add_cquad4(1, 1, [1, 2, 3, 4])
    m.add_cquad4(2, 1, [2, 5, 6, 3])
    return gen


@pytest.fixture
def hex_gen():
    """Single unit hex + one tet for solid-quality tests."""
    gen = NastranModelGenerator()
    m = gen.model
    for nid, xyz in zip(
        range(1, 9),
        ([0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
         [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]),
    ):
        m.add_grid(nid, xyz)
    for nid, xyz in zip(
        (11, 12, 13, 14),
        ([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]),
    ):
        m.add_grid(nid, xyz)
    m.add_mat1(1, E=2e11, G=None, nu=0.3, rho=7800.0)
    m.add_psolid(1, mid=1)
    m.add_chexa(1, 1, [1, 2, 3, 4, 5, 6, 7, 8])
    m.add_ctetra(2, 1, [11, 12, 13, 14])
    return gen
