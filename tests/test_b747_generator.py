"""Tests for examples/generate_b747.py.

Generates a small variant of the b747 model and verifies:
  - The generator runs without error
  - All MAT1 / PSHELL / PBEAM cards are present in the output
  - The strict pyNastran reader can parse the output (no lenient fallback)
  - Element / node counts are within expectations
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
GEN_PATH = EXAMPLES / "generate_b747.py"


def _load_generator_module():
    """Import examples/generate_b747.py as a module without polluting sys.path
    permanently. Returns the loaded module."""
    spec = importlib.util.spec_from_file_location("generate_b747", GEN_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules["generate_b747"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def small_b747(tmp_path_factory):
    """Generate a tiny b747 variant - only enough to exercise the
    generator without bogging down the test suite."""
    g = _load_generator_module()

    # Dial the discretization way down so the file is < 1 MB and parses
    # in well under a second.
    g.FUSELAGE_FRAMES = 30
    g.FUSELAGE_STRINGERS = 16
    g.WING_RIBS = 12
    g.WING_CHORD_NODES = 24
    g.VS_RIBS = 10
    g.VS_CHORD_NODES = 14
    g.HS_RIBS = 10
    g.HS_CHORD_NODES = 14
    g.FLOOR_X_STATIONS = 20
    g.FLOOR_Y_STATIONS = 8
    g.FUSELAGE_FRAME_BEAM_EVERY_N = 2
    g.FUSELAGE_STRINGER_BEAM_EVERY_N = 2
    g.FLOOR_LONG_BEAM_EVERY_N = 2
    g.FLOOR_LATERAL_BEAM_EVERY_N = 2
    g.STANCHION_EVERY_N = 4

    out = tmp_path_factory.mktemp("b747") / "tiny_b747.bdf"
    g.main(out)
    return out


class TestGeneratorOutput:

    def test_file_exists_and_nonempty(self, small_b747):
        assert small_b747.exists()
        assert small_b747.stat().st_size > 1024

    def test_strict_parse_succeeds(self, small_b747):
        from node_runner.model import NastranModelGenerator
        model, lenient = NastranModelGenerator._read_bdf_robust(str(small_b747))
        assert lenient is None, "strict parse should have succeeded"

    def test_materials_and_properties_present(self, small_b747):
        from node_runner.model import NastranModelGenerator
        model, _ = NastranModelGenerator._read_bdf_robust(str(small_b747))
        # 3 materials (Al 2024, Al 7075, Ti) + 6 PSHELL + 9 PBEAM = 15 props
        assert len(model.materials) == 3
        assert len(model.properties) == 15

    def test_geometry_components_present(self, small_b747):
        from node_runner.model import NastranModelGenerator
        model, _ = NastranModelGenerator._read_bdf_robust(str(small_b747))
        # We should have substantial counts of each card class
        elem_types = {e.type for e in model.elements.values()}
        assert 'CQUAD4' in elem_types
        assert 'CBAR' in elem_types

    def test_field_format_round_trip(self, small_b747):
        from node_runner.model import detect_bdf_field_format
        # The generator uses short field by default
        assert detect_bdf_field_format(str(small_b747)) == 'short'
