"""Engine-level tests: BDF round-trip, format detection, parallel parser, expression evaluator."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

from node_runner.model import (
    NastranModelGenerator,
    detect_bdf_field_format,
    _convert_bdf_to_free,
    _convert_fixed_line_to_free,
    FIELD_FORMATS,
)
from node_runner.expression import evaluate, ExpressionError


# ---------------------------------------------------------------------------
# Field-format detection and conversion
# ---------------------------------------------------------------------------

class TestFieldFormatDetection:

    def test_short_field_detection(self, tiny_gen, tmp_path):
        path = tmp_path / "m.bdf"
        tiny_gen._write_bdf(str(path), field_format='short')
        assert detect_bdf_field_format(str(path)) == 'short'

    def test_long_field_detection(self, tiny_gen, tmp_path):
        path = tmp_path / "m.bdf"
        tiny_gen._write_bdf(str(path), field_format='long')
        assert detect_bdf_field_format(str(path)) == 'long'

    def test_free_field_detection(self, tiny_gen, tmp_path):
        path = tmp_path / "m.bdf"
        tiny_gen._write_bdf(str(path), field_format='free')
        assert detect_bdf_field_format(str(path)) == 'free'

    def test_missing_file_returns_short_default(self):
        assert detect_bdf_field_format("definitely_does_not_exist.bdf") == 'short'

    def test_format_constants(self):
        assert FIELD_FORMATS == ('short', 'long', 'free')


class TestFreeFieldConversion:

    def test_short_line_to_free(self):
        # Short field: 8-char card name + 8-char fields
        line = "GRID           1               0.      0.      0.\n"
        out = _convert_fixed_line_to_free(line)
        assert out.startswith("GRID,")
        assert "1" in out
        # Trailing empty fields are dropped
        assert not out.rstrip().endswith(',')

    def test_long_line_to_free(self):
        line = "GRID*                  1                              0.              0.\n"
        out = _convert_fixed_line_to_free(line)
        assert out.startswith("GRID*,")

    def test_blank_line_unchanged(self):
        assert _convert_fixed_line_to_free("\n") == "\n"
        assert _convert_fixed_line_to_free("") == ""

    def test_full_punch_file_round_trip(self, tiny_gen, tmp_path):
        path = tmp_path / "m.bdf"
        tiny_gen._write_bdf(str(path), field_format='short')
        # Convert short to free in place, then re-detect
        _convert_bdf_to_free(str(path))
        assert detect_bdf_field_format(str(path)) == 'free'


# ---------------------------------------------------------------------------
# BDF round-trip via _read_bdf_robust
# ---------------------------------------------------------------------------

class TestBdfRoundTrip:

    @pytest.mark.parametrize("fmt", ['short', 'long', 'free'])
    def test_round_trip_each_format(self, tiny_gen, tmp_path, fmt):
        path = tmp_path / f"rt_{fmt}.bdf"
        tiny_gen._write_bdf(str(path), field_format=fmt)
        model, lenient = NastranModelGenerator._read_bdf_robust(str(path))
        assert len(model.nodes) == 6
        assert len(model.elements) == 2
        assert len(model.materials) == 1
        assert len(model.properties) == 2
        # Strict path should succeed for all three formats
        assert lenient is None

    def test_parallel_falls_back_for_small_files(self, tiny_gen, tmp_path):
        """Parallel parser must fall back to single-threaded on tiny files."""
        path = tmp_path / "small.bdf"
        tiny_gen._write_bdf(str(path), field_format='short')
        model, lenient = NastranModelGenerator._read_bdf_parallel(str(path))
        # Files under 50k lines bypass the parallel path entirely - we
        # still get a valid model from the fallback.
        assert len(model.nodes) == 6


# ---------------------------------------------------------------------------
# Expression evaluator
# ---------------------------------------------------------------------------

class TestExpressionEvaluator:

    def test_basic_arithmetic(self):
        result = evaluate("a + b * 2", {"a": np.array([1.0, 2.0]), "b": np.array([3.0, 4.0])})
        np.testing.assert_array_almost_equal(result, [7.0, 10.0])

    def test_numpy_functions(self):
        result = evaluate(
            "sqrt(stress_vm**2 + 0.3*stress_max**2)",
            {"stress_vm": np.array([100.0, 200.0]),
             "stress_max": np.array([10.0, 20.0])},
        )
        assert result.shape == (2,)
        np.testing.assert_array_less(0.0, result)

    def test_safety_factor_pattern(self):
        result = evaluate(
            "1 - vm/450e6",
            {"vm": np.array([100e6, 300e6, 500e6])},
        )
        np.testing.assert_array_almost_equal(
            result, [1 - 100e6/450e6, 1 - 300e6/450e6, 1 - 500e6/450e6],
        )

    def test_constants(self):
        result = evaluate("pi", {})
        assert abs(float(result) - np.pi) < 1e-12

    def test_ifexp_allowed_with_scalar(self):
        # Python's `a if b else c` ternary needs a scalar condition.
        # Element-wise selection should use np.where (covered separately).
        result = evaluate("a if b > 0 else -a",
                          {"a": np.array([1.0, 2.0]), "b": 1.0})
        np.testing.assert_array_almost_equal(result, [1.0, 2.0])

    def test_where_for_elementwise_selection(self):
        # np.where is the right tool for vectorized branching.
        result = evaluate("where(b > 0, a, -a)",
                          {"a": np.array([1.0, 2.0]), "b": np.array([1.0, -1.0])})
        np.testing.assert_array_almost_equal(result, [1.0, -2.0])

    @pytest.mark.parametrize("malicious", [
        "__import__('os').system('echo hi')",
        "open('/etc/passwd').read()",
        "x.__class__.__bases__",
        "lambda x: x",  # we only allow Expression, not Lambda
    ])
    def test_sandbox_rejects_malicious(self, malicious):
        with pytest.raises(ExpressionError):
            evaluate(malicious, {"x": np.array([1.0])})

    def test_undefined_name_rejected(self):
        with pytest.raises(ExpressionError):
            evaluate("undefined_var + 1", {})

    def test_disallowed_function_rejected(self):
        with pytest.raises(ExpressionError):
            evaluate("eval('1+1')", {})

    def test_syntax_error_wrapped(self):
        with pytest.raises(ExpressionError):
            evaluate("a + + +", {"a": np.array([1.0])})


# ---------------------------------------------------------------------------
# Quality metrics for shells and solids
# ---------------------------------------------------------------------------

class TestQualityMetrics:

    def test_unit_quad_metrics(self, tiny_gen):
        q = tiny_gen.calculate_element_quality()
        assert 101 in q
        m = q[101]
        # Unit square: warp=0, aspect=1, skew=0, jacobian=1
        assert abs(m['warp']) < 1e-9
        assert abs(m['aspect'] - 1.0) < 1e-9
        assert abs(m['skew']) < 1e-9
        assert abs(m['jacobian'] - 1.0) < 1e-9

    def test_unit_hex_volume_is_one(self, hex_gen):
        q = hex_gen.calculate_element_quality()
        assert abs(q[1]['volume'] - 1.0) < 1e-9
        assert q[1]['aspect'] == pytest.approx(1.0, abs=1e-9)
        assert q[1]['jacobian'] == pytest.approx(1.0, abs=1e-9)

    def test_corner_tet_volume(self, hex_gen):
        q = hex_gen.calculate_element_quality()
        # Right-angle unit tet has volume 1/6
        assert abs(q[2]['volume'] - 1.0 / 6.0) < 1e-9
