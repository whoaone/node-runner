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


# ---------------------------------------------------------------------------
# INCLUDE handling - real-world multi-file aerospace decks
# ---------------------------------------------------------------------------

class TestIncludeStatements:
    """Covers .dat decks with INCLUDE statements (.dat -> .bdf relative).

    Regression net for v3.0.1: relative paths like ``..\\BULK\\foo.bdf``
    must resolve correctly even when _read_bdf_robust falls back to
    temp-file strategies, and lenient parsing must see cards inside
    included files (previously it scanned only the top-level file).
    """

    def _build_deck(self, root):
        (root / 'exec').mkdir()
        (root / 'bulk').mkdir()
        (root / 'bulk' / 'grids.bdf').write_text(
            'GRID,1,,0.0,0.0,0.0\n'
            'GRID,2,,1.0,0.0,0.0\n'
            'GRID,3,,0.0,1.0,0.0\n'
            'GRID,4,,1.0,1.0,0.0\n'
        )
        (root / 'bulk' / 'elems.bdf').write_text(
            'CQUAD4,1,1,1,2,4,3\n'
            'PSHELL,1,1,0.1\n'
            'MAT1,1,1.0E7,,0.3\n'
        )
        (root / 'exec' / 'main.dat').write_text(
            "SOL 101\nCEND\nTITLE = Test\nBEGIN BULK\n"
            "INCLUDE '..\\bulk\\grids.bdf'\n"
            "INCLUDE '..\\bulk\\elems.bdf'\n"
            "ENDDATA\n"
        )
        return root / 'exec' / 'main.dat'

    def test_strict_include_resolution(self, tmp_path):
        main = self._build_deck(tmp_path)
        m, lenient = NastranModelGenerator._read_bdf_robust(str(main))
        assert len(m.nodes) == 4
        assert len(m.elements) == 1
        assert lenient is None  # clean strict parse

    def test_lenient_path_sees_includes(self, tmp_path):
        # Inject a duplicate-ID GRID into the included file to force
        # fallback into the lenient card-by-card path. The included
        # cards must still arrive in the final model.
        main = self._build_deck(tmp_path)
        grids = tmp_path / 'bulk' / 'grids.bdf'
        grids.write_text(grids.read_text() + 'GRID,3,,9.0,9.0,9.0\n')
        m, lenient = NastranModelGenerator._read_bdf_robust(str(main))
        assert len(m.nodes) == 4
        assert len(m.elements) == 1
        assert lenient is not None
        assert len(lenient.skipped) >= 1

    def test_file_has_includes_detection(self, tmp_path):
        main = self._build_deck(tmp_path)
        assert NastranModelGenerator._file_has_includes(str(main))
        plain = tmp_path / 'plain.bdf'
        plain.write_text('GRID,1,,0.0,0.0,0.0\n')
        assert not NastranModelGenerator._file_has_includes(str(plain))

    def test_inline_includes_flattens(self, tmp_path):
        main = self._build_deck(tmp_path)
        text, missing = NastranModelGenerator._inline_includes(str(main))
        assert missing == []
        # The flattened deck must contain GRIDs and the CQUAD4 inline.
        assert 'GRID,1' in text
        assert 'CQUAD4,1' in text
        # And must NOT still contain raw INCLUDE statements.
        for line in text.splitlines():
            assert not line.lstrip().upper().startswith('INCLUDE ')

    def test_inline_includes_reports_missing(self, tmp_path):
        (tmp_path / 'exec').mkdir()
        (tmp_path / 'exec' / 'main.dat').write_text(
            "BEGIN BULK\nINCLUDE 'nope.bdf'\nENDDATA\n"
        )
        text, missing = NastranModelGenerator._inline_includes(
            str(tmp_path / 'exec' / 'main.dat'))
        assert len(missing) == 1
        assert 'nope.bdf' in missing[0].replace('\\', '/')


# ---------------------------------------------------------------------------
# Multi-line INCLUDE statements
# ---------------------------------------------------------------------------

class TestMultiLineInclude:
    """v3.0.2: real-world decks split long INCLUDE paths across lines.

    Three Nastran continuation styles must work:
      * comma-trailer (the QRG-documented form)
      * unclosed-quote continuation
      * leading + or * marker on the continuation line
    """

    def _make_target(self, root):
        (root / 'sub').mkdir()
        target = root / 'sub' / 'piece.bdf'
        target.write_text('GRID,1,,0.0,0.0,0.0\n')
        return target

    def test_multi_line_path_with_leading_whitespace(self, tmp_path):
        # Most common multi-line form: open quote, path continues on
        # subsequent indented line(s), close quote at the end. Matches
        # pyNastran's documented multi-line example.
        self._make_target(tmp_path)
        main = tmp_path / 'main.dat'
        main.write_text(
            "BEGIN BULK\n"
            "INCLUDE 'sub\n"
            "         /piece.bdf'\n"  # <- indented continuation
            "ENDDATA\n"
        )
        text, missing = NastranModelGenerator._inline_includes(str(main))
        assert missing == [], f'unexpected missing: {missing}'
        assert 'GRID,1' in text

    def test_unclosed_quote_continuation(self, tmp_path):
        self._make_target(tmp_path)
        main = tmp_path / 'main.dat'
        main.write_text(
            "BEGIN BULK\n"
            "INCLUDE 'sub/\n"        # <- no close quote, no comma
            "piece.bdf'\n"
            "ENDDATA\n"
        )
        text, missing = NastranModelGenerator._inline_includes(str(main))
        assert missing == [], f'unexpected missing: {missing}'
        assert 'GRID,1' in text

    def test_continuation_line_with_plus_marker(self, tmp_path):
        self._make_target(tmp_path)
        main = tmp_path / 'main.dat'
        main.write_text(
            "BEGIN BULK\n"
            "INCLUDE 'sub/\n"
            "+piece.bdf'\n"   # <- leading + (Nastran cont marker)
            "ENDDATA\n"
        )
        text, missing = NastranModelGenerator._inline_includes(str(main))
        assert missing == []
        assert 'GRID,1' in text

    def test_two_line_path_resolves_through_robust(self, tmp_path):
        self._make_target(tmp_path)
        main = tmp_path / 'main.dat'
        main.write_text(
            "SOL 101\nCEND\nBEGIN BULK\n"
            "INCLUDE 'sub/\n"
            "piece.bdf'\n"
            "ENDDATA\n"
        )
        m, lenient = NastranModelGenerator._read_bdf_robust(str(main))
        assert len(m.nodes) == 1

    def test_alllods_bdf_style_2line_include(self, tmp_path):
        """The exact 2-line form seen in user's LD2.5.1.8/LOAD/AllLoads.bdf:

           INCLUDE '..\\LOAD\\EXT\\
           LD2-DLL-Entry_E0att009cA1XA_..bdf'

        Line 1 opens the quote and ends mid-path (after a backslash
        separator). Line 2 starts at column 1 with the filename and
        closes the quote. Many real aerospace decks use this form so
        the path stays under MSC Nastran's 72-column field limit.
        """
        (tmp_path / 'LOAD' / 'EXT').mkdir(parents=True)
        (tmp_path / 'EXEC').mkdir()
        target = (tmp_path / 'LOAD' / 'EXT' /
                  'LD2-DLL-Entry_E0att009cA1XA.bdf')
        target.write_text('GRID,999,,1.0,2.0,3.0\n')
        main = tmp_path / 'EXEC' / 'main.bdf'
        main.write_text(
            "BEGIN BULK\n"
            "INCLUDE '..\\LOAD\\EXT\\\n"
            "LD2-DLL-Entry_E0att009cA1XA.bdf'\n"
            "ENDDATA\n"
        )
        m, lenient = NastranModelGenerator._read_bdf_robust(str(main))
        assert 999 in m.nodes


# ---------------------------------------------------------------------------
# Streaming parser with progress callback
# ---------------------------------------------------------------------------

class TestStreamingParser:
    """v3.0.2: _read_bdf_streaming gives per-file progress for big decks."""

    def _build_deck(self, root, n_includes=3):
        (root / 'bulk').mkdir()
        for k in range(n_includes):
            (root / 'bulk' / f'g{k}.bdf').write_text(
                f'GRID,{k+1},,{k}.0,0.0,0.0\n')
        main = root / 'main.dat'
        body = "SOL 101\nCEND\nBEGIN BULK\n"
        for k in range(n_includes):
            body += f"INCLUDE 'bulk/g{k}.bdf'\n"
        body += "ENDDATA\n"
        main.write_text(body)
        return main

    def test_streaming_imports_all_nodes(self, tmp_path):
        main = self._build_deck(tmp_path, n_includes=3)
        m, lenient = NastranModelGenerator._read_bdf_streaming(str(main))
        assert len(m.nodes) == 3

    def test_streaming_emits_progress(self, tmp_path):
        main = self._build_deck(tmp_path, n_includes=4)
        events = []
        def prog(stage, msg, frac):
            events.append((stage, msg, frac))
        m, _ = NastranModelGenerator._read_bdf_streaming(
            str(main), progress=prog)
        assert len(m.nodes) == 4
        # We expect at least one 'inline' event and one 'done' event.
        stages = [e[0] for e in events]
        assert 'inline' in stages
        assert 'done' in stages
        # The final event's fraction should be 1.0.
        assert events[-1][2] == 1.0

    def test_gather_include_chain_lists_files(self, tmp_path):
        main = self._build_deck(tmp_path, n_includes=2)
        files, total, missing = (
            NastranModelGenerator._gather_include_chain(str(main)))
        assert missing == []
        # Root + 2 includes = 3 files.
        assert len(files) == 3
        assert total > 0
