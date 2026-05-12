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


# ---------------------------------------------------------------------------
# Chunked-strict parser preserves cross-chunk references
# ---------------------------------------------------------------------------

class TestChunkedCrossRefs:
    """Regression: a card in chunk N that references an ID defined in
    chunk N-1 must still resolve correctly in the merged master model.

    This works because:
      * chunks read with xref=False, validate=False: pyNastran doesn't
        check that referenced IDs exist at parse time. Cards store the
        IDs as plain integers.
      * after each chunk parses, we dict-merge its entities into the
        master BDF. The master ends up with every card in original
        order.
      * referenced IDs in later chunks resolve against the merged
        master, not against the chunk that contained the reference.
    """

    def _force_chunked_path(self, deck_text):
        """Inject a duplicate-ID card at the top to force strict-fail."""
        return deck_text.replace(
            'BEGIN BULK\n',
            'BEGIN BULK\nGRID,1,,9.9,9.9,9.9\n',
        )

    def test_element_in_later_chunk_refs_grid_in_earlier_chunk(self, tmp_path):
        # Build a deck guaranteed to span >= 2 chunks of 20k cards each.
        deck = ['BEGIN BULK\n']
        n_grids = 30_000
        for k in range(1, n_grids + 1):
            x = (k % 100) * 0.1
            y = ((k // 100) % 100) * 0.1
            deck.append(f'GRID,{k},,{x:.3f},{y:.3f},0.0\n')
        deck.append('MAT1,100,1.0E7,,0.3\n')
        deck.append('PSHELL,1,100,0.1\n')
        # Elements in chunk 2+ that reference GRIDs from chunk 1.
        n_elems = 5000
        for k in range(1, n_elems + 1):
            deck.append(f'CQUAD4,{k},1,{k},{k+1},{k+2},{k+3}\n')
        deck.append('ENDDATA\n')
        text = self._force_chunked_path(''.join(deck))

        path = tmp_path / 'deck.bdf'
        path.write_text(text)

        m, lenient = NastranModelGenerator._read_bdf_streaming(str(path))
        # All grids + all elements survived the chunk boundary.
        assert len(m.nodes) == n_grids, \
            f'expected {n_grids:,} nodes, got {len(m.nodes):,}'
        assert len(m.elements) == n_elems, \
            f'expected {n_elems:,} elements, got {len(m.elements):,}'

        # Pick a mid-deck element. Its node IDs lived in chunk 1, the
        # element itself lived in chunk 2+. Both must be in the merged
        # master.
        eid = n_elems // 2
        elem = m.elements[eid]
        for nid in elem.node_ids:
            assert nid in m.nodes, \
                f'GRID {nid} (referenced by element {eid}) missing from master'

        # The merged master must also cross-reference cleanly.
        m.cross_reference()

    def test_loads_in_later_chunk_resolve_to_earlier_grids(self, tmp_path):
        deck = ['BEGIN BULK\n']
        n_grids = 25_000
        for k in range(1, n_grids + 1):
            deck.append(f'GRID,{k},,0.0,0.0,0.0\n')
        # FORCEs and SPCs referencing GRIDs from chunk 1.
        for k in range(1, 200):
            deck.append(f'FORCE,1,{k},,1.0,0.0,0.0,1.0\n')
        for k in range(1, 50):
            deck.append(f'SPC1,2,123456,{k}\n')
        deck.append('ENDDATA\n')
        text = self._force_chunked_path(''.join(deck))

        path = tmp_path / 'deck.bdf'
        path.write_text(text)
        m, _ = NastranModelGenerator._read_bdf_streaming(str(path))

        assert len(m.nodes) == n_grids
        total_loads = sum(len(v) for v in m.loads.values())
        assert total_loads >= 199, f'lost loads across chunk boundary: {total_loads}'
        # And cross_reference still resolves.
        m.cross_reference()

    def test_vectorized_node_coords_match_get_position(self, tmp_path):
        """v3.2.0 regression: build_node_coords_vectorized must produce
        the same xyz array as the per-node get_position() loop, both
        for basic-coord nodes and for non-basic-coord nodes."""
        from node_runner.scene_build import build_node_coords_vectorized
        deck = (
            "SOL 101\nCEND\nBEGIN BULK\n"
            "CORD2R,10,0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n"
            "GRID,1,0,1.0,0.0,0.0\n"
            "GRID,2,0,0.0,2.0,0.0\n"
            "GRID,3,10,1.0,0.0,0.0\n"
            "GRID,4,10,0.0,1.0,0.0\n"
            "ENDDATA\n"
        )
        path = tmp_path / 'deck.bdf'
        path.write_text(deck)
        m, _ = NastranModelGenerator._read_bdf_robust(str(path))
        NastranModelGenerator._finalize_for_viewer(m)

        sorted_ids = sorted(m.nodes.keys())
        # Per-node reference
        ref = np.array([m.nodes[n].get_position() for n in sorted_ids])
        vec = build_node_coords_vectorized(m, sorted_ids)
        assert vec.shape == ref.shape
        assert np.allclose(vec, ref, atol=1e-9), \
            f'vectorized != per-node:\n{vec}\nvs\n{ref}'

    def test_vectorized_element_arrays(self, tmp_path):
        """build_element_arrays_vectorized must produce the right cell
        counts + matching EID arrays."""
        from node_runner.scene_build import build_element_arrays_vectorized
        deck = (
            "BEGIN BULK\n"
            "GRID,1,,0.0,0.0,0.0\n"
            "GRID,2,,1.0,0.0,0.0\n"
            "GRID,3,,0.0,1.0,0.0\n"
            "GRID,4,,1.0,1.0,0.0\n"
            "PSHELL,1,1,0.1\n"
            "MAT1,1,1.0E7,,0.3\n"
            "CQUAD4,1,1,1,2,4,3\n"
            "CTRIA3,2,1,1,2,3\n"
            "CBAR,3,1,1,2,1.0,0.0,0.0\n"
            "ENDDATA\n"
        )
        path = tmp_path / 'deck.bdf'
        path.write_text(deck)
        m, _ = NastranModelGenerator._read_bdf_robust(str(path))
        NastranModelGenerator._finalize_for_viewer(m)
        sorted_ids = sorted(m.nodes.keys())
        node_map = {n: i for i, n in enumerate(sorted_ids)}
        out = build_element_arrays_vectorized(m, sorted_ids, node_map)
        assert out['eid'].size == 3
        assert set(out['eid'].tolist()) == {1, 2, 3}
        # Cells array shape: CQUAD4 has 5 ints (4+1), CTRIA3 has 4,
        # CBAR has 3. Total = 12.
        assert out['cells'].size == 5 + 4 + 3

    def test_apply_display_lod_passes_through_small_grids(self, tmp_path):
        """Under-threshold grids must be returned untouched."""
        from node_runner.scene_build import apply_display_lod
        import pyvista as pv
        grid = pv.UnstructuredGrid(
            np.array([4, 0, 1, 2, 3], dtype=np.int64),
            np.array([pv.CellType.QUAD], dtype=np.uint8),
            np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float),
        )
        display, info = apply_display_lod(grid, threshold=10)
        assert not info['lod_active']
        assert display is grid

    def test_finalize_for_viewer_topo_sorts_coords(self, tmp_path):
        """v3.3.0: coords whose rid references another coord must
        be resolved in dependency order, even when added out-of-order.

        Before v3.3.0 the loop walked coords.items() in dict order and
        could leave a child's rid_ref pointing at an un-setup parent,
        which crashed later when pyNastran tried to compute local
        unit vectors from the parent's transform.
        """
        from pyNastran.bdf.bdf import BDF
        m = BDF(debug=False)
        # Add cid=2 BEFORE its parent cid=1. _finalize must resolve
        # parent first.
        m.add_cord2r(cid=2, origin=[10, 0, 0],
                     zaxis=[10, 0, 1], xzplane=[11, 0, 0], rid=1)
        m.add_cord2r(cid=1, origin=[5, 0, 0],
                     zaxis=[5, 0, 1], xzplane=[6, 0, 0], rid=0)
        NastranModelGenerator._finalize_for_viewer(m)
        # Both coords have their unit vectors set.
        for cid in (1, 2):
            cs = m.coords[cid]
            assert cs.i is not None
            assert cs.j is not None
            assert cs.k is not None
            assert np.allclose(np.linalg.norm(cs.i), 1.0)
        # No spurious warnings.
        assert m._coord_resolution_warnings == []

    def test_finalize_for_viewer_handles_dangling_rid(self, tmp_path):
        """v3.3.0: a coord whose rid references a coord not in the deck
        must not raise. The user hit this on FUSE-XWFF_04A: CORD2R
        cid=920000 rid=200000, but 200000 was never defined.

        Expected: rid_ref falls back to basic, a warning is recorded
        on model._coord_resolution_warnings, and node.get_position()
        works on GRIDs that reference that coord.
        """
        from pyNastran.bdf.bdf import BDF
        m = BDF(debug=False)
        # rid=200000 doesn't exist in this model.
        m.add_cord2r(cid=920000, origin=[10, 0, 0],
                     zaxis=[10, 0, 1], xzplane=[11, 0, 0], rid=200000)
        m.add_grid(1, [1.0, 2.0, 3.0], cp=920000)
        # Must not raise.
        NastranModelGenerator._finalize_for_viewer(m)
        assert m._coord_resolution_warnings
        assert any(w[0] == 920000 for w in m._coord_resolution_warnings)
        # And node.get_position works on a node that references it.
        pos = m.nodes[1].get_position()
        assert pos is not None and len(pos) == 3

    def test_vectorized_element_arrays_returns_group_counts(self, tmp_path):
        """v3.3.0: build_element_arrays_vectorized piggy-backs the
        By-Type / By-Shape group counts onto its existing single pass
        over all elements so the model tree doesn't have to do its own
        1.2M-iteration Counter loop."""
        from node_runner.scene_build import build_element_arrays_vectorized
        deck = (
            "BEGIN BULK\n"
            "GRID,1,,0.0,0.0,0.0\n"
            "GRID,2,,1.0,0.0,0.0\n"
            "GRID,3,,0.0,1.0,0.0\n"
            "GRID,4,,1.0,1.0,0.0\n"
            "GRID,5,,0.0,0.0,1.0\n"
            "GRID,6,,1.0,0.0,1.0\n"
            "GRID,7,,1.0,1.0,1.0\n"
            "GRID,8,,0.0,1.0,1.0\n"
            "PSHELL,1,1,0.1\n"
            "PSOLID,2,1\n"
            "MAT1,1,1.0E7,,0.3\n"
            "CQUAD4,1,1,1,2,3,4\n"
            "CQUAD4,2,1,1,2,3,4\n"
            "CHEXA,3,2,1,2,3,4,5,6,7,8\n"
            "RBE2,4,1,123,2,3\n"
            "ENDDATA\n"
        )
        path = tmp_path / 'deck.bdf'
        path.write_text(deck)
        m, _ = NastranModelGenerator._read_bdf_robust(str(path))
        NastranModelGenerator._finalize_for_viewer(m)
        sorted_ids = sorted(m.nodes.keys())
        node_map = {n: i for i, n in enumerate(sorted_ids)}
        out = build_element_arrays_vectorized(m, sorted_ids, node_map)
        assert 'by_type_counts' in out
        assert 'by_shape_counts' in out
        assert out['by_type_counts'].get('Plates') == 2
        assert out['by_type_counts'].get('Solids') == 1
        assert out['by_type_counts'].get('Rigid') == 1
        assert out['by_shape_counts'].get('Quad') == 2
        assert out['by_shape_counts'].get('Hex') == 1
        assert out['by_shape_counts'].get('Rigid') == 1

    def test_finalize_for_viewer_resolves_cp_ref(self, tmp_path):
        """v3.1.5 regression: a deck whose GRIDs reference a non-basic
        coord system must be usable for get_position() even without
        a full pyNastran cross_reference pass.

        Before v3.1.5: get_position() crashed with
            'NoneType' object has no attribute 'transform_node_to_global'
        because we read with xref=False for speed and never resolved
        node.cp_ref. The post-parse _finalize_for_viewer fixup wires
        every node's cp_ref so the viewer's scene build works.
        """
        deck = (
            "SOL 101\nCEND\nBEGIN BULK\n"
            "CORD2R,10,0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n"
            "GRID,1,10,1.0,0.0,0.0\n"
            "GRID,2,0,5.0,0.0,0.0\n"
            "GRID,3,99,0.0,1.0,0.0\n"   # cp=99 doesn't exist; falls back to basic
            "ENDDATA\n"
        )
        path = tmp_path / 'deck.bdf'
        path.write_text(deck)
        m, _ = NastranModelGenerator._read_bdf_robust(str(path))
        NastranModelGenerator._finalize_for_viewer(m)
        # Every node now has cp_ref set; get_position works.
        for nid in (1, 2, 3):
            pos = m.nodes[nid].get_position()
            assert pos is not None
            assert len(pos) == 3

    def test_dmig_matrix_survives_chunked_path(self, tmp_path):
        """A DMIG matrix with many entries (more than chunk_card_count)
        must keep all its entries when the chunked path is taken.

        v3.1.4 pre-extracts DMIG-family cards into a dedicated punch
        sub-parse so chunking can't split a matrix and the chunk's
        lenient fallback can't drop entries during add_card. Before
        v3.1.4 this test showed KAAX with 0 entries.
        """
        # Build a deck with > 8000 DMIG entries so chunking would
        # normally split it. Also include a duplicate GRID to force
        # the strict path to fail and trigger the chunked fallback.
        n_dmig = 12_000  # > default chunk size of 8000
        deck = ['BEGIN BULK\n']
        deck.append('GRID,1,,9.9,9.9,9.9\n')  # duplicate forces chunked
        for k in range(1, 101):
            deck.append(f'GRID,{k},,0.0,0.0,0.0\n')
        deck.append('DMIG,KAAX,0,6,1,0\n')  # matrix header
        for k in range(1, n_dmig + 1):
            gj = (k % 100) + 1
            cj = (k % 6) + 1
            deck.append(f'DMIG,KAAX,{gj},{cj},,{gj},{cj},{1.0 + k:.4f}\n')
        deck.append('ENDDATA\n')

        path = tmp_path / 'deck.bdf'
        path.write_text(''.join(deck))
        m, _ = NastranModelGenerator._read_bdf_streaming(str(path))

        assert 'KAAX' in m.dmig, 'KAAX matrix missing from chunked import'
        kaax = m.dmig['KAAX']
        # GCi / GCj / Real must each have all n_dmig entries.
        assert len(kaax.GCi) == n_dmig, \
            f'lost DMIG entries: {len(kaax.GCi)} != {n_dmig}'
        assert len(kaax.GCj) == n_dmig
        assert len(kaax.Real) == n_dmig
