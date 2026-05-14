"""v4.0.0 (Phase D) parallel-parse tests.

The per-INCLUDE ProcessPoolExecutor parser must produce the same model
counts as the serial path, and must gracefully fall back when the
engagement heuristics say "don't bother" (small decks, DMIG, etc).
"""
import os

import pytest


def _write_multi_include_deck(tmp_path, n_includes=4, grids_per=200):
    """Write a master BDF with n_includes INCLUDE files.

    Each include holds grids + a CQUAD4 shell mesh + materials/properties
    so the per-file workers have realistic merge work.
    """
    master = tmp_path / 'master.bdf'
    lines = ['BEGIN BULK\n']
    for i in range(n_includes):
        inc_name = f'inc_{i}.bdf'
        inc_path = tmp_path / inc_name
        inc_lines = []
        base_nid = 1 + i * grids_per
        base_eid = 1 + i * grids_per
        for k in range(grids_per):
            nid = base_nid + k
            x = float(k % 50)
            y = float(k // 50)
            z = float(i)
            inc_lines.append(f'GRID,{nid},,{x},{y},{z}\n')
        # Make a tiny QUAD strip from the first 4 nodes for completeness.
        if grids_per >= 4:
            inc_lines.append(
                f'CQUAD4,{base_eid},{1000+i},{base_nid},'
                f'{base_nid+1},{base_nid+2},{base_nid+3}\n')
        inc_lines.append(f'PSHELL,{1000+i},{2000+i},0.1\n')
        inc_lines.append(f'MAT1,{2000+i},2.0e11,,0.3,7850.0\n')
        # FORCE on a node so loads / load_combinations exercise the
        # list-extend merge path.
        inc_lines.append(f'FORCE,{500},{base_nid},,1.0,1.0,0.0,0.0\n')
        inc_path.write_text(''.join(inc_lines))
        lines.append(f"INCLUDE '{inc_name}'\n")
    lines.append('ENDDATA\n')
    master.write_text(''.join(lines))
    return master


def test_should_use_parallel_thresholds(tmp_path):
    """Small deck (1 file, tiny) → skip parallel."""
    from node_runner.parallel_parse import should_use_parallel
    p = tmp_path / 'small.bdf'
    p.write_text('GRID,1,,0.0,0.0,0.0\n')
    assert should_use_parallel(str(p)) is False


def test_parallel_equals_serial_on_synthetic_deck(tmp_path):
    """Force parallel run on a tiny multi-INCLUDE deck and compare
    merge result to the serial streaming path."""
    from node_runner.parallel_parse import parse_bdf_parallel
    from node_runner.model import NastranModelGenerator

    master = _write_multi_include_deck(tmp_path, n_includes=3,
                                       grids_per=50)

    # Parallel
    pmodel, _ = parse_bdf_parallel(str(master), n_workers=2)
    assert pmodel is not None, "parallel parse should not fall back"
    p_nodes = len(pmodel.nodes)
    p_elements = len(pmodel.elements)
    p_props = len(pmodel.properties)
    p_mats = len(pmodel.materials)
    # FORCE sid 500 appears in every INCLUDE → list extended.
    p_loads_500 = len(pmodel.loads.get(500, []))

    # Serial
    smodel, _ = NastranModelGenerator._read_bdf_streaming(str(master))
    s_nodes = len(smodel.nodes)
    s_elements = len(smodel.elements)
    s_props = len(smodel.properties)
    s_mats = len(smodel.materials)
    s_loads_500 = len(smodel.loads.get(500, []))

    assert p_nodes == s_nodes
    assert p_elements == s_elements
    assert p_props == s_props
    assert p_mats == s_mats
    assert p_loads_500 == s_loads_500


def test_parallel_returns_none_when_no_files(tmp_path):
    """Unparseable / missing master returns (None, None) cleanly."""
    from node_runner.parallel_parse import parse_bdf_parallel
    missing = tmp_path / 'nope.bdf'
    model, lenient = parse_bdf_parallel(str(missing))
    assert model is None
    assert lenient is None


def test_dmig_deck_skipped_by_heuristic(tmp_path):
    """Decks containing DMIG cards must fall back to serial because
    per-file workers can't resolve cross-INCLUDE refs."""
    from node_runner.parallel_parse import should_use_parallel
    master = _write_multi_include_deck(tmp_path, n_includes=5,
                                       grids_per=200)
    # Inject a DMIG sentinel into one of the includes.
    inc = tmp_path / 'inc_2.bdf'
    text = inc.read_text()
    inc.write_text('DMIG,K2GG,0,6,1,0\n' + text)
    assert should_use_parallel(str(master)) is False
