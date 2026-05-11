"""Tests for node_runner/skip_cards.py - the analysis-only card stripper.

v3.2.2 pre-strips cards Node Runner doesn't render (TEMP, TLOAD, tables,
optimization, etc.) before pyNastran sees them. This dropped a real
user's deck from 'never finishes' to under a minute. These tests pin
the contract: strip the right cards, preserve everything else, and
make sure the round-trip via _write_bdf actually writes the cards
back.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from node_runner.skip_cards import (
    ANALYSIS_ONLY_CARDS, strip_analysis_only_cards, _card_name,
)
from node_runner.model import NastranModelGenerator


class TestSkipCards:

    def test_skip_set_includes_temp_and_tload(self):
        """Pin the cards that must be in the skip list."""
        for name in ('TEMP', 'TEMPD', 'TLOAD1', 'RLOAD1',
                     'TABLED1', 'DCONSTR', 'SET1'):
            assert name in ANALYSIS_ONLY_CARDS, \
                f'{name} must be in ANALYSIS_ONLY_CARDS'

    def test_skip_set_excludes_geometry_cards(self):
        """Pin the cards that must NEVER be in the skip list."""
        for name in ('GRID', 'CQUAD4', 'CTRIA3', 'CHEXA',
                     'CBAR', 'CBEAM', 'CBUSH', 'RBE2', 'RBE3',
                     'PSHELL', 'PBEAM', 'MAT1', 'CORD2R',
                     'FORCE', 'MOMENT', 'PLOAD4', 'SPC1'):
            assert name not in ANALYSIS_ONLY_CARDS, \
                f'{name} must NOT be in ANALYSIS_ONLY_CARDS'

    def test_card_name_parses_fixed_format(self):
        assert _card_name('GRID    1       0       0.0     0.0     0.0') == 'GRID'
        assert _card_name('TEMP    100     1       72.5') == 'TEMP'
        assert _card_name('CQUAD4  10      1       1       2       3       4') == 'CQUAD4'

    def test_card_name_parses_free_format(self):
        assert _card_name('GRID,1,,0.0,0.0,0.0\n') == 'GRID'
        assert _card_name('TEMP,100,1,72.5\n') == 'TEMP'

    def test_card_name_handles_continuation_lines(self):
        assert _card_name('+CONT   1.0     2.0     3.0') == ''
        assert _card_name('*       1.0     2.0     3.0') == ''
        assert _card_name('        2       3       4') == ''  # leading spaces

    def test_card_name_handles_comments_and_blanks(self):
        assert _card_name('$ this is a comment') == ''
        assert _card_name('') == ''
        assert _card_name('\n') == ''

    def test_card_name_handles_long_format_suffix(self):
        assert _card_name('GRID*   1') == 'GRID'  # asterisk stripped

    def test_strip_keeps_geometry_drops_temp(self):
        sample = (
            'GRID,1,,0.0,0.0,0.0\n'
            'GRID,2,,1.0,0.0,0.0\n'
            'CQUAD4,1,1,1,2,1,2\n'
            'TEMP,100,1,72.0\n'
            'TEMP,100,2,73.0\n'
            'FORCE,1,1,,1.0,0.0,0.0,1.0\n'
        )
        kept, skipped, counts = strip_analysis_only_cards(sample)
        # GRIDs / CQUAD4 / FORCE survive
        assert 'GRID,1' in kept
        assert 'CQUAD4,1' in kept
        assert 'FORCE,1' in kept
        # TEMPs are out of kept and in the stash
        assert 'TEMP,100,1' not in kept
        assert 'TEMP,100,1' in skipped
        assert 'TEMP,100,2' in skipped
        assert counts['TEMP'] == 2

    def test_strip_includes_continuation_lines_with_their_card(self):
        """A continuation line that follows a skipped card should also
        end up in the stash, not in kept."""
        sample = (
            'GRID,1,,0.0,0.0,0.0\n'
            'TABLED1,42,0.0\n'
            '+,1.0,1.5,2.0,2.5\n'        # continuation of TABLED1
            '+,3.0,3.5,4.0,4.5\n'        # continuation of TABLED1
            'GRID,2,,1.0,0.0,0.0\n'
        )
        kept, skipped, counts = strip_analysis_only_cards(sample)
        assert '+,1.0' in skipped
        assert '+,3.0' in skipped
        assert '+,1.0' not in kept
        assert '+,3.0' not in kept
        # The GRIDs straddling the table both survive
        assert 'GRID,1' in kept
        assert 'GRID,2' in kept

    def test_strip_empty_when_no_analysis_cards(self):
        """If there's nothing to skip, kept_text == flat_text exactly."""
        sample = (
            'GRID,1,,0.0,0.0,0.0\n'
            'GRID,2,,1.0,0.0,0.0\n'
            'CQUAD4,1,1,1,2,1,2\n'
        )
        kept, skipped, counts = strip_analysis_only_cards(sample)
        assert kept == sample
        assert skipped == ''
        assert dict(counts) == {}

    def test_strip_counts_match(self):
        """The counts match the number of stripped card-start lines."""
        sample = ''.join(
            f'TEMP,100,{i},{72.0 + i * 0.1:.1f}\n' for i in range(1, 1001)
        )
        kept, skipped, counts = strip_analysis_only_cards(sample)
        assert counts['TEMP'] == 1000
        assert kept == ''
        assert skipped.count('TEMP,100') == 1000


class TestSkipCardsRoundTrip:
    """End-to-end: import a TEMP-heavy deck, save, re-import, verify."""

    def test_temp_cards_preserved_through_export(self, tmp_path):
        deck = ['SOL 101\nCEND\nBEGIN BULK\n']
        deck.append('GRID,1,,0.0,0.0,0.0\n')
        deck.append('GRID,2,,1.0,0.0,0.0\n')
        deck.append('GRID,3,,0.0,1.0,0.0\n')
        # 50 TEMP cards
        for i in range(1, 51):
            for nid in (1, 2, 3):
                deck.append(f'TEMP,{i},{nid},{72.0 + i * 0.1:.1f}\n')
        deck.append('ENDDATA\n')
        original = tmp_path / 'in.bdf'
        original.write_text(''.join(deck))

        # Import via streaming reader (the only path that strips).
        # Use _read_bdf_streaming directly so we don't depend on Qt.
        gen = NastranModelGenerator()
        # streaming requires INCLUDE detection - it picks the path
        # depending on file_has_includes. This deck has no INCLUDEs,
        # so we exercise the strip directly:
        from node_runner.skip_cards import strip_analysis_only_cards
        flat = original.read_text()
        kept, skipped, counts = strip_analysis_only_cards(flat)
        assert counts['TEMP'] == 150, \
            f"expected 150 TEMP, got {counts.get('TEMP', 0)}"

        # Read the kept-only buffer through pyNastran
        kept_path = tmp_path / 'kept.bdf'
        kept_path.write_text(kept)
        m, _ = NastranModelGenerator._read_bdf_robust(str(kept_path))
        # The kept model should NOT have any TEMP loads (we stripped them).
        all_loads = [
            l for sid_list in m.loads.values()
            for l in sid_list
            if hasattr(l, 'type') and l.type == 'TEMP'
        ]
        assert len(all_loads) == 0, "TEMP cards leaked through the strip"

        # Stash the raw lines on the model and write it out.
        m._skipped_raw_lines = skipped
        m._skipped_card_counts = dict(counts)
        gen.model = m
        out_path = tmp_path / 'out.bdf'
        gen._write_bdf(str(out_path), field_format='short')

        # The output should contain the TEMP cards we stashed.
        out_text = out_path.read_text()
        # Count TEMP lines (look for card-start lines, not substring).
        temp_lines = [ln for ln in out_text.splitlines()
                      if _card_name(ln) == 'TEMP']
        assert len(temp_lines) == 150, \
            f"expected 150 TEMP lines on export, got {len(temp_lines)}"
        # And ENDDATA still terminates the deck.
        assert 'ENDDATA' in out_text
