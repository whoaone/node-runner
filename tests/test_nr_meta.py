"""v5.1.0 item 28: tests for the $ NR-META v1 line protocol."""

from __future__ import annotations

import pytest

from node_runner.nr_meta import (
    dump_meta, parse_meta, strip_meta_block, ParsedMeta,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _sample_groups() -> dict:
    return {
        "Skin": {
            "gid": 1,
            "nodes": [1, 2, 3, 4],
            "elements": [101, 102, 103],
            "properties": [1, 2],
            "materials": [1],
            "coords": [],
        },
        "Frames": {
            "gid": 2,
            "nodes": [5, 6, 7],
            "elements": [201, 202],
            "properties": [3],
            "materials": [1],
            "coords": [10],
        },
    }


# ---------------------------------------------------------------------------
# Round trip
# ---------------------------------------------------------------------------

class TestNrMetaRoundTrip:

    def test_dump_then_parse_preserves_groups(self):
        groups = _sample_groups()
        text = dump_meta(groups, hidden_groups={"Frames"},
                        nr_version="5.1.0-dev")
        meta = parse_meta(text)
        assert isinstance(meta, ParsedMeta)
        assert meta.version == "v1"
        assert set(meta.groups.keys()) == {"Skin", "Frames"}
        assert meta.groups["Skin"]["elements"] == [101, 102, 103]
        assert meta.groups["Skin"]["properties"] == [1, 2]
        assert meta.groups["Frames"]["coords"] == [10]
        assert meta.hidden_groups == ["Frames"]
        assert meta.isolate_mode is None
        assert not meta.bad_lines

    def test_empty_input_returns_empty_block(self):
        text = dump_meta({}, hidden_groups=set(), isolate_mode=None)
        assert text == ""

    def test_hidden_only_emits_tree_state_section(self):
        text = dump_meta({}, hidden_groups={"A", "B"},
                        nr_version="5.1.0-dev")
        assert "section=tree_state" in text
        meta = parse_meta(text)
        assert meta.hidden_groups == ["A", "B"]

    def test_isolate_mode_round_trips(self):
        text = dump_meta(_sample_groups(),
                        hidden_groups=set(),
                        isolate_mode="Skin",
                        nr_version="5.1.0-dev")
        meta = parse_meta(text)
        assert meta.isolate_mode == "Skin"


# ---------------------------------------------------------------------------
# Strip helper -- the import-pipeline entry point
# ---------------------------------------------------------------------------

class TestStripMetaBlock:

    def test_strip_in_middle_of_bdf(self):
        meta = dump_meta(_sample_groups(), hidden_groups={"Frames"},
                         nr_version="5.1.0-dev")
        bdf = (
            "$ Bulk data\n"
            "GRID,1,0,0.0,0.0,0.0\n"
            + meta
            + "GRID,2,0,1.0,0.0,0.0\n"
            "ENDDATA\n"
        )
        cleaned, block = strip_meta_block(bdf)
        # Cleaned text has no NR-META lines.
        assert "NR-META" not in cleaned
        # pyNastran-visible content survived.
        assert "GRID,1,0,0.0,0.0,0.0" in cleaned
        assert "GRID,2,0,1.0,0.0,0.0" in cleaned
        assert "ENDDATA" in cleaned
        # Block round-trips.
        meta_parsed = parse_meta(block)
        assert set(meta_parsed.groups.keys()) == {"Skin", "Frames"}

    def test_strip_when_no_block_present(self):
        bdf = "GRID,1,0,0.0,0.0,0.0\nENDDATA\n"
        cleaned, block = strip_meta_block(bdf)
        assert cleaned == bdf
        assert block == ""


# ---------------------------------------------------------------------------
# Robustness
# ---------------------------------------------------------------------------

class TestNrMetaRobustness:

    def test_bad_json_line_is_reported_but_doesnt_abort(self):
        bdf_block = (
            "$ NR-META v1 BEGIN\n"
            "$ NR-META source=node_runner version=5.1.0 written=2026-05-14T00:00:00\n"
            "$ NR-META section=groups\n"
            "$ NR-META {\"gid\": 1, \"name\": \"OK\", \"nodes\": [1]}\n"
            "$ NR-META {malformed json here\n"
            "$ NR-META {\"gid\": 2, \"name\": \"AlsoOK\", \"nodes\": [2]}\n"
            "$ NR-META v1 END\n"
        )
        meta = parse_meta(bdf_block)
        assert len(meta.bad_lines) == 1
        # Good lines on either side of the bad one still landed.
        assert set(meta.groups.keys()) == {"OK", "AlsoOK"}

    def test_unknown_section_is_ignored(self):
        bdf_block = (
            "$ NR-META v1 BEGIN\n"
            "$ NR-META section=unknown_future_section\n"
            "$ NR-META {\"some\":\"future_data\"}\n"
            "$ NR-META section=groups\n"
            "$ NR-META {\"gid\": 1, \"name\": \"X\", \"nodes\": []}\n"
            "$ NR-META v1 END\n"
        )
        meta = parse_meta(bdf_block)
        assert "X" in meta.groups

    def test_future_version_block_is_skipped(self):
        bdf_block = (
            "$ NR-META v99 BEGIN\n"
            "$ NR-META section=groups\n"
            "$ NR-META {\"gid\": 1, \"name\": \"FromFuture\", \"nodes\": []}\n"
            "$ NR-META v99 END\n"
        )
        meta = parse_meta(bdf_block)
        # The future block is recognised version-wise but its content
        # is not consumed.
        assert meta.version == "v99"
        assert meta.groups == {}

    def test_block_without_begin_marker_is_no_op(self):
        # Lines that look NR-META-ish but aren't inside a BEGIN..END
        # block should be ignored entirely.
        text = "$ NR-META {\"gid\":1,\"name\":\"ghost\",\"nodes\":[]}\n"
        meta = parse_meta(text)
        assert meta.groups == {}


# ---------------------------------------------------------------------------
# Bulk-data compatibility -- the block must look like comments to pyNastran
# ---------------------------------------------------------------------------

class TestBdfCompatibility:

    def test_every_emitted_line_starts_with_dollar(self):
        """Per the v5.1.0 plan: every NR-META line must be a syntactically
        valid Nastran comment so MSC / NX / MYSTRAN ignore the block.
        Smoke check by walking the dumped output."""
        text = dump_meta(_sample_groups(), hidden_groups={"Frames"},
                        nr_version="5.1.0-dev")
        for line in text.splitlines():
            assert line.startswith("$"), (
                f"NR-META line must start with $ to be a Nastran comment: "
                f"{line!r}")
