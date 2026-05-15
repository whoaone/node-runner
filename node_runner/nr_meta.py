"""Node Runner metadata embedded in Nastran BDF comments.

v5.1.0 item 28 introduces the ``$ NR-META v1`` line protocol so groups
and tree state survive a BDF round-trip without needing a sidecar file.
Every line is a syntactically valid Nastran comment (``$`` at column 1)
which means MSC, NX, and MYSTRAN parse them as no-ops; the Node Runner
import path strips the marker prefix and re-hydrates the metadata.

Block shape::

    $ NR-META v1 BEGIN
    $ NR-META source=node_runner version=5.1.0 written=2026-05-14T16:48:53
    $ NR-META section=groups
    $ NR-META {"gid": 1, "name": "Skin", "nodes": [...], ...}
    $ NR-META {"gid": 2, "name": "Frames", ...}
    $ NR-META section=tree_state
    $ NR-META {"hidden_groups": ["Frames"], "isolate_mode": null}
    $ NR-META v1 END

Design notes:

- One JSON object per line so a single corrupt line doesn't poison the
  rest of the block.
- Every line carries the ``$ NR-META`` prefix so pyNastran's bulk-data
  parser sees ordinary comments and no state is required there.
- Section headers (``section=groups``, ``section=tree_state``) make the
  parser dispatcher trivial: read until next ``section=`` or ``END``,
  parse each line as JSON.
- Versioned (``v1``); newer Node Runner versions can refuse unknown
  versions cleanly rather than misinterpret.

Per the v5.1.0 plan, entity-level color overrides are **deliberately not
persisted**. Properties and materials themselves already round-trip via
the standard ``PROP*`` / ``MAT*`` cards; the group-level ID lists below
are pure references.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from datetime import datetime
from typing import Iterable, Optional


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VERSION = "v1"
_PREFIX = "$ NR-META "
_BEGIN_MARKER = f"$ NR-META {_VERSION} BEGIN"
_END_MARKER = f"$ NR-META {_VERSION} END"
_SECTION_GROUPS = "section=groups"
_SECTION_TREE_STATE = "section=tree_state"
_HEADER_LINE_RE = re.compile(
    r"^\$\s*NR-META\s+source=([\w_-]+)\s+version=([\w.\-+]+)"
    r"\s+written=(\S+)\s*$"
)


@dataclass
class ParsedMeta:
    """Result of :func:`parse_meta` over a BDF text blob."""
    version: Optional[str] = None
    source: Optional[str] = None
    written: Optional[str] = None
    groups: dict[str, dict] = field(default_factory=dict)
    hidden_groups: list[str] = field(default_factory=list)
    isolate_mode: Optional[str] = None
    bad_lines: list[tuple[int, str, str]] = field(default_factory=list)

    @property
    def has_content(self) -> bool:
        return bool(self.groups) or bool(self.hidden_groups)


# ---------------------------------------------------------------------------
# dump
# ---------------------------------------------------------------------------

def dump_meta(
        groups: dict[str, dict],
        hidden_groups: Optional[Iterable[str]] = None,
        isolate_mode: Optional[str] = None,
        nr_version: str = "",
        ) -> str:
    """Return the embedded-block text (newline-terminated).

    ``groups`` is the mainwindow ``self.groups`` dict (name -> dict
    with keys gid / nodes / elements / properties / materials /
    coords). ``hidden_groups`` is a set / list of group names that are
    currently hidden. ``isolate_mode`` is the name of an isolated
    group, or None.

    Returns an empty string if there's no real content to emit, so
    callers can skip writing the block entirely.
    """
    hidden = sorted(set(hidden_groups or []))
    if not groups and not hidden and isolate_mode is None:
        return ""

    lines: list[str] = []
    lines.append(_BEGIN_MARKER)
    lines.append(
        f"{_PREFIX}source=node_runner version={nr_version or 'unknown'} "
        f"written={datetime.now().replace(microsecond=0).isoformat()}"
    )
    # ---- groups section ----
    if groups:
        lines.append(f"{_PREFIX}{_SECTION_GROUPS}")
        for name, data in sorted(groups.items(),
                                 key=lambda kv: kv[1].get("gid", 0)):
            obj = {
                "gid": int(data.get("gid", 0)),
                "name": str(name),
                "nodes": [int(n) for n in data.get("nodes") or []],
                "elements": [int(e) for e in data.get("elements") or []],
                "properties": [int(p) for p in data.get("properties") or []],
                "materials": [int(m) for m in data.get("materials") or []],
                "coords": [int(c) for c in data.get("coords") or []],
            }
            lines.append(_PREFIX + json.dumps(obj, separators=(",", ":")))
    # ---- tree_state section ----
    if hidden or isolate_mode is not None:
        lines.append(f"{_PREFIX}{_SECTION_TREE_STATE}")
        obj = {
            "hidden_groups": hidden,
            "isolate_mode": isolate_mode,
        }
        lines.append(_PREFIX + json.dumps(obj, separators=(",", ":")))
    lines.append(_END_MARKER)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# parse
# ---------------------------------------------------------------------------

def parse_meta(text: str) -> ParsedMeta:
    """Scan a BDF text blob for embedded ``$ NR-META`` block(s) and
    return a :class:`ParsedMeta`.

    Tolerant of unknown sections, bad JSON on individual lines (those
    are reported via ``ParsedMeta.bad_lines`` but don't abort), and
    blocks anywhere in the file. Only the first ``v1`` block is
    consumed; subsequent blocks are ignored (the export path only ever
    writes one).
    """
    meta = ParsedMeta()
    inside = False
    # Tracks the version of the currently-open block. When this is
    # something other than _VERSION (e.g. a v99 block from a newer
    # Node Runner), we consume the BEGIN/END markers but skip all
    # data lines so we don't misinterpret a future format.
    cur_block_version: Optional[str] = None
    cur_section: Optional[str] = None

    for lineno, raw in enumerate(text.splitlines(), 1):
        stripped = raw.rstrip()
        # Begin marker (any version).
        m_begin = re.match(r"^\$\s*NR-META\s+(v\d+)\s+BEGIN\s*$", stripped)
        if m_begin and not inside:
            inside = True
            cur_block_version = m_begin.group(1)
            meta.version = cur_block_version
            cur_section = None
            continue
        # End marker (any version).
        if inside and re.match(
                r"^\$\s*NR-META\s+v\d+\s+END\s*$", stripped):
            inside = False
            cur_block_version = None
            cur_section = None
            continue
        if not inside:
            continue
        # Inside an unknown-version block: ignore everything until END.
        if cur_block_version != _VERSION:
            continue

        if not stripped.startswith(_PREFIX):
            # Inside block, but no marker -- ignore (allows interleaved
            # blank/garbage lines without aborting).
            continue
        body = stripped[len(_PREFIX):].strip()

        # Header line (key=value pairs)?
        m = _HEADER_LINE_RE.match(stripped)
        if m:
            meta.source = m.group(1)
            meta.written = m.group(3)
            continue

        # Section switch?
        if body == _SECTION_GROUPS:
            cur_section = "groups"
            continue
        if body == _SECTION_TREE_STATE:
            cur_section = "tree_state"
            continue
        if body.startswith("section="):
            cur_section = body[len("section="):].strip() or None
            continue

        # JSON data line.
        if cur_section in ("groups", "tree_state"):
            try:
                obj = json.loads(body)
            except Exception as exc:
                meta.bad_lines.append((lineno, body[:80], str(exc)[:80]))
                continue
            if cur_section == "groups":
                _consume_group_obj(meta, obj)
            else:
                _consume_tree_state_obj(meta, obj)

    return meta


def _consume_group_obj(meta: ParsedMeta, obj: dict) -> None:
    if not isinstance(obj, dict):
        return
    name = obj.get("name")
    if not name:
        return
    meta.groups[str(name)] = {
        "gid": int(obj.get("gid", 0)),
        "nodes": [int(n) for n in obj.get("nodes") or []],
        "elements": [int(e) for e in obj.get("elements") or []],
        "properties": [int(p) for p in obj.get("properties") or []],
        "materials": [int(m) for m in obj.get("materials") or []],
        "coords": [int(c) for c in obj.get("coords") or []],
    }


def _consume_tree_state_obj(meta: ParsedMeta, obj: dict) -> None:
    if not isinstance(obj, dict):
        return
    hg = obj.get("hidden_groups")
    if isinstance(hg, list):
        meta.hidden_groups = [str(x) for x in hg]
    iso = obj.get("isolate_mode")
    if iso is None or isinstance(iso, str):
        meta.isolate_mode = iso


# ---------------------------------------------------------------------------
# strip helper -- used by the import path before pyNastran sees the deck
# ---------------------------------------------------------------------------

def strip_meta_block(text: str) -> tuple[str, str]:
    """Return ``(cleaned_text, meta_block_text)``.

    Removes the ``$ NR-META`` BEGIN..END block from ``text`` and returns
    it separately. pyNastran sees ``cleaned_text``; the caller passes
    ``meta_block_text`` to :func:`parse_meta` if it wants to restore
    Node Runner state.

    Tolerates blocks anywhere in the input. Only the first complete
    ``v1`` block is stripped; any further blocks are left in place
    (they'll continue to look like comments to pyNastran).
    """
    out_lines: list[str] = []
    meta_lines: list[str] = []
    inside = False
    for raw in text.splitlines(keepends=True):
        stripped = raw.rstrip("\r\n")
        if not inside and re.match(r"^\$\s*NR-META\s+v\d+\s+BEGIN\s*$",
                                   stripped):
            inside = True
            meta_lines.append(raw)
            continue
        if inside:
            meta_lines.append(raw)
            if re.match(r"^\$\s*NR-META\s+v\d+\s+END\s*$", stripped):
                inside = False
            continue
        out_lines.append(raw)
    return "".join(out_lines), "".join(meta_lines)
