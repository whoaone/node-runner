"""v4.0.0 (Phase D): per-INCLUDE parallel BDF parser.

Runs pyNastran's ``BDF.read_bdf`` on each INCLUDE file in its own
subprocess (``ProcessPoolExecutor``) so the per-card Python overhead is
GIL-immune. The master process merges the per-include dicts into a
single BDF object.

This replaces the v3.x ``_read_bdf_parallel`` experiment at
``model.py:_read_bdf_parallel`` which used ``ThreadPoolExecutor``; that
path was retired because Python threads cannot beat the GIL for
per-card Python work.

Engagement is gated by ``should_use_parallel(filepath)`` so small decks
keep the existing fast serial path. Heuristic:

  - total INCLUDE-file size > ``_PARALLEL_MIN_BYTES``
  - n_includes >= ``_PARALLEL_MIN_INCLUDES``
  - DMIG / superelement cards absent (cross-INCLUDE refs)

On any error the caller is expected to fall back to
``_read_bdf_streaming``; ``parse_bdf_parallel`` returns
``(None, None)`` on unrecoverable failure so the caller can detect.
"""
from __future__ import annotations

import concurrent.futures
import multiprocessing
import os
from typing import Callable, Optional, Tuple

_PARALLEL_MIN_BYTES = 50 * 1024 * 1024     # 50 MB
_PARALLEL_MIN_INCLUDES = 4
_PARALLEL_MAX_WORKERS = 8
_DMIG_SENTINELS = ('DMIG', 'DMIAX', 'DMI', 'DMIK', 'SUPER', 'SEBULK')


def should_use_parallel(filepath: str) -> bool:
    """Return True iff the deck is worth a parallel parse.

    Reads include chain via the existing
    ``NastranModelGenerator._gather_include_chain`` and short-circuits
    on DMIG / superelement detection (cross-INCLUDE refs that the
    per-file workers can't resolve without xref).
    """
    try:
        from node_runner.model import NastranModelGenerator
    except ImportError:
        return False
    files, total_bytes, _ = NastranModelGenerator._gather_include_chain(
        filepath)
    if len(files) < _PARALLEL_MIN_INCLUDES + 1:  # +1 for master file
        return False
    if total_bytes < _PARALLEL_MIN_BYTES:
        return False
    # Cheap DMIG sniff: scan first 4KB of every reachable file.
    for f in files:
        try:
            with open(f, 'rb') as fh:
                head = fh.read(4096).decode('latin-1', errors='replace')
            up = head.upper()
            if any(s in up for s in _DMIG_SENTINELS):
                return False
        except OSError:
            return False
    return True


def _strip_includes_and_directives(text):
    """Drop INCLUDE statements and executive/case-control sections from
    a single file's text so the worker can parse it standalone in punch
    mode without pyNastran trying to recurse into siblings (which would
    fail in the worker because the working dir is wrong) or rejecting
    the file as not-a-deck.

    The INCLUDE chain is already enumerated by ``_gather_include_chain``
    in the master process; per-file workers only need this file's own
    bulk cards.
    """
    import io
    out = []
    in_continuation = False
    for raw in text.splitlines(keepends=False):
        up = raw.lstrip().upper()
        # Strip INCLUDE statements (single-line and continuation forms).
        if up.startswith('INCLUDE'):
            in_continuation = raw.rstrip().endswith(',')
            continue
        if in_continuation:
            # Continuation marker (+, *) or a line ending with comma.
            stripped = raw.lstrip()
            if stripped[:1] in ('+', '*') or raw.rstrip().endswith(','):
                in_continuation = raw.rstrip().endswith(',')
                continue
            in_continuation = False
        # Skip executive / case-control markers so punch parses cleanly.
        if up.startswith(('SOL ', 'CEND', 'BEGIN BULK', 'ENDDATA',
                          'TIME ', 'TITLE', 'SUBCASE', 'LABEL',
                          'ECHO', 'OUTPUT', 'DIAG')):
            continue
        out.append(raw)
    return '\n'.join(out) + '\n'


def _worker_parse_include(args):
    """Top-level worker function (module-level so it's picklable on
    Windows).

    Reads one BDF file in punch mode (so executive/case decks are
    optional) and returns a minimal pickle-friendly dict of the cards
    the viewer cares about. INCLUDE statements are stripped before
    parse so the worker never tries to recurse into sibling files
    (the master process enumerated them already and they get their
    own workers).
    """
    filepath, apply_strip = args
    import tempfile
    from pyNastran.bdf.bdf import BDF

    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
            raw = fh.read()
    except OSError as exc:
        return {'error': f"open failed: {exc}", 'filepath': filepath}

    text = _strip_includes_and_directives(raw)
    if apply_strip:
        try:
            from node_runner.skip_cards import strip_analysis_only_cards
            text, _skipped, _counts = strip_analysis_only_cards(text)
        except Exception:
            pass

    with tempfile.NamedTemporaryFile(
            mode='w', suffix='.bdf', delete=False,
            encoding='utf-8') as tmp:
        tmp.write(text)
        tmp_path = tmp.name

    bdf = BDF(debug=False)
    try:
        bdf.read_bdf(tmp_path, xref=False, validate=False, punch=True)
    except Exception as exc:
        return {'error': f"{type(exc).__name__}: {exc}", 'filepath': filepath}
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    return {
        'filepath': filepath,
        'nodes': dict(bdf.nodes),
        'elements': dict(bdf.elements),
        'rigid_elements': dict(bdf.rigid_elements),
        'properties': dict(bdf.properties),
        'materials': dict(bdf.materials),
        'coords': dict(bdf.coords),
        'loads': {k: list(v) for k, v in bdf.loads.items()},
        'load_combinations': {k: list(v) for k, v in
                              bdf.load_combinations.items()},
        'spcs': {k: list(v) for k, v in bdf.spcs.items()},
        'spcadds': {k: list(v) for k, v in bdf.spcadds.items()},
        'mpcs': {k: list(v) for k, v in bdf.mpcs.items()},
        'mpcadds': {k: list(v) for k, v in bdf.mpcadds.items()},
        'tempds': dict(bdf.tempds),
        'masses': dict(bdf.masses),
        'plotels': dict(bdf.plotels),
    }


def _merge_into_master(master, payload):
    """Merge one worker payload into the master BDF. Dict-shaped
    slots (nodes/elements/etc.) overwrite on collision (last writer
    wins). List-shaped slots (loads/spcs/etc.) extend per-SID.

    Collisions are silently ignored to match pyNastran's own behavior
    on duplicate cards across INCLUDEs (it warns but doesn't fail).
    """
    for slot in ('nodes', 'elements', 'rigid_elements', 'properties',
                 'materials', 'coords', 'tempds', 'masses', 'plotels'):
        getattr(master, slot).update(payload[slot])

    for slot in ('loads', 'load_combinations', 'spcs', 'spcadds',
                 'mpcs', 'mpcadds'):
        target = getattr(master, slot)
        for sid, cards in payload[slot].items():
            target.setdefault(sid, []).extend(cards)


def parse_bdf_parallel(
        filepath: str,
        n_workers: Optional[int] = None,
        progress: Optional[Callable] = None) -> Tuple[object, object]:
    """Parse a BDF deck using one worker process per INCLUDE file.

    Returns ``(model, lenient_or_None)`` matching the shape of the
    other ``_read_bdf_*`` entry points. On unrecoverable failure
    returns ``(None, None)`` so the caller can fall back to serial.
    """
    from node_runner.model import NastranModelGenerator, ImportProgress
    from pyNastran.bdf.bdf import BDF

    files, total_bytes, missing = (
        NastranModelGenerator._gather_include_chain(filepath))
    if not files:
        return (None, None)

    if n_workers is None:
        n_workers = min(_PARALLEL_MAX_WORKERS,
                        max(2, multiprocessing.cpu_count() // 2),
                        len(files))

    def _emit(stage, label, frac, source=''):
        if progress is None:
            return
        try:
            progress(ImportProgress(
                stage=stage, label=label, fraction=frac,
                source_file=source))
        except Exception:
            pass

    _emit('parse',
          f'Parallel parse: {len(files)} files across {n_workers} workers',
          0.05)

    work_items = [(f, True) for f in files]
    payloads = []
    completed = 0
    try:
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=n_workers) as pool:
            futures = {pool.submit(_worker_parse_include, item): item[0]
                       for item in work_items}
            for fut in concurrent.futures.as_completed(futures):
                payload = fut.result()
                if 'error' in payload:
                    # One worker failed -> bail out to serial.
                    return (None, None)
                payloads.append(payload)
                completed += 1
                _emit('parse',
                      f'Parsed {completed}/{len(files)} INCLUDE files',
                      0.1 + 0.7 * (completed / len(files)),
                      source=os.path.basename(payload['filepath']))
    except Exception:
        return (None, None)

    _emit('parse', 'Merging worker results...',
          0.85, source=os.path.basename(filepath))

    master = BDF(debug=False)
    for p in payloads:
        try:
            _merge_into_master(master, p)
        except Exception:
            return (None, None)

    _emit('parse', 'Parse complete', 1.0)
    return (master, None)
