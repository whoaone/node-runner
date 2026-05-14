"""v4.0.1 (Stage 1): diagnostic instrumentation.

Centralized perf-stage timer. Enabled by ``NR_PROFILE=1``. Stays
zero-overhead when disabled.

Three sinks (all active when enabled):
  - stdout (interleaves with the existing import-dialog progress).
  - status-bar callback (registered by MainWindow); mid-stage progress
    is visible to the user.
  - profile log file at ``~/.node_runner/profile_<timestamp>.log``.

Usage::

    from node_runner.profiling import perf_stage, perf_event

    with perf_stage('viewer', 'build_node_coords'):
        ...

    perf_event('shading', 'actor_count',
               renderer_actors=42, plotter_actors=7)
"""
from __future__ import annotations

import contextlib
import datetime
import os
import sys
import threading
import time
from pathlib import Path
from typing import Callable, Optional


_ENABLED: Optional[bool] = None
_LOG_PATH: Optional[Path] = None
_LOG_FH = None  # opened lazily
_STATUS_CB: Optional[Callable[[str], None]] = None
_LOCK = threading.Lock()


def enabled() -> bool:
    """Cache + return the env-var gate. Once set per-process."""
    global _ENABLED
    if _ENABLED is None:
        _ENABLED = os.environ.get('NR_PROFILE', '').strip() in ('1', 'true', 'True', 'TRUE')
    return _ENABLED


def _ensure_log() -> Optional[Path]:
    """Open the profile log file on first use. Returns the path or None."""
    global _LOG_PATH, _LOG_FH
    if not enabled():
        return None
    if _LOG_FH is not None:
        return _LOG_PATH
    try:
        d = Path.home() / '.node_runner'
        d.mkdir(parents=True, exist_ok=True)
        stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        _LOG_PATH = d / f'profile_{stamp}.log'
        _LOG_FH = open(_LOG_PATH, 'w', encoding='utf-8', buffering=1)
        # Header
        import platform
        _LOG_FH.write(f"# Node Runner profile log\n")
        _LOG_FH.write(f"# started: {datetime.datetime.now().isoformat()}\n")
        _LOG_FH.write(f"# platform: {platform.platform()}\n")
        _LOG_FH.write(f"# python: {sys.version.split()[0]}\n")
        _LOG_FH.write(f"# cores: {os.cpu_count()}\n")
        _LOG_FH.write(f"#\n")
    except Exception:
        _LOG_FH = None
        _LOG_PATH = None
    return _LOG_PATH


def set_status_callback(cb: Optional[Callable[[str], None]]) -> None:
    """Register the status-bar sink. MainWindow calls this in its ctor."""
    global _STATUS_CB
    _STATUS_CB = cb


def log_path() -> Optional[Path]:
    """Return the profile log path (or None if profiling is off)."""
    return _LOG_PATH


def _emit(line: str, *, to_status: bool = False) -> None:
    if not enabled():
        return
    _ensure_log()
    with _LOCK:
        try:
            print(line, flush=True)
        except Exception:
            pass
        if _LOG_FH is not None:
            try:
                _LOG_FH.write(line + '\n')
            except Exception:
                pass
    if to_status and _STATUS_CB is not None:
        try:
            _STATUS_CB(line)
        except Exception:
            pass


def perf_event(category: str, stage: str, **extra) -> None:
    """One-shot event (no timing). Used for actor counts, diagnostic
    snapshots, branch decisions."""
    if not enabled():
        return
    parts = [f"[perf] {category}.{stage}"]
    for k, v in extra.items():
        parts.append(f"{k}={v}")
    _emit('  '.join(parts))


@contextlib.contextmanager
def perf_stage(category: str, stage: str, *, status_msg: Optional[str] = None,
               status: bool = False, **extra):
    """Time a block. Prints ``[perf] <cat>.<stage>  wall=<s>``.

    Parameters
    ----------
    category : str
        Top-level bucket: 'viewer', 'visibility', 'tree', 'constraints',
        'shading', 'lod', etc.
    stage : str
        Sub-step name.
    status_msg : str, optional
        Human-readable status line for the import dialog / status bar.
        Emitted at the START of the stage (so the user sees what's
        running, not just what finished).
    status : bool
        If True, the completion line is also forwarded to the status sink.
    extra : kwargs
        Extra fields appended to the log line (e.g. n_cells=2400000).
    """
    if not enabled():
        yield
        return
    if status_msg and _STATUS_CB is not None:
        try:
            _STATUS_CB(status_msg)
        except Exception:
            pass
    t0 = time.perf_counter()
    try:
        yield
    finally:
        wall = time.perf_counter() - t0
        parts = [f"[perf] {category}.{stage}", f"wall={wall:.3f}s"]
        for k, v in extra.items():
            parts.append(f"{k}={v}")
        _emit('  '.join(parts), to_status=status)


# Convenience: top-of-file "deck header" event for the profile log
def log_deck_header(deck_path: str, n_nodes: int = 0, n_elements: int = 0,
                    n_includes: int = 0) -> None:
    if not enabled():
        return
    try:
        size_mb = os.path.getsize(deck_path) / 1024 / 1024
    except OSError:
        size_mb = 0.0
    perf_event(
        'deck', 'header',
        path=os.path.basename(deck_path),
        size_mb=f"{size_mb:.1f}",
        n_nodes=n_nodes,
        n_elements=n_elements,
        n_includes=n_includes,
    )
