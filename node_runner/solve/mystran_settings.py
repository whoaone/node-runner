"""QSettings helpers for the MYSTRAN integration.

Centralizes the ``mystran/*`` key prefix so direct QSettings calls don't
get scattered through the solve package. Mirrors the pattern used by
:mod:`node_runner.dialogs.preferences` for v5.0.0 brightness +
materials-library settings.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

try:
    from PySide6.QtCore import QSettings
except Exception:  # pragma: no cover -- PySide6 always present at runtime
    QSettings = None  # type: ignore


# Defaults are applied when a key is absent. They mirror what the
# Preferences MYSTRAN tab displays as the "factory" state.
DEFAULTS = {
    "executable_path":      "",
    "scratch_root":         str(Path.home() / ".node_runner" / "mystran_runs"),
    "keep_intermediates":   True,
    # MYSTRAN's SOLLIB only accepts BANDED or SPARSE (verified against
    # MYSTRAN v18.0.0). SPARSE is the right default for typical decks
    # because it scales to large DOF counts; BANDED is faster for tiny
    # tightly-banded structural decks but blows up on most production
    # models.
    "default_sollib":       "SPARSE",   # SPARSE | BANDED
    "default_quad4typ":     "MIN4T",    # MIN4T | MIN4
    "default_wtmass":       1.0,
    "auto_load_results":    True,
    "auto_open_browser":    True,
    "cleanup_after_days":   30,         # 0 = never
}


def _qs():
    if QSettings is None:  # pragma: no cover
        return None
    return QSettings("NodeRunner", "NodeRunner")


def get(key: str, default=None):
    """Read ``mystran/<key>`` with type-correct fallback to DEFAULTS.

    Includes a one-shot migration for the v5.0.0 first-build SOLLIB
    default which was incorrectly "IntMKL" -- MYSTRAN only accepts
    "SPARSE" or "BANDED". Any persisted IntMKL/LAPACK value is
    silently rewritten to the new default.
    """
    qs = _qs()
    if qs is None:
        return DEFAULTS.get(key, default)
    fallback = DEFAULTS.get(key, default)
    if key == "default_sollib":
        raw = qs.value("mystran/default_sollib", fallback, type=str)
        if str(raw).upper() not in ("SPARSE", "BANDED"):
            qs.setValue("mystran/default_sollib", "SPARSE")
            return "SPARSE"
        return str(raw).upper()
    if isinstance(fallback, bool):
        # QSettings.value returns strings on some platforms; coerce.
        v = qs.value(f"mystran/{key}", fallback)
        if isinstance(v, str):
            return v.lower() in ("true", "1", "yes")
        return bool(v)
    if isinstance(fallback, int):
        try:
            return int(qs.value(f"mystran/{key}", fallback))
        except (TypeError, ValueError):
            return fallback
    if isinstance(fallback, float):
        try:
            return float(qs.value(f"mystran/{key}", fallback))
        except (TypeError, ValueError):
            return fallback
    return qs.value(f"mystran/{key}", fallback, type=str)


def set_value(key: str, value):
    qs = _qs()
    if qs is None:
        return
    qs.setValue(f"mystran/{key}", value)


def get_all() -> dict:
    """Return every MYSTRAN setting (used by Preferences load)."""
    return {k: get(k) for k in DEFAULTS}


def save_all(payload: dict):
    """Persist every key in payload that we recognize."""
    qs = _qs()
    if qs is None:
        return
    for key in DEFAULTS:
        if key in payload:
            qs.setValue(f"mystran/{key}", payload[key])
    qs.sync()


def get_executable_path() -> Optional[Path]:
    """Returns the user-configured MYSTRAN exe path, or None."""
    p = get("executable_path", "")
    if p:
        return Path(p)
    # Env var fallback (runner does its own discovery beyond this).
    env = os.environ.get("MYSTRAN_EXE", "").strip()
    return Path(env) if env else None


def get_scratch_root() -> Path:
    p = Path(get("scratch_root", DEFAULTS["scratch_root"]))
    return p
