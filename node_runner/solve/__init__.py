"""MYSTRAN solver integration package (v5.0.0).

Public surface is the dataclasses below + the four submodules:

  - :mod:`mystran_runner`   -- QProcess wrapper + executable discovery
  - :mod:`mystran_export`   -- MSC/NX BDF -> MYSTRAN-compatible BDF translator
  - :mod:`mystran_preflight` -- card / PARAM compatibility scanner
  - :mod:`mystran_results`  -- OP2-first / F06-fallback adapter
  - :mod:`mystran_settings` -- QSettings helpers + Preferences-tab key prefix

The dataclasses live here (not in the submodules) so tests can import them
without pulling in QtCore.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Literal, Optional


SUPPORTED_SOL_NUMBERS: tuple[int, ...] = (1, 3, 5, 101, 103, 105)


@dataclass
class MystranRun:
    """Bundles all the artifacts and metadata produced by one run."""

    bdf_path: Path
    f06_path: Optional[Path] = None
    op2_path: Optional[Path] = None
    log_path: Optional[Path] = None
    return_code: Optional[int] = None
    wall_time: Optional[float] = None
    stage_timings: dict[str, float] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    started_at: datetime = field(default_factory=datetime.now)
    sol: int = 101
    n_dof: Optional[int] = None
    analysis_set_name: str = ""

    # Whether the result-adapter took the OP2 path or fell back to F06.
    # Set by mystran_results.load_mystran_results; None until results
    # have been loaded.
    results_source: Optional[Literal["op2", "f06"]] = None

    def to_dict(self) -> dict:
        """JSON-serializable view for ``run_meta.json``."""
        d = asdict(self)
        d['bdf_path'] = str(self.bdf_path) if self.bdf_path else None
        d['f06_path'] = str(self.f06_path) if self.f06_path else None
        d['op2_path'] = str(self.op2_path) if self.op2_path else None
        d['log_path'] = str(self.log_path) if self.log_path else None
        d['started_at'] = self.started_at.isoformat()
        return d

    def save_meta(self, dir_path: Path) -> Path:
        """Write run_meta.json into dir_path. Returns the path written."""
        meta = dir_path / "run_meta.json"
        meta.write_text(json.dumps(self.to_dict(), indent=2), encoding='utf-8')
        return meta


@dataclass
class PreflightIssue:
    """One finding from the MYSTRAN compatibility scan."""

    severity: Literal["blocking", "warning", "info"]
    code: str             # stable short id, e.g. "AERO_CARD_PRESENT"
    message: str          # one-line human description
    detail: str = ""      # multi-line context
    card_type: Optional[str] = None
    card_count: Optional[int] = None


@dataclass
class PreflightReport:
    """Result of one ``scan_for_mystran(model)`` call."""

    issues: list[PreflightIssue] = field(default_factory=list)
    solution_supported: bool = False
    estimated_dof: int = 0
    estimated_wall_time_seconds: Optional[float] = None
    scan_wall_s: float = 0.0

    @property
    def blocking_count(self) -> int:
        return sum(1 for i in self.issues if i.severity == "blocking")

    @property
    def warning_count(self) -> int:
        return sum(1 for i in self.issues if i.severity == "warning")

    @property
    def info_count(self) -> int:
        return sum(1 for i in self.issues if i.severity == "info")

    @property
    def has_blocking(self) -> bool:
        return self.blocking_count > 0

    def issues_by_severity(
            self, severity: str) -> list[PreflightIssue]:
        return [i for i in self.issues if i.severity == severity]


__all__ = [
    "SUPPORTED_SOL_NUMBERS",
    "MystranRun",
    "PreflightIssue",
    "PreflightReport",
]
