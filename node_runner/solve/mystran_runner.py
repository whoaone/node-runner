"""QProcess-based MYSTRAN runner.

Public surface:

  - :func:`discover_mystran_executable` -> Path | None
  - :func:`detect_mystran_version`      -> str
  - :class:`SolverRunWorker`            QObject that wraps a QProcess
                                        and emits progress signals
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    from PySide6.QtCore import (
        QObject, QProcess, QTimer, Signal, QByteArray,
    )
except Exception:  # pragma: no cover -- PySide6 always present at runtime
    QObject = object
    QProcess = None
    QTimer = None
    Signal = lambda *a, **kw: None  # type: ignore
    QByteArray = None

try:
    from node_runner.profiling import perf_event
except Exception:  # pragma: no cover
    def perf_event(*_a, **_kw):
        pass

from node_runner.solve import MystranRun
from node_runner.solve import mystran_settings


# ---------------------------------------------------------------------------
# Executable discovery
# ---------------------------------------------------------------------------

_COMMON_PATHS: list[Path] = [
    Path(r"C:\MYSTRAN\mystran.exe"),
    Path(r"C:\Program Files\MYSTRAN\mystran.exe"),
    Path(r"C:\Program Files (x86)\MYSTRAN\mystran.exe"),
    Path("/usr/local/bin/mystran"),
    Path("/opt/mystran/mystran"),
    Path.home() / "MYSTRAN" / "mystran.exe",
    Path.home() / "mystran" / "mystran",
]


def discover_mystran_executable() -> Optional[Path]:
    """Try, in order:

    1. ``mystran/executable_path`` QSetting
    2. ``MYSTRAN_EXE`` env var
    3. ``shutil.which("mystran")`` / ``shutil.which("mystran.exe")``
    4. Platform-specific common install locations

    Returns the first existing path or ``None``.
    """
    saved = mystran_settings.get_executable_path()
    if saved and saved.is_file():
        perf_event('mystran.exe', 'discovered',
                   path=str(saved), source='qsettings')
        return saved

    env = os.environ.get("MYSTRAN_EXE", "").strip()
    if env:
        p = Path(env)
        if p.is_file():
            perf_event('mystran.exe', 'discovered',
                       path=str(p), source='env')
            return p

    for name in ("mystran", "mystran.exe", "MYSTRAN"):
        found = shutil.which(name)
        if found:
            p = Path(found)
            perf_event('mystran.exe', 'discovered',
                       path=str(p), source='path')
            return p

    for candidate in _COMMON_PATHS:
        if candidate.is_file():
            perf_event('mystran.exe', 'discovered',
                       path=str(candidate), source='common_location')
            return candidate

    perf_event('mystran.exe', 'not_found')
    return None


_VERSION_RE = re.compile(
    r"MYSTRAN\s+(?:VERSION|VER\.?|v)?\s*([0-9][0-9.]*[A-Za-z0-9]*)",
    re.IGNORECASE,
)


def detect_mystran_version(exe_path) -> str:
    """Run ``<exe> --version`` (1s timeout) and parse the banner.

    Returns a short string like ``"MYSTRAN 19.3"`` or, if the version
    flag isn't supported by this build, the first stdout line. Failures
    return an empty string -- the caller treats that as "found, but
    couldn't probe".
    """
    p = Path(exe_path)
    if not p.is_file():
        return ""
    for args in (["--version"], ["-v"], []):
        try:
            cp = subprocess.run(
                [str(p)] + args,
                capture_output=True, text=True, timeout=1.5,
                # MYSTRAN's no-args mode prompts for stdin -- close it
                stdin=subprocess.DEVNULL,
            )
            text = (cp.stdout or "") + "\n" + (cp.stderr or "")
            m = _VERSION_RE.search(text)
            if m:
                return f"MYSTRAN {m.group(1)}"
            # Fall back to first non-empty banner line.
            for line in text.splitlines():
                line = line.strip()
                if "MYSTRAN" in line.upper():
                    return line[:80]
        except (subprocess.TimeoutExpired, OSError):
            continue
    return ""


# ---------------------------------------------------------------------------
# Stage detection from MYSTRAN stdout
# ---------------------------------------------------------------------------

# (compiled regex, user-facing stage label, fraction estimate)
_STAGE_TABLE: list[tuple[re.Pattern, str, float]] = [
    (re.compile(r"READING\s+BULK\s+DATA", re.IGNORECASE),
     "Reading bulk data", 0.10),
    (re.compile(r"PROCESSING\s+CASE\s+CONTROL", re.IGNORECASE),
     "Processing case control", 0.15),
    (re.compile(r"ASSEMBLING\s+STIFFNESS", re.IGNORECASE),
     "Assembling stiffness matrix", 0.30),
    (re.compile(r"ASSEMBLING\s+MASS", re.IGNORECASE),
     "Assembling mass matrix", 0.45),
    (re.compile(r"(DECOMPOSING|FACTOR(IZ|IS)ING)", re.IGNORECASE),
     "Factorizing", 0.60),
    (re.compile(r"SOLVING\b", re.IGNORECASE),
     "Solving", 0.80),
    (re.compile(r"RECOVERING\s+(STRESSES|STRAINS|FORCES)", re.IGNORECASE),
     "Recovering element results", 0.90),
    (re.compile(r"EIGENVALUE\s+EXTRACTION", re.IGNORECASE),
     "Extracting eigenvalues", 0.70),
    (re.compile(r"WRITING\s+(F06|OP2|OUTPUT)", re.IGNORECASE),
     "Writing output", 0.95),
    (re.compile(r"END\s+OF\s+MYSTRAN", re.IGNORECASE),
     "Finishing", 0.98),
]


def _classify_stage(line: str) -> Optional[tuple[str, float]]:
    """Return ``(label, fraction)`` if the line matches a known phase."""
    for pat, label, frac in _STAGE_TABLE:
        if pat.search(line):
            return label, frac
    return None


# ---------------------------------------------------------------------------
# QProcess worker
# ---------------------------------------------------------------------------

class SolverRunWorker(QObject):
    """Drives one MYSTRAN run via QProcess.

    Signals:
        progress(stage: str, detail: str, fraction: float)
        finished(MystranRun)
        failed(str)
        cancelled()

    Typical lifecycle:
        worker = SolverRunWorker(parent)
        worker.progress.connect(dialog.update)
        worker.finished.connect(on_done)
        worker.failed.connect(on_fail)
        worker.start(bdf_path=..., exe_path=..., output_dir=...)
    """

    progress = Signal(str, str, float)
    finished = Signal(object)   # MystranRun
    failed = Signal(str)
    cancelled = Signal()

    # ----- construction -----
    def __init__(self, parent=None):
        super().__init__(parent)
        self._proc: Optional[QProcess] = None
        self._run: Optional[MystranRun] = None
        self._started_wall: Optional[float] = None
        self._latest_stage = "Starting"
        self._latest_detail = ""
        self._latest_fraction = 0.0
        self._cancelled = False
        # 5-s SIGTERM -> SIGKILL escalation timer.
        self._kill_timer: Optional[QTimer] = None

    # ----- public API -----
    def start(self, *, bdf_path: Path, exe_path: Path, output_dir: Path,
              sol: int = 101, analysis_set_name: str = "") -> None:
        """Spawn MYSTRAN as a subprocess in ``output_dir`` against
        ``bdf_path``. Emits ``progress``, then ``finished`` (success) or
        ``failed`` / ``cancelled``.
        """
        import time as _time
        self._cancelled = False
        self._started_wall = _time.perf_counter()
        bdf_path = Path(bdf_path)
        exe_path = Path(exe_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        log_path = output_dir / "input.log"
        f06_guess = bdf_path.with_suffix(".F06")
        op2_guess = bdf_path.with_suffix(".OP2")

        self._run = MystranRun(
            bdf_path=bdf_path,
            f06_path=f06_guess if f06_guess.exists() else None,
            op2_path=op2_guess if op2_guess.exists() else None,
            log_path=log_path,
            started_at=datetime.now(),
            sol=sol,
            analysis_set_name=analysis_set_name,
        )

        perf_event('mystran.run', 'starting',
                   exe=str(exe_path), bdf=str(bdf_path), sol=sol)

        if QProcess is None:  # pragma: no cover
            self.failed.emit("Qt is not available; cannot start subprocess.")
            return

        self._proc = QProcess(self)
        self._proc.setWorkingDirectory(str(output_dir))
        # Merge stdout + stderr so MYSTRAN's stderr doesn't race ahead.
        self._proc.setProcessChannelMode(QProcess.MergedChannels)
        self._proc.readyReadStandardOutput.connect(self._on_ready_read)
        self._proc.finished.connect(self._on_finished)
        self._proc.errorOccurred.connect(self._on_error)

        # MYSTRAN's command-line is typically: mystran <input.bdf>
        # Some builds also accept --silent or similar; we leave it
        # default so the stdout stream stays informative.
        self._proc.start(str(exe_path), [str(bdf_path)])
        self.progress.emit("Starting MYSTRAN", str(bdf_path.name), 0.02)

    def cancel(self) -> None:
        """Terminate the child process; escalate to kill() after 5 s."""
        if self._proc is None or self._cancelled:
            return
        self._cancelled = True
        try:
            self._proc.terminate()
        except Exception:
            pass
        if QTimer is not None:
            self._kill_timer = QTimer(self)
            self._kill_timer.setSingleShot(True)
            self._kill_timer.timeout.connect(self._escalate_kill)
            self._kill_timer.start(5000)

    # ----- internal slots -----
    def _escalate_kill(self):
        if self._proc is None:
            return
        try:
            if self._proc.state() != QProcess.NotRunning:
                self._proc.kill()
        except Exception:
            pass

    def _on_ready_read(self):
        if self._proc is None or self._run is None:
            return
        try:
            data: QByteArray = self._proc.readAllStandardOutput()
            text = bytes(data).decode("utf-8", errors="replace")
        except Exception:
            return
        # Strip a UTF-8 BOM if the Windows build emitted one.
        if text.startswith("﻿"):
            text = text.lstrip("﻿")
        # Append to log file (best-effort).
        try:
            with open(self._run.log_path, "a", encoding="utf-8") as fh:
                fh.write(text)
        except Exception:
            pass
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped:
                continue
            self._latest_detail = stripped[:200]
            stage = _classify_stage(stripped)
            if stage is not None:
                label, frac = stage
                if frac > self._latest_fraction:
                    self._latest_fraction = frac
                self._latest_stage = label
                # Stage transition catcher.
                import time as _time
                wall = (_time.perf_counter() - self._started_wall
                        if self._started_wall else 0.0)
                self._run.stage_timings[label] = round(wall, 3)
                perf_event('mystran.run', 'stage',
                           stage=label, wall_s=round(wall, 3))
            self.progress.emit(self._latest_stage, self._latest_detail,
                               self._latest_fraction)

    def _on_finished(self, exit_code: int, exit_status):
        import time as _time
        if self._run is None:
            return
        wall = (_time.perf_counter() - self._started_wall
                if self._started_wall else 0.0)
        self._run.wall_time = round(wall, 3)
        self._run.return_code = int(exit_code)

        # Re-resolve F06 / OP2 -- MYSTRAN writes them next to the BDF.
        for attr, ext in (("f06_path", ".F06"), ("op2_path", ".OP2")):
            guess = self._run.bdf_path.with_suffix(ext)
            if guess.exists():
                setattr(self._run, attr, guess)
            else:
                lower = self._run.bdf_path.with_suffix(ext.lower())
                if lower.exists():
                    setattr(self._run, attr, lower)

        if self._cancelled:
            perf_event('mystran.run', 'cancelled',
                       after_wall_s=round(wall, 3),
                       stage=self._latest_stage)
            self.cancelled.emit()
            return

        if exit_code != 0:
            err = (f"MYSTRAN exited with code {exit_code} after "
                   f"{wall:.1f}s. See log: {self._run.log_path}")
            self._run.errors.append(err)
            perf_event('mystran.run', 'failed',
                       wall_s=round(wall, 3), return_code=int(exit_code))
            self.failed.emit(err)
            return

        perf_event('mystran.run', 'total_wall',
                   wall_s=round(wall, 3), sol=self._run.sol,
                   return_code=int(exit_code))
        # Persist run_meta.json next to the BDF for Analysis History.
        try:
            self._run.save_meta(self._run.bdf_path.parent)
        except Exception:
            pass
        self.finished.emit(self._run)

    def _on_error(self, err):
        if self._run is None:
            return
        msg = (f"QProcess error: {err}. Check that {self._run.bdf_path} "
               f"and the MYSTRAN executable are accessible.")
        self._run.errors.append(msg)
        perf_event('mystran.run', 'process_error', err=str(err))
        # finished() will fire shortly with non-zero exit code in most
        # cases; if not, surface failure now.
        if self._proc is not None and self._proc.state() == QProcess.NotRunning:
            self.failed.emit(msg)
