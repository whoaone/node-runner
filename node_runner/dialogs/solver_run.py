"""MYSTRAN run-flow dialogs (v5.0.0 items 15 + 17).

Four dialogs:

  - :class:`MystranPreflightDialog` -- modal report of pre-flight issues
    before launching MYSTRAN. Disables [Run] when blocking issues exist.
  - :class:`RunMystranDialog` -- per-run overrides (SOL, output requests,
    working dir). Picks an existing AnalysisSet via the AnalysisSet
    manager rather than re-implementing that editor.
  - :class:`SolverProgressDialog` -- non-blocking progress dialog driven
    by SolverRunWorker signals.
  - :class:`AnalysisHistoryDialog` -- browses ``<scratch_root>/runs/``
    subfolders and re-loads results from a prior run.
"""

from __future__ import annotations

import os
import json
from datetime import datetime
from pathlib import Path
from typing import Optional

from PySide6 import QtCore
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QPushButton,
    QFrame, QTableWidget, QTableWidgetItem, QHeaderView, QCheckBox,
    QRadioButton, QButtonGroup, QGroupBox, QComboBox, QLineEdit,
    QFileDialog, QSpinBox, QDoubleSpinBox, QProgressBar, QMessageBox,
    QDialogButtonBox,
)

from node_runner.solve import (
    MystranRun, PreflightReport, PreflightIssue, SUPPORTED_SOL_NUMBERS,
)
from node_runner.solve import mystran_tips as tips


# ---------------------------------------------------------------------------
# 1. Pre-flight report dialog
# ---------------------------------------------------------------------------

_SEVERITY_BANNER = {
    "blocking": ("background-color: #f38ba8; color: #1e1e2e; "
                 "padding: 10px; border-radius: 4px; font-size: 13px;"),
    "warning":  ("background-color: #f9e2af; color: #1e1e2e; "
                 "padding: 10px; border-radius: 4px; font-size: 13px;"),
    "ok":       ("background-color: #a6e3a1; color: #1e1e2e; "
                 "padding: 10px; border-radius: 4px; font-size: 13px;"),
}
_SEVERITY_ROW = {
    "blocking": QColor("#f38ba8"),
    "warning":  QColor("#f9e2af"),
    "info":     QColor("#89b4fa"),
}


class MystranPreflightDialog(QDialog):
    """Modal pre-flight report. ``exec()`` returns Accepted iff the user
    chose to proceed."""

    def __init__(self, report: PreflightReport, parent=None):
        super().__init__(parent)
        self._report = report
        self.setWindowTitle("MYSTRAN Pre-flight Report")
        # v5.1.1 item 32: compact when there are zero issues; the
        # green-banner case used to leave a big empty space below the
        # banner. Sized for the worst-case (many issues) only when we
        # actually have a table to fill.
        if report.issues:
            self.setMinimumSize(720, 460)
            self.resize(820, 540)
        else:
            self.setMinimumSize(560, 260)
            self.resize(620, 300)

        layout = QVBoxLayout(self)

        # --- Banner --- (v5.1.1 item 32: copy now explains what
        # pre-flight does so first-time users aren't guessing.)
        scan_purpose = (
            "Scanned for cards MYSTRAN can't run "
            "(aero, nonlinear, contact, optimization, unsupported "
            "elements)."
        )
        if report.has_blocking:
            banner_text = (
                f"<b>Pre-flight cannot proceed.</b><br>"
                f"{scan_purpose} Found "
                f"{report.blocking_count} blocking issue"
                f"{'s' if report.blocking_count != 1 else ''}"
                f"{f' and {report.warning_count} warning' if report.warning_count else ''}"
                f"{'s' if report.warning_count > 1 else ''}. "
                f"Resolve the blocking items and try again."
            )
            banner_style = _SEVERITY_BANNER['blocking']
        elif report.warning_count:
            banner_text = (
                f"<b>Pre-flight passed with "
                f"{report.warning_count} warning"
                f"{'s' if report.warning_count != 1 else ''}.</b><br>"
                f"{scan_purpose} The MYSTRAN export will translate or "
                f"drop the items below; the run can proceed."
            )
            banner_style = _SEVERITY_BANNER['warning']
        else:
            banner_text = (
                f"<b>Pre-flight passed.</b><br>"
                f"{scan_purpose} None found. "
                f"~{report.estimated_dof:,} DOF; "
                f"estimated wall time "
                f"{_fmt_secs(report.estimated_wall_time_seconds)}."
            )
            banner_style = _SEVERITY_BANNER['ok']
        banner = QLabel(banner_text)
        banner.setStyleSheet(banner_style)
        banner.setWordWrap(True)
        banner.setToolTip(tips.PREFLIGHT_BANNER)
        layout.addWidget(banner)

        # --- Summary line ---
        summary_text = (
            f"Estimated DOF: <b>{report.estimated_dof:,}</b> | "
            f"Estimated wall: <b>"
            f"{_fmt_secs(report.estimated_wall_time_seconds)}</b> | "
            f"Scan: {report.scan_wall_s * 1000:.0f} ms"
        )
        summary = QLabel(summary_text)
        summary.setStyleSheet(
            "color: #94e2d5; font-family: Consolas, monospace; "
            "padding: 4px 0;")
        layout.addWidget(summary)

        # --- Issue table ---
        if report.issues:
            layout.addWidget(QLabel("<b>Issues:</b>"))
            table = QTableWidget(len(report.issues), 3)
            table.setHorizontalHeaderLabels(["Severity", "Code", "Message"])
            table.setEditTriggers(QTableWidget.NoEditTriggers)
            table.setAlternatingRowColors(True)
            table.verticalHeader().setVisible(False)
            table.horizontalHeader().setSectionResizeMode(
                0, QHeaderView.ResizeToContents)
            table.horizontalHeader().setSectionResizeMode(
                1, QHeaderView.ResizeToContents)
            table.horizontalHeader().setStretchLastSection(True)
            # Order: blocking, then warning, then info.
            sort_key = {"blocking": 0, "warning": 1, "info": 2}
            for row, issue in enumerate(
                    sorted(report.issues, key=lambda i: sort_key.get(
                        i.severity, 9))):
                sev_item = QTableWidgetItem(issue.severity.title())
                sev_color = _SEVERITY_ROW.get(issue.severity)
                if sev_color is not None:
                    sev_item.setForeground(sev_color)
                    f = sev_item.font()
                    f.setBold(True)
                    sev_item.setFont(f)
                table.setItem(row, 0, sev_item)
                table.setItem(row, 1, QTableWidgetItem(issue.code))
                msg_item = QTableWidgetItem(issue.message)
                msg_item.setToolTip(
                    f"{issue.message}\n\n{issue.detail}"
                    if issue.detail else issue.message)
                table.setItem(row, 2, msg_item)
            table.resizeColumnToContents(0)
            layout.addWidget(table)
        else:
            ok = QLabel("No issues found. Ready to run.")
            ok.setStyleSheet("color: #a6e3a1; padding: 8px 0; "
                             "font-style: italic;")
            layout.addWidget(ok)

        # --- Buttons ---
        button_row = QHBoxLayout()
        copy_btn = QPushButton("Copy to Clipboard")
        copy_btn.clicked.connect(self._copy_to_clipboard)
        button_row.addWidget(copy_btn)
        button_row.addStretch(1)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_row.addWidget(cancel_btn)
        self._run_btn = QPushButton(
            "Run Anyway" if report.warning_count else "Run")
        self._run_btn.setDefault(True)
        self._run_btn.setEnabled(not report.has_blocking)
        if report.has_blocking:
            self._run_btn.setToolTip(
                "Resolve the blocking issues above before running.")
        self._run_btn.clicked.connect(self.accept)
        button_row.addWidget(self._run_btn)
        layout.addLayout(button_row)

    def _copy_to_clipboard(self):
        try:
            from PySide6.QtWidgets import QApplication
            lines = [
                f"MYSTRAN Pre-flight Report",
                f"Estimated DOF: {self._report.estimated_dof:,}",
                f"Blocking: {self._report.blocking_count}, "
                f"Warnings: {self._report.warning_count}, "
                f"Info: {self._report.info_count}",
                "",
            ]
            for issue in self._report.issues:
                lines.append(
                    f"[{issue.severity.upper()}] {issue.code}: "
                    f"{issue.message}")
                if issue.detail:
                    lines.append(f"    {issue.detail}")
            QApplication.clipboard().setText("\n".join(lines))
        except Exception:
            pass


def _fmt_secs(s: Optional[float]) -> str:
    if s is None:
        return "?"
    if s < 1.0:
        return f"<1 s"
    if s < 60:
        return f"~{s:.0f} s"
    if s < 3600:
        return f"~{s / 60:.1f} min"
    return f"~{s / 3600:.1f} h"


# ---------------------------------------------------------------------------
# 2. Run-options dialog
# ---------------------------------------------------------------------------

class RunMystranDialog(QDialog):
    """Per-run override dialog. Returns selected options via
    :meth:`get_options`.
    """

    def __init__(self, main_window, parent=None):
        super().__init__(parent or main_window)
        self.main_window = main_window
        self.setWindowTitle("Run Analysis (MYSTRAN)")
        self.setMinimumWidth(540)

        # Load defaults from settings.
        from node_runner.solve import mystran_settings as _ms
        self._defaults = _ms.get_all()

        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 10, 10, 8)

        intro = QLabel(
            "Picks the AnalysisSet to use and lets you override SOL and "
            "output requests for this run only -- the AnalysisSet itself "
            "is not modified.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #cdd6f4; font-size: 11px;")
        outer.addWidget(intro)

        form = QFormLayout()

        # --- AnalysisSet picker ---
        self._set_combo = QComboBox()
        self._set_combo.setToolTip(tips.ANALYSIS_SET)
        self._refresh_set_combo()
        as_label = QLabel("AnalysisSet:")
        as_label.setToolTip(tips.ANALYSIS_SET)
        form.addRow(as_label, self._set_combo)
        edit_sets_btn = QPushButton("Edit AnalysisSets…")
        edit_sets_btn.setToolTip(
            "Open the AnalysisSet Manager. Each AnalysisSet bundles "
            "a SOL number, subcases, output requests, and PARAMs.")
        edit_sets_btn.clicked.connect(self._edit_sets)
        form.addRow(QLabel(""), edit_sets_btn)

        # --- SOL radio ---
        sol_box = QGroupBox("Solution")
        sol_box.setToolTip(
            "Pick the MYSTRAN solution sequence. SOL 1/101 = linear "
            "static, SOL 3/103 = normal modes, SOL 5/105 = buckling. "
            "Same numbering as MSC Nastran. MYSTRAN does NOT support "
            "SOL 106 (nonlinear), SOL 107-112, SOL 144, or SOL 200.")
        sol_lay = QHBoxLayout(sol_box)
        self._sol_group = QButtonGroup(self)
        self._sol_buttons: dict[int, QRadioButton] = {}
        sol_tips = {101: tips.SOL_101, 103: tips.SOL_103, 105: tips.SOL_105}
        for sol_num, label in [
                (101, "SOL 1 / 101 (Linear Static)"),
                (103, "SOL 3 / 103 (Normal Modes)"),
                (105, "SOL 5 / 105 (Buckling)")]:
            rb = QRadioButton(label)
            rb.setToolTip(sol_tips[sol_num])
            self._sol_group.addButton(rb)
            sol_lay.addWidget(rb)
            self._sol_buttons[sol_num] = rb
        self._sol_buttons[101].setChecked(True)
        form.addRow(sol_box)

        # --- Output requests ---
        out_box = QGroupBox("Output Requests")
        out_box.setToolTip(
            "Which result tables MYSTRAN should write. DISP is always "
            "useful; STRESS/STRAIN are needed for Color-By-Results "
            "contours; FORCE is for hand-checks + free-body diagrams.")
        out_lay = QHBoxLayout(out_box)
        self._out_disp = QCheckBox("DISP"); self._out_disp.setChecked(True)
        self._out_disp.setToolTip(tips.OUTPUT_DISP)
        self._out_stress = QCheckBox("STRESS"); self._out_stress.setChecked(True)
        self._out_stress.setToolTip(tips.OUTPUT_STRESS)
        self._out_strain = QCheckBox("STRAIN"); self._out_strain.setChecked(True)
        self._out_strain.setToolTip(tips.OUTPUT_STRAIN)
        self._out_force = QCheckBox("FORCE"); self._out_force.setChecked(True)
        self._out_force.setToolTip(tips.OUTPUT_FORCE)
        for cb in (self._out_disp, self._out_stress,
                   self._out_strain, self._out_force):
            out_lay.addWidget(cb)
        form.addRow(out_box)

        # --- Working dir ---
        wd_row = QHBoxLayout()
        scratch_root = Path(self._defaults.get(
            'scratch_root', str(Path.home() / '.node_runner' / 'mystran_runs')))
        bdf_stem = self._infer_bdf_stem()
        default_wd = (scratch_root / 'runs'
                      / f"{datetime.now():%Y%m%d_%H%M%S}_{bdf_stem}")
        self._wd_edit = QLineEdit(str(default_wd))
        self._wd_edit.setToolTip(tips.WORKING_DIR)
        wd_row.addWidget(self._wd_edit, 1)
        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._pick_wd)
        wd_row.addWidget(browse_btn)
        wd_label = QLabel("Working directory:")
        wd_label.setToolTip(tips.WORKING_DIR)
        form.addRow(wd_label, wd_row)

        outer.addLayout(form)

        # --- Buttons ---
        sep = QFrame(); sep.setFrameShape(QFrame.HLine)
        outer.addWidget(sep)
        btns = QDialogButtonBox()
        run_btn = btns.addButton("Run Pre-flight + Solve",
                                  QDialogButtonBox.AcceptRole)
        run_btn.setDefault(True)
        cancel_btn = btns.addButton(QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        outer.addWidget(btns)

    def _infer_bdf_stem(self) -> str:
        gen = getattr(self.main_window, 'current_generator', None)
        filepath = getattr(self.main_window, '_loaded_filepath', None)
        if filepath:
            return Path(filepath).stem
        return "untitled"

    def _refresh_set_combo(self):
        """Populate from main_window analysis sets if available.

        MainWindow stores them on ``self.analysis_sets`` (a dict keyed
        by set id; see _open_analysis_set_manager).
        """
        self._set_combo.clear()
        self._set_combo.addItem("(Default — model SOL + Output Requests)", None)
        try:
            sets = getattr(self.main_window, 'analysis_sets', None) or {}
            for sid, aset in sets.items():
                name = getattr(aset, 'name', f"Set {sid}")
                self._set_combo.addItem(f"{name} (id={sid})", sid)
        except Exception:
            pass

    def _edit_sets(self):
        try:
            self.main_window._open_analysis_set_manager()
        except Exception:
            QMessageBox.information(
                self, "AnalysisSet Manager",
                "Couldn't open the AnalysisSet manager from here. "
                "Use Analysis -> Analysis Set Manager... directly.")
        self._refresh_set_combo()

    def _pick_wd(self):
        from PySide6.QtCore import QStandardPaths
        start = (self._wd_edit.text()
                 or QStandardPaths.writableLocation(
                     QStandardPaths.HomeLocation))
        path = QFileDialog.getExistingDirectory(
            self, "Pick working directory", start)
        if path:
            self._wd_edit.setText(path)

    def get_options(self) -> dict:
        sol = 101
        for sol_num, rb in self._sol_buttons.items():
            if rb.isChecked():
                sol = sol_num
                break
        return {
            'analysis_set_id': self._set_combo.currentData(),
            'analysis_set_name': self._set_combo.currentText(),
            'sol': sol,
            'output_disp': bool(self._out_disp.isChecked()),
            'output_stress': bool(self._out_stress.isChecked()),
            'output_strain': bool(self._out_strain.isChecked()),
            'output_force': bool(self._out_force.isChecked()),
            'working_dir': self._wd_edit.text().strip(),
        }


# ---------------------------------------------------------------------------
# 3. Progress dialog
# ---------------------------------------------------------------------------

class SolverProgressDialog(QDialog):
    """Non-blocking progress dialog driven by SolverRunWorker signals."""

    cancel_requested = QtCore.Signal()

    def __init__(self, bdf_name: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Running MYSTRAN")
        self.setMinimumWidth(520)
        self.setModal(False)
        self._wall_start: Optional[float] = None

        layout = QVBoxLayout(self)
        self._header = QLabel(f"<b>Solving:</b> {bdf_name}")
        layout.addWidget(self._header)

        self._stage_label = QLabel("Starting MYSTRAN…")
        self._stage_label.setStyleSheet("font-weight: 500; font-size: 12px;")
        layout.addWidget(self._stage_label)

        self._detail_label = QLabel("")
        self._detail_label.setStyleSheet(
            "color: #6c7086; font-family: Consolas, monospace; "
            "font-size: 11px;")
        self._detail_label.setWordWrap(True)
        layout.addWidget(self._detail_label)

        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setValue(2)
        layout.addWidget(self._progress)

        self._wall_label = QLabel("Elapsed: 0.0 s")
        self._wall_label.setStyleSheet(
            "color: #94e2d5; font-family: Consolas, monospace;")
        layout.addWidget(self._wall_label)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self._cancel_btn = QPushButton("Cancel")
        self._cancel_btn.clicked.connect(self._on_cancel)
        button_row.addWidget(self._cancel_btn)
        layout.addLayout(button_row)

        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(200)
        self._timer.timeout.connect(self._tick)

    def start_timer(self):
        import time as _time
        self._wall_start = _time.perf_counter()
        self._timer.start()

    def stop_timer(self):
        self._timer.stop()

    # Slot connected to SolverRunWorker.progress
    def update_progress(self, stage: str, detail: str, fraction: float):
        self._stage_label.setText(stage)
        self._detail_label.setText(detail)
        try:
            self._progress.setValue(int(max(0.0, min(0.99, fraction)) * 100))
        except Exception:
            pass

    def _tick(self):
        import time as _time
        if self._wall_start is None:
            return
        wall = _time.perf_counter() - self._wall_start
        if wall < 60:
            self._wall_label.setText(f"Elapsed: {wall:.1f} s")
        else:
            self._wall_label.setText(
                f"Elapsed: {wall / 60:.1f} min ({wall:.0f} s)")

    def _on_cancel(self):
        self._cancel_btn.setEnabled(False)
        self._cancel_btn.setText("Cancelling…")
        self._stage_label.setText("Cancel requested — terminating MYSTRAN…")
        self.cancel_requested.emit()


# ---------------------------------------------------------------------------
# 4. Analysis History
# ---------------------------------------------------------------------------

class AnalysisHistoryDialog(QDialog):
    """Browse <scratch_root>/runs/ subfolders; reload a prior run."""

    reload_requested = QtCore.Signal(object)  # MystranRun

    def __init__(self, scratch_root: Path, parent=None):
        super().__init__(parent)
        self.scratch_root = Path(scratch_root)
        self.setWindowTitle("Analysis History (MYSTRAN)")
        self.setMinimumSize(820, 460)

        layout = QVBoxLayout(self)
        intro = QLabel(
            f"Runs in <code>{self.scratch_root}</code>. Double-click a "
            f"row to reload its results.")
        intro.setStyleSheet("color: #cdd6f4; font-size: 11px;")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        self._table = QTableWidget(0, 6)
        self._table.setHorizontalHeaderLabels(
            ["Started", "BDF", "SOL", "Exit", "Wall", "DOF"])
        self._table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectRows)
        self._table.setAlternatingRowColors(True)
        self._table.verticalHeader().setVisible(False)
        self._table.horizontalHeader().setStretchLastSection(False)
        self._table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.Stretch)
        self._table.doubleClicked.connect(self._on_double_clicked)
        layout.addWidget(self._table)

        button_row = QHBoxLayout()
        refresh_btn = QPushButton("Refresh")
        refresh_btn.clicked.connect(self._populate)
        button_row.addWidget(refresh_btn)
        open_folder_btn = QPushButton("Open Folder")
        open_folder_btn.clicked.connect(self._open_folder)
        button_row.addWidget(open_folder_btn)
        button_row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        button_row.addWidget(close_btn)
        layout.addLayout(button_row)

        self._populate()

    def _populate(self):
        self._table.setRowCount(0)
        runs_root = self.scratch_root / 'runs'
        if not runs_root.is_dir():
            return
        rows = []
        for run_dir in sorted(runs_root.iterdir(), reverse=True):
            meta = run_dir / 'run_meta.json'
            if not meta.is_file():
                continue
            try:
                d = json.loads(meta.read_text(encoding='utf-8'))
            except Exception:
                continue
            rows.append((run_dir, d))
        self._table.setRowCount(len(rows))
        for i, (run_dir, d) in enumerate(rows):
            started = d.get('started_at') or ''
            try:
                started = datetime.fromisoformat(started).strftime(
                    "%Y-%m-%d %H:%M:%S")
            except Exception:
                pass
            self._table.setItem(i, 0, QTableWidgetItem(started))
            bdf = d.get('bdf_path') or ''
            self._table.setItem(i, 1, QTableWidgetItem(Path(bdf).name))
            self._table.setItem(i, 2,
                                QTableWidgetItem(str(d.get('sol', ''))))
            rc = d.get('return_code')
            rc_item = QTableWidgetItem(str(rc) if rc is not None else '?')
            if rc == 0:
                rc_item.setForeground(QColor("#a6e3a1"))
            elif rc is not None:
                rc_item.setForeground(QColor("#f38ba8"))
            self._table.setItem(i, 3, rc_item)
            wall = d.get('wall_time') or 0
            self._table.setItem(i, 4, QTableWidgetItem(f"{wall:.1f} s"))
            self._table.setItem(i, 5,
                                QTableWidgetItem(str(d.get('n_dof', '?'))))
            self._table.setRowHeight(i, 22)
            # Stash the dir on row 0's UserRole for retrieval.
            self._table.item(i, 0).setData(Qt.UserRole, str(run_dir))

    def _on_double_clicked(self, idx):
        row = idx.row()
        item = self._table.item(row, 0)
        if not item:
            return
        run_dir = item.data(Qt.UserRole)
        if not run_dir:
            return
        try:
            d = json.loads((Path(run_dir) / 'run_meta.json').read_text(
                encoding='utf-8'))
        except Exception:
            return
        run = MystranRun(
            bdf_path=Path(d['bdf_path']) if d.get('bdf_path') else Path(
                run_dir) / 'input.bdf',
            f06_path=Path(d['f06_path']) if d.get('f06_path') else None,
            op2_path=Path(d['op2_path']) if d.get('op2_path') else None,
            log_path=Path(d['log_path']) if d.get('log_path') else None,
            return_code=d.get('return_code'),
            wall_time=d.get('wall_time'),
            sol=d.get('sol', 101),
            analysis_set_name=d.get('analysis_set_name', ''),
        )
        self.reload_requested.emit(run)
        self.accept()

    def _open_folder(self):
        import subprocess, sys
        path = str(self.scratch_root)
        try:
            if sys.platform == 'win32':
                os.startfile(path)
            elif sys.platform == 'darwin':
                subprocess.Popen(['open', path])
            else:
                subprocess.Popen(['xdg-open', path])
        except Exception:
            pass
