"""BDF import progress UI.

Earlier versions tried to run pyNastran's parser on a QThread so the UI
stayed responsive during big imports. In practice the interaction
between QProgressDialog's busy-indicator (min=max=0 spinner) and
pyNastran's pure-Python parsing path caused the worker thread to be
effectively starved of CPU on Windows: the parse that completes in 6 s
synchronously took 120 s+ via the worker, with the UI showing
"Parsing BDF..." indefinitely.

The pragmatic fix is to drop the worker thread and parse synchronously
on the main thread inside a busy-cursor + simple status dialog. The UI
is briefly frozen during the parse (about 6 s for a 60k-node BDF),
which is preferable to an indefinite hang. The previous QThread / Worker
classes are kept for back-compat but the public entry point now calls
the sync path.
"""

import os
from PySide6.QtCore import QObject, QThread, QTimer, Signal, Slot, Qt
from PySide6.QtGui import QCursor
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QProgressBar,
    QProgressDialog, QApplication,
)


class BdfReadWorker(QObject):
    """Run NastranModelGenerator's robust read in a background thread.

    Signals:
        finished(generator, lenient_result_or_none, detected_format)
            Emitted on successful parse. The MainWindow consumes these
            three values just like the old synchronous return tuple.
        failed(error_message)
            Emitted if every parsing strategy fails.
        progress(text)
            Optional status updates for the progress dialog (currently
            just stage labels: 'detecting format', 'parsing', etc.).
    """

    finished = Signal(object, object, str)
    failed = Signal(str)
    progress = Signal(str)

    def __init__(self, filepath, parent=None):
        super().__init__(parent)
        self._filepath = filepath
        self._cancelled = False

    def cancel(self):
        # pyNastran's read_bdf is a single blocking C-extension call which
        # we cannot interrupt mid-parse. The cancel flag means: "ignore the
        # result if the work finishes after the user already cancelled."
        self._cancelled = True

    @Slot()
    def run(self):
        try:
            from node_runner.model import (
                NastranModelGenerator, detect_bdf_field_format,
            )
            self.progress.emit("Detecting field format...")
            detected_format = detect_bdf_field_format(self._filepath)
            if self._cancelled:
                return

            # NOTE: the parallel parser path (model._read_bdf_parallel) is
            # currently SLOWER than the single-threaded one for typical
            # files - pyNastran's card parsing is mostly pure Python and
            # the GIL serializes the workers, so 4 threads end up
            # competing for the same lock and fighting cache locality.
            # The path is kept in model.py for future improvement, but
            # the worker does not invoke it. Re-enabling it requires
            # replacing pyNastran's Python parsers with something that
            # actually releases the GIL during the bulk-data scan.
            self.progress.emit("Parsing BDF...")
            generator = NastranModelGenerator()
            model, lenient_result = NastranModelGenerator._read_bdf_robust(
                self._filepath,
            )
            generator.model = model
            if self._cancelled:
                return

            self.finished.emit(generator, lenient_result, detected_format)
        except Exception as exc:
            if not self._cancelled:
                self.failed.emit(str(exc))


class BdfImportProgressDialog(QProgressDialog):
    """Modal busy dialog with a Cancel button while a BDF parses.

    Indeterminate progress (range 0,0) since pyNastran does not expose
    parse progress. Cancel sets the worker's flag; the worker will discard
    its result when it completes.
    """

    def __init__(self, parent=None, label="Reading BDF..."):
        super().__init__(label, "Cancel", 0, 0, parent)
        self.setWindowTitle("Import")
        self.setWindowModality(Qt.WindowModal)
        self.setMinimumDuration(0)  # show immediately
        self.setAutoReset(False)
        self.setAutoClose(False)


class _SyncImportDialog(QDialog):
    """Status dialog shown while a BDF imports.

    v3.2.0 layout, top to bottom:
      1. **Header**  - "Importing: <filename>" - persistent, set once.
      2. **Stage**   - "Stage X/Y: <short label>" - what we're doing now.
      3. **Bar**     - determinate, 0..100%, filled by the parser.
      4. **Detail**  - the helpful line: source include file, current
                       card type, line# (if known), most recent error,
                       cards/sec rate, ETA.
      5. **Cancel**  - sets a flag the parser checks between chunks;
                       parser exits cleanly, viewer keeps previous state.

    A QTimer watches the time since the last progress event and flips
    the detail line to "Waiting for parser... N s since last update"
    if it goes more than 30 s without activity. This keeps the dialog
    informative even if a single pyNastran call genuinely takes minutes.
    """

    def __init__(self, filepath: str = '', parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import")
        self.setMinimumWidth(520)
        self.setModal(False)
        # Track cancel state; the parser polls via was_cancelled().
        self._cancel_requested = False
        self._closed_by_user = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(6)

        # 1. Persistent header - the file the user actually selected.
        self._header = QLabel("")
        self._header.setStyleSheet(
            "font-weight: 600; font-size: 12px; color: #cdd6f4;")
        self._header.setWordWrap(True)
        layout.addWidget(self._header)

        # 2. Stage label.
        self._stage = QLabel("Preparing...")
        self._stage.setWordWrap(True)
        layout.addWidget(self._stage)

        # 3. Determinate bar.
        self._bar = QProgressBar(self)
        self._bar.setRange(0, 100)
        self._bar.setValue(0)
        self._bar.setTextVisible(True)
        layout.addWidget(self._bar)

        # 4. Detail line - the actually-helpful info.
        self._detail = QLabel("")
        self._detail.setWordWrap(True)
        self._detail.setStyleSheet("color: #b4befe; font-size: 11px;")
        # Reserve two lines of vertical space so the dialog doesn't
        # jump around when text wraps.
        self._detail.setMinimumHeight(32)
        layout.addWidget(self._detail)

        hint = QLabel(
            "The window may freeze briefly during pyNastran's bulk-data "
            "scan. Multi-file decks update file-by-file."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #6c7086; font-size: 10px;")
        layout.addWidget(hint)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self._cancel_btn = QPushButton("Cancel")
        self._cancel_btn.clicked.connect(self._on_cancel)
        button_row.addWidget(self._cancel_btn)
        layout.addLayout(button_row)

        if filepath:
            self.set_filepath(filepath)

        # Watchdog: if we don't see a progress event for 30 s, the
        # detail line gets an emergency stamp so the user knows we're
        # waiting on pyNastran (and that Cancel is available).
        import time as _time
        self._last_progress_time = _time.time()
        self._watchdog = QTimer(self)
        self._watchdog.setInterval(2000)  # check every 2 s
        self._watchdog.timeout.connect(self._tick_watchdog)
        self._watchdog.start()

    # ----- public API used by run_bdf_import_threaded -----

    def set_filepath(self, filepath: str) -> None:
        import os as _os
        self._header.setText(
            f"Importing: <b>{_os.path.basename(filepath)}</b>")

    def apply_progress(self, record) -> None:
        """Render an ImportProgress record. ``record`` is a dataclass
        instance from node_runner.model.ImportProgress."""
        import time as _time
        self._last_progress_time = _time.time()

        if record.label:
            self._stage.setText(record.label)
        if record.fraction is not None:
            try:
                v = max(0, min(100, int(float(record.fraction) * 100)))
                self._bar.setValue(v)
            except (TypeError, ValueError):
                pass

        # Build the detail line. Compose only the fields that are set
        # so we don't leave " | " separators dangling.
        parts = []
        if record.error:
            # Errors get center stage with file/card prefix when available.
            prefix = []
            if record.source_file:
                prefix.append(record.source_file)
            if record.line_number:
                prefix.append(f"line {record.line_number}")
            if record.card_type:
                prefix.append(record.card_type)
            prefix_str = " : ".join(prefix)
            if prefix_str:
                parts.append(
                    f"<span style='color:#f38ba8'>{prefix_str}  -  "
                    f"{record.error}</span>")
            else:
                parts.append(
                    f"<span style='color:#f38ba8'>{record.error}</span>")
        else:
            if record.source_file:
                parts.append(f"from <b>{record.source_file}</b>")
            if record.card_type:
                parts.append(record.card_type)
            if record.line_number:
                parts.append(f"line {record.line_number}")
            if record.counter_total:
                parts.append(
                    f"{record.counter_done:,} / {record.counter_total:,} cards")
            if record.rate:
                parts.append(f"{int(record.rate):,} cards/s")
            if record.eta_seconds:
                eta = record.eta_seconds
                if eta > 90:
                    parts.append(f"ETA {eta / 60:.1f} min")
                else:
                    parts.append(f"ETA {int(eta)} s")
        self._detail.setText("  |  ".join(parts) if parts else "")

    def was_cancelled(self) -> bool:
        return self._cancel_requested

    # ----- backwards-compat shim (legacy callers) -----

    def set_message(self, text: str):
        self._stage.setText(text)

    def set_detail(self, text: str):
        self._detail.setText(text)

    def set_fraction(self, frac):
        if frac is None:
            return
        try:
            v = max(0, min(100, int(float(frac) * 100)))
            self._bar.setValue(v)
        except (TypeError, ValueError):
            pass

    # ----- internals -----

    def _on_cancel(self):
        self._cancel_requested = True
        self._cancel_btn.setEnabled(False)
        self._cancel_btn.setText("Cancelling...")
        self._stage.setText("Cancellation requested - waiting for safe stop point...")

    def _tick_watchdog(self):
        import time as _time
        idle = _time.time() - self._last_progress_time
        if idle > 30 and not self._cancel_requested:
            # Don't overwrite the detail line; append a hint.
            current = self._detail.text()
            stamp = (f"<span style='color:#fab387'>"
                     f"(no parser activity for {int(idle)}s - "
                     f"pyNastran may be in a single big call)"
                     f"</span>")
            if "no parser activity" not in current:
                self._detail.setText(
                    current + ("<br/>" if current else "") + stamp)
            else:
                # Update the existing stamp's elapsed counter
                import re
                self._detail.setText(re.sub(
                    r'\(no parser activity for \d+s',
                    f'(no parser activity for {int(idle)}s',
                    current,
                ))


def run_bdf_import_threaded(parent, filepath, on_success, on_failure):
    """Run a BDF import. Despite the legacy name this is now SYNCHRONOUS.

    Walks: detect_bdf_field_format -> _read_bdf_robust (or streaming
    for INCLUDE-heavy decks) -> on_success. Shows a small status dialog
    with a determinate progress bar; the main thread blocks for the
    parse (in the same way the v3.0.0 sync-import fix works, without
    the spinner that starved the worker thread of CPU).

    Return shape kept the same `(thread, worker, dialog)` tuple so the
    callers in MainWindow don't need to be rewritten.
    """
    from node_runner.model import (
        NastranModelGenerator, detect_bdf_field_format, ImportProgress,
    )

    dialog = _SyncImportDialog(filepath, parent)
    dialog.show()
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    QApplication.processEvents()  # one paint pass so the dialog appears

    class _Cancelled(Exception):
        pass

    def _progress(record):
        """Receives ImportProgress dataclass instances from the model
        layer. Pushes them into the dialog and checks for cancel."""
        if dialog.was_cancelled():
            raise _Cancelled()
        dialog.apply_progress(record)
        QApplication.processEvents()

    try:
        dialog.apply_progress(ImportProgress(
            stage='inline',
            label='Stage 1/5: Detecting field format',
            source_file=os.path.basename(filepath) if filepath else '',
            fraction=0.0,
        ))
        QApplication.processEvents()
        detected_format = detect_bdf_field_format(filepath)

        # If the deck has INCLUDE statements, take the streaming path.
        # That path enumerates files first so the progress bar can be
        # determinate (file N / M), pre-inlines into a flat buffer, and
        # then parses once with no per-include round trips inside
        # pyNastran's reader.
        has_includes = NastranModelGenerator._file_has_includes(filepath)
        generator = NastranModelGenerator()

        if has_includes:
            dialog.apply_progress(ImportProgress(
                stage='inline',
                label='Stage 2/5: Scanning INCLUDE chain',
                source_file='Counting files + sizes...',
                fraction=0.0,
            ))
            QApplication.processEvents()
            model, lenient_result = (
                NastranModelGenerator._read_bdf_streaming(
                    filepath, progress=_progress))
        else:
            dialog.apply_progress(ImportProgress(
                stage='parse',
                label='Stage 2/5: Parsing BDF',
                source_file='No INCLUDE statements - reading directly',
                fraction=0.1,
            ))
            QApplication.processEvents()
            model, lenient_result = (
                NastranModelGenerator._read_bdf_robust(filepath))
            dialog.apply_progress(ImportProgress(
                stage='done', label='Parse complete', fraction=1.0))

        generator.model = model

        # Resolve cp_ref / cd_ref on every node so get_position() works.
        # Without this, scene-building crashes with
        #   'NoneType' object has no attribute 'transform_node_to_global'
        # because we read with xref=False for speed and pyNastran's
        # node.get_position() requires cp_ref to be set.
        if dialog.was_cancelled():
            raise _Cancelled()
        dialog.apply_progress(ImportProgress(
            stage='render',
            label='Stage 5/5: Resolving coordinate references',
            source_file='Wiring up grid coord systems for the viewer...',
            fraction=1.0,
        ))
        QApplication.processEvents()
        try:
            NastranModelGenerator._finalize_for_viewer(model)
        except Exception:
            # Best effort; visualization may fail later but we don't
            # want to lose the whole import over an unexpected hiccup.
            pass
    except _Cancelled:
        QApplication.restoreOverrideCursor()
        dialog.close()
        on_failure("Import cancelled by user.")
        return None, None, None
    except Exception as exc:
        QApplication.restoreOverrideCursor()
        dialog.close()
        on_failure(str(exc))
        return None, None, None

    QApplication.restoreOverrideCursor()
    dialog.close()
    QApplication.processEvents()  # let the close paint before on_success runs
    on_success(generator, lenient_result, detected_format)
    return None, None, None
