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
        # Track the last file the parser successfully read. Persists
        # across stage transitions so the dialog never "forgets" what
        # was being processed last (v3.2.1 fix).
        self._last_file_seen = ''

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(6)

        # 1. Persistent header - the file the user actually selected.
        self._header = QLabel("")
        self._header.setStyleSheet(
            "font-weight: 600; font-size: 12px; color: #cdd6f4;")
        self._header.setWordWrap(True)
        layout.addWidget(self._header)

        # 1b. "Last file read" - persistent across stage transitions
        # so the user can always see which include we got furthest with.
        self._last_file_label = QLabel("")
        self._last_file_label.setStyleSheet(
            "color: #94e2d5; font-size: 11px;")
        self._last_file_label.setWordWrap(True)
        layout.addWidget(self._last_file_label)

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

        # Update the persistent "Last file read" line whenever a record
        # carries a non-empty source_file. We DO NOT clear it on records
        # that lack source_file - this is the v3.2.1 behaviour. v3.2.2
        # additionally guarantees that source_file is STRICTLY a
        # filename (not a sentence) because the model layer's _emit no
        # longer stuffs message tails into source_file.
        if record.source_file:
            self._last_file_seen = record.source_file
            self._last_file_label.setText(
                f"Last file read: <b>{record.source_file}</b>")

        # v3.2.2: structured ``detail`` field on the record carries any
        # human-readable detail the model layer wants displayed under
        # the bar (e.g. "Skipped 4,638,690 analysis-only cards (TEMP:
        # 4,638,690). Will be written back on export."). If the record
        # also has card_type / error / counter / rate, those get
        # appended in the format below.
        parts = []
        if record.detail:
            parts.append(record.detail)
        if record.error:
            prefix = []
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


class _BdfImportThread(QThread):
    """Background QThread that runs the BDF import off the UI thread.

    v3.2.1: the v3.1.x sync import meant pyNastran held the GIL for
    the entire parse and Windows declared the app "Not Responding" on
    decks where parsing took more than ~5 s of unsponsored time. v3.0
    originally went sync because the v2.x progress dialog had a
    marching-ants spinner that hammered the paint pipeline and starved
    the worker. Since v3.0.2 the dialog has been determinate (static
    bar + label updates only on signal), so the original concern no
    longer applies and threading is safe.

    Emits:
        progress(ImportProgress) - forwarded from the model layer
        done(generator, lenient_result, detected_format) - happy path
        failed(error_message) - any exception (incl. user cancel)
    """

    progress = Signal(object)
    done = Signal(object, object, str)
    failed = Signal(str)

    def __init__(self, filepath: str, parent=None):
        super().__init__(parent)
        self._filepath = filepath
        self._cancel = False

    def cancel(self):
        """Set the cancel flag. The next progress emit raises and the
        run loop returns. Called from the main thread when the user
        hits the dialog's Cancel button."""
        self._cancel = True

    def run(self):
        from node_runner.model import (
            NastranModelGenerator, detect_bdf_field_format, ImportProgress,
        )

        class _Cancelled(Exception):
            pass

        def _bridge(record):
            if self._cancel:
                raise _Cancelled()
            # Cross-thread signal: Qt marshals this onto the main
            # thread's event loop automatically.
            self.progress.emit(record)

        try:
            self.progress.emit(ImportProgress(
                stage='inline',
                label='Stage 1/5: Detecting field format',
                source_file=os.path.basename(self._filepath)
                if self._filepath else '',
                fraction=0.0,
            ))
            detected_format = detect_bdf_field_format(self._filepath)
            has_includes = NastranModelGenerator._file_has_includes(
                self._filepath)

            if has_includes:
                self.progress.emit(ImportProgress(
                    stage='inline',
                    label='Stage 2/5: Scanning INCLUDE chain',
                    source_file='Counting files + sizes...',
                    fraction=0.0,
                ))
                model, lenient_result = (
                    NastranModelGenerator._read_bdf_streaming(
                        self._filepath, progress=_bridge))
            else:
                self.progress.emit(ImportProgress(
                    stage='parse',
                    label='Stage 2/5: Parsing BDF',
                    source_file='No INCLUDE statements - reading directly',
                    fraction=0.1,
                ))
                model, lenient_result = (
                    NastranModelGenerator._read_bdf_robust(self._filepath))
                self.progress.emit(ImportProgress(
                    stage='done', label='Parse complete', fraction=1.0))

            generator = NastranModelGenerator()
            generator.model = model

            if self._cancel:
                raise _Cancelled()
            self.progress.emit(ImportProgress(
                stage='render',
                label='Stage 5/5: Resolving coordinate references',
                source_file='Wiring up grid coord systems for the viewer...',
                fraction=1.0,
            ))
            try:
                NastranModelGenerator._finalize_for_viewer(model)
            except Exception:
                pass

            self.done.emit(generator, lenient_result, detected_format)
        except _Cancelled:
            self.failed.emit("Import cancelled by user.")
        except Exception as exc:
            self.failed.emit(str(exc))


def run_bdf_import_threaded(parent, filepath, on_success, on_failure):
    """Run a BDF import in a background QThread.

    v3.2.1: actually threaded again (the v3.1.x sync version held the
    GIL through pyNastran's parse and Windows would declare the app
    "Not Responding" on big decks). The dialog uses a determinate
    progress bar that updates only on signal emission - this is
    GIL-friendly and gives the main thread enough time to repaint.

    Returns ``(thread, None, dialog)`` so old callers that hold a
    reference to ``thread`` for lifetime management still work.
    """
    dialog = _SyncImportDialog(filepath, parent)
    dialog.show()
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    QApplication.processEvents()  # paint the dialog before work begins

    thread = _BdfImportThread(filepath, parent=parent)

    def _on_progress(record):
        try:
            dialog.apply_progress(record)
        except Exception:
            pass

    def _on_done(generator, lenient_result, detected_format):
        QApplication.restoreOverrideCursor()
        dialog.close()
        try:
            on_success(generator, lenient_result, detected_format)
        finally:
            thread.deleteLater()

    def _on_failed(msg):
        QApplication.restoreOverrideCursor()
        dialog.close()
        try:
            on_failure(msg)
        finally:
            thread.deleteLater()

    def _on_cancel_clicked():
        # The dialog's existing _on_cancel sets _cancel_requested AND
        # disables the button. Wire it through to the worker.
        thread.cancel()

    thread.progress.connect(_on_progress)
    thread.done.connect(_on_done)
    thread.failed.connect(_on_failed)
    dialog._cancel_btn.clicked.connect(_on_cancel_clicked)

    thread.start()
    return thread, None, dialog
