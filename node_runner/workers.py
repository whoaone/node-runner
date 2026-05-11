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
    """Plain non-modal status dialog shown while we parse synchronously.

    Deliberately NOT a QProgressDialog with busy indicator - that
    animation poisoned the GIL and starved the parse worker on Windows
    (see v3.0.0 history). The progress bar here is DETERMINATE - a
    static rectangle that fills as the parser walks INCLUDEs, not a
    marching-ants spinner. That avoids the GIL-starvation issue while
    still giving the user feedback that work is happening.

    The Hide button dismisses the dialog without interrupting the
    parse (pyNastran's bulk-data read is a single blocking call we
    can't preempt).
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import")
        self.setMinimumWidth(440)
        self.setModal(False)

        layout = QVBoxLayout(self)
        self._label = QLabel("Reading BDF...")
        layout.addWidget(self._label)

        # Determinate progress bar. The streaming reader updates the
        # value as it walks INCLUDE files; legacy callers that don't
        # emit progress just leave it at zero.
        self._bar = QProgressBar(self)
        self._bar.setRange(0, 100)
        self._bar.setValue(0)
        self._bar.setTextVisible(True)
        layout.addWidget(self._bar)

        self._detail = QLabel("")
        self._detail.setWordWrap(True)
        self._detail.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(self._detail)

        hint = QLabel(
            "The window may freeze briefly while pyNastran parses the "
            "bulk-data section. Multi-file decks update file-by-file."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(hint)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self._close_btn = QPushButton("Hide")
        self._close_btn.clicked.connect(self.close)
        button_row.addWidget(self._close_btn)
        layout.addLayout(button_row)

    def set_message(self, text: str):
        self._label.setText(text)

    def set_detail(self, text: str):
        self._detail.setText(text)

    def set_fraction(self, frac):
        """Update the bar; ``frac`` is 0.0..1.0 or None for no change."""
        if frac is None:
            return
        try:
            v = max(0, min(100, int(float(frac) * 100)))
            self._bar.setValue(v)
        except (TypeError, ValueError):
            pass


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
        NastranModelGenerator, detect_bdf_field_format,
    )

    dialog = _SyncImportDialog(parent)
    dialog.show()
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    QApplication.processEvents()  # one paint pass so the dialog appears

    try:
        dialog.set_message("Detecting field format...")
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
            dialog.set_message("Stage 1/4: Scanning INCLUDE chain")
            dialog.set_detail("Counting files + sizes...")
            QApplication.processEvents()

            def _progress(stage, message, frac):
                # The model layer emits messages of the form
                #   "Stage X/Y: <stage label> - <details>"
                # We split on the FIRST ' - ' so the top of the dialog
                # gets the high-level stage label and the line under
                # the progress bar gets the concrete details (file,
                # card, error, etc.). If there's no separator, the
                # whole message is treated as stage and details are
                # blanked.
                if ' - ' in message:
                    stage_label, _, detail_msg = message.partition(' - ')
                    stage_label = stage_label.strip()
                    detail_msg = detail_msg.strip()
                else:
                    stage_label = message.strip()
                    detail_msg = ''
                dialog.set_message(stage_label)
                dialog.set_detail(detail_msg)
                dialog.set_fraction(frac)
                QApplication.processEvents()

            model, lenient_result = (
                NastranModelGenerator._read_bdf_streaming(
                    filepath, progress=_progress))
        else:
            dialog.set_message("Parsing BDF (this may take a moment)...")
            dialog.set_fraction(0.1)
            QApplication.processEvents()
            model, lenient_result = (
                NastranModelGenerator._read_bdf_robust(filepath))
            dialog.set_fraction(1.0)

        generator.model = model

        # Resolve cp_ref / cd_ref on every node so get_position() works.
        # Without this, scene-building crashes with
        #   'NoneType' object has no attribute 'transform_node_to_global'
        # because we read with xref=False for speed and pyNastran's
        # node.get_position() requires cp_ref to be set.
        dialog.set_message("Stage 5/5: Resolving coordinate references")
        dialog.set_detail("Wiring up grid coord systems for the viewer...")
        QApplication.processEvents()
        try:
            NastranModelGenerator._finalize_for_viewer(model)
        except Exception:
            # Best effort; visualization may fail later but we don't
            # want to lose the whole import over an unexpected hiccup.
            pass
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
