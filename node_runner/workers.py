"""Background workers and progress UI for long-running engine operations.

Phase 2 currently provides BDF import off the UI thread. Future phases
will extend this with parallel parsing.
"""

from PySide6.QtCore import QObject, QThread, Signal, Slot, Qt
from PySide6.QtWidgets import QProgressDialog, QApplication


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
            from PySide6.QtCore import QSettings
            from node_runner.model import (
                NastranModelGenerator, detect_bdf_field_format,
            )
            self.progress.emit("Detecting field format...")
            detected_format = detect_bdf_field_format(self._filepath)
            if self._cancelled:
                return

            settings = QSettings("NodeRunner", "NodeRunner")
            use_parallel = bool(settings.value(
                "import/parallel_parsing", False, type=bool,
            ))

            generator = NastranModelGenerator()
            if use_parallel:
                self.progress.emit("Parsing BDF (parallel, experimental)...")
                model, lenient_result = NastranModelGenerator._read_bdf_parallel(
                    self._filepath,
                )
            else:
                self.progress.emit("Parsing BDF...")
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


def run_bdf_import_threaded(parent, filepath, on_success, on_failure):
    """Helper: spawn a BdfReadWorker, show progress, wire callbacks.

    on_success(generator, lenient_result, detected_format) - called on
        the UI thread when parsing completes successfully and the user
        has not cancelled.
    on_failure(error_message) - called on the UI thread when parsing
        fails for any reason other than cancellation.

    Returns a tuple (thread, worker, dialog) so the caller can keep them
    alive (Qt requires the QThread instance to outlive its work).
    """
    thread = QThread(parent)
    worker = BdfReadWorker(filepath)
    worker.moveToThread(thread)

    dialog = BdfImportProgressDialog(parent)

    # Wire signals.
    thread.started.connect(worker.run)
    worker.progress.connect(dialog.setLabelText)

    cancelled = {"flag": False}

    def _on_finished(generator, lenient_result, detected_format):
        dialog.close()
        thread.quit()
        thread.wait()
        if not cancelled["flag"]:
            on_success(generator, lenient_result, detected_format)

    def _on_failed(msg):
        dialog.close()
        thread.quit()
        thread.wait()
        if not cancelled["flag"]:
            on_failure(msg)

    def _on_canceled():
        cancelled["flag"] = True
        worker.cancel()
        dialog.close()
        # Don't call thread.quit() here - the underlying parse cannot be
        # interrupted; let it finish in the background. The worker drops
        # the result because of the cancel flag. The QThread cleans up
        # when its run() returns.

    worker.finished.connect(_on_finished)
    worker.failed.connect(_on_failed)
    dialog.canceled.connect(_on_canceled)

    thread.start()
    # Process events so the dialog actually paints.
    QApplication.processEvents()
    return thread, worker, dialog
