"""Background workers and progress UI for long-running engine operations.

Phase 2 currently provides BDF import off the UI thread. Future phases
will extend this with parallel parsing.
"""

from PySide6.QtCore import QObject, QThread, QTimer, Signal, Slot, Qt
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

    def _hide_dialog():
        # CRITICAL: QProgressDialog.close() emits the `canceled` signal
        # (because its cancel button is the default close action). If we
        # leave _on_canceled connected, programmatic close triggers the
        # cancel cascade and we drop the parsed result. Disconnect first.
        try:
            dialog.canceled.disconnect(_on_canceled)
        except (TypeError, RuntimeError):
            # Already disconnected (e.g., second call). Safe to ignore.
            pass
        try:
            dialog.reset()
        except Exception:
            pass
        dialog.hide()
        dialog.close()

    def _cleanup_thread():
        # quit() asks the worker thread's event loop to exit. We do NOT
        # call thread.wait() here - that would block the UI thread. The
        # worker's run() has already returned (it just emitted finished/
        # failed), so the event loop will exit on its own.
        thread.quit()

    def _on_finished(generator, lenient_result, detected_format):
        was_cancelled = cancelled["flag"]
        _hide_dialog()
        _cleanup_thread()
        if was_cancelled:
            return
        on_success(generator, lenient_result, detected_format)

    def _on_failed(msg):
        was_cancelled = cancelled["flag"]
        _hide_dialog()
        _cleanup_thread()
        if was_cancelled:
            return
        on_failure(msg)

    def _on_canceled():
        cancelled["flag"] = True
        worker.cancel()
        _hide_dialog()
        # Underlying parse cannot be interrupted; let the worker finish in
        # the background. It drops its result because cancelled is set.

    worker.finished.connect(_on_finished)
    worker.failed.connect(_on_failed)
    dialog.canceled.connect(_on_canceled)

    # Keep strong references to the Python closures alive for as long as
    # the worker / dialog exist. PySide6 can drop weakly-held closure
    # connections after the enclosing function returns, which would silently
    # skip on_success / on_failure when the signal fires.
    worker._kept_alive = (_on_finished, _on_failed, _hide_dialog, _cleanup_thread)
    dialog._kept_alive = (_on_canceled, _hide_dialog)

    thread.start()
    # Don't call QApplication.processEvents() here. It used to be a "make
    # the dialog paint" hint, but it can re-enter _on_finished if the worker
    # finishes parsing very quickly, blocking indefinitely on Windows. The
    # dialog paints fine through Qt's normal event loop once we return.
    return thread, worker, dialog
