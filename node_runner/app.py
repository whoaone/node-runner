import sys
from PySide6.QtWidgets import QApplication
from node_runner.theme import dark_palette, DARK_STYLESHEET


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setPalette(dark_palette)
    app.setStyleSheet(DARK_STYLESHEET)

    from node_runner.mainwindow import MainWindow
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
