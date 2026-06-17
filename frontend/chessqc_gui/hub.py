"""CHESS-QC desktop landing page (hub).

Lists the available applications grouped by ACES functional area; clicking one
opens its calculator. The hub stays open so the user can switch between apps.
"""
from __future__ import annotations

from .qt import QtCore, QtWidgets, Qt
from .app_shell import CalculatorWindow
from . import settings
from .theme import load_qss

# canonical functional-area order
AREA_ORDER = [
    "Wave Prediction", "Wave Theory", "Wave Transformation", "Structural Design",
    "Wave Runup, Transmission, and Overtopping", "Littoral Processes",
    "Inlet Processes", "Harbor Design", "Storm Surge", "Miscellaneous Routines",
]


class HubWindow(QtWidgets.QMainWindow):
    def __init__(self, apps: dict, parent=None):
        super().__init__(parent)
        self.apps = apps                 # {aces_id: module}
        self._open_windows: list[CalculatorWindow] = []
        self.setWindowTitle("CHESS-QC (Quick Compute)")
        self._build()
        self.resize(620, 560)

    def _build(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        v = QtWidgets.QVBoxLayout(central)

        head = QtWidgets.QHBoxLayout()
        col = QtWidgets.QVBoxLayout()
        title = QtWidgets.QLabel("CHESS-QC")
        title.setStyleSheet("font-size:22px; font-weight:800;")
        sub = QtWidgets.QLabel("Coastal Hazards, Engineering, and Structures System (CHESS) — Quick Compute (QC)")
        sub.setObjectName("fieldLabel")
        col.addWidget(title)
        col.addWidget(sub)
        head.addLayout(col)
        head.addStretch(1)
        self.lbl_mode = QtWidgets.QLabel()
        self.lbl_mode.setObjectName("fieldLabel")
        self.sw_mode = QtWidgets.QCheckBox()        # styled as a switch via QSS indicator images
        self.sw_mode.setObjectName("ModeSwitch")
        self.sw_mode.setCursor(Qt.PointingHandCursor)
        self.sw_mode.setChecked(settings.get_theme() == "Dark")
        self.sw_mode.toggled.connect(self._on_mode_toggle)
        self._sync_mode_label()
        # flat header (no nested stretch layouts) -> the switch keeps its place on restyle
        head.addWidget(self.lbl_mode, 0, Qt.AlignVCenter)
        head.addWidget(self.sw_mode, 0, Qt.AlignVCenter)
        v.addLayout(head)
        v.addSpacing(6)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QtWidgets.QWidget()
        iv = QtWidgets.QVBoxLayout(inner)

        by_area: dict[str, list] = {}
        for mod in self.apps.values():
            by_area.setdefault(mod.APP_META.area, []).append(mod)
        ordered = [a for a in AREA_ORDER if a in by_area] + \
                  [a for a in by_area if a not in AREA_ORDER]

        if not self.apps:
            iv.addWidget(QtWidgets.QLabel("No applications found in applications/."))
        for area in ordered:
            box = QtWidgets.QGroupBox(area)
            grid = QtWidgets.QGridLayout(box)
            mods = sorted(by_area[area], key=lambda m: m.APP_META.aces_id)
            for i, mod in enumerate(mods):
                grid.addWidget(self._app_button(mod), i // 2, i % 2)
            iv.addWidget(box)
        iv.addStretch(1)
        scroll.setWidget(inner)
        v.addWidget(scroll, 1)

    def _app_button(self, mod) -> QtWidgets.QWidget:
        m = mod.APP_META
        btn = QtWidgets.QPushButton(f"{m.aces_id}    {m.name}")
        btn.setObjectName("AppCard")               # styled (incl. area-colored border) by the QSS
        btn.setProperty("area", m.area)
        btn.setMinimumHeight(42)
        btn.setCursor(Qt.PointingHandCursor)
        btn.setToolTip(m.cite)
        btn.clicked.connect(lambda _=False, mm=mod: self._open(mm))
        return btn

    def _open(self, mod):
        win = CalculatorWindow(mod, on_home=self._raise_hub,
                               apps=self.apps, on_switch=self._open_by_id)
        self._open_windows.append(win)
        win.show()
        win.raise_()
        win.activateWindow()

    def _open_by_id(self, aces_id):
        mod = self.apps.get(aces_id)
        if mod is not None:
            self._open(mod)

    def _sync_mode_label(self):
        self.lbl_mode.setText("Dark mode" if settings.get_theme() == "Dark" else "Light mode")

    def _on_mode_toggle(self, checked: bool):
        new = "Dark" if checked else "Light"
        settings.set_theme(new)                      # persists as the new default
        self._sync_mode_label()
        # Defer the global restyle out of the QCheckBox's own toggled signal.
        QtCore.QTimer.singleShot(0, self._apply_mode_style)

    def _apply_mode_style(self):
        app = QtWidgets.QApplication.instance()
        if app is not None:                          # restyle hub + any open calculators
            app.setStyleSheet(load_qss(settings.get_theme(), settings.get_vibe(), settings.get_badge()))
        # matplotlib canvases aren't QSS-styled -> redraw each open calculator's plot
        for win in list(self._open_windows):
            try:
                win.apply_theme()
            except RuntimeError:                     # window already closed/destroyed
                self._open_windows.remove(win)

    def _raise_hub(self):
        self.show()
        self.raise_()
        self.activateWindow()
