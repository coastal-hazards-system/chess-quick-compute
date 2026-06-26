"""Generic CHESS-QC calculator window.

`CalculatorWindow(app)` renders any chessqc_* application from its contract
(APP_META, INPUTS, OUTPUTS, compute). Layout: inputs left, Compute center,
results right (value rows + Plot/Table tabs); SI/US unit toggle + example/reset.
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

from .qt import QtCore, QtGui, QtWidgets, Qt
from common import units
from . import settings
from .theme import get_plot_palette, vibe_label_opts

# area-10 workflow hand-off: the upstream window stows its full-resolution series
# here; the next window injects it into its first CSV field on open, then clears it.
_HANDOFF: dict = {}

# bundled station CSVs live at the repo root: data/<data_dir>/<id>.csv
_DATA_ROOT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data")

import io

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

# LaTeX equation -> QPixmap via matplotlib mathtext (the same `tex` strings the web
# front-end feeds to KaTeX). Cached by (tex, color, dpr); transparent background so it
# sits on any themed panel. Returns None if mathtext can't parse the string.
_MATH_CACHE: dict = {}


def _latex_pixmap(tex: str, color: str, dpr: float = 2.0, fontsize: int = 13):
    key = (tex, color, round(dpr, 2), fontsize)
    if key in _MATH_CACHE:
        return _MATH_CACHE[key]
    fig = Figure(dpi=100.0 * dpr)
    fig.patch.set_alpha(0.0)
    FigureCanvasAgg(fig)
    fig.text(0.0, 0.0, f"${tex}$", color=color, fontsize=fontsize)
    buf = io.BytesIO()
    try:
        fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0.04, transparent=True)
    except Exception:
        _MATH_CACHE[key] = None
        return None
    pm = QtGui.QPixmap()
    pm.loadFromData(buf.getvalue(), "PNG")
    pm.setDevicePixelRatio(dpr)
    _MATH_CACHE[key] = pm
    return pm


_GREEK = {"alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota",
          "kappa", "lambda", "mu", "nu", "xi", "pi", "rho", "sigma", "tau", "upsilon",
          "phi", "chi", "psi", "omega"}


def _sym_tok(t: str) -> str:
    base, sub = t, ""
    u = t.find("_")
    if u >= 0:
        base, sub = t[:u], t[u + 1:]
    if base.lower() in _GREEK:
        base = "\\" + base
    if sub == "":
        return base
    return base + "_" + ("{" + sub + "}" if len(sub) > 1 else sub)


def _sym_to_tex(s: str) -> str:
    """ABOUT symbol shorthand ("U_obs", "u_*", "phi", "Delta T", "bar t") -> LaTeX. Mirrors
    the web driver's symToTex so the symbol column typesets identically in both front-ends."""
    s = str(s).strip()
    if "," in s:
        return ",\\; ".join(_sym_to_tex(p) for p in s.split(","))
    toks = s.split()
    out, i = [], 0
    while i < len(toks):
        if toks[i] == "bar" and i + 1 < len(toks):
            out.append("\\bar{" + _sym_tok(toks[i + 1]) + "}"); i += 2
        else:
            out.append(_sym_tok(toks[i])); i += 1
    return " ".join(out)

_BIG = 1.0e9  # finite bound for "unrestricted" spin boxes
# fidelity classification -> compact badge letter and QSS [badge="..."] style
_CLASS_LETTER = {"exact": "I", "standard": "II", "provisional": "III"}
_CLASS_BADGE = {"exact": "success", "standard": "warning", "provisional": "danger"}


def _clamp_bound(x: float) -> float:
    if x == math.inf:
        return _BIG
    if x == -math.inf:
        return -_BIG
    return float(x)


class CalculatorWindow(QtWidgets.QMainWindow):
    def __init__(self, app, on_home=None, apps=None, on_switch=None, parent=None):
        super().__init__(parent)
        self.app = app
        self._on_home = on_home
        self._apps = apps or {}              # {aces_id: module} for the in-window app switcher
        self._on_switch = on_switch          # callback(aces_id) -> open that app
        self.meta = app.APP_META
        self.inputs = list(app.INPUTS)
        self.outputs = list(app.OUTPUTS)
        # launcher-set units (DEFAULT_UNITS) take precedence; else the app's own default
        self.system = settings.get_units(getattr(self.meta, "default_system", "SI"))
        self.decimals = settings.get_decimals()
        self._tc_speed = settings.get_tcwind_speed()   # SI unit for wind/TC speeds (km/h|m/s)
        self._vibe_opts = vibe_label_opts(settings.get_vibe())
        self._widgets: dict[str, QtWidgets.QWidget] = {}
        self._value_labels: dict[str, QtWidgets.QLabel] = {}
        self._last_result = None

        self.setWindowTitle(f"CHESS-QC · {self.meta.aces_id} {self.meta.name}")
        self._build_ui()
        self.resize(960, 700)
        self._apply_handoff()   # inject a series carried from the previous workflow step
        self._on_compute()  # populate inputs/outputs from defaults (no CSV fetch)
        self._focus_first_input()

    def _apply_handoff(self):
        """If a previous workflow step stowed a series, inject it into this app's
        first CSV field (e.g. detrended WL -> 10-2/10-3; NTR -> 10-3)."""
        h = _HANDOFF.pop("series", None) if _HANDOFF else None
        if not h:
            return
        for f in self.inputs:
            if f.kind == "csv":
                w = self._widgets.get(f.key)
                if w is not None:
                    w._csv_text = h.get("csv", "")
                    w._combo.blockSignals(True); w._combo.setCurrentIndex(0); w._combo.blockSignals(False)
                    w._last_idx = 0
                    w._status.setText(f"← {h.get('label', 'from previous step')}")
                break

    def _focus_first_input(self):
        """Focus the first input widget on open (parity with the web front-end)."""
        for f in self.inputs:
            w = self._widgets.get(f.key)
            if w is not None:
                w.setFocus()
                break

    def keyPressEvent(self, event):
        """Enter/Return runs Compute (parity with the web front-end). A focused table
        cell editor consumes the key first, so this does not disrupt table editing."""
        if event.key() in (Qt.Key_Return, Qt.Key_Enter):
            self._on_compute()
            return
        super().keyPressEvent(event)

    def _on_switch_app(self, idx):
        """Switch to another application from the in-window combo (parity with the web
        app switcher): open the selected app, then close this window."""
        aces_id = self.app_combo.itemData(idx)
        if aces_id and aces_id != self.meta.aces_id and self._on_switch is not None:
            self._on_switch(aces_id)
            QtCore.QTimer.singleShot(0, self.close)

    # ---- field unit helpers ----
    def _si_unit(self, u: str) -> str:
        """Wind/TC speeds (SI unit km/h) follow the launcher choice (km/h|m/s);
        water/wave m/s and all other units are returned unchanged."""
        return self._tc_speed if u == "km/h" else u

    def _unit(self, f) -> str:
        return self._si_unit(f.unit_si) if self.system == "SI" else f.unit_us

    def _out_unit(self, o) -> str:
        return self._si_unit(o.unit_si) if self.system == "SI" else o.unit_us

    def _fmt(self, x: float) -> str:
        """Fixed-decimal display, with a significant-figure fallback so a small
        nonzero value (e.g. 0.003) is never shown as 0 (kills the '-0.00' artifact
        and works for any field without per-field decimal settings)."""
        x = float(x)
        d = self.decimals
        if x == 0 or not math.isfinite(x):
            return f"{0.0:.{d}f}" if x == 0 else str(x)
        s = f"{x:.{d}f}"
        if float(s) != 0:
            return s
        return f"{x:.{max(d, 2)}g}"

    def _field_label(self, text: str) -> QtWidgets.QLabel:
        """Field-name label styled per the active vibe (Forge: uppercase mono, spaced)."""
        if self._vibe_opts["label_upper"]:
            text = text.upper()
        lbl = QtWidgets.QLabel(text)
        lbl.setObjectName("fieldLabel")
        if self._vibe_opts["mono_labels"]:
            f = lbl.font()
            f.setFamilies(["JetBrains Mono", "Cascadia Mono", "Consolas", "monospace"])
            f.setLetterSpacing(QtGui.QFont.SpacingType.PercentageSpacing, 106)
            lbl.setFont(f)
        return lbl

    # ---- UI construction ----
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        outer = QtWidgets.QVBoxLayout(central)

        # top bar: [home] units + example + reset
        top = QtWidgets.QHBoxLayout()
        if self._on_home is not None:
            btn_home = QtWidgets.QPushButton("← Applications")
            btn_home.clicked.connect(self._on_home)
            top.addWidget(btn_home)
        if self._apps and self._on_switch is not None:
            self.app_combo = QtWidgets.QComboBox()   # switch apps without returning to the hub
            self.app_combo.setMaxVisibleItems(20)
            for aid in sorted(self._apps):
                m = self._apps[aid].APP_META
                self.app_combo.addItem(f"{m.aces_id}  {m.name}", aid)
            ix = self.app_combo.findData(self.meta.aces_id)
            if ix >= 0:
                self.app_combo.setCurrentIndex(ix)
            self.app_combo.activated.connect(self._on_switch_app)
            top.addWidget(self.app_combo)
        top.addWidget(QtWidgets.QLabel("Units:"))
        self.rb_si = QtWidgets.QRadioButton("SI")
        self.rb_us = QtWidgets.QRadioButton("US")
        (self.rb_us if self.system == "US" else self.rb_si).setChecked(True)
        self.rb_si.toggled.connect(self._on_units_changed)
        top.addWidget(self.rb_si)
        top.addWidget(self.rb_us)
        top.addStretch(1)
        badge = QtWidgets.QLabel(_CLASS_LETTER.get(self.meta.classification, self.meta.classification))
        badge.setToolTip(self.meta.classification)
        badge.setProperty("badge", _CLASS_BADGE.get(self.meta.classification, "accent"))
        top.addWidget(QtWidgets.QLabel(self.meta.area))
        top.addWidget(badge)
        btn_reset = QtWidgets.QPushButton("Reset")
        btn_reset.clicked.connect(self._on_reset)
        top.addWidget(btn_reset)
        outer.addLayout(top)

        # middle: inputs | compute | results
        split = QtWidgets.QSplitter(Qt.Horizontal)
        in_scroll = QtWidgets.QScrollArea()        # tall input forms (e.g. a table) scroll
        in_scroll.setWidgetResizable(True)
        in_scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        in_scroll.setWidget(self._build_inputs_box())
        split.addWidget(in_scroll)
        split.addWidget(self._build_center())
        split.addWidget(self._build_results_box())
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 1)
        split.setStretchFactor(2, 4)

        # "Method & equations" panel (per-app ABOUT). Built after the splitter so the
        # input widgets exist for active-method highlighting, but shown above it.
        method_panel = self._build_method_panel()
        if method_panel is not None:
            outer.addWidget(method_panel)
        outer.addWidget(split, 1)

        self.status = self.statusBar()
        self._set_status("ready", ok=True)

    def _build_inputs_box(self) -> QtWidgets.QGroupBox:
        box = QtWidgets.QGroupBox("Inputs")
        form = QtWidgets.QFormLayout(box)
        form.setLabelAlignment(Qt.AlignRight)
        self._inputs_form = form
        self._input_rows = {}
        for f in self.inputs:
            w = self._make_widget(f)
            self._widgets[f.key] = w
            lbl = self._field_label(f.label)
            lbl.setToolTip(f.note or "")
            form.addRow(lbl, w)
            self._input_rows[f.key] = (lbl, w)
        # reactive show_if / enable_if: re-evaluate when a controlling input changes
        controllers = {c[0] for f in self.inputs if (c := getattr(f, "show_if", ()))}
        controllers |= {c[0] for f in self.inputs if (c := getattr(f, "enable_if", ()))}
        for k in controllers:
            cw = self._widgets.get(k)
            if isinstance(cw, QtWidgets.QComboBox):
                cw.currentIndexChanged.connect(lambda *_: (self._apply_show_if(), self._apply_enable_if()))
            elif isinstance(cw, QtWidgets.QCheckBox):
                cw.toggled.connect(lambda *_: (self._apply_show_if(), self._apply_enable_if()))
        self._apply_show_if()
        self._apply_enable_if()
        return box

    def _apply_enable_if(self):
        """Gray out (disable) inputs whose `enable_if = (other_key, value)` is not met
        (e.g. the Van der Meer parameters while the Hudson method is selected). The row
        stays visible, unlike show_if, so the user sees it exists but isn't used."""
        for f in self.inputs:
            cond = getattr(f, "enable_if", ())
            if not cond:
                continue
            cw = self._widgets.get(cond[0])
            if isinstance(cw, QtWidgets.QComboBox):
                cur = cw.currentText()
            elif isinstance(cw, QtWidgets.QCheckBox):
                cur = cw.isChecked()
            else:
                cur = None
            on = str(cur) == str(cond[1])
            lbl, w = self._input_rows[f.key]
            w.setEnabled(on)
            lbl.setEnabled(on)

    def _apply_show_if(self):
        """Hide inputs whose `show_if = (other_key, value)` condition is not met
        (e.g. 'Specified slope' shows only when Slope source = Specified slope)."""
        for f in self.inputs:
            cond = getattr(f, "show_if", ())
            if not cond:
                continue
            cw = self._widgets.get(cond[0])
            if isinstance(cw, QtWidgets.QComboBox):
                cur = cw.currentText()
            elif isinstance(cw, QtWidgets.QCheckBox):
                cur = cw.isChecked()
            else:
                cur = None
            vis = str(cur) == str(cond[1])
            lbl, w = self._input_rows[f.key]
            try:
                self._inputs_form.setRowVisible(w, vis)
            except (AttributeError, TypeError):
                lbl.setVisible(vis); w.setVisible(vis)

    # --- "Method & equations" panel ------------------------------------------------
    def _build_method_panel(self):
        about = getattr(self.app, "ABOUT", None)
        if not about:
            return None
        self._about = about
        self._method_sections = {}
        container = QtWidgets.QWidget()
        cv = QtWidgets.QVBoxLayout(container)
        cv.setContentsMargins(0, 0, 0, 0)
        cv.setSpacing(0)
        self._method_toggle = QtWidgets.QToolButton()
        self._method_toggle.setText("Method & equations")
        self._method_toggle.setCheckable(True)
        self._method_toggle.setChecked(False)
        self._method_toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self._method_toggle.setArrowType(Qt.RightArrow)
        self._method_toggle.setAutoRaise(True)
        f = self._method_toggle.font(); f.setBold(True); self._method_toggle.setFont(f)
        self._method_toggle.toggled.connect(self._on_method_toggle)
        cv.addWidget(self._method_toggle)
        self._method_body = QtWidgets.QFrame()
        self._method_body.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self._method_body.setVisible(False)
        cv.addWidget(self._method_body)
        self._populate_method_body()
        # re-highlight the active formulation when the method-choice input changes
        key = about.get("method_key")
        if key and len(about.get("methods", ())) > 1:
            w = self._widgets.get(key)
            if isinstance(w, QtWidgets.QComboBox):
                w.currentIndexChanged.connect(lambda *_: self._highlight_active_method())
        return container

    def _on_method_toggle(self, on: bool):
        self._method_toggle.setArrowType(Qt.DownArrow if on else Qt.RightArrow)
        self._method_body.setVisible(on)

    def _populate_method_body(self):
        about = self._about
        pal = get_plot_palette(settings.get_theme())
        fg, muted = pal["fg"], pal["text"]
        lay = QtWidgets.QVBoxLayout(self._method_body)
        lay.setContentsMargins(16, 10, 16, 12)
        lay.setSpacing(8)
        if about.get("summary"):
            s = QtWidgets.QLabel(about["summary"])
            s.setWordWrap(True)
            lay.addWidget(s)
        for m in about.get("methods", ()):
            sec = QtWidgets.QFrame()
            sv = QtWidgets.QVBoxLayout(sec)
            sv.setContentsMargins(0, 6, 0, 6)
            sv.setSpacing(4)
            head = QtWidgets.QLabel(m.get("name", ""))
            hf = head.font(); hf.setBold(True); head.setFont(hf)
            tag = m.get("tag")
            if tag:
                head.setText(f"{m.get('name', '')}    [{tag}]")
            sv.addWidget(head)
            if m.get("note"):
                n = QtWidgets.QLabel(m["note"]); n.setWordWrap(True)
                n.setStyleSheet(f"color:{muted};")
                sv.addWidget(n)
            for eq in m.get("equations", ()):
                row = QtWidgets.QHBoxLayout()
                row.setSpacing(12)
                lbl = QtWidgets.QLabel()
                pm = _latex_pixmap(eq.get("tex", ""), fg, self.devicePixelRatioF() or 2.0)
                if pm is not None:
                    lbl.setPixmap(pm)
                else:
                    lbl.setText(eq.get("tex", ""))
                row.addWidget(lbl)
                if eq.get("desc"):
                    d = QtWidgets.QLabel(eq["desc"])
                    d.setStyleSheet(f"color:{muted};")
                    row.addWidget(d)
                row.addStretch(1)
                sv.addLayout(row)
            lay.addWidget(sec)
            self._method_sections[m.get("when", "")] = sec
        if about.get("symbols"):
            cap = QtWidgets.QLabel("Symbols")
            cf = cap.font(); cf.setBold(True); cap.setFont(cf)
            lay.addWidget(cap)
            grid = QtWidgets.QGridLayout()
            grid.setHorizontalSpacing(14); grid.setVerticalSpacing(2)
            for i, (sym, desc) in enumerate(about["symbols"]):
                sy = QtWidgets.QLabel()
                pm = _latex_pixmap(_sym_to_tex(sym), fg, self.devicePixelRatioF() or 2.0, 11)
                if pm is not None:
                    sy.setPixmap(pm)
                else:
                    sy.setText(sym)
                ds = QtWidgets.QLabel(desc); ds.setStyleSheet(f"color:{muted};")
                grid.addWidget(sy, i, 0, Qt.AlignVCenter); grid.addWidget(ds, i, 1)
            grid.setColumnStretch(1, 1)
            lay.addLayout(grid)
        if about.get("references"):
            refs = QtWidgets.QLabel("References: " + "  ·  ".join(about["references"]))
            refs.setWordWrap(True)
            refs.setStyleSheet(f"color:{muted};")
            lay.addWidget(refs)
        self._highlight_active_method()

    def _highlight_active_method(self):
        about = getattr(self, "_about", None)
        if not about:
            return
        key = about.get("method_key")
        methods = about.get("methods", ())
        single = (not key) or len(methods) < 2
        active = None
        if not single:
            w = self._widgets.get(key)
            if isinstance(w, QtWidgets.QComboBox):
                active = w.currentText()
        for when, sec in self._method_sections.items():
            on = single or (when == str(active))
            eff = sec.graphicsEffect()
            if not isinstance(eff, QtWidgets.QGraphicsOpacityEffect):
                eff = QtWidgets.QGraphicsOpacityEffect(sec)
                sec.setGraphicsEffect(eff)
            eff.setOpacity(1.0 if on else 0.4)

    def _build_center(self) -> QtWidgets.QWidget:
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        v.addStretch(1)
        self.btn_compute = QtWidgets.QPushButton("▶  Compute")
        self.btn_compute.setProperty("accent", "true")   # -> spec QPushButton[accent="true"]
        self.btn_compute.setDefault(True)
        self.btn_compute.clicked.connect(self._on_compute)
        v.addWidget(self.btn_compute)
        v.addStretch(1)
        return w

    def _build_results_box(self) -> QtWidgets.QGroupBox:
        box = QtWidgets.QGroupBox("Outputs")
        v = QtWidgets.QVBoxLayout(box)
        # scalar + point value rows
        rows = QtWidgets.QWidget()
        self._rows_form = QtWidgets.QFormLayout(rows)
        for o in self.outputs:
            if o.kind in ("profile", "grid", "vline", "scatter", "scatter_x", "data"):
                continue
            lab = QtWidgets.QLabel("-")
            lab.setObjectName("resultValue")
            lab.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self._value_labels[o.key] = lab
            name = self._field_label(o.label)
            self._rows_form.addRow(name, lab)
        v.addWidget(rows)
        # area-10 workflow: "send to next app" buttons (carry the series forward)
        nexts = getattr(self.meta, "next_apps", ()) or ()
        if nexts and self._on_switch is not None:
            wf = QtWidgets.QHBoxLayout()
            for tid, tlabel in nexts:
                b = QtWidgets.QPushButton(f"Send to {tid} {tlabel} →")
                b.setObjectName("WorkflowBtn")
                b.clicked.connect(lambda _=False, t=tid: self._next_step(t))
                wf.addWidget(b)
            wf.addStretch(1)
            v.addLayout(wf)
        # tabs: Plot | Table  (only if the app has profile or grid outputs)
        self._has_profiles = any(o.kind == "profile" for o in self.outputs)
        self._has_grid = any(o.kind == "grid" for o in self.outputs)
        if self._has_profiles or self._has_grid:
            tabs = QtWidgets.QTabWidget()
            self.figure = Figure(figsize=(4.6, 3.6), dpi=110, constrained_layout=True)
            self.canvas = FigureCanvas(self.figure)
            self.canvas.setMinimumHeight(240)
            # plot tab = navigation toolbar (pan / zoom / home / save) over the canvas
            plot_tab = QtWidgets.QWidget()
            pv = QtWidgets.QVBoxLayout(plot_tab)
            pv.setContentsMargins(0, 0, 0, 0)
            pv.setSpacing(0)
            self._nav = NavigationToolbar(self.canvas, plot_tab)
            pv.addWidget(self._nav)
            pv.addWidget(self.canvas, 1)
            tabs.addTab(plot_tab, "Plot")
            self.table = QtWidgets.QTableWidget()
            tabs.addTab(self.table, "Table")
            tabs.setMinimumHeight(270)
            v.addWidget(tabs, 1)
        # export bar: Report (all tools); Copy/CSV/PNG for graph tools
        btns = QtWidgets.QHBoxLayout()
        btns.addStretch(1)
        b_report = QtWidgets.QPushButton("Save Report")
        b_report.clicked.connect(self._save_report)
        btns.addWidget(b_report)
        if self._has_profiles or self._has_grid:
            b_copy = QtWidgets.QPushButton("Copy")
            b_copy.clicked.connect(self._copy_table)
            b_csv = QtWidgets.QPushButton("Export CSV")
            b_csv.clicked.connect(self._export_csv)
            b_png = QtWidgets.QPushButton("Save Plot")
            b_png.clicked.connect(self._save_png)
            b_pop = QtWidgets.QPushButton("Pop out")
            b_pop.clicked.connect(self._pop_out_plot)
            btns.addWidget(b_copy)
            btns.addWidget(b_csv)
            btns.addWidget(b_png)
            btns.addWidget(b_pop)
        v.addLayout(btns)
        return box

    def _make_widget(self, f) -> QtWidgets.QWidget:
        if f.kind == "choice":
            cb = QtWidgets.QComboBox()
            cb.addItems([str(c) for c in f.choices])
            if f.default in f.choices:
                cb.setCurrentIndex(list(f.choices).index(f.default))
            return cb
        if f.kind == "bool":
            return QtWidgets.QCheckBox()
        if f.kind == "table":
            return self._make_table(f)
        if f.kind == "csv":
            return self._make_csv(f)
        if f.kind in ("list", "matrix"):             # structured numeric input edited as raw JSON
            ed = QtWidgets.QPlainTextEdit()
            ed.setPlaceholderText("JSON (SI units); empty = use the app's built-in geometry")
            if f.default is not None:
                ed.setPlainText(json.dumps(f.default))
            ed.setFixedHeight(96 if f.kind == "matrix" else 56)
            return ed
        # float / int / angle -> double spin box (values in current display unit)
        sb = QtWidgets.QDoubleSpinBox()
        unit = self._unit(f)
        dec = self._input_decimals(f, unit)          # enough precision for small values
        sb.setDecimals(dec)
        lo = units.from_si(_clamp_bound(f.lo), unit)
        hi = units.from_si(_clamp_bound(f.hi), unit)
        sb.setRange(min(lo, hi), max(lo, hi))
        sb.setValue(units.from_si(float(f.default), unit))
        if unit:
            sb.setSuffix(f" {unit}")
        sb.setSingleStep(0.1 if dec <= 2 else 10 ** -(dec - 1))
        return sb

    # ---- CSV-record field (kind == "csv") ----
    def _make_csv(self, f) -> QtWidgets.QWidget:
        """A bundled-station dropdown + file upload; the loaded CSV text is held on
        the container (`_csv_text`) and read back by _gather_si."""
        box = QtWidgets.QWidget()
        lay = QtWidgets.QHBoxLayout(box)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)
        combo = QtWidgets.QComboBox()
        combo.addItem("— sample (built-in) —")
        for c in f.choices:
            sid, _, name = str(c).partition("|")
            combo.addItem(f"{sid} — {name}" if name else sid)
            combo.setItemData(combo.count() - 1, sid)
        combo.addItem("Upload a file…")
        combo.setItemData(combo.count() - 1, "__upload__")
        status = QtWidgets.QLabel("built-in sample")
        status.setObjectName("CsvStatus")
        box._combo = combo
        box._status = status
        box._csv_text = "" if f.default is None else str(f.default)
        box._default = box._csv_text
        box._last_idx = 0
        combo.currentIndexChanged.connect(
            lambda idx, w=box, fld=f: self._on_csv_changed(w, fld, idx))
        lay.addWidget(combo, 1)
        lay.addWidget(status, 2)
        return box

    def _on_csv_changed(self, w, f, idx: int) -> None:
        if idx == 0:                                  # built-in sample
            w._csv_text = "" if f.default is None else str(f.default)
            w._status.setText("built-in sample")
            w._last_idx = 0
            self._on_compute()
            return
        data = w._combo.itemData(idx)
        if data == "__upload__":
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Select water-level CSV", "", "CSV files (*.csv);;All files (*)")
            if not path:                              # cancelled -> revert selection
                w._combo.blockSignals(True)
                w._combo.setCurrentIndex(w._last_idx)
                w._combo.blockSignals(False)
                return
            try:
                with open(path, encoding="utf-8", errors="replace") as fh:
                    w._csv_text = fh.read()
            except OSError as e:
                w._status.setText(f"read failed: {e}")
                return
            w._status.setText(f"uploaded {os.path.basename(path)}")
        else:                                         # bundled station id
            sub = getattr(f, "data_dir", "water_levels") or "water_levels"
            try:
                with open(os.path.join(_DATA_ROOT, sub, f"{data}.csv"),
                          encoding="utf-8", errors="replace") as fh:
                    w._csv_text = fh.read()
            except OSError as e:
                w._status.setText(f"load failed: {e}")
                return
            w._status.setText(f"loaded {w._combo.itemText(idx)}")
        w._last_idx = idx
        self._on_compute()

    def _input_decimals(self, f, unit) -> int:
        """Decimal places for an input spin box: at least the global display count,
        but enough that the default / bounds for this field are not rounded to zero
        (the global decimals controls *output* formatting, not input precision)."""
        vals = [units.from_si(float(f.default), unit)]
        for b in (f.lo, f.hi):
            if math.isfinite(b):
                vals.append(units.from_si(b, unit))
        need = self.decimals
        for v in vals:
            v = abs(v)
            if 0 < v < 1:
                d = self.decimals
                while d < 8 and round(v, d) <= 0:    # first place where it shows
                    d += 1
                need = max(need, min(d + 1, 8))       # +1 for a little resolution
        return need

    # ---- table-input field (kind == "table") ----
    def _table_cols(self, f) -> list:
        """Column specs (label, unit_si, unit_us). A table field may declare its own
        `columns`; otherwise it is a single column taking the field's own units."""
        cols = getattr(f, "columns", ()) or ()
        return [tuple(c) for c in cols] if cols else [(f.label, f.unit_si, f.unit_us)]

    def _table_headers(self, tbl, cols, system):
        heads = []
        for lab, usi, uus in cols:
            u = self._si_unit(usi) if system == "SI" else uus
            heads.append(f"{lab} ({u})" if u else lab)
        tbl.setHorizontalHeaderLabels(heads)
        tbl.horizontalHeader().setStretchLastSection(True)

    def _fill_table(self, tbl, cols, rows_si, system):
        tbl.setRowCount(len(rows_si))
        for r, row in enumerate(rows_si):
            for c, (lab, usi, uus) in enumerate(cols):
                u = self._si_unit(usi) if system == "SI" else uus
                si = float(row[c]) if c < len(row) else 0.0
                tbl.setItem(r, c, QtWidgets.QTableWidgetItem(self._fmt(units.from_si(si, u))))

    def _make_table(self, f) -> QtWidgets.QWidget:
        cols = self._table_cols(f)
        wrap = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(wrap)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        tbl = QtWidgets.QTableWidget(0, len(cols))
        tbl.setMinimumHeight(160)
        tbl.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        tbl.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self._table_headers(tbl, cols, self.system)
        self._fill_table(tbl, cols, list(f.default or []), self.system)
        bar = QtWidgets.QHBoxLayout()
        add = QtWidgets.QPushButton("+ row")
        rem = QtWidgets.QPushButton("− row")
        add.clicked.connect(lambda: tbl.insertRow(tbl.rowCount()))
        rem.clicked.connect(lambda: tbl.rowCount() and tbl.removeRow(tbl.rowCount() - 1))
        bar.addStretch(1)
        bar.addWidget(add)
        bar.addWidget(rem)
        lay.addWidget(tbl)
        lay.addLayout(bar)
        wrap._table = tbl           # accessed by gather / unit-toggle / reset
        wrap._cols = cols
        return wrap

    def _gather_table(self, wrap) -> list:
        """Read a table widget -> list of rows (each a list of SI cell values).
        Blank rows are skipped; blank cells in a kept row become 0.0."""
        tbl, cols, rows = wrap._table, wrap._cols, []
        for r in range(tbl.rowCount()):
            vals, blank = [], True
            for c, (lab, usi, uus) in enumerate(cols):
                it = tbl.item(r, c)
                txt = it.text().strip() if it and it.text() else ""
                if txt == "":
                    vals.append(0.0)
                else:
                    blank = False
                    u = self._si_unit(usi) if self.system == "SI" else uus
                    vals.append(units.to_si(float(txt), u))
            if not blank:
                rows.append(vals)
        return rows

    # ---- compute / render ----
    def _gather_si(self) -> dict:
        out = {}
        for f in self.inputs:
            w = self._widgets[f.key]
            if f.kind == "choice":
                out[f.key] = w.currentText()
            elif f.kind == "bool":
                out[f.key] = w.isChecked()
            elif f.kind == "table":
                out[f.key] = self._gather_table(w)
            elif f.kind == "csv":
                out[f.key] = getattr(w, "_csv_text", "")
            elif f.kind in ("list", "matrix"):       # raw JSON (SI); empty -> None (built-in geometry)
                txt = w.toPlainText().strip()
                out[f.key] = json.loads(txt) if txt else None
            else:
                out[f.key] = units.to_si(w.value(), self._unit(f))
        return out

    def _on_compute(self):
        try:
            inp = self._gather_si()
            self._last_result = self.app.compute(inp)
        except Exception as exc:  # validation or numeric error
            self._set_status(f"error: {exc}", ok=False)
            return
        self._render(self._last_result)

    def _render(self, r):
        # scalar/point value rows (convert SI -> current unit)
        for o in self.outputs:
            if o.kind in ("profile", "grid", "vline", "scatter", "scatter_x", "data"):
                continue
            val_si = getattr(r, o.key, None)
            if val_si is None:
                continue
            u = self._out_unit(o)
            try:                                   # numeric -> convert to display units + format
                txt = f"{self._fmt(units.from_si(float(val_si), u))} {u}".rstrip()
            except (TypeError, ValueError):        # string/sentinel output (e.g. a regime label) -> show as-is
                txt = f"{val_si} {u}".rstrip() if u else str(val_si)
            self._value_labels[o.key].setText(txt)
        if self._has_grid:
            self._render_heatmap(r)
            self._render_grid_table(r)
        elif self._has_profiles:
            self._render_plot(r)
            self._render_table(r)
        note = getattr(r, "notes", "")
        self._set_status(f"units: {self.system} · valid · {note} · {self.meta.cite}", ok=True)

    def apply_theme(self):
        """Re-render the matplotlib plot for the current Light/Dark mode. Qt widgets
        are restyled by the application stylesheet, but plot colors are baked in at
        draw time, so the canvas must be redrawn when the mode toggles."""
        if self._last_result is not None:
            if self._has_grid:
                self._render_heatmap(self._last_result)
            elif self._has_profiles:
                self._render_plot(self._last_result)

    def _style_ax(self, ax, pal):
        ax.set_facecolor(pal["bg"])
        for s in ax.spines.values():
            s.set_color(pal["axis"])
        ax.tick_params(labelsize=7, colors=pal["text"])
        ax.grid(True, color=pal["grid"], alpha=.9, linewidth=.6)
        ax.locator_params(axis="y", nbins=4)

    def _profile_series(self, r):
        """Discover the profile outputs and return (x, groups), shared by the plot and
        the table. `x` = (label, unit, array) for the independent axis (`profile_X`);
        `groups` = y-series grouped by display unit, each {unit, series:[(label, array)]}.
        Labels and units come from the Out descriptors, so any app is self-describing."""
        profs = [o for o in self.outputs if o.kind == "profile"]
        xo = next((o for o in profs if o.key == "profile_X"), profs[0])
        short = lambda o: o.label.split(":", 1)[-1].strip()   # drop the "Profile: " prefix
        x = (short(xo), self._out_unit(xo),
             units.from_si(np.asarray(getattr(r, xo.key), dtype=float), self._out_unit(xo)))
        ys = [o for o in profs if o is not xo]
        groups = []
        if any(getattr(o, "group", "") for o in ys):
            # explicit grouping: profiles sharing a `group` plot on the same panel,
            # colored by within-panel index so overlays (data + trend) contrast
            pal = get_plot_palette(settings.get_theme())
            ecolors = [pal["eta"], pal["w"], pal["fg"], pal["text"]]
            by = {}
            for o in ys:
                gk = getattr(o, "group", "") or o.key
                u = self._out_unit(o)
                arr = units.from_si(np.asarray(getattr(r, o.key), dtype=float), u)
                if gk not in by:
                    by[gk] = {"unit": u, "gid": gk, "series": []}
                    groups.append(by[gk])
                idx = len(by[gk]["series"])
                by[gk]["series"].append((short(o), arr, ecolors[idx % len(ecolors)]))
            return x, groups
        for i, o in enumerate(ys):
            u = self._out_unit(o)
            arr = units.from_si(np.asarray(getattr(r, o.key), dtype=float), u)
            # the first y-series gets its own panel (prominence + scale); the rest group
            # by unit so same-quantity series (e.g. u & w, or two depths) share an axis
            for g in (groups[1:] if i else []):
                if g["unit"] == u:
                    g["series"].append((short(o), arr))
                    break
            else:
                groups.append({"unit": u, "series": [(short(o), arr)]})
        return x, groups

    def _render_heatmap(self, r, fig=None, canvas=None):
        """Render the first grid output as a heatmap (parity with the web drawHeatmap):
        grid_x / grid_y give the axes; the 2-D field has shape (len(grid_y), len(grid_x))."""
        fig = fig if fig is not None else self.figure
        canvas = canvas if canvas is not None else self.canvas
        pal = get_plot_palette(settings.get_theme())
        fig.clear()
        fig.set_facecolor(pal["bg"])
        profs = [o for o in self.outputs if o.kind == "profile"]
        grids = [o for o in self.outputs if o.kind == "grid"]
        if not profs or not grids:
            canvas.draw_idle(); return
        xo = next((o for o in profs if o.key == "grid_x"), profs[0])
        yo = next((o for o in profs if o.key == "grid_y"), profs[1] if len(profs) > 1 else profs[0])
        zo = grids[0]
        short = lambda o: o.label.split(":", 1)[-1].strip()
        xu, yu, zu = self._out_unit(xo), self._out_unit(yo), self._out_unit(zo)
        xs = np.asarray(units.from_si(getattr(r, xo.key), xu), dtype=float)
        ys = np.asarray(units.from_si(getattr(r, yo.key), yu), dtype=float)
        Z = np.asarray(units.from_si(getattr(r, zo.key), zu), dtype=float)
        ax = fig.add_subplot(1, 1, 1)
        mesh = ax.pcolormesh(xs, ys, Z, shading="nearest", cmap="viridis")
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.ax.tick_params(labelsize=7, colors=pal["text"])
        cbar.set_label(f"{short(zo)}{(' (' + zu + ')') if zu else ''}", fontsize=8, color=pal["fg"])
        ax.set_xlabel(f"{short(xo)} ({xu})", fontsize=8, color=pal["fg"])
        ax.set_ylabel(f"{short(yo)} ({yu})", fontsize=8, color=pal["fg"])
        self._style_ax(ax, pal)
        ax.grid(False)
        canvas.draw_idle()

    def _render_grid_table(self, r):
        """Fill the data table with the first grid output (rows = grid_y, cols = grid_x)."""
        profs = [o for o in self.outputs if o.kind == "profile"]
        grids = [o for o in self.outputs if o.kind == "grid"]
        if not profs or not grids:
            return
        xo = next((o for o in profs if o.key == "grid_x"), profs[0])
        yo = next((o for o in profs if o.key == "grid_y"), profs[1] if len(profs) > 1 else profs[0])
        zo = grids[0]
        xs = np.asarray(units.from_si(getattr(r, xo.key), self._out_unit(xo)), dtype=float)
        ys = np.asarray(units.from_si(getattr(r, yo.key), self._out_unit(yo)), dtype=float)
        Z = np.asarray(units.from_si(getattr(r, zo.key), self._out_unit(zo)), dtype=float)
        self.table.clear()
        self.table.setRowCount(len(ys))
        self.table.setColumnCount(len(xs))
        self.table.setHorizontalHeaderLabels([self._fmt(v) for v in xs])
        self.table.setVerticalHeaderLabels([self._fmt(v) for v in ys])
        for j in range(len(ys)):
            for i in range(len(xs)):
                self.table.setItem(j, i, QtWidgets.QTableWidgetItem(self._fmt(float(Z[j][i]))))

    def _render_plot(self, r, fig=None, canvas=None):
        fig = fig if fig is not None else self.figure
        canvas = canvas if canvas is not None else self.canvas
        pal = get_plot_palette(settings.get_theme())
        fig.clear()
        fig.set_facecolor(pal["bg"])
        (xlab, xu, X), groups = self._profile_series(r)
        # vertical markers (kind "vline"), e.g. the NTDE midpoint; NaN -> not drawn
        vlines = []
        for o in self.outputs:
            if o.kind != "vline":
                continue
            val = getattr(r, o.key, None)
            try:
                fv = units.from_si(float(val), xu)
            except (TypeError, ValueError):
                continue
            if math.isfinite(fv):
                vlines.append((o.label.split(":", 1)[-1].strip(), fv))
        # scatter overlays (kind "scatter") -> markers at their own (x_key, key), by group
        scatters = []
        for o in self.outputs:
            if o.kind != "scatter" or not getattr(o, "x_key", ""):
                continue
            xs = units.from_si(np.asarray(getattr(r, o.x_key), dtype=float), xu)
            ys = units.from_si(np.asarray(getattr(r, o.key), dtype=float), self._out_unit(o))
            scatters.append((getattr(o, "group", ""), o.label.split(":", 1)[-1].strip(), xs, ys))
        colors = [pal["eta"], pal["u"], pal["w"], pal["fg"], pal["text"]]   # series color cycle
        n, ci, axes = len(groups), 0, []
        for gi, g in enumerate(groups):
            ax = fig.add_subplot(n, 1, gi + 1, sharex=axes[0] if axes else None)
            for s in g["series"]:
                lab, arr = s[0], s[1]
                color = s[2] if len(s) > 2 and s[2] else colors[ci % len(colors)]
                if not (len(s) > 2 and s[2]):
                    ci += 1
                ax.plot(X, arr, label=lab, color=color)
            for sgid, slab, sxs, sys_ in scatters:
                if sgid == g.get("gid"):
                    ax.plot(sxs, sys_, "o", color="#d62728", markersize=3.5,
                            markeredgecolor="white", markeredgewidth=0.3, label=slab, zorder=5)
            ax.axhline(0, color=pal["axis"], lw=.6)
            for vlab, vx in vlines:
                ax.axvline(vx, color=pal["fg"], ls="--", lw=1.0, alpha=.75)
                if gi == 0:
                    ax.text(vx, 1.0, f" {vlab} {vx:g}", transform=ax.get_xaxis_transform(),
                            fontsize=7, color=pal["fg"], va="top", ha="left", clip_on=True)
            ylab = g["series"][0][0] + " " if len(g["series"]) == 1 else ""
            ax.set_ylabel(f"{ylab}({g['unit']})".strip(), fontsize=8, color=pal["fg"])
            self._style_ax(ax, pal)
            if len(g["series"]) > 1:
                leg = ax.legend(loc="upper right", fontsize=7, ncol=len(g["series"]),
                                facecolor=pal["bg"], edgecolor=pal["axis"], labelcolor=pal["fg"])
                leg.get_frame().set_alpha(.85)
            axes.append(ax)
        for ax in axes[:-1]:
            ax.tick_params(labelbottom=False)
        axes[-1].set_xlabel(f"{xlab} ({xu})", fontsize=8, color=pal["fg"])
        canvas.draw_idle()

    def _next_step(self, target_id: str):
        """Re-run this app for its full-resolution series, stow it, and open the
        target app, which injects the series into its first CSV field on open."""
        try:
            inp = self._gather_si()
            inp["handoff"] = "1"
            r = self.app.compute(inp)
        except Exception as exc:
            self._set_status(f"cannot hand off: {exc}", ok=False)
            return
        _HANDOFF["series"] = {"csv": getattr(r, "handoff_csv", ""),
                              "label": f"from {self.meta.aces_id} {self.meta.name}"}
        if self._on_switch is not None:
            self._on_switch(target_id)

    def _pop_out_plot(self):
        """Open the current plot in a separate, larger, resizable window with its own
        navigation toolbar (pan / zoom / save)."""
        if self._last_result is None:
            return
        win = QtWidgets.QMainWindow(self)
        win.setWindowTitle(f"CHESS-QC {self.meta.aces_id} {self.meta.name} — plot")
        win.resize(940, 680)
        fig = Figure(figsize=(8.5, 6.0), dpi=110, constrained_layout=True)
        canvas = FigureCanvas(fig)
        central = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(central)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(NavigationToolbar(canvas, win))
        lay.addWidget(canvas, 1)
        win.setCentralWidget(central)
        if self._has_grid:
            self._render_heatmap(self._last_result, fig, canvas)
        else:
            self._render_plot(self._last_result, fig, canvas)
        self._popouts = getattr(self, "_popouts", [])
        self._popouts.append(win)            # keep a reference so it is not GC'd
        win.show()

    def _render_table(self, r):
        (xlab, xu, X), groups = self._profile_series(r)
        cols = [(f"{xlab} ({xu})", X)]
        for g in groups:
            for s in g["series"]:
                cols.append((f"{s[0]} ({g['unit']})", s[1]))
        self.table.setRowCount(len(X))
        self.table.setColumnCount(len(cols))
        self.table.setHorizontalHeaderLabels([h for h, _ in cols])
        for i in range(len(X)):
            for j, (_, arr) in enumerate(cols):
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(self._fmt(arr[i])))

    # ---- unit toggle ----
    def _on_units_changed(self):
        new = "SI" if self.rb_si.isChecked() else "US"
        if new == self.system:
            return
        old = self.system
        for f in self.inputs:
            if f.kind == "table":                    # reconvert cells + headers
                wrap = self._widgets[f.key]
                tbl, cols = wrap._table, wrap._cols
                for r in range(tbl.rowCount()):
                    for c, (lab, usi, uus) in enumerate(cols):
                        it = tbl.item(r, c)
                        if it and it.text().strip():
                            old_u = self._si_unit(usi) if old == "SI" else uus
                            new_u = self._si_unit(usi) if new == "SI" else uus
                            si = units.to_si(float(it.text()), old_u)
                            it.setText(self._fmt(units.from_si(si, new_u)))
                self._table_headers(tbl, cols, new)
                continue
            if f.kind not in ("float", "int", "angle"):
                continue
            w = self._widgets[f.key]
            old_u = self._si_unit(f.unit_si) if old == "SI" else f.unit_us
            new_u = self._si_unit(f.unit_si) if new == "SI" else f.unit_us
            si_val = units.to_si(w.value(), old_u)
            w.blockSignals(True)
            dec = self._input_decimals(f, new_u)
            w.setDecimals(dec)
            w.setSingleStep(0.1 if dec <= 2 else 10 ** -(dec - 1))
            lo = units.from_si(_clamp_bound(f.lo), new_u)
            hi = units.from_si(_clamp_bound(f.hi), new_u)
            w.setRange(min(lo, hi), max(lo, hi))
            w.setSuffix(f" {new_u}" if new_u else "")
            w.setValue(units.from_si(si_val, new_u))
            w.blockSignals(False)
        self.system = new
        if self._last_result is not None:
            self._render(self._last_result)

    def _on_reset(self):
        for f in self.inputs:
            w = self._widgets[f.key]
            if f.kind == "choice":
                if f.default in f.choices:
                    w.setCurrentIndex(list(f.choices).index(f.default))
            elif f.kind == "bool":
                w.setChecked(bool(f.default))
            elif f.kind == "table":
                self._fill_table(w._table, w._cols, list(f.default or []), self.system)
            elif f.kind == "csv":
                w._combo.blockSignals(True)
                w._combo.setCurrentIndex(0)
                w._combo.blockSignals(False)
                w._csv_text = "" if f.default is None else str(f.default)
                w._last_idx = 0
                w._status.setText("built-in sample")
            elif f.kind in ("list", "matrix"):
                w.setPlainText("" if f.default is None else json.dumps(f.default))
            else:
                w.blockSignals(True)
                w.setValue(units.from_si(float(f.default), self._unit(f)))
                w.blockSignals(False)
        self._on_compute()

    # ---- table copy/export ----
    def _table_tsv(self) -> str:
        cols = self.table.columnCount()
        heads = [self.table.horizontalHeaderItem(c).text() for c in range(cols)]
        lines = ["\t".join(heads)]
        for rrow in range(self.table.rowCount()):
            lines.append("\t".join(self.table.item(rrow, c).text() for c in range(cols)))
        return "\n".join(lines)

    def _copy_table(self):
        QtWidgets.QApplication.clipboard().setText(self._table_tsv())
        self._set_status("table copied to clipboard", ok=True)

    def _export_csv(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export CSV", "profile.csv", "CSV (*.csv)")
        if not path:
            return
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self._table_tsv().replace("\t", ","))
        self._set_status(f"saved {os.path.basename(path)}", ok=True)

    def _report_text(self) -> str:
        """Plain-text summary report (inputs + outputs + notes) in the current units."""
        m = self.meta
        lines = [f"CHESS-QC {m.aces_id}  {m.name}  [{_CLASS_LETTER.get(m.classification, m.classification)} {m.classification}]",
                 f"{m.cite}",
                 f"Functional area: {m.area}",
                 f"Units: {self.system}",
                 "", "Inputs:"]
        for f in self.inputs:
            w = self._widgets[f.key]
            if f.kind == "choice":
                val = w.currentText()
            elif f.kind == "bool":
                val = "yes" if w.isChecked() else "no"
            elif f.kind == "table":
                rows = self._gather_table(w)
                val = f"{len(rows)} rows"
            elif f.kind == "csv":
                txt = getattr(w, "_csv_text", "")
                val = f"{txt.count(chr(10)) + 1 if txt else 0} lines loaded"
            elif f.kind in ("list", "matrix"):
                txt = w.toPlainText().strip()
                val = "(custom JSON)" if txt else "(app default)"
            else:
                val = f"{w.value():g} {self._unit(f)}".rstrip()
            lines.append(f"  {f.label}: {val}")
        lines.append("")
        lines.append("Outputs:")
        for o in self.outputs:
            if o.kind == "profile":
                continue
            lab = self._value_labels.get(o.key)
            if lab is not None:
                lines.append(f"  {o.label}: {lab.text()}")
        note = getattr(self._last_result, "notes", "") if self._last_result else ""
        if note:
            lines += ["", f"Notes: {note}"]
        return "\n".join(lines)

    def _save_report(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save Report", f"chessqc_{self.meta.aces_id}_report.txt", "Text (*.txt)")
        if not path:
            return
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self._report_text())
        self._set_status(f"saved {os.path.basename(path)}", ok=True)

    def _save_png(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save Plot", f"chessqc_{self.meta.aces_id}_plot.png", "PNG (*.png)")
        if not path:
            return
        self.figure.savefig(path, dpi=200, facecolor=self.figure.get_facecolor())
        self._set_status(f"saved {os.path.basename(path)}", ok=True)

    def _set_status(self, msg: str, ok: bool):
        self.status.setStyleSheet("" if ok else "color:#c53030;")
        self.status.showMessage(("✓ " if ok else "✗ ") + msg)
