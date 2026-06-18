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

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

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
        self._speed_si = settings.get_speed_si()   # SI metric unit for m/s quantities (m/s|km/h)
        self._vibe_opts = vibe_label_opts(settings.get_vibe())
        self._widgets: dict[str, QtWidgets.QWidget] = {}
        self._value_labels: dict[str, QtWidgets.QLabel] = {}
        self._last_result = None

        self.setWindowTitle(f"CHESS-QC · {self.meta.aces_id} {self.meta.name}")
        self._build_ui()
        self.resize(960, 700)
        self._on_compute()  # populate with defaults
        self._focus_first_input()

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
        """Substitute the launcher-chosen SI speed unit (m/s|km/h) for m/s quantities."""
        return self._speed_si if u == "m/s" else u

    def _unit(self, f) -> str:
        return self._si_unit(f.unit_si) if self.system == "SI" else f.unit_us

    def _out_unit(self, o) -> str:
        return self._si_unit(o.unit_si) if self.system == "SI" else o.unit_us

    def _fmt(self, x: float) -> str:
        """Fixed-decimal display (kills the '-0.00' artifact)."""
        s = f"{float(x):.{self.decimals}f}"
        if float(s) == 0:
            s = f"{0.0:.{self.decimals}f}"
        return s

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
        outer.addWidget(split, 1)

        self.status = self.statusBar()
        self._set_status("ready", ok=True)

    def _build_inputs_box(self) -> QtWidgets.QGroupBox:
        box = QtWidgets.QGroupBox("Inputs")
        form = QtWidgets.QFormLayout(box)
        form.setLabelAlignment(Qt.AlignRight)
        for f in self.inputs:
            w = self._make_widget(f)
            self._widgets[f.key] = w
            lbl = self._field_label(f.label)
            lbl.setToolTip(f.note or "")
            form.addRow(lbl, w)
        return box

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
            if o.kind in ("profile", "grid"):
                continue
            lab = QtWidgets.QLabel("-")
            lab.setObjectName("resultValue")
            lab.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self._value_labels[o.key] = lab
            name = self._field_label(o.label)
            self._rows_form.addRow(name, lab)
        v.addWidget(rows)
        # tabs: Plot | Table  (only if the app has profile or grid outputs)
        self._has_profiles = any(o.kind == "profile" for o in self.outputs)
        self._has_grid = any(o.kind == "grid" for o in self.outputs)
        if self._has_profiles or self._has_grid:
            tabs = QtWidgets.QTabWidget()
            self.figure = Figure(figsize=(4.6, 3.6), constrained_layout=True)
            self.canvas = FigureCanvas(self.figure)
            self.canvas.setMinimumHeight(240)
            tabs.addTab(self.canvas, "Plot")
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
            btns.addWidget(b_copy)
            btns.addWidget(b_csv)
            btns.addWidget(b_png)
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
            if o.kind in ("profile", "grid"):
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
        self._set_status(f"valid · {note} · {self.meta.cite}", ok=True)

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

    def _render_heatmap(self, r):
        """Render the first grid output as a heatmap (parity with the web drawHeatmap):
        grid_x / grid_y give the axes; the 2-D field has shape (len(grid_y), len(grid_x))."""
        pal = get_plot_palette(settings.get_theme())
        self.figure.clear()
        self.figure.set_facecolor(pal["bg"])
        profs = [o for o in self.outputs if o.kind == "profile"]
        grids = [o for o in self.outputs if o.kind == "grid"]
        if not profs or not grids:
            self.canvas.draw_idle(); return
        xo = next((o for o in profs if o.key == "grid_x"), profs[0])
        yo = next((o for o in profs if o.key == "grid_y"), profs[1] if len(profs) > 1 else profs[0])
        zo = grids[0]
        short = lambda o: o.label.split(":", 1)[-1].strip()
        xu, yu, zu = self._out_unit(xo), self._out_unit(yo), self._out_unit(zo)
        xs = np.asarray(units.from_si(getattr(r, xo.key), xu), dtype=float)
        ys = np.asarray(units.from_si(getattr(r, yo.key), yu), dtype=float)
        Z = np.asarray(units.from_si(getattr(r, zo.key), zu), dtype=float)
        ax = self.figure.add_subplot(1, 1, 1)
        mesh = ax.pcolormesh(xs, ys, Z, shading="nearest", cmap="viridis")
        cbar = self.figure.colorbar(mesh, ax=ax)
        cbar.ax.tick_params(labelsize=7, colors=pal["text"])
        cbar.set_label(f"{short(zo)}{(' (' + zu + ')') if zu else ''}", fontsize=8, color=pal["fg"])
        ax.set_xlabel(f"{short(xo)} ({xu})", fontsize=8, color=pal["fg"])
        ax.set_ylabel(f"{short(yo)} ({yu})", fontsize=8, color=pal["fg"])
        self._style_ax(ax, pal)
        ax.grid(False)
        self.canvas.draw_idle()

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

    def _render_plot(self, r):
        pal = get_plot_palette(settings.get_theme())
        self.figure.clear()
        self.figure.set_facecolor(pal["bg"])
        (xlab, xu, X), groups = self._profile_series(r)
        colors = [pal["eta"], pal["u"], pal["w"], pal["fg"], pal["text"]]   # series color cycle
        n, ci, axes = len(groups), 0, []
        for gi, g in enumerate(groups):
            ax = self.figure.add_subplot(n, 1, gi + 1, sharex=axes[0] if axes else None)
            for lab, arr in g["series"]:
                ax.plot(X, arr, label=lab, color=colors[ci % len(colors)]); ci += 1
            ax.axhline(0, color=pal["axis"], lw=.6)
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
        self.canvas.draw_idle()

    def _render_table(self, r):
        (xlab, xu, X), groups = self._profile_series(r)
        cols = [(f"{xlab} ({xu})", X)]
        for g in groups:
            for lab, arr in g["series"]:
                cols.append((f"{lab} ({g['unit']})", arr))
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
            old_u = f.unit_si if old == "SI" else f.unit_us
            new_u = f.unit_si if new == "SI" else f.unit_us
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
