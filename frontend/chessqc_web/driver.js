/* CHESS-QC serverless web driver.
 * Boots Pyodide, loads the SAME chessqc_* application + bridge.py, reads the contract,
 * builds the calculator, runs compute() in-browser, and renders results + a canvas
 * plot of the profile arrays. No backend; works from a static host (or bundled offline).
 */
"use strict";

// Application manifest comes from apps.js (CHESSQC_APPS); add applications there.

// Fidelity classification -> compact badge letter (A exact, B standard, C provisional).
const CLASS_LETTER = { exact: "I", standard: "II", provisional: "III" };
const CLASS_KEY = { exact: "a", standard: "b", provisional: "c" };   // stable CSS/colour key
const classLetter = (c) => CLASS_LETTER[c] || c || "";
const classKey = (c) => CLASS_KEY[c] || "a";

// SI conversion: value_in_unit * TO_SI[unit] = value_SI  (mirrors chessqc_gui/units.py)
const TO_SI = {
  "": 1, "m": 1, "ft": 0.3048, "km": 1000, "mi": 1609.344, "nm": 1852, "hr": 3600,
  "mph": 0.44704, "kt": 0.514444, "C": 1, "hPa": 100, "mb": 100, "inHg": 3386.389,
  "s": 1, "m/s": 1, "km/h": 1000 / 3600, "ft/s": 0.3048,
  "m/s^2": 1, "ft/s^2": 0.3048, "Pa": 1, "psf": 47.880259,
  "N/m": 1, "lb/ft": 14.593903, "N/s": 1, "lb/s": 4.4482216, "deg": 1, "rad": 1,
  "%": 1, "N/m^3": 1, "kN/m^3": 1000, "lb/ft^3": 157.08746,
  "N": 1, "kN": 1000, "lb": 4.4482216, "tons": 8896.4432,
  "N-m/m": 1, "lb-ft/ft": 4.4482216,
  "m^3/yr": 1, "yd^3/yr": 0.764554858,
  "m/yr": 1, "ft/yr": 0.3048, "mm/yr": 0.001, "in/yr": 0.0254,
  "m^3": 1, "yd^3": 0.764554858, "phi": 1,
  "m^3/s/m": 1, "ft^3/s/ft": 0.09290304,
  // areas, flows, viscosity, inverse length, nautical mile
  "m^2": 1, "ft^2": 0.09290304,
  "m^3/s": 1, "ft^3/s": 0.028316846592,
  "m^2/s": 1, "ft^2/s": 0.09290304,
  "m^2/s^2": 1, "ft^2/s^2": 0.09290304,
  "N*s/m^2": 1, "lb*s/ft^2": 47.880259,
  "1/m": 1, "1/ft": 3.280839895,
  "nmi": 1852,
  // unit strings that are identical in SI and US (no conversion)
  "kg/m^3": 1, "deg C": 1, "1/s": 1, "1/yr": 1, "yr": 1, "h": 1, "min": 1, "kg": 1,
  "tonne": 1, "x1000": 1, "m^(1/3)": 1, "m^B": 1, "s or 1/s": 1, "g": 1, "mm": 1,
  "log10 yr": 1,
};
const f = (u) => (u in TO_SI ? TO_SI[u] : 1);
const toSI = (v, u) => v * f(u);
const fromSI = (v, u) => v / f(u);

const $ = (id) => document.getElementById(id);
// Focus the first editable input (skips the SI/US unit radios); used on load and reset.
const focusFirstInput = () => {
  const el = document.querySelector("#inputs input:not([type=radio]), #inputs select, #inputs textarea");
  if (el) el.focus();
};
let py, bridge, appmod, contract, system = "SI", lastRes = null;

// fixed-decimal display (user-chosen, default 2); kills the "-0.00" artifact
const decimals = () => (window.CHESSQCUI ? window.CHESSQCUI.getDecimals() : 2);
const fmtNum = (x) => {
  if (!Number.isFinite(x)) return String(x);
  const d = decimals();
  if (x === 0) return (0).toFixed(d);
  const s = x.toFixed(d);
  if (parseFloat(s) !== 0) return s;            // normal fixed-decimal display
  // nonzero but would round to 0 at d places -> show significant digits instead,
  // so small values (e.g. 0.003) are never displayed as 0 (no per-field config)
  return (+x.toPrecision(Math.max(d, 2))).toString();
};

// --- Parquet station records / uploads -----------------------------------------
// The bundled records are Parquet written by the PyStorm engines: the columns plus a
// `pystorm` footer carrying the station, product, units and vertical datum. Reading
// it here means the datum travels with the data instead of being retyped, and the
// calculator keeps taking plain (date, value) text the way it always has.
const PARQUET_MAGIC = "PAR1";

const isParquet = (buf) =>
  buf.byteLength > 4 && new TextDecoder().decode(new Uint8Array(buf, 0, 4)) === PARQUET_MAGIC;

// decimal year -> YYYY-MM-15, so a monthly-mean table lands on its own calendar month
const _monthStamp = (y, m) => `${y}-${String(m).padStart(2, "0")}-15`;

async function readParquetSeries(buf) {
  const api = window.CHESSQC_PARQUET;
  if (!api) throw new Error("Parquet reader not loaded");
  const meta = await api.parquetMetadataAsync(buf);
  const kv = {};
  for (const e of meta.key_value_metadata || []) kv[e.key] = e.value;
  let info = {};
  try { info = kv.pystorm ? JSON.parse(kv.pystorm) : {}; } catch (e) { info = {}; }
  const rows = await api.parquetReadObjects({ file: buf });
  if (!rows.length) throw new Error("Parquet file holds no rows");
  const cols = Object.keys(rows[0]);
  const has = (c) => cols.includes(c);
  const lines = ["date,value"];
  if (has("year") && has("month") && has("mean")) {
    // monthly-mean table: keep the engine's own completeness verdict, and write the
    // months it rejected as gaps so the fit excludes them rather than re-deriving
    for (const r of rows) {
      const ok = has("complete") ? r.complete : r.mean != null;
      lines.push(`${_monthStamp(r.year, r.month)},${ok && r.mean != null ? r.mean : ""}`);
    }
  } else {
    const tc = cols.find((c) => /date|time/i.test(c)) || cols[0];
    const vc = cols.find((c) => c !== tc && /level|value|mean|dwl|wl/i.test(c))
      || cols.find((c) => c !== tc);
    for (const r of rows) {
      const t = r[tc], v = r[vc];
      const stamp = t instanceof Date ? t.toISOString().slice(0, 16).replace("T", " ") : String(t);
      lines.push(`${stamp},${v == null ? "" : v}`);
    }
  }
  return { csv: lines.join("\n"), datum: info.datum_label || "", unit: info.unit_label || "",
           station: info.station_id || "", rows: rows.length };
}

// Fill the vertical-datum control from a record that carries one. A CSV cannot, so
// there the operator's own selection stands.
function applyDatum(datum) {
  if (!datum) return;
  const sel = document.querySelector('[data-key="vdatum"]');
  if (sel && [...sel.options].some((o) => o.value === datum)) {
    sel.value = datum;
    applyShowIf(); applyEnableIf();
  }
}

// The bridge sends non-finite scalars as string tokens (JSON has no Infinity/NaN).
// Show them as the symbols an engineer expects rather than the raw token: a neutrally
// stratified boundary layer really does have an infinite Monin-Obukhov length.
const _TOKENS = { inf: "∞", "-inf": "−∞", nan: "n/a" };
const fmtToken = (s) => _TOKENS[s] || String(s);

// INPUT-field formatter: preserve precision (up to 6 significant figures, no trailing
// zeros). The output `decimals` setting must NOT round inputs -- a small coefficient like
// 0.0025 would collapse to "0.00" and fail the app's range validation on load.
const fmtIn = (x) => {
  if (!Number.isFinite(x)) return String(x);
  if (x === 0) return "0";
  return String(+x.toPrecision(6));
};

// Wind / TC-translation speeds use km/h by default; the launcher can switch them to m/s.
// Only fields whose SI unit is km/h (the wind/TC speeds) are affected; water/wave m/s stays m/s.
const tcSpeed = () => ((window.CHESSQC_PREFS || {}).tcwind_speed === "m/s" ? "m/s" : "km/h");
const siUnit = (u) => (u === "km/h" ? tcSpeed() : u);
const fieldUnit = (fld) => (system === "SI" ? siUnit(fld.unit_si) : fld.unit_us);
const outUnit = (o) => (system === "SI" ? siUnit(o.unit_si) : o.unit_us);

async function boot() {
  try {
    py = await loadPyodide();
    // bridge.py is written and imported ONCE; app modules are (re)loaded per selection in
    // loadApp(), so Pyodide and numpy stay warm when switching apps (no per-app re-boot).
    // `no-cache` revalidates instead of serving from the browser cache: these two
    // sources ARE the calculator, so a stale copy silently computes the old answers
    // after a deploy (unchanged files still come back as a cheap 304).
    const brSrc = await (await fetch("../../common/bridge.py", { cache: "no-cache" })).text();
    py.FS.writeFile("/bridge.py", brSrc);
    py.runPython(`
import sys, importlib.util
_s = importlib.util.spec_from_file_location("bridge", "/bridge.py")
_m = importlib.util.module_from_spec(_s); sys.modules["bridge"] = _m; _s.loader.exec_module(_m)
`);
    bridge = py.pyimport("bridge");
    buildAppSelect();
    const id = new URLSearchParams(location.search).get("app") || (CHESSQC_APPS[0] && CHESSQC_APPS[0].id);
    await loadApp(id);
    $("overlay").style.display = "none";
  } catch (e) { fail(String(e)); }
}

// Populate the in-header application switcher (grouped by functional area).
function buildAppSelect() {
  const sel = $("appSelect");
  if (!sel) return;
  const byArea = {};
  for (const a of CHESSQC_APPS) { if (a.comingSoon) continue; (byArea[a.area] ||= []).push(a); }
  const order = CHESSQC_AREAS.filter((x) => x in byArea)
    .concat(Object.keys(byArea).filter((x) => !CHESSQC_AREAS.includes(x)));
  sel.innerHTML = "";
  for (const area of order) {
    const og = document.createElement("optgroup"); og.label = area;
    for (const a of byArea[area].sort((x, y) => x.id.localeCompare(y.id))) {
      const o = document.createElement("option");
      o.value = a.id; o.textContent = `${a.id} · ${a.name}`;
      og.appendChild(o);
    }
    sel.appendChild(og);
  }
  sel.addEventListener("change", () => {
    const id = sel.value;
    history.pushState({ app: id }, "", `?app=${encodeURIComponent(id)}`);
    loadApp(id);
  });
}

// Load (or switch to) an application in the already-running interpreter. Only the app
// module is re-imported; Pyodide, numpy and bridge.py stay warm.
async function loadApp(id) {
  const entry = CHESSQC_APPS.find((a) => a.id === id);
  if (!entry) { fail(`unknown app '${id}'`); return; }
  if (entry.comingSoon) { fail(`${id} ${entry.name} is coming soon`); return; }
  setStatus(`loading ${id}…`, true);
  const appSrc = await (await fetch(entry.src, { cache: "no-cache" })).text();
  // Load whatever the module actually imports (numpy, scipy, ...) and nothing more, so
  // stdlib-only apps stay fast. Pyodide scans the source itself: hard-coding the package
  // list here silently breaks any app that later picks up a new dependency.
  const pkgs = [...(entry.packages || [])];
  if (pkgs.length) await py.loadPackage(pkgs);
  await py.loadPackagesFromImports(appSrc);
  py.FS.writeFile("/chessqc_app.py", appSrc);
  py.runPython(`
import sys, importlib.util
sys.modules.pop("chessqc_app", None)   # drop the previous app so the re-import is fresh
_s = importlib.util.spec_from_file_location("chessqc_app", "/chessqc_app.py")
_m = importlib.util.module_from_spec(_s); sys.modules["chessqc_app"] = _m; _s.loader.exec_module(_m)
`);
  appmod = py.pyimport("chessqc_app");
  contract = JSON.parse(bridge.contract(appmod));
  // launcher-set units (CHESSQC_PREFS.units) take precedence; else the app's own default
  const prefU = (window.CHESSQC_PREFS || {}).units;
  const appSys = contract.meta && contract.meta.default_system === "US" ? "US" : "SI";
  system = prefU === "SI" || prefU === "US" ? prefU : appSys;
  const rb = document.querySelector(`input[name="units"][value="${system}"]`);
  if (rb) rb.checked = true;
  const sel = $("appSelect"); if (sel) sel.value = id;
  buildForm();
  clearOutputs();       // inputs only on open; the outputs follow from Compute
  focusFirstInput();
}

// Empty-output state: the app opens with its defaults loaded and nothing computed,
// so what is on screen is always the result of the inputs currently in the form.
function clearOutputs() {
  lastRes = null;
  const vbox = $("values");
  if (vbox) {
    vbox.innerHTML = "";
    const hint = document.createElement("div");
    hint.className = "await-compute";
    hint.textContent = "Press Compute to evaluate the inputs.";
    vbox.appendChild(hint);
  }
  const cv = $("plot");
  if (cv) {
    cv._res = null; cv._plotBox = null; cv._viewX = null;
    const ctx = cv.getContext("2d");
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, cv.width, cv.height);
  }
  const tp = $("tablePane"); if (tp) tp.innerHTML = "";
  setStatus(`units: ${system} · ready · ${contract.meta.cite}`, true);
}

function buildForm() {
  const pc = $("plot"); if (pc) pc._viewX = null;   // reset any prior zoom on app switch
  const m = contract.meta;
  $("appTitle").textContent = `CHESS-QC · ${m.aces_id} ${m.name}`;
  $("appArea").textContent = m.area;
  const badge = $("classBadge");
  const letter = classLetter(m.classification);
  badge.textContent = letter; badge.style.display = "";
  badge.title = m.classification || "";
  badge.className = "badge cls-" + classKey(m.classification);
  document.title = `CHESS-QC — ${m.name}`;

  // Scalar-only apps have no profile/grid outputs, so the plot/table area would just be an
  // empty placeholder. Hide it (and the plot-only export buttons) for those apps.
  const hasPlottable = contract.outputs.some((o) => o.kind === "profile" || o.kind === "grid");
  const tabs = document.querySelector(".tabs");
  if (tabs) tabs.style.display = hasPlottable ? "" : "none";
  if ($("png")) $("png").style.display = hasPlottable ? "" : "none";
  if ($("csv")) $("csv").style.display = hasPlottable ? "" : "none";
  if ($("popout")) $("popout").style.display = hasPlottable ? "" : "none";

  const box = $("inputs"); box.innerHTML = "";
  for (const fld of contract.inputs) {
    if (fld.kind === "table") {
      const block = document.createElement("div"); block.style.gridColumn = "1 / -1";
      block.dataset.fieldwrap = fld.key;
      const lab = document.createElement("label");
      lab.textContent = fld.label + ":"; lab.title = fld.note || "";
      lab.style.color = "var(--label)"; lab.style.display = "block"; lab.style.margin = "4px 0";
      block.append(lab, buildTable(fld));
      box.appendChild(block);
      continue;
    }
    if (fld.kind === "list" || fld.kind === "matrix") {
      // structured numeric input edited as raw JSON (no unit conversion); the default
      // is pre-loaded, and an empty box sends null so the app uses its built-in geometry.
      const block = document.createElement("div"); block.style.gridColumn = "1 / -1";
      block.dataset.fieldwrap = fld.key;
      const u = siUnit(fld.unit_si) || "";
      const lab = document.createElement("label");
      lab.textContent = `${fld.label}, JSON${u ? ` (${u})` : ""}:`; lab.title = fld.note || "";
      lab.style.color = "var(--label)"; lab.style.display = "block"; lab.style.margin = "4px 0";
      const ta = document.createElement("textarea");
      ta.dataset.key = fld.key; ta.dataset.struct = "1";
      ta.rows = fld.kind === "matrix" ? 6 : 2; ta.style.width = "100%";
      ta.style.fontFamily = "var(--mono, monospace)"; ta.style.fontSize = "11px";
      ta.value = fld.default == null ? "" : JSON.stringify(fld.default);
      block.append(lab, ta); box.appendChild(block);
      continue;
    }
    if (fld.kind === "csv") {
      // CSV record input: a bundled-station dropdown + a file upload. The loaded
      // text is held on the block (not in the DOM) and read back by gatherSI.
      const block = document.createElement("div"); block.style.gridColumn = "1 / -1";
      block.dataset.fieldwrap = fld.key;
      block.dataset.key = fld.key; block.dataset.csv = "1";
      block._csvText = fld.default == null ? "" : String(fld.default);
      const lab = document.createElement("label");
      lab.textContent = fld.label + ":"; lab.title = fld.note || "";
      lab.style.color = "var(--label)"; lab.style.display = "block"; lab.style.margin = "4px 0";
      const ctrls = document.createElement("div");
      ctrls.style.display = "flex"; ctrls.style.gap = "8px";
      ctrls.style.flexWrap = "wrap"; ctrls.style.alignItems = "center";
      const sel = document.createElement("select");
      const opt0 = document.createElement("option");
      opt0.value = ""; opt0.textContent = "— sample (built-in) —"; sel.appendChild(opt0);
      for (const c of (fld.choices || [])) {
        const [id, name] = String(c).split("|");
        const o = document.createElement("option");
        o.value = id; o.textContent = name ? `${id} — ${name}` : id; sel.appendChild(o);
      }
      const file = document.createElement("input");
      file.type = "file"; file.accept = ".parquet,.csv,text/csv";
      const status = document.createElement("span"); status.className = "unit";
      status.textContent = "built-in sample";
      const orlab = document.createElement("span");
      orlab.textContent = "or upload:"; orlab.style.color = "var(--muted, #888)";
      sel.addEventListener("change", async () => {
        file.value = "";
        if (!sel.value) {
          block._csvText = String(fld.default || ""); status.textContent = "built-in sample";
          clearOutputs(); return;      // different data: drop the previous result
        }
        status.textContent = "loading…";
        try {
          const dir = fld.data_dir || "water_levels";
          const resp = await fetch(`../../data/${dir}/${sel.value}.parquet`);
          if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
          const rec = await readParquetSeries(await resp.arrayBuffer());
          block._csvText = rec.csv;
          applyDatum(rec.datum);       // the station record carries its own datum
          status.textContent = `loaded ${sel.options[sel.selectedIndex].textContent}`
            + (rec.datum ? ` · ${rec.datum}` : "");
          clearOutputs();              // different data: drop the previous result
        } catch (e) { status.textContent = `load failed: ${e.message}`; }
      });
      file.addEventListener("change", async () => {
        const f = file.files && file.files[0]; if (!f) return;
        sel.value = ""; status.textContent = "loading…";
        try {
          const buf = await f.arrayBuffer();
          if (isParquet(buf)) {
            // a native Parquet record states its own datum; a CSV cannot, so there
            // the operator's selection is what stands
            const rec = await readParquetSeries(buf);
            block._csvText = rec.csv;
            applyDatum(rec.datum);
            status.textContent = `loaded ${f.name}` + (rec.datum ? ` · ${rec.datum}` : "")
              + (rec.datum ? "" : " · set the datum below");
          } else {
            block._csvText = new TextDecoder().decode(buf);
            status.textContent = `loaded ${f.name} · CSV: set the datum below`;
          }
          clearOutputs();
        } catch (e) { status.textContent = `read failed: ${e.message}`; }
      });
      ctrls.append(sel, orlab, file, status);
      block.append(lab, ctrls); box.appendChild(block);
      // Default to the built-in sample (no network); a station CSV is fetched only
      // when the user selects one. Nothing is computed until Compute is pressed.
      continue;
    }
    const row = document.createElement("div"); row.className = "row"; row.dataset.fieldwrap = fld.key;
    if (fld.note) { row.title = fld.note; row.classList.add("has-tip"); }   // hover anywhere in the row
    const lab = document.createElement("label"); lab.innerHTML = fmtLabel(fld.label);
    let ctl;
    if (fld.kind === "choice") {
      ctl = document.createElement("select");
      for (const c of fld.choices) { const o = document.createElement("option"); o.value = o.textContent = c; ctl.appendChild(o); }
      ctl.value = fld.default;
    } else if (fld.kind === "bool") {
      ctl = document.createElement("input"); ctl.type = "checkbox"; ctl.checked = !!fld.default;
    } else {
      ctl = document.createElement("input"); ctl.type = "number"; ctl.step = "any";
      const u = fieldUnit(fld);
      if (fld.lo != null) ctl.min = sig(fromSI(fld.lo, u));
      if (fld.hi != null) ctl.max = sig(fromSI(fld.hi, u));
      ctl.value = fmtIn(fromSI(Number(fld.default), u));
    }
    ctl.dataset.key = fld.key;
    const unit = document.createElement("span"); unit.className = "unit"; unit.dataset.key = fld.key;
    unit.innerHTML = fld.kind === "choice" || fld.kind === "bool" ? "" : fmtLabel(fieldUnit(fld));
    row.append(lab, ctl, unit);
    box.appendChild(row);
  }
  applyShowIf();
  applyEnableIf();
  buildWorkflowButtons();
  renderMethodPanel();
  applyHandoff();
}

// Gray out (disable) inputs whose `enable_if: [otherKey, value]` is not currently met
// — e.g. the Van der Meer parameters when the Hudson method is selected. Unlike show_if,
// the field stays visible so the user sees it exists but isn't used by this method.
function applyEnableIf() {
  for (const fld of contract.inputs) {
    const cond = fld.enable_if;
    if (!cond || !cond.length) continue;
    const ctrl = document.querySelector(`[data-key="${cond[0]}"]`);
    const cur = ctrl ? (ctrl.type === "checkbox" ? ctrl.checked : ctrl.value) : undefined;
    const on = String(cur) === String(cond[1]);
    const wrap = document.querySelector(`[data-fieldwrap="${fld.key}"]`);
    if (wrap) wrap.classList.toggle("disabled", !on);
    const own = document.querySelector(`[data-key="${fld.key}"]`);
    if (own && "disabled" in own) own.disabled = !on;
  }
}

// --- "Method & equations" panel (per-app ABOUT, KaTeX-typeset) -------------------
const _escHTML = (s) => String(s).replace(/[&<>]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));

function renderTeX(tex, el, display) {
  if (window.katex) {
    try { katex.render(tex, el, { displayMode: !!display, throwOnError: false }); return; }
    catch (e) { /* fall through to plain text */ }
  }
  el.textContent = tex;   // KaTeX not yet loaded / failed: show the raw source
}

// Greek-letter names that the ABOUT symbol shorthand spells out (lower-case keys).
const _GREEK = new Set(["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
  "iota", "kappa", "lambda", "mu", "nu", "xi", "pi", "rho", "sigma", "tau", "upsilon", "phi",
  "chi", "psi", "omega"]);

// sub/sup token grammar (shared by every symbol/unit/prose typesetter): a SUBSCRIPT may be a
// digit fraction (R_1/3), include % (R_2%), or be a word/*; a SUPERSCRIPT is an exponent only
// (so the unit "m^3/s" keeps the "/s" outside the exponent).
const _SUBC = "(?:\\([^)]+\\)|\\d+/\\d+|\\d+\\.\\d+|[A-Za-z0-9'%]+|\\*)";
const _SUPC = "(?:\\([^)]+\\)|\\d+\\.\\d+|[A-Za-z0-9']+|\\*)";
const _SSGRP = `(?:_${_SUBC}|\\^${_SUPC})`;
const _SYMTOK_RE = new RegExp(`^([A-Za-z]+)((?:${_SSGRP})*)$`);
const _ssGroups = (s) => [...s.matchAll(/([_^])(\([^)]+\)|\d+\/\d+|\d+\.\d+|[A-Za-z0-9'%]+|\*)/g)].map((m) => [m[1], m[2]]);

// Turn an ABOUT symbol shorthand (e.g. "U_obs", "u_*", "phi", "Delta T", "bar t", "R_1/3") into
// LaTeX so it typesets with real subscripts/Greek. Both front-ends share this convention.
function _symTok(t) {
  const m = _SYMTOK_RE.exec(t);
  if (!m) return t;                                  // not a clean symbol token
  let out = _GREEK.has(m[1].toLowerCase()) ? "\\" + m[1] : m[1];
  for (const [op, txt] of _ssGroups(m[2])) {
    const t = txt.replace(/%/g, "\\%");        // % is a comment char in (KaTeX/mathtext) LaTeX
    out += op + (t.length > 1 ? "{" + t + "}" : t);
  }
  return out;
}
function symToTex(s) {
  s = String(s).trim();
  if (s.includes(",")) return s.split(",").map((p) => symToTex(p)).join(",\\; ");
  const toks = s.split(/\s+/).filter(Boolean);
  const out = [];
  for (let i = 0; i < toks.length; i++) {
    if (toks[i] === "bar" && i + 1 < toks.length) out.push("\\bar{" + _symTok(toks[++i]) + "}");
    else out.push(_symTok(toks[i]));
  }
  return out.join(" ");
}

// Inline symbols inside prose (descriptions, notes): a token is "math" if it has an
// underscore subscript (z_obs, u_*) or is a spelled-out Greek letter. Typeset just those,
// leave the rest as text. Returns an HTML string (math via KaTeX renderToString).
const _MATHTOK = new RegExp(
  `([A-Za-z]+(?:${_SSGRP})+|(?<![A-Za-z\\\\])(?:` +
  [..._GREEK].join("|") + "|Delta|Theta|Lambda|Xi|Pi|Sigma|Phi|Psi|Omega|Gamma" +
  ")(?![A-Za-z]))", "g");

function renderProse(text) {
  if (!text) return "";
  const parts = String(text).split(_MATHTOK);
  let html = "";
  for (let i = 0; i < parts.length; i++) {
    if (i % 2 === 1 && window.katex) {
      try { html += katex.renderToString(symToTex(parts[i]), { throwOnError: false }); continue; }
      catch (e) { /* fall through */ }
    }
    html += _escHTML(parts[i]);
  }
  return html;
}

// Lightweight, NON-italic typesetting for labels/units/status (kN/m^3, ft^2, Hs/(d_n50),
// K_D=...): super/subscripts via <sub>/<sup>, Greek names to Unicode. Returns an HTML
// string. (The equation panel uses KaTeX italic math; labels/units read better upright.)
const _GREEK_UNI = {
  alpha: "α", beta: "β", gamma: "γ", delta: "δ", epsilon: "ε", zeta: "ζ", eta: "η",
  theta: "θ", iota: "ι", kappa: "κ", lambda: "λ", mu: "μ", nu: "ν", xi: "ξ", pi: "π",
  rho: "ρ", sigma: "σ", tau: "τ", upsilon: "υ", phi: "φ", chi: "χ", psi: "ψ", omega: "ω",
  Delta: "Δ", Theta: "Θ", Lambda: "Λ", Xi: "Ξ", Pi: "Π", Sigma: "Σ", Phi: "Φ", Psi: "Ψ",
  Omega: "Ω", Gamma: "Γ",
};
const _greekUni = (w) => _GREEK_UNI[w] || _GREEK_UNI[w.toLowerCase()];
function _tokHTML(tok) {
  const m = _SYMTOK_RE.exec(tok);
  if (!m) return _greekUni(tok) || _escHTML(tok);
  let out = _greekUni(m[1]) || _escHTML(m[1]);
  for (const [op, txt] of _ssGroups(m[2])) { const tag = op === "_" ? "sub" : "sup"; out += `<${tag}>${_escHTML(txt)}</${tag}>`; }
  return out;
}
function fmtLabel(text) {
  if (text == null || text === "") return "";
  return String(text).split(_MATHTOK).map((p, i) => (i % 2 ? _tokHTML(p) : _escHTML(p))).join("");
}

function renderMethodPanel() {
  const box = $("methodBox"), body = $("methodBody");
  if (!box || !body) return;
  const about = contract.about;
  if (!about) { box.style.display = "none"; return; }
  box.style.display = "";
  body.innerHTML = "";

  if (about.summary) {
    const p = document.createElement("p"); p.className = "method-summary";
    p.innerHTML = renderProse(about.summary); body.appendChild(p);
  }
  for (const m of (about.methods || [])) {
    const sec = document.createElement("section");
    sec.className = "method"; sec.dataset.when = m.when || "";
    const h = document.createElement("div"); h.className = "method-head";
    const nm = document.createElement("span"); nm.className = "method-name"; nm.textContent = m.name || "";
    h.appendChild(nm);
    if (m.tag) { const t = document.createElement("span"); t.className = "chip chip-" + m.tag; t.textContent = m.tag; h.appendChild(t); }
    sec.appendChild(h);
    if (m.note) { const n = document.createElement("p"); n.className = "method-note"; n.innerHTML = renderProse(m.note); sec.appendChild(n); }
    for (const eq of (m.equations || [])) {
      const wrap = document.createElement("div"); wrap.className = "eq";
      const e = document.createElement("div"); e.className = "eq-tex"; renderTeX(eq.tex, e, true); wrap.appendChild(e);
      if (eq.desc) { const d = document.createElement("span"); d.className = "eq-desc"; d.innerHTML = renderProse(eq.desc); wrap.appendChild(d); }
      sec.appendChild(wrap);
    }
    body.appendChild(sec);
  }
  if (about.symbols && about.symbols.length) {
    const dl = document.createElement("dl"); dl.className = "symbols";
    for (const [sym, desc] of about.symbols) {
      const dt = document.createElement("dt"); renderTeX(symToTex(sym), dt, false);
      const dd = document.createElement("dd"); dd.innerHTML = renderProse(desc);
      dl.append(dt, dd);
    }
    const cap = document.createElement("div"); cap.className = "method-subhead"; cap.textContent = "Symbols";
    body.append(cap, dl);
  }
  if (about.references && about.references.length) {
    const r = document.createElement("p"); r.className = "method-refs";
    r.innerHTML = "<span class='method-subhead'>References</span> " + about.references.map(_escHTML).join(" · ");
    body.appendChild(r);
  }
  highlightActiveMethod();
}

// Dim the methods that are not the currently-selected formulation (method_key choice).
function highlightActiveMethod() {
  const about = contract.about; if (!about) return;
  const tagEl = $("methodTag"); if (tagEl) tagEl.textContent = "";
  const key = about.method_key;
  const single = !key || (about.methods || []).length < 2;
  let active = null;
  if (!single) {
    const ctrl = document.querySelector(`[data-key="${key}"]`);
    active = ctrl ? ctrl.value : null;
  }
  for (const sec of document.querySelectorAll("#methodBody .method")) {
    const on = single || sec.dataset.when === String(active);
    sec.classList.toggle("inactive", !on);
    if (on && !single) {
      const nm = sec.querySelector(".method-name");
      if (nm && tagEl) tagEl.textContent = "· " + nm.textContent;
    }
  }
}

// --- area-10 workflow: "send to next app" buttons + in-memory series hand-off ---
// The web app-switcher (loadApp) does not reload the page, so a JS variable carries
// the full-resolution series between steps without storage limits.
function buildWorkflowButtons() {
  const wf = $("workflow");
  if (!wf) return;
  wf.innerHTML = "";
  const nexts = contract.meta.next_apps || [];
  for (const [tid, tlabel] of nexts) {
    const btn = document.createElement("button");
    btn.className = "wf-btn";
    btn.textContent = `Send to ${tid} ${tlabel} →`;
    btn.addEventListener("click", () => nextStep(tid));
    wf.appendChild(btn);
  }
  wf.style.display = nexts.length ? "" : "none";
}

function nextStep(targetId) {
  // re-run the current app asking for its full-resolution series, then carry it forward
  let res;
  try {
    const inp = gatherSI(); inp.handoff = "1";
    res = JSON.parse(bridge.run(appmod, JSON.stringify(inp)));
  } catch (e) { setStatus(String(e), false); return; }
  if (res._error) { setStatus(res._error, false); return; }
  window.CHESSQC_HANDOFF = { csv: res.handoff_csv || "",
                            label: `from ${contract.meta.aces_id} ${contract.meta.name}` };
  history.pushState({ app: targetId }, "", `?app=${encodeURIComponent(targetId)}`);
  loadApp(targetId);
}

function applyHandoff() {
  const h = window.CHESSQC_HANDOFF;
  if (!h || !h.csv) return;
  const block = document.querySelector('[data-csv="1"]');   // first CSV field of the target
  if (block) {
    block._csvText = h.csv;
    const st = block.querySelector("span.unit");
    if (st) st.textContent = `← ${h.label}`;
  }
  window.CHESSQC_HANDOFF = null;
}

// Hide inputs whose `show_if: [otherKey, value]` condition is not currently met
// (e.g. "Specified slope" only shows when Slope source = Specified slope).
function applyShowIf() {
  for (const fld of contract.inputs) {
    const cond = fld.show_if;
    if (!cond || !cond.length) continue;
    const ctrl = document.querySelector(`[data-key="${cond[0]}"]`);
    const cur = ctrl ? (ctrl.type === "checkbox" ? ctrl.checked : ctrl.value) : undefined;
    const wrap = document.querySelector(`[data-fieldwrap="${fld.key}"]`);
    if (wrap) wrap.style.display = String(cur) === String(cond[1]) ? "" : "none";
  }
}

// --- table-input field (kind === "table") ---
const tableCols = (fld) => (fld.columns && fld.columns.length)
  ? fld.columns.map((c) => ({ label: c[0], si: c[1], us: c[2] }))
  : [{ label: fld.label, si: fld.unit_si, us: fld.unit_us }];
const colUnit = (c) => (system === "SI" ? siUnit(c.si) : c.us);

function addTableRow(wrap, rowSI) {
  const tb = wrap.querySelector("tbody"), tr = document.createElement("tr");
  wrap._cols.forEach((c, i) => {
    const td = document.createElement("td"), inp = document.createElement("input");
    inp.type = "number"; inp.step = "any"; inp.style.width = "84px";
    if (rowSI && rowSI[i] != null) inp.value = fmtIn(fromSI(Number(rowSI[i]), colUnit(c)));
    td.appendChild(inp); tr.appendChild(td);
  });
  tb.appendChild(tr);
}

function buildTable(fld) {
  const cols = tableCols(fld);
  const wrap = document.createElement("div"); wrap.dataset.tableKey = fld.key; wrap._cols = cols;
  wrap.className = "tablewrap";
  const t = document.createElement("table"); t.className = "intable";
  const htr = document.createElement("tr");
  for (const c of cols) {
    const th = document.createElement("th"); const u = colUnit(c);
    th.innerHTML = u ? `${fmtLabel(c.label)} (${fmtLabel(u)})` : fmtLabel(c.label); htr.appendChild(th);
  }
  const thead = document.createElement("thead"); thead.appendChild(htr);
  t.append(thead, document.createElement("tbody")); wrap.appendChild(t);
  for (const r of (fld.default || [])) addTableRow(wrap, r);
  const bar = document.createElement("div"); bar.style.marginTop = "4px";
  const add = document.createElement("button"); add.type = "button"; add.textContent = "+ row";
  const rem = document.createElement("button"); rem.type = "button"; rem.textContent = "− row";
  add.addEventListener("click", () => addTableRow(wrap, null));
  rem.addEventListener("click", () => { const b = t.querySelector("tbody"); if (b.rows.length) b.deleteRow(b.rows.length - 1); });
  bar.append(add, rem); wrap.appendChild(bar);
  return wrap;
}

function gatherTable(fld) {
  const wrap = document.querySelector(`[data-table-key="${fld.key}"]`), rows = [];
  for (const tr of wrap.querySelectorAll("tbody tr")) {
    const vals = []; let blank = true;
    [...tr.querySelectorAll("input")].forEach((inp, i) => {
      const txt = inp.value.trim();
      if (txt === "") vals.push(0.0);
      else { blank = false; vals.push(toSI(parseFloat(txt), colUnit(wrap._cols[i]))); }
    });
    if (!blank) rows.push(vals);
  }
  return rows;
}

function gatherSI() {
  const inp = {};
  for (const fld of contract.inputs) {
    if (fld.kind === "table") { inp[fld.key] = gatherTable(fld); continue; }
    if (fld.kind === "csv") {
      const el = document.querySelector(`[data-key="${fld.key}"][data-csv="1"]`);
      inp[fld.key] = el && el._csvText != null ? el._csvText : "";
      continue;
    }
    if (fld.kind === "list" || fld.kind === "matrix") {
      const ta = document.querySelector(`[data-key="${fld.key}"]`);
      const txt = ta.value.trim();
      if (txt === "") { inp[fld.key] = null; continue; }
      try { inp[fld.key] = JSON.parse(txt); }
      catch (e) { throw new Error(`${fld.label}: invalid JSON, ${e.message}`); }
      continue;
    }
    const ctl = document.querySelector(`[data-key="${fld.key}"]`);
    if (fld.kind === "choice") inp[fld.key] = ctl.value;
    else if (fld.kind === "bool") inp[fld.key] = ctl.checked;
    else inp[fld.key] = toSI(parseFloat(ctl.value), fieldUnit(fld));
  }
  return inp;
}

function doCompute() {
  let res;
  try { res = JSON.parse(bridge.run(appmod, JSON.stringify(gatherSI()))); }
  catch (e) { setStatus(String(e), false); return; }
  if (res._error) { setStatus(res._error, false); return; }
  lastRes = res;
  render(res);
}

function render(res) {
  // scalar/point value rows
  const vbox = $("values"); vbox.innerHTML = "";
  for (const o of contract.outputs) {
    if (["profile", "grid", "vline", "scatter", "scatter_x", "data"].includes(o.kind) || !(o.key in res)) continue;
    const u = outUnit(o), raw = res[o.key];
    // numeric -> convert+format; string (a label like "Beta-Rayleigh", or an inf/nan
    // sentinel from the bridge) -> show as-is with the unit appended when present
    const numTxt = typeof raw === "number" ? fmtNum(fromSI(raw, u)) : fmtToken(raw);
    const row = document.createElement("div"); row.className = "vrow";
    const n = document.createElement("span"); n.innerHTML = fmtLabel(o.label);
    if (o.note) { n.title = o.note; n.classList.add("has-tip"); }   // hover definition
    const v = document.createElement("span"); v.className = "val";
    v.innerHTML = _escHTML(numTxt) + (u ? " " + fmtLabel(u) : "");
    row.append(n, v); vbox.appendChild(row);
  }
  drawPlot(res); fillTable(res);
  setStatus(`units: ${system} · valid · ${res._notes || ""} · ${contract.meta.cite}`, true);
}

// --- profile discovery (shared by the canvas plot and the data table) ---
// Reads the output descriptors so any app is self-describing: x = profile_X; the
// remaining profile series are grouped by display unit (first series alone for scale,
// the rest grouped) exactly as the desktop renderer does.
const _shortLabel = (o) => o.label.split(":").slice(-1)[0].trim();
function profileSeries(res) {
  const profs = contract.outputs.filter((o) => o.kind === "profile" && o.key in res);
  const xo = profs.find((o) => o.key === "profile_X") || profs[0];
  const x = { label: _shortLabel(xo), unit: outUnit(xo),
              data: res[xo.key].map((v) => fromSI(v, outUnit(xo))) };
  const colors = [cssVar("--plot-eta", "#008cb0"), cssVar("--plot-u", "#008cb0"),
                  cssVar("--plot-w", "#a86f2f"), cssVar("--plot-fg", "#333"),
                  cssVar("--plot-text", "#777")];
  // within-panel color cycle (distinct adjacent colors for overlays like data + trend)
  const gcolors = [cssVar("--plot-eta", "#008cb0"), cssVar("--plot-w", "#a86f2f"),
                   cssVar("--plot-fg", "#333"), cssVar("--plot-text", "#777")];
  const ys = profs.filter((o) => o !== xo);
  const groups = [];
  if (ys.some((o) => o.group)) {
    // explicit grouping: profiles sharing a `group` plot on the same panel
    const by = {};
    for (const o of ys) {
      const gk = o.group || o.key, u = outUnit(o);
      if (!(gk in by)) { by[gk] = { unit: u, gid: gk, series: [] }; groups.push(by[gk]); }
      const idx = by[gk].series.length;
      by[gk].series.push({ name: _shortLabel(o), data: res[o.key].map((v) => fromSI(v, u)),
                           color: gcolors[idx % gcolors.length] });
    }
  } else {
    ys.forEach((o, i) => {
      const u = outUnit(o);
      const s = { name: _shortLabel(o), data: res[o.key].map((v) => fromSI(v, u)),
                  color: colors[i % colors.length] };
      const g = (i ? groups.slice(1) : []).find((g) => g.unit === u);
      if (g) g.series.push(s); else groups.push({ unit: u, series: [s] });
    });
  }
  return { x, groups };
}
// Axis tick label. When the tick step is known the decimal count follows the step
// rather than the magnitude, so adjacent ticks can never collapse to the same text
// (a 30-day record on a decimal-year axis reads 2000.000, 2000.016, ... instead of
// "2000.0" six times over, and a near-flat series still shows its true range).
const fmtTick = (v, step) => {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (step > 0 && a < 1e5) {
    return v.toFixed(Math.max(0, Math.min(10, Math.ceil(-Math.log10(step)) + 1)));
  }
  if (a >= 1e5 || a < 1e-3) return v.toExponential(1);
  // keep the integer part intact for large values (e.g. years: 1992.5, not 1990)
  if (a >= 100) return String(Math.round(v * 10) / 10);
  return String(+v.toPrecision(4));
};

// profile series with at most this many samples are drawn with point markers
const _MARK_MAX = 20;

// Common tick offset, as matplotlib's axes use on the desktop side: when the span is
// tiny next to the values themselves (a 30-day record on a decimal-year axis, say),
// factor a shared base out of the labels and annotate it once on the axis. Returns 0
// when the plain labels already read well.
const tickBase = (lo, hi) => {
  const span = hi - lo, mag = Math.max(Math.abs(lo), Math.abs(hi));
  if (!(span > 0) || !(mag / span >= 1000)) return 0;
  const decade = Math.pow(10, Math.floor(Math.log10(span)));
  return Math.floor(lo / decade) * decade;
};
// "+2000" / "−0.0031" annotation for a factored-out tick base (exact, no padding zeros)
const fmtBase = (b) => (b > 0 ? "+" : "−") + String(+Math.abs(b).toPrecision(12));

const cssVar = (name, fallback) => {
  const v = getComputedStyle(document.body).getPropertyValue(name).trim();
  return v || fallback;
};

function drawPlot(res) { drawPlotInto(res, $("plot")); }

// X-axis zoom (wheel) + pan (drag) + reset (double-click); Y auto-scales to the
// visible window. Attached once per canvas; no-op on heatmaps (no _plotBox).
function attachPlotZoom(cv) {
  if (cv._zoomAttached) return;
  cv._zoomAttached = true;
  const redraw = () => { if (cv._res) drawPlotInto(cv._res, cv); };
  cv.addEventListener("wheel", (e) => {
    if (!cv._plotBox) return;
    e.preventDefault();
    const { L, R, xmin, xmax } = cv._plotBox, [fx0, fx1] = cv._xfull;
    const rect = cv.getBoundingClientRect();
    const frac = Math.min(1, Math.max(0, ((e.clientX - rect.left) - L) / (R - L)));
    const xc = xmin + frac * (xmax - xmin), f = e.deltaY < 0 ? 0.8 : 1.25;
    let a = xc - (xc - xmin) * f, b = xc + (xmax - xc) * f;
    a = Math.max(fx0, a); b = Math.min(fx1, b);
    if (b - a < (fx1 - fx0) * 1e-4) return;
    cv._viewX = (a <= fx0 && b >= fx1) ? null : [a, b];
    redraw();
  }, { passive: false });
  let drag = false, lastX = 0;
  cv.addEventListener("mousedown", (e) => { if (cv._plotBox) { drag = true; lastX = e.clientX; cv.style.cursor = "grabbing"; } });
  window.addEventListener("mousemove", (e) => {
    if (!drag || !cv._plotBox) return;
    const { L, R, xmin, xmax } = cv._plotBox, [fx0, fx1] = cv._xfull, span = xmax - xmin;
    const dd = -(e.clientX - lastX) / (R - L) * span; lastX = e.clientX;
    let a = xmin + dd, b = xmax + dd;
    if (a < fx0) { a = fx0; b = fx0 + span; }
    if (b > fx1) { b = fx1; a = fx1 - span; }
    cv._viewX = [a, b]; redraw();
  });
  window.addEventListener("mouseup", () => { if (drag) { drag = false; cv.style.cursor = ""; } });
  cv.addEventListener("dblclick", () => { cv._viewX = null; redraw(); });
}

function drawPlotInto(res, cv) {
  const ctx = cv.getContext("2d");
  // high-DPI: back the canvas at devicePixelRatio for crisp lines/text
  const dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth || 520, H = Math.max(cv.clientHeight || 280, 280);
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  cv._res = res;
  const C = {
    bg: cssVar("--plot-bg", "#ffffff"), fg: cssVar("--plot-fg", "#33475b"),
    grid: cssVar("--plot-grid", "#eceef1"), axis: cssVar("--plot-axis", "#8a93a0"),
    text: cssVar("--plot-text", "#6b7785"),
  };
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = C.bg; ctx.fillRect(0, 0, W, H);
  const gridOuts = contract.outputs.filter((o) => o.kind === "grid" && o.key in res);
  if (gridOuts.length) { cv._plotBox = null; drawHeatmap(ctx, W, H, C, res, gridOuts[0]); return; }
  const haveProfiles = contract.outputs.some((o) => o.kind === "profile" && o.key in res);
  if (!haveProfiles) { cv._plotBox = null; return; }
  const { x, groups } = profileSeries(res);
  const X = x.data, Xf = X.filter(Number.isFinite);
  const fullMin = Math.min(...Xf), fullMax = Math.max(...Xf);
  cv._xfull = [fullMin, fullMax];
  const [xmin, xmax] = cv._viewX && cv._viewX[1] > cv._viewX[0] ? cv._viewX : [fullMin, fullMax];
  // vertical markers (kind "vline") on the x axis, e.g. the NTDE midpoint; "nan" -> omitted
  const vlines = contract.outputs
    .filter((o) => o.kind === "vline" && o.key in res)
    .map((o) => ({ x: fromSI(Number(res[o.key]), x.unit), label: _shortLabel(o) }))
    .filter((v) => Number.isFinite(v.x));
  // scatter overlays (kind "scatter") draw markers at their own (x_key, key); grouped onto a panel by `group`
  const pkColor = cssVar("--plot-peaks", "#d62728");
  const scatters = contract.outputs
    .filter((o) => o.kind === "scatter" && o.key in res && o.x_key && o.x_key in res)
    .map((o) => ({ gid: o.group,
                   xs: res[o.x_key].map((v) => fromSI(Number(v), x.unit)),
                   ys: res[o.key].map((v) => fromSI(Number(v), outUnit(o))) }));

  const panel = (y0, y1, series, label, gid) => {
    const L = 48, R = W - 12, T = y0 + 18, B = y1 - 26;       // plot box
    let lo = Infinity, hi = -Infinity;
    for (const s of series) X.forEach((xx, i) => {
      if (xx < xmin || xx > xmax) return;
      const v = s.data[i];
      if (!Number.isFinite(v)) return;
      if (v < lo) lo = v; if (v > hi) hi = v;
    });
    if (!isFinite(lo) || !isFinite(hi)) { lo = 0; hi = 1; }
    if (lo === hi) { lo -= 1; hi += 1; }
    const padv = (hi - lo) * 0.08 || 1; lo -= padv; hi += padv;
    const sx = (xx) => L + (xx - xmin) / (xmax - xmin) * (R - L);
    const sy = (v) => B - (v - lo) / (hi - lo) * (B - T);

    ctx.font = "10px sans-serif"; ctx.lineWidth = 1;
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    const ystep = (hi - lo) / 4, xstep = (xmax - xmin) / 5;
    const ybase = tickBase(lo, hi), xbase = tickBase(xmin, xmax);
    for (let i = 0; i <= 4; i++) {
      const val = lo + ystep * i, y = sy(val);
      ctx.strokeStyle = C.grid; ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(R, y); ctx.stroke();
      ctx.fillStyle = C.text; ctx.fillText(fmtTick(val - ybase, ystep), L - 5, y);
    }
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    for (let i = 0; i <= 5; i++) {
      const val = xmin + xstep * i, xp = sx(val);
      ctx.strokeStyle = C.grid; ctx.beginPath(); ctx.moveTo(xp, T); ctx.lineTo(xp, B); ctx.stroke();
      ctx.fillStyle = C.text; ctx.fillText(fmtTick(val - xbase, xstep), xp, B + 4);
    }
    // the factored-out y base, annotated above the axis (as matplotlib does); the x
    // base is shared by every panel and goes on the x-axis label below the plot
    if (ybase) {
      ctx.fillStyle = C.text; ctx.font = "9px sans-serif";
      ctx.textAlign = "left"; ctx.textBaseline = "bottom";
      ctx.fillText(fmtBase(ybase), L - 44, T - 2);
      ctx.font = "10px sans-serif";
    }
    ctx.strokeStyle = C.axis; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(L, T); ctx.lineTo(L, B); ctx.lineTo(R, B); ctx.stroke();
    if (lo < 0 && hi > 0) {
      ctx.strokeStyle = C.grid; const yz = sy(0);
      ctx.beginPath(); ctx.moveTo(L, yz); ctx.lineTo(R, yz); ctx.stroke();
    }
    // data: clipped to the box, with the line broken across NaN gaps. Series of only
    // a few points (a modal spectrum, a directional band table) get their samples
    // marked, so a discrete series does not read as a continuous curve.
    const marks = X.length <= _MARK_MAX;
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, R - L, B - T); ctx.clip();
    for (const s of series) {
      ctx.strokeStyle = s.color; ctx.lineWidth = 1.8; ctx.beginPath();
      let pen = false;
      X.forEach((xx, i) => {
        const v = s.data[i];
        if (!Number.isFinite(v)) { pen = false; return; }
        const px = sx(xx), py = sy(v);
        if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; }
      });
      ctx.stroke();
      if (marks) {
        ctx.fillStyle = s.color;
        X.forEach((xx, i) => {
          const v = s.data[i];
          if (!Number.isFinite(v)) return;
          ctx.beginPath(); ctx.arc(sx(xx), sy(v), 2.2, 0, 2 * Math.PI); ctx.fill();
        });
      }
    }
    // vertical markers (e.g. NTDE midpoint), dashed, clipped to the box
    if (vlines.length) {
      ctx.strokeStyle = C.fg; ctx.lineWidth = 1.1; ctx.setLineDash([5, 4]);
      ctx.font = "10px sans-serif"; ctx.textBaseline = "top"; ctx.textAlign = "left";
      for (const v of vlines) {
        if (v.x < xmin || v.x > xmax) continue;
        const px = sx(v.x);
        ctx.beginPath(); ctx.moveTo(px, T); ctx.lineTo(px, B); ctx.stroke();
        if (y0 < 1) { ctx.fillStyle = C.fg; ctx.fillText(`${v.label} ${fmtTick(v.x, xstep)}`, px + 4, T + 2); }
      }
      ctx.setLineDash([]);
    }
    // scatter markers (e.g. POT peaks) for this panel's group
    for (const sc of scatters) {
      if (sc.gid !== gid) continue;
      ctx.fillStyle = pkColor;
      for (let i = 0; i < sc.xs.length; i++) {
        const xx = sc.xs[i], v = sc.ys[i];
        if (!Number.isFinite(xx) || !Number.isFinite(v) || xx < xmin || xx > xmax) continue;
        ctx.beginPath(); ctx.arc(sx(xx), sy(v), 2.6, 0, 2 * Math.PI); ctx.fill();
      }
    }
    ctx.restore();
    ctx.fillStyle = C.fg; ctx.textAlign = "left"; ctx.textBaseline = "top"; ctx.font = "11px sans-serif";
    ctx.fillText(label, L + 2, y0 + 2);
    if (series.length > 1) {
      ctx.font = "10px sans-serif";
      const items = series.map((s) => ({ s, w: ctx.measureText(s.name).width + 22 }));
      let lx = R - items.reduce((a, it) => a + it.w, 0);
      for (const { s, w } of items) {
        ctx.fillStyle = s.color; ctx.fillRect(lx, y0 + 7, 12, 3);
        ctx.fillStyle = C.fg; ctx.textBaseline = "top"; ctx.fillText(s.name, lx + 16, y0 + 2); lx += w;
      }
    }
    return { L, R };
  };

  const n = groups.length;
  let box = null;
  groups.forEach((g, i) => {
    const label = g.series.length === 1 ? `${g.series[0].name} (${g.unit})` : `(${g.unit})`;
    box = panel((H * i) / n, (H * (i + 1)) / n, g.series, label, g.gid);
  });
  // x-axis label, carrying the factored-out tick base when there is one
  const xb = tickBase(xmin, xmax);
  ctx.fillStyle = C.text; ctx.textAlign = "right"; ctx.textBaseline = "bottom"; ctx.font = "10px sans-serif";
  ctx.fillText(`${x.label} (${x.unit})${xb ? ` ${fmtBase(xb)}` : ""}`, W - 12, H - 2);
  cv._plotBox = box ? { L: box.L, R: box.R, xmin, xmax } : null;
  attachPlotZoom(cv);
}

// --- 2-D field heatmap (kind === "grid"; e.g. 3-4 wedge grid) ---
const _heat = (t) => {
  t = Math.max(0, Math.min(1, t));                      // blue -> cyan -> yellow -> red
  const r = Math.round(255 * Math.max(0, Math.min(1, 1.5 - Math.abs(4 * t - 3))));
  const g = Math.round(255 * Math.max(0, Math.min(1, 1.5 - Math.abs(4 * t - 2))));
  const b = Math.round(255 * Math.max(0, Math.min(1, 1.5 - Math.abs(4 * t - 1))));
  return `rgb(${r},${g},${b})`;
};

function drawHeatmap(ctx, W, H, C, res, zout) {
  const profs = contract.outputs.filter((o) => o.kind === "profile" && o.key in res);
  const xo = profs.find((o) => o.key === "grid_x") || profs[0];
  const yo = profs.find((o) => o.key === "grid_y") || profs[1] || profs[0];
  if (!xo || !yo) return;
  const xs = res[xo.key].map((v) => fromSI(v, outUnit(xo)));
  const ys = res[yo.key].map((v) => fromSI(v, outUnit(yo)));
  const zu = outUnit(zout);
  const Z = res[zout.key].map((row) => row.map((v) => fromSI(v, zu)));
  const ny = Z.length, nx = (Z[0] || []).length;
  if (!nx || !ny) return;
  let lo = Infinity, hi = -Infinity;
  for (const r of Z) for (const v of r) { if (v < lo) lo = v; if (v > hi) hi = v; }
  if (lo === hi) { lo -= 1; hi += 1; }
  const L = 52, R = W - 64, T = 22, B = H - 30;
  const sx = (i) => L + (nx <= 1 ? 0.5 : i / (nx - 1)) * (R - L);
  const sy = (j) => B - (ny <= 1 ? 0.5 : j / (ny - 1)) * (B - T);
  const cw = (R - L) / Math.max(nx - 1, 1), ch = (B - T) / Math.max(ny - 1, 1);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    ctx.fillStyle = _heat((Z[j][i] - lo) / (hi - lo));
    ctx.fillRect(sx(i) - cw / 2, sy(j) - ch / 2, cw + 1, ch + 1);
  }
  // frame + axis ticks
  ctx.strokeStyle = C.axis; ctx.lineWidth = 1;
  ctx.strokeRect(L, T, R - L, B - T);
  ctx.font = "10px sans-serif"; ctx.fillStyle = C.text;
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  const xstep = (xs[xs.length - 1] - xs[0]) / 4, ystep = (ys[ys.length - 1] - ys[0]) / 4;
  for (let i = 0; i <= 4; i++) {
    ctx.fillText(fmtTick(xs[0] + xstep * i, xstep), L + (R - L) * i / 4, B + 4);
  }
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (let j = 0; j <= 4; j++) {
    ctx.fillText(fmtTick(ys[0] + ystep * j, ystep), L - 5, B - (B - T) * j / 4);
  }
  // colorbar
  const cbx = R + 14, cbw = 12;
  for (let p = 0; p <= 50; p++) {
    ctx.fillStyle = _heat(p / 50);
    ctx.fillRect(cbx, B - (B - T) * p / 50, cbw, (B - T) / 50 + 1);
  }
  ctx.strokeStyle = C.axis; ctx.strokeRect(cbx, T, cbw, B - T);
  ctx.fillStyle = C.text; ctx.textAlign = "left"; ctx.textBaseline = "middle"; ctx.font = "9px sans-serif";
  ctx.fillText(fmtTick(hi, (hi - lo) / 4), cbx + cbw + 2, T);
  ctx.fillText(fmtTick(lo, (hi - lo) / 4), cbx + cbw + 2, B);
  // titles
  ctx.fillStyle = C.fg; ctx.textAlign = "left"; ctx.textBaseline = "top"; ctx.font = "11px sans-serif";
  ctx.fillText(`${_shortLabel(zout)}${zu ? ` (${zu})` : ""}`, L + 2, 4);
  ctx.fillStyle = C.text; ctx.textAlign = "right"; ctx.textBaseline = "bottom"; ctx.font = "10px sans-serif";
  ctx.fillText(`${_shortLabel(xo)} (${outUnit(xo)})  vs  ${_shortLabel(yo)} (${outUnit(yo)})`, W - 12, H - 2);
}

function fillTable(res) {
  const gridOuts = contract.outputs.filter((o) => o.kind === "grid" && o.key in res);
  if (gridOuts.length) { fillGridTable(res, gridOuts[0]); return; }
  const haveProfiles = contract.outputs.some((o) => o.kind === "profile" && o.key in res);
  if (!haveProfiles) { $("tablePane").innerHTML = ""; return; }
  const { x, groups } = profileSeries(res);
  const cols = [[`${x.label} (${x.unit})`, x.data]];
  for (const g of groups) for (const s of g.series) cols.push([`${s.name} (${g.unit})`, s.data]);
  let html = "<table><thead><tr>" + cols.map(([h]) => `<th>${h}</th>`).join("") + "</tr></thead><tbody>";
  for (let i = 0; i < x.data.length; i++)
    html += "<tr>" + cols.map(([, a]) => `<td>${fmtNum(a[i])}</td>`).join("") + "</tr>";
  $("tablePane").innerHTML = html + "</tbody></table>";
}

function fillGridTable(res, zout) {
  const profs = contract.outputs.filter((o) => o.kind === "profile" && o.key in res);
  const xo = profs.find((o) => o.key === "grid_x") || profs[0];
  const yo = profs.find((o) => o.key === "grid_y") || profs[1] || profs[0];
  if (!xo || !yo) { $("tablePane").innerHTML = ""; return; }
  const xs = res[xo.key].map((v) => fromSI(v, outUnit(xo)));
  const ys = res[yo.key].map((v) => fromSI(v, outUnit(yo)));
  const zu = outUnit(zout);
  const Z = res[zout.key].map((row) => row.map((v) => fromSI(v, zu)));
  let html = `<table><thead><tr><th>${_shortLabel(yo)}\\${_shortLabel(xo)}</th>`;
  html += xs.map((v) => `<th>${fmtNum(v)}</th>`).join("") + "</tr></thead><tbody>";
  for (let j = 0; j < ys.length; j++) {
    html += `<tr><th>${fmtNum(ys[j])}</th>`;
    html += (Z[j] || []).map((v) => `<td>${fmtNum(v)}</td>`).join("") + "</tr>";
  }
  $("tablePane").innerHTML = html + "</tbody></table>";
}

function tableCSV(sep) {
  const t = $("tablePane").querySelector("table"); if (!t) return "";
  return [...t.querySelectorAll("tr")].map((tr) =>
    [...tr.children].map((c) => c.textContent).join(sep)).join("\n");
}

// --- plain-text summary report (inputs + outputs + notes), current units ---
function reportText() {
  const m = contract.meta;
  let s = `CHESS-QC ${m.aces_id}  ${m.name}  [${classLetter(m.classification)} ${m.classification}]\n${m.cite}\n`;
  s += `Functional area: ${m.area}\nUnits: ${system}\n\nInputs:\n`;
  for (const fld of contract.inputs) {
    let val;
    if (fld.kind === "table") { val = `${gatherTable(fld).length} rows`; s += `  ${fld.label}: ${val}\n`; continue; }
    if (fld.kind === "csv") {
      const el = document.querySelector(`[data-key="${fld.key}"][data-csv="1"]`);
      const nrows = el && el._csvText ? el._csvText.split("\n").length : 0;
      s += `  ${fld.label}: ${nrows} lines loaded\n`; continue;
    }
    if (fld.kind === "list" || fld.kind === "matrix") {
      const ta = document.querySelector(`[data-key="${fld.key}"]`);
      val = ta.value.trim() === "" ? "(app default)" : "(custom JSON)";
      s += `  ${fld.label}: ${val}\n`; continue;
    }
    const ctl = document.querySelector(`[data-key="${fld.key}"]`);
    if (fld.kind === "choice") val = ctl.value;
    else if (fld.kind === "bool") val = ctl.checked ? "yes" : "no";
    else val = `${ctl.value} ${fieldUnit(fld)}`.trim();
    s += `  ${fld.label}: ${val}\n`;
  }
  s += "\nOutputs:\n";
  for (const o of contract.outputs) {
    if (["profile", "grid", "vline", "scatter", "scatter_x", "data"].includes(o.kind) || !lastRes || !(o.key in lastRes)) continue;
    const u = outUnit(o), raw = lastRes[o.key];
    const txt = typeof raw === "number"
      ? `${fmtNum(fromSI(raw, u))} ${u}`.trim()
      : `${fmtToken(raw)}${u ? " " + u : ""}`;
    s += `  ${o.label}: ${txt}\n`;
  }
  if (lastRes && lastRes.notes) s += `\nNotes: ${lastRes.notes}\n`;
  return s;
}

// --- unit toggle ---
function onUnits(newSys) {
  if (newSys === system) return;
  for (const fld of contract.inputs) {
    if (["choice", "bool", "table", "list", "matrix", "csv"].includes(fld.kind)) continue;
    const ctl = document.querySelector(`[data-key="${fld.key}"]`);
    const oldU = fieldUnit(fld);
    const newU = newSys === "SI" ? siUnit(fld.unit_si) : fld.unit_us;
    const si = toSI(parseFloat(ctl.value), oldU);
    if (fld.lo != null) ctl.min = sig(fromSI(fld.lo, newU)); else ctl.removeAttribute("min");
    if (fld.hi != null) ctl.max = sig(fromSI(fld.hi, newU)); else ctl.removeAttribute("max");
    ctl.value = fmtIn(fromSI(si, newU));
    document.querySelector(`span.unit[data-key="${fld.key}"]`).innerHTML = fmtLabel(newU);
  }
  for (const fld of contract.inputs) {           // tables: reconvert cells + headers
    if (fld.kind !== "table") continue;
    const wrap = document.querySelector(`[data-table-key="${fld.key}"]`), cols = wrap._cols;
    for (const tr of wrap.querySelectorAll("tbody tr"))
      [...tr.querySelectorAll("input")].forEach((inp, i) => {
        if (inp.value.trim() === "") return;
        const oldU = system === "SI" ? siUnit(cols[i].si) : cols[i].us;
        const newU = newSys === "SI" ? siUnit(cols[i].si) : cols[i].us;
        inp.value = fmtIn(fromSI(toSI(parseFloat(inp.value), oldU), newU));
      });
    [...wrap.querySelectorAll("thead th")].forEach((th, i) => {
      const u = newSys === "SI" ? siUnit(cols[i].si) : cols[i].us;
      th.innerHTML = u ? `${fmtLabel(cols[i].label)} (${fmtLabel(u)})` : fmtLabel(cols[i].label);
    });
  }
  system = newSys;
  // Results belong to the inputs they were computed from: re-render them in the new
  // units when there are any, otherwise just restate the ready line.
  if (lastRes) render(lastRes);
  else setStatus(`units: ${system} · ready · ${contract.meta.cite}`, true);
}

const sig = (x) => (Number.isFinite(x) ? +x.toPrecision(6) : x);
function setStatus(msg, ok) { const s = $("status"); s.innerHTML = (ok ? "✓ " : "✗ ") + fmtLabel(msg); s.className = ok ? "" : "err"; }
function fail(msg) { $("overlay").textContent = "Error: " + msg; setStatus(msg, false); }

// --- wire up controls ---
document.addEventListener("DOMContentLoaded", () => {
  $("compute").addEventListener("click", doCompute);
  $("reset").addEventListener("click", () => { buildForm(); clearOutputs(); focusFirstInput(); });
  // Enter anywhere in the inputs (except multi-line JSON/table textareas) runs Compute.
  $("inputs").addEventListener("keydown", (e) => {
    if (e.key === "Enter" && e.target.tagName !== "TEXTAREA") { e.preventDefault(); doCompute(); }
  });
  $("inputs").addEventListener("change", () => { applyShowIf(); applyEnableIf(); highlightActiveMethod(); });
  // KaTeX is deferred; if the panel painted before it loaded, re-typeset once it's ready.
  window.addEventListener("load", () => {
    if (typeof contract !== "undefined" && contract && contract.about) renderMethodPanel();
  });
  // Browser back/forward switches the active app (history entries set in buildAppSelect).
  window.addEventListener("popstate", () => {
    const id = new URLSearchParams(location.search).get("app");
    if (id) loadApp(id);
  });
  document.querySelectorAll('input[name="units"]').forEach((r) =>
    r.addEventListener("change", (e) => onUnits(e.target.value)));
  document.querySelectorAll(".tabbar button").forEach((b) =>
    b.addEventListener("click", () => {
      document.querySelectorAll(".tabbar button").forEach((x) => x.classList.remove("active"));
      b.classList.add("active");
      $("plotPane").style.display = b.dataset.tab === "plot" ? "" : "none";
      $("tablePane").style.display = b.dataset.tab === "table" ? "" : "none";
    }));
  // every export needs something computed first; say so rather than emitting a blank
  const needResult = () => {
    if (lastRes) return true;
    setStatus("nothing to export yet - press Compute first", false);
    return false;
  };
  $("copy").addEventListener("click", () => {
    if (needResult()) navigator.clipboard.writeText(tableCSV("\t"));
  });
  $("csv").addEventListener("click", () => {
    if (!needResult()) return;
    const blob = new Blob([tableCSV(",")], { type: "text/csv" });
    const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = "profile.csv"; a.click();
  });
  $("report").addEventListener("click", () => {
    if (!needResult()) return;
    const blob = new Blob([reportText()], { type: "text/plain" });
    const a = document.createElement("a"); a.href = URL.createObjectURL(blob);
    a.download = `chessqc_${contract.meta.aces_id}_report.txt`; a.click();
  });
  $("png").addEventListener("click", () => {
    if (!needResult()) return;
    const a = document.createElement("a"); a.href = $("plot").toDataURL("image/png");
    a.download = `chessqc_${contract.meta.aces_id}_plot.png`; a.click();
  });
  // pop-out: render the current plot into a large, zoomable full-window canvas
  const big = $("plotBig");
  $("popout").addEventListener("click", () => {
    if (!lastRes) return;
    $("plotModal").style.display = "flex";
    big._viewX = ($("plot") && $("plot")._viewX) || null;   // inherit current zoom
    requestAnimationFrame(() => drawPlotInto(lastRes, big));  // size known after layout
  });
  $("plotClose").addEventListener("click", () => { $("plotModal").style.display = "none"; });
  $("bigReset").addEventListener("click", () => {
    if (!lastRes) return;
    big._viewX = null; drawPlotInto(lastRes, big);
  });
  $("bigPng").addEventListener("click", () => {
    const a = document.createElement("a"); a.href = big.toDataURL("image/png");
    a.download = `chessqc_${contract.meta.aces_id}_plot.png`; a.click();
  });
  // Method & Equations pop-out: clone the (already-typeset) panel body into a wide modal.
  $("methodPop").addEventListener("click", (e) => {
    e.preventDefault(); e.stopPropagation();          // don't toggle the <details>
    $("methodModalBody").innerHTML = $("methodBody").innerHTML;
    $("methodModalTag").textContent = $("methodTag").textContent;
    $("methodModal").style.display = "flex";
  });
  $("methodModalClose").addEventListener("click", () => { $("methodModal").style.display = "none"; });
  window.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && $("plotModal").style.display !== "none") $("plotModal").style.display = "none";
    if (e.key === "Escape" && $("methodModal").style.display !== "none") $("methodModal").style.display = "none";
  });
  window.addEventListener("resize", () => {
    if ($("plotModal").style.display !== "none" && lastRes) drawPlotInto(lastRes, big);
  });
  boot();
});
