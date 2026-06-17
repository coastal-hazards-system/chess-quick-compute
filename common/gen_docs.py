"""Auto-generate CHESS-QC documentation from the application contracts.

Single source of truth = each app's embedded contract (APP_META / INPUTS / OUTPUTS)
plus its module docstring. Emits two markdown docs that therefore never drift from
the code:

  docs/USER_MANUAL.md    : per-app pages: overview, inputs (units/range/default/notes),
                           outputs, citation, and a status banner.
  docs/OUTDATED_APPS.md  : every app whose status != "Current" (or that has a newer
                           method in `superseded_by`), extracted from the same metadata.

The hand-written Technical Reference (from docs/EQUATIONS.md) is the third doc and
the source of the Status & Caveats wording; this tool covers the two auto-generated ones.

Run:  python common/gen_docs.py
"""
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # common/ -> repo root
APP_DIR = os.path.join(ROOT, "backend", "applications")
DOC_DIR = os.path.join(ROOT, "docs")


def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# common/units.py is pure (no Qt); load it directly for SI->display conversion
_UNITS = _load_module(os.path.join(ROOT, "common", "units.py"), "doc_units")

# canonical functional-area order (matches the hub)
AREA_ORDER = [
    "Wave Prediction", "Wave Theory", "Wave Transformation", "Structural Design",
    "Wave Runup, Transmission, and Overtopping", "Littoral Processes",
    "Inlet Processes", "Harbor Design", "Storm Surge", "Miscellaneous Routines",
]

# Fidelity classification: A exact, B standard, C provisional.
_CLASS_LETTER = {"exact": "A", "standard": "B", "provisional": "C"}


def _letter(meta) -> str:
    c = getattr(meta, "classification", "")
    return _CLASS_LETTER.get(c, c)


def _load(path):
    name = "doc_" + os.path.basename(path)[:-3]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod                       # register before exec (dataclass annotations)
    spec.loader.exec_module(mod)
    return mod


def _overview(mod):
    """First paragraph of the module docstring (skip the leading title line)."""
    doc = (mod.__doc__ or "").strip()
    if not doc:
        return ""
    lines = doc.splitlines()
    body = lines[1:] if len(lines) > 1 else lines    # drop the "CHESS-QC application ..." title
    para = []
    for ln in body:
        if ln.strip() == "" and para:
            break
        if ln.strip():
            para.append(ln.strip())
    return " ".join(para)


def _units(unit_si, unit_us):
    if unit_si and unit_us and unit_si != unit_us:
        return f"{unit_us} / {unit_si}"
    return unit_si or unit_us or "(none)"


def _fmt_num(x):
    if isinstance(x, float):
        if x in (float("inf"), float("-inf")):
            return ""
        return f"{x:g}"
    return str(x)


def _input_rows(mod, system):
    """Rows for the inputs table, numeric values converted SI -> the display unit
    (US if the app opens in US, else SI)."""
    rows = []
    for f in mod.INPUTS:
        unit = (f.unit_si if system == "SI" else f.unit_us) or f.unit_si or f.unit_us
        conv = lambda v: _UNITS.from_si(v, unit) if isinstance(v, (int, float)) else v
        if f.kind == "choice":
            rng = "choices: " + ", ".join(str(c) for c in f.choices)
            default = str(f.default)
        elif f.kind == "bool":
            rng = "yes / no"
            default = "yes" if f.default else "no"
        elif f.kind == "table":
            cols = getattr(f, "columns", ()) or ()
            rng = "table: " + ", ".join(c[0] for c in cols) if cols else "table"
            default = f"{len(f.default or [])} default rows"
        else:
            lo, hi = _fmt_num(conv(f.lo)), _fmt_num(conv(f.hi))
            rng = f"{lo} to {hi}".strip() if (lo or hi) else "any"
            default = _fmt_num(conv(f.default))
        rows.append((f.label, f.key, _units(f.unit_si, f.unit_us), rng, default, f.note or ""))
    return rows


def _md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "| " + " | ".join("---" for _ in headers) + " |"]
    for r in rows:
        out.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(out)


def discover():
    apps = []
    for fn in sorted(os.listdir(APP_DIR)):
        if fn.startswith("chessqc_") and fn.endswith(".py"):
            mod = _load(os.path.join(APP_DIR, fn))
            apps.append((mod.APP_META, mod, fn))
    return apps


def _area_key(area):
    return (AREA_ORDER.index(area) if area in AREA_ORDER else len(AREA_ORDER), area)


def gen_user_manual(apps):
    lines = ["# CHESS-QC — User Manual",
             "*Coastal Hazards, Engineering, and Structures System (CHESS) — Quick Compute (QC)*",
             "",
             "Each tool computes in SI internally and displays in the selected unit system; "
             "ranges/defaults below are shown in display units. *Auto-generated from the "
             "application contracts by `common/gen_docs.py`; do not edit by hand.*",
             "",
             "Each app header carries a fidelity class: **A** exact · **B** standard · "
             "**C** provisional.",
             ""]
    by_area = {}
    for meta, mod, fn in apps:
        by_area.setdefault(meta.area, []).append((meta, mod, fn))
    # table of contents
    lines.append("## Contents")
    for area in sorted(by_area, key=_area_key):
        lines.append(f"- **{area}**")
        for meta, mod, fn in sorted(by_area[area], key=lambda t: t[0].aces_id):
            lines.append(f"  - {meta.aces_id} {meta.name}")
    lines.append("")
    for area in sorted(by_area, key=_area_key):
        lines.append(f"\n## {area}\n")
        for meta, mod, fn in sorted(by_area[area], key=lambda t: t[0].aces_id):
            status = getattr(meta, "status", "Current")
            superseded = getattr(meta, "superseded_by", "")
            lines.append(f"### {meta.aces_id} — {meta.name}  `[{_letter(meta)}]`\n")
            banner = f"**Status:** {status}."
            if superseded:
                banner += f" Newer method: {superseded}."
            lines.append(banner + "\n")
            ov = _overview(mod)
            if ov:
                lines.append(ov + "\n")
            sys_disp = getattr(meta, "default_system", "SI")
            lines.append(f"**Inputs** (values in {sys_disp} units)\n")
            lines.append(_md_table(["Input", "key", "units (US/SI)", "range", "default", "notes"],
                                   _input_rows(mod, sys_disp)))
            lines.append("\n**Outputs**\n")
            orows = [(o.label, o.key, _units(o.unit_si, o.unit_us), o.kind) for o in mod.OUTPUTS]
            lines.append(_md_table(["Output", "key", "units (US/SI)", "kind"], orows))
            lines.append(f"\n*Reference:* {meta.cite}\n")
            lines.append(f"*Module:* `backend/applications/{fn}`\n")
    return "\n".join(lines)


def gen_outdated(apps):
    flagged = [(m, mod, fn) for (m, mod, fn) in apps
               if getattr(m, "status", "Current") != "Current" or getattr(m, "superseded_by", "")]
    lines = ["# CHESS-QC — Outdated / Superseded Applications",
             "*Coastal Hazards, Engineering, and Structures System (CHESS) — Quick Compute (QC)*",
             "",
             "*Auto-generated by `common/gen_docs.py` from each app's `AppMeta.status` / "
             "`superseded_by`.* These tools remain available (often valuable for screening, "
             "teaching, or consistency with legacy studies) but a newer method is preferred "
             "for routine design. See the Technical Reference for the full Status & Caveats.",
             ""]
    if not flagged:
        lines.append("_No applications are currently flagged._")
        return "\n".join(lines)
    lines.append(_md_table(
        ["App", "Name", "Area", "Status", "Preferred / newer method"],
        [(m.aces_id, m.name, m.area, getattr(m, "status", "Current"),
          getattr(m, "superseded_by", "") or "(none)") for (m, mod, fn) in
         sorted(flagged, key=lambda t: _area_key(t[0].area))]))
    lines.append("")
    return "\n".join(lines)


def main():
    os.makedirs(DOC_DIR, exist_ok=True)
    apps = discover()
    man = os.path.join(DOC_DIR, "USER_MANUAL.md")
    out = os.path.join(DOC_DIR, "OUTDATED_APPS.md")
    with open(man, "w", encoding="utf-8") as fh:
        fh.write(gen_user_manual(apps))
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(gen_outdated(apps))
    n_flagged = sum(1 for (m, _, _) in apps
                    if getattr(m, "status", "Current") != "Current" or getattr(m, "superseded_by", ""))
    print(f"gen_docs: {len(apps)} apps -> {os.path.relpath(man, ROOT)}, "
          f"{os.path.relpath(out, ROOT)} ({n_flagged} flagged)")


if __name__ == "__main__":
    main()
