"""Generate docs/APPLICATION_STATUS.md: what each application is, what source material
exists for it, and what it has been checked against.

Run:  python common/gen_app_status.py   ->  docs/APPLICATION_STATUS.md

Everything that can be read from the repository is read, so the document cannot drift
from it: the application list, names, areas, classifications and citations come from each
module's APP_META, the ACES provenance from its docstring, the DOS-workbook coverage from
the recorded corpus summary, and the compiled-source drivers from the build script. The
three things that cannot be derived -- which ACES FORTRAN file implements an application,
which of the Hawaii ports covers it, and the one-line verification summary -- are tables
below, each entry checked against the deliverable.
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APP_DIR = os.path.join(ROOT, "backend", "applications")
DOC = os.path.join(ROOT, "docs", "APPLICATION_STATUS.md")
CORPUS = os.path.join(ROOT, "tests", "aces_oracle", "corpus", "SUMMARY.json")
FORTRAN = os.path.join(ROOT, "tests", "aces_oracle", "fortran")
MANUAL = os.path.join(ROOT, "tests", "test_manual_oracle.py")

# --- ACES FORTRAN source, from the dispatch table in INTDRV.FOR ------------------
# CHESS-QC id -> (ACES application, driver routine, source file). "-" where CHESS-QC
# has no ACES counterpart. Note 6-3 and 6-5: ACES numbers them differently.
ACES_SOURCE = {
    "1-1": ("1-1", "WAWG",     "WAWG.FOR"),
    "1-2": ("1-2", "BETRDRV",  "BETAR.FOR"),
    "1-3": ("1-3", "EXTWDRV",  "EXTGODA.FOR"),
    "1-4": ("1-4", "TIDRV",    "TIDES.FOR, CNS2LIB.FOR"),
    "1-5": ("-",   "-",        "-"),
    "1-6": ("-",   "-",        "-"),
    "2-1": ("2-1", "LWT",      "LWT.FOR"),
    "2-2": ("2-2", "CWT",      "CWT.FOR"),
    "2-3": ("2-3", "FWT",      "FWT.FOR"),
    "2-4": ("-",   "-",        "-"),
    "2-5": ("-",   "-",        "-"),
    "3-1": ("3-1", "LWTS",     "LWTS.FOR"),
    "3-2": ("3-2", "WSU",      "WSU.FOR"),
    "3-3": ("3-3", "DFRAC",    "DFRAC.FOR"),
    "3-4": ("3-3", "DFRAC",    "DFRAC.FOR (grid evaluation of the same kernel)"),
    "4-1": ("4-1", "HUDSON",   "HUDSON.FOR"),
    "4-2": ("4-2", "TOEPRO",   "TOEPRO.FOR"),
    "4-3": ("4-3", "WFVW",     "WFVW.FOR"),
    "4-4": ("4-4", "RUBBLE",   "RUBBLE.FOR"),
    "5-1": ("5-1", "RUBCH",    "RUBCH.FOR"),
    "5-2": ("5-2", "RUNOVT",   "RUNOVT.FOR"),
    "5-3": ("5-3", "RUNTRN",   "RUNTRN.FOR"),
    "5-4": ("5-4", "MADTRANS", "MADTRANS.FOR"),
    "5-5": ("-",   "-",        "-"),
    "6-1": ("6-1", "LSTDRV",   "Lstran.for (2/98) and lstran.forg (3/92)"),
    "6-2": ("6-2", "KDDRV1",   "KDMAIN.FOR, KD2-KD18"),
    "6-3": ("6-1", "LSTDRV",   "Lstran.for, CEDRS capability (lscband, lscomp, LSCTRN)"),
    "6-4": ("6-4", "BFDRV",    "BEACH.FOR"),
    "6-5": ("6-3", "CGS",      "CGS.FOR, CGSANL.FOR"),
    "7-1": ("7-1", "IHDRV",    "INLMAIN.FOR, IH1-IH14"),
    "7-2": ("-",   "-",        "-"),
    "8-1": ("-",   "-",        "-"),
    "8-2": ("-",   "-",        "-"),
    "8-3": ("-",   "-",        "-"),
    "9-1": ("-",   "-",        "-"),
    "10-1": ("-",  "-",        "-"),
    "10-2": ("-",  "-",        "-"),
    "10-3": ("-",  "-",        "-"),
    "10-4": ("-",  "-",        "-"),
    "M-1": ("-",   "HLIMIT",   "HLIMIT.FOR (shared breaker utilities, not a menu entry)"),
}

# --- Hawaii 2017 ports, by driver file name --------------------------------------
HAWAII_MATLAB = {
    "1-1": "wind_adj.m", "1-2": "beta_rayleigh.m", "1-3": "ext_Hs_analysis.m",
    "1-4": "tide_generation.m", "2-1": "linear_wave_theory.m",
    "2-2": "cnoidal_wave_theory.m", "2-3": "fourier.m", "3-1": "snells_law.m",
    "3-2": "irr_wave_trans.m", "3-3": "refdiff_vert_wedge.m",
    "4-1": "breakwater_Hudson.m", "4-2": "toe_design.m", "4-3": "wave_forces.m",
    "4-4": "rubble_mound.m", "5-1": "irregular_runup.m", "5-2": "runup_overtopping.m",
    "5-3": "wavetrans_imperm.m", "5-4": "wavetrans_perm.m",
    "6-1": "longshore_trans.m", "6-4": "beach_nourishment.m",
    "6-5": "compositeGrain.m",
}
HAWAII_PYTHON = {k: v.replace(".m", ".py") for k, v in HAWAII_MATLAB.items()
                 if k not in ("2-3", "6-1", "6-5")}

# --- compiled-source drivers, from tests/aces_oracle/fortran/build.sh ------------
# 5-5 is deliberately absent: WSU (the `goda` driver) produces the setup that ACES
# reports for 3-2, not the closed-form model CHESS-QC 5-5 implements. Treating it as an
# oracle for 5-5 was an earlier mapping error and is recorded in FINDINGS.
COMPILED = {
    "1-4": "tides", "3-2": "goda", "6-5": "cgs", "7-1": "inlet",
    "6-1": "lst / lst92", "6-3": "lst / lst92", "6-2": "kd", "2-3": "fwt",
}

# --- one-line verification summary ----------------------------------------------
# Each is the headline result recorded in the module docstring, the manual-oracle
# test, or tests/aces_oracle/FINDINGS.md.
VERIFIED = {
    "1-1": "1,017 recorded DOS cases; wind adjustment residual 0.79% tracks the 0.60% "
           "in the wind itself (FINDINGS D5). Growth outputs open at 1.4x tolerance",
    "1-2": "315 DOS cases, all 10 outputs clean. The H1/10 6.55-vs-6.30 dispute was "
           "settled from BETAR.FOR: ACES's own 100-bin quadrature (A1)",
    "1-3": "two recorded runs of the original program, all fitted distributions",
    "1-4": "the FORTRAN recompiled on its own deck: the full 481-point record to "
           "0.0002 ft, and the node factors and phase constants of all 37 "
           "constituents to 5e-7 and 0.01 deg (E14)",
    "1-5": "analytic: the Ekman-layer relations and their limits",
    "1-6": "analytic: Holland's profile, its gradient-wind balance and known limits",
    "2-1": "516 DOS cases over 4 sheets; 48 of 56 fields clean",
    "2-2": "DOS cases at both first and second order (the workbook stacks them, C5)",
    "2-3": "the ACES FORTRAN recompiled and run: 8 cases spanning deep and shallow "
           "water, near-breaking steepness, both celerity definitions and both current "
           "directions, to 0.008% on wavelength and celerity (E13)",
    "2-4": "analytic identities against 2-1, of which it is a presentation",
    "2-5": "analytic: the solitary-wave relations and their conservation properties",
    "3-1": "4,579 DOS cases. Three workbook columns were found to be mislabelled and "
           "were re-identified on 100% of rows (C1)",
    "3-2": "the FORTRAN recompiled: all 11 outputs of a 133-case sweep to print "
           "rounding. Two Hawaii mistranscriptions corrected from it (E1, E2)",
    "3-3": "149 DOS cases; the far-field metric rows sit at ~2% where ACES sums up to "
           "200 alternating single-precision terms (B1, B2 fixed)",
    "3-4": "the same kernel as 3-3 on a grid; checked against it pointwise",
    "4-1": "DOS cases across both unit systems. The armour-weight column was found to "
           "switch units mid-column and is handled (C2)",
    "4-2": "DOS cases; the wide-apron stability number is chaotic by construction and "
           "is documented rather than forced (B7)",
    "4-3": "DOS cases, Sainflou and Miche-Rundgren, crest and trough",
    "4-4": "DOS cases; runup exact",
    "5-1": "DOS cases against User's Guide Example 5-1-4",
    "5-2": "DOS cases, User's Guide Examples 1-7",
    "5-3": "DOS cases, User's Guide Examples 1-4",
    "5-4": "DOS cases over three breakwater geometries; five source corrections made "
           "from MADTRANS.FOR. Three isolated rows remain open (B8)",
    "5-5": "closed-form Longuet-Higgins/Weggel. NOT an ACES application: the setup ACES "
           "reports comes from WSU, which serves 3-2. Earlier mapping was wrong",
    "6-1": "the FORTRAN recompiled, both shipped vintages: the published examples to "
           "6 figures on the 3/92 source, and 5e-5 on both with quartz (E8)",
    "6-2": "the FORTRAN recompiled: 4 of the 5 shipped XSHORE decks, volume and all "
           "four contour changes to better than half the source's printed digit (E11)",
    "6-3": "the FORTRAN recompiled over the shipped CEDRS station file: all 9 "
           "contributing bands and the net to within 0.6 yd^3/yr (E12)",
    "6-4": "13 usable DOS rows - all the workbook holds; ACES has no multi-case mode "
           "here. ACES's normal CDF is invalid for negative arguments (B6)",
    "6-5": "the FORTRAN recompiled: all 9 statistics. It found a half-class bias in "
           "CHESS-QC's moment mean that the old self-referential test could not (E5)",
    "7-1": "the FORTRAN recompiled on the shipped deck: the full equal-discharge flow "
           "net, peak discharges to 0.01% and exchange volumes to 0.04% (E6)",
    "7-2": "analytic: dispersion and action identities, blocking, no-current unity",
    "8-1": "analytic: Merian, 2-D and Helmholtz closed forms and their reductions",
    "8-2": "analytic: Schijf continuity and energy, the 35.26 deg deep-water crest",
    "8-3": "analytic: single-line period, symmetry, projection and scaling laws",
    "9-1": "analytic against Bodine (1971) TM-35; beyond the ACES set",
    "10-1": "analytic: least-squares trend removal and its invariants",
    "10-2": "analytic; chained from 10-1",
    "10-3": "analytic; chained from 10-2",
    "10-4": "analytic",
    "M-1": "the breaker and steepness limits against their published forms, and for "
           "consistency with 2-1 and 3-1 which use them",
}


# --- one- or two-sentence description of each application -----------------------
DESCRIPTION = {
    "1-1": "Adjusts an observed wind to a 10 m, 1-hour, overwater equivalent, then "
           "grows a wave field from it over open water or a restricted fetch.",
    "1-2": "Gives the full shallow-water wave-height distribution from a significant "
           "height, period and depth, using the beta-Rayleigh form rather than "
           "assuming Rayleigh.",
    "1-3": "Fits extremal distributions to a series of annual maxima and returns "
           "design wave heights with confidence intervals for chosen return periods.",
    "1-4": "Predicts a water-level time series from tidal harmonic constituents, with "
           "the Schureman node factors and equilibrium arguments for the epoch.",
    "1-5": "Converts a geostrophic wind to the near-surface wind over water through "
           "the Ekman-layer relations, accounting for air-sea temperature difference.",
    "1-6": "Builds the radial wind and pressure profile of a tropical cyclone from "
           "Holland's two-parameter model.",
    "2-1": "First-order (Airy) wave theory: celerity, length, kinematics, pressure and "
           "energy for a wave of given height, period and depth.",
    "2-2": "Finite-amplitude periodic long-wave theory in terms of the cnoidal "
           "functions, for waves too long and steep for linear theory.",
    "2-3": "Steady progressive wave of permanent form solved as an N-term "
           "stream-function Fourier series, valid to near breaking where the analytic "
           "theories fail.",
    "2-4": "The classical wave-parameter presentation of linear theory: the "
           "dimensionless ratios and the deep, transitional and shallow-water limits.",
    "2-5": "Solitary wave theory: a single wave of translation lying entirely above the "
           "still-water level, with its celerity, profile and kinematics.",
    "3-1": "Transforms a wave from one depth to another, and to deep water, by linear "
           "theory with Snell's law refraction and energy-flux shoaling.",
    "3-2": "Transforms an irregular wave train shoreward by Goda's method, marching the "
           "spectrum through shoaling, refraction, breaking and wave setup.",
    "3-3": "Wave height, phase and pressure in the lee of a semi-infinite vertical "
           "wedge, combining diffraction and reflection at a single point.",
    "3-4": "The same wedge kernel evaluated over a uniform grid, giving a field rather "
           "than a single point.",
    "4-1": "Sizes primary armour units for a rubble-mound breakwater or revetment from "
           "the Hudson stability equation, with crest width, layer thickness and "
           "placement density.",
    "4-2": "Designs the toe apron width and toe-stone weight for a vertical wall, "
           "bulkhead or revetment.",
    "4-3": "Standing-wave (clapotis) forces and overturning moments on a vertical wall, "
           "by both the Sainflou and Miche-Rundgren methods.",
    "4-4": "Sizes the armour and filter layers of a rubble-mound revetment, with runup "
           "and the resulting crest elevation.",
    "5-1": "Estimates irregular-wave runup on a beach or a rough impermeable slope.",
    "5-2": "Runup and overtopping rate on an impermeable structure, including the "
           "effect of an onshore wind.",
    "5-3": "Wave transmission over an impermeable structure by overtopping, and the "
           "transmitted height behind it.",
    "5-4": "Wave transmission through and over a permeable rubble-mound structure, by "
           "Madsen's method with the internal flow resistance.",
    "5-5": "Wave-induced change in mean water level across the surf zone: the set-down "
           "seaward of breaking and the set-up shoreward of it.",
    "6-1": "Potential longshore sand transport rate from a single wave condition, by "
           "the CERC energy-flux method, from either deepwater or breaking inputs.",
    "6-2": "Time-marches beach and dune profile change through a storm by the XSHORE "
           "explicit finite-difference scheme, reporting eroded volume and contour "
           "recession.",
    "6-3": "Net and gross longshore transport from a full CEDRS directional wave "
           "climate, refracting and shoaling each height/period/direction combination "
           "in to breaking.",
    "6-4": "Tells a nourishment designer how much borrow material a given beach fill "
           "needs, from the grain-size distributions of the native and borrow sand.",
    "6-5": "Grain-size statistics of a composited sediment sample: the method of "
           "moments and the Folk graphic measures.",
    "7-1": "Time-marches inlet discharge and bay water level through a tidal cycle for "
           "a multi-section inlet, resolving the channel flow net across each "
           "cross-section.",
    "7-2": "Wave-current interaction in a channel: how a following or opposing current "
           "changes wave height and length, and when the wave is blocked.",
    "8-1": "Natural resonant oscillation periods of a rectangular harbour basin, open "
           "or closed, including the Helmholtz mode.",
    "8-2": "Waves generated by a vessel moving in a channel, and the resulting drawdown "
           "and return current.",
    "8-3": "Surge motion of a vessel held at a berth by mooring lines, and the "
           "resulting line loads.",
    "9-1": "Open-coast hurricane surge along a single cross-shelf traverse by the "
           "quasi-1D bathystrophic method.",
    "10-1": "Removes the long-term linear sea-level trend from a water-level record by "
            "least-squares regression.",
    "10-2": "Separates the non-tidal residual from a water-level record by removing the "
            "predicted astronomical tide.",
    "10-3": "Extracts independent storm peaks from a residual series by a "
            "peaks-over-threshold declustering.",
    "10-4": "Builds a synthetic storm population from a fitted distribution and "
            "propagates it to a hazard curve.",
    "M-1": "The shared breaker and steepness relations the other applications call: the "
           "Miche steepness limit, the McCowan and Weggel breaking criteria, and the "
           "breaker height at a structure.",
}

_AREA_ORDER = ["Wave Prediction", "Wave Theory", "Wave Transformation",
               "Structural Design", "Wave Runup, Transmission, and Overtopping",
               "Littoral Processes", "Inlet Processes", "Harbor Design",
               "Storm Surge", "Coastal Hazards", "Miscellaneous Routines"]


def _load(path):
    spec = importlib.util.spec_from_file_location("m_" + os.path.basename(path)[:-3], path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _brief(doc: str) -> str:
    """One or two sentences of description, from the module docstring."""
    body = doc.split("\n", 1)[1] if "\n" in doc else doc
    # drop the "Originating ACES ..." lead-in and keep what follows
    m = re.search(r"\(functional area:[^)]*\)\.\s*(.+?)(?:\n\n|Classification:)", body, re.S)
    if not m:
        m = re.search(r"\n\n(.+?)(?:\n\n|Classification:)", body, re.S)
    text = " ".join(m.group(1).split()) if m else ""
    # trim to the first two sentences
    parts = re.split(r"(?<=[.])\s+", text)
    out = " ".join(parts[:2]).strip()
    return out or "(see the module docstring)"


def main():
    apps = []
    for fn in sorted(os.listdir(APP_DIR)):
        if not (fn.startswith("chessqc_") and fn.endswith(".py")):
            continue
        mod = _load(os.path.join(APP_DIR, fn))
        meta = getattr(mod, "APP_META", None)
        if meta is None:
            continue
        apps.append((meta, mod.__doc__ or "", fn))

    summary = {}
    if os.path.exists(CORPUS):
        summary = json.load(open(CORPUS, encoding="utf-8"))
    manual = set()
    if os.path.exists(MANUAL):
        txt = open(MANUAL, encoding="utf-8").read()
        for m in re.finditer(r"^def test_(\w+?)\(", txt, re.M):
            manual.add(m.group(1))

    def has_manual(aid):
        key = aid.replace("-", "_").lower()
        return any(t.startswith(key + "_") or t == key for t in manual)

    def sortkey(item):
        meta = item[0]
        try:
            a = _AREA_ORDER.index(meta.area)
        except ValueError:
            a = len(_AREA_ORDER)
        bits = meta.aces_id.replace("M-", "99-").split("-")
        return (a, int(bits[0]), int(bits[1]))

    apps.sort(key=sortkey)

    out = []
    w = out.append
    w("# CHESS-QC — application status")
    w("")
    w("What each application is, what source material exists for it, and what it has")
    w("been checked against. Generated from the repository by `common/gen_app_status.py`;")
    w("re-run it after any change so this cannot drift from the code.")
    w("")
    w("")
    mapped = [m.aces_id for m, _, _ in apps
              if ACES_SOURCE.get(m.aces_id, ("-",))[0] != "-"]
    distinct = {ACES_SOURCE[a][0] for a in mapped}
    w(f"{len(apps)} applications. {len(mapped)} of them implement the "
      f"{len(distinct)} applications of the original")
    w("ACES suite: ACES 3-3 and 6-1 each serve two here, because CHESS-QC splits the")
    w("wedge kernel into point and grid forms, and the longshore method into")
    w("single-condition and wave-climate forms. The remaining "
      f"{len(apps) - len(mapped)} are Quick Compute")
    w("additions in the same style, with no ACES counterpart.")
    w("")
    w("**Classification** is the module's own: *exact* means the method is reproduced")
    w("with nothing fitted or guessed; *standard* means a documented approximation or an")
    w("unresolved residual remains. Both are stated per application below.")
    w("")

    # --- summary table ---
    w("## Summary")
    w("")
    w("| ID | Application | Class | ACES | FORTRAN | MATLAB | Python | DOS cases | Source run |")
    w("|---|---|---|---|---|---|---|---|---|")
    for meta, doc, fn in apps:
        aid = meta.aces_id
        aces_app, _routine, srcf = ACES_SOURCE.get(aid, ("-", "-", "-"))
        cases = summary.get(aid, {}).get("cases")
        w(f"| {aid} | {meta.name} | {meta.classification} | "
          f"{aces_app} | "
          f"{'yes' if srcf != '-' else '—'} | "
          f"{'yes' if aid in HAWAII_MATLAB else '—'} | "
          f"{'yes' if aid in HAWAII_PYTHON else '—'} | "
          f"{cases if cases else '—'} | "
          f"{COMPILED.get(aid, '—')} |")
    w("")
    w("*ACES* is the application's number in the original suite, which differs from")
    w("CHESS-QC's in two places: CHESS-QC 6-3 is a capability of ACES 6-1, and CHESS-QC")
    w("6-5 is ACES 6-3. *DOS cases* counts recorded cases from the Hawaii comparison")
    w("workbooks. *Source run* names the driver that compiles and runs the original")
    w("FORTRAN (`tests/aces_oracle/fortran`).")
    w("")

    # --- per-application detail ---
    last_area = None
    for meta, doc, fn in apps:
        aid = meta.aces_id
        if meta.area != last_area:
            w(f"## {meta.area}")
            w("")
            last_area = meta.area
        aces_app, routine, srcfile = ACES_SOURCE.get(aid, ("-", "-", "-"))
        w(f"### {aid} — {meta.name}")
        w("")
        w(f"*{meta.classification}*. "
          f"{DESCRIPTION.get(aid) or _brief(doc)}")
        w("")
        w("**Sources available**")
        w("")
        if aces_app != "-":
            w(f"- ACES Technical Reference and User's Guide, chapter {aces_app} "
              f"(`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)")
        else:
            w("- No numbered ACES application; built from the cited literature below")
        if srcfile != "-":
            w(f"- ACES FORTRAN source: `{srcfile}`"
              + (f", driver `{routine}`" if routine != "-" else ""))
        if aid in HAWAII_MATLAB:
            w(f"- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/{HAWAII_MATLAB[aid]}`")
        if aid in HAWAII_PYTHON:
            w(f"- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/{HAWAII_PYTHON[aid]}`")
        if aid in summary:
            s = summary[aid]
            w(f"- Recorded DOS output: {s['cases']} cases over {s['sheets']} sheet(s), "
              f"{s['clean_fields']} of {s['fields']} fields comparing clean")
        if aid in COMPILED:
            w(f"- **The original FORTRAN, compiled and runnable**: "
              f"`sh tests/aces_oracle/fortran/build.sh {COMPILED[aid].split(' /')[0]}`")
        w(f"- Literature: {meta.cite}")
        w("")
        w("**Verification and validation**")
        w("")
        w(f"- {VERIFIED.get(aid, '(not recorded)')}")
        w(f"- Self-tests in the module; {'a dedicated case in `tests/test_manual_oracle.py`'
                                        if has_manual(aid) else 'covered by the shared suites'}")
        if aid in summary:
            w("- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)")
        w("")

    w("---")
    w("")
    w("## How to read the verification column")
    w("")
    w("Three kinds of evidence appear, and they are not equally strong:")
    w("")
    w("1. **The original FORTRAN, recompiled and run.** The strongest available: the")
    w("   ACES computational routines separate cleanly from their screen layer and build")
    w("   with the compiler shipped alongside them, so results can be compared at full")
    w("   precision on any input. "
      f"{len({a for a in COMPILED if a in {m.aces_id for m, _, _ in apps}})} "
      "applications now have this.")
    w("2. **Recorded DOS output.** 28,742 cases from the Hawaii workbooks, compared with")
    w("   absolute tolerances derived from the printed field widths rather than")
    w(f"   percentages. {len([a for a in summary])} applications have this.")
    w("3. **Analytic checks.** For applications with no ACES counterpart: closed-form")
    w("   identities, limiting cases and conservation properties. Weaker than either of")
    w("   the above, and labelled as such.")
    w("")
    w("Where a published example and the source disagree, the source wins, and the")
    w("disagreement is recorded. Where the source and a recorded DOS listing disagree,")
    w("the listing wins. Open items and the full evidence are in")
    w("`tests/aces_oracle/FINDINGS.md`.")
    w("")

    os.makedirs(os.path.dirname(DOC), exist_ok=True)
    open(DOC, "w", encoding="utf-8", newline="\n").write("\n".join(out))
    print(f"gen_app_status: {len(apps)} apps -> docs/APPLICATION_STATUS.md")


if __name__ == "__main__":
    main()
