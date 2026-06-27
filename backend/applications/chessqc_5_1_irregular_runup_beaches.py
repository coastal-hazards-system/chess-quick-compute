"""CHESS-QC application 5-1 — Irregular Wave Runup on Beaches.

Originating ACES grouping: 5-1 "Irregular Wave Runup on Beaches" (functional area:
Wave Runup, Transmission, and Overtopping). Estimates how high irregular (real,
many-period) waves run up a smooth, uniform, impermeable beach slope, reporting five
runup statistics: the maximum, the level exceeded by 2 percent of runups, and the
averages of the highest tenth, highest third, and all runups.

Classification: exact (Mase 1989 runup power-law with the published (a,b) coefficient pairs
known for each statistic, nothing guessed; reproduces the User's Guide Example 5-1-4 to
0.01 ft).
Theory and references: each runup statistic is a power law in the Iribarren
(surf-similarity) number (TR chapter 5-1, eqs 1-2 in docs/EQUATIONS.md):

    R_p = H_s0 * a_p * xi^(b_p)            (1)
    xi  = tan(theta) / sqrt(H_s0 / L_0)    (2)   L_0 = g T_p^2 / (2 pi)

with H_s0 the deepwater significant wave height, T_p the peak period, tan(theta) the
foreshore slope, and (a_p, b_p) the statistic-specific coefficients of Mase (1989).
The coefficients are cited but not tabulated in the ACES TR; the published Mase (1989)
set is used here and reproduces the ACES worked example to 0.01 ft (see self-tests).

Self-containment: zero sibling imports; embeds its own Field/AppMeta/Out/Result
dataclasses. stdlib + numpy only. Runnable standalone:
    python chessqc_5_1_irregular_runup_beaches.py
which runs the ACES-oracle and consistency self-tests, then prints an ACES-style
tabulation.

Validation strategy: the five statistics reproduce the ACES User's Guide worked
example (H_s0 = 4.60 ft, T_p = 9.50 s, cot theta = 13.0 -> R_max = 8.74, R_2% = 7.11,
R_1/10 = 6.50, R_1/3 = 5.29, R_bar = 3.38 ft) to 0.01 ft; ordering and the linear
scaling in H_s0 are checked as closed-form consistency tests.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# --- standard physical constants (overridable; SI internal) ---------------------
G_SI = 9.80665           # m/s^2


# --- embedded contract dataclasses (self-contained; identical across all apps) --
@dataclass(frozen=True)
class AppMeta:
    aces_id: str
    name: str
    area: str
    classification: str          # "exact" | "standard" | "provisional"
    cite: str
    default_system: str = "SI"   # unit system the GUI opens in ("SI" | "US")
    status: str = "Current"      # operational currency: Current | Screening only | ...
    superseded_by: str = ""


@dataclass(frozen=True)
class Field:
    """One GUI/contract input field."""
    key: str
    label: str
    kind: str = "float"            # float | int | choice | bool | angle | file
    unit_si: str = ""
    unit_us: str = ""
    default: object = 0.0
    lo: float = -math.inf
    hi: float = math.inf
    choices: tuple = ()
    note: str = ""


@dataclass(frozen=True)
class Out:
    """One output descriptor (metadata; values live in Result)."""
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"           # scalar | point | profile
    note: str = ""                 # hover definition shown on the output label


# --- application metadata --------------------------------------------------------
APP_META = AppMeta(
    aces_id="5-1",
    name="Irregular Wave Runup on Beaches",
    area="Wave Runup, Transmission, and Overtopping",
    classification="exact",
    cite="Mase (1989); Hunt (1959); Walton & Ahrens (1989)",
    default_system="US",
)

_FT = 0.3048
# Default example = ACES User's Guide worked example (5-1-4): deepwater swell on a
# 1:13 foreshore.  Inputs stored SI-internal (heights in metres).
INPUTS = (
    Field("Hs0", "Deepwater significant wave height", "float", "m", "ft",
          default=4.60 * _FT, lo=1e-4, hi=1e3, note="> 0 (H_s0)"),
    Field("Tp", "Peak energy wave period", "float", "s", "s",
          default=9.50, lo=1e-2, hi=1e3, note="> 0"),
    Field("cot_theta", "Cotangent of foreshore slope", "float", "", "",
          default=13.0, lo=1e-3, hi=1e4, note="cot(theta) > 0; tan(theta) = 1/cot(theta)"),
)

OUTPUTS = (
    Out("L0",     "Deepwater wave length",                  "m", "ft", "scalar",
        note="Deepwater wavelength L_0 = g T_p^2 / (2 pi) from the peak wave period."),
    Out("xi",     "Iribarren (surf-similarity) number",     "",  "",   "scalar",
        note="Iribarren (surf-similarity) number xi = tan(theta) / sqrt(H_s0 / L_0); larger xi means steeper slope or longer/lower waves and higher runup."),
    Out("R_max",  "Maximum runup",                          "m", "ft", "scalar",
        note="Maximum wave runup, the highest vertical excursion above still-water level, R = H_s0 a xi^b with Mase (1989) coefficients (2.32, 0.77)."),
    Out("R_2",    "Runup exceeded by 2% of runups",         "m", "ft", "scalar",
        note="Runup height exceeded by only 2 percent of runups, above still-water level, with Mase (1989) coefficients (1.86, 0.71)."),
    Out("R_1_10", "Average of highest 1/10 of runups",      "m", "ft", "scalar",
        note="Average vertical runup of the highest one-tenth of runups above still-water level, with Mase (1989) coefficients (1.70, 0.71)."),
    Out("R_1_3",  "Average of highest 1/3 of runups",       "m", "ft", "scalar",
        note="Significant runup, the average of the highest one-third of runups above still-water level, with Mase (1989) coefficients (1.38, 0.70)."),
    Out("R_mean", "Average runup",                          "m", "ft", "scalar",
        note="Mean vertical runup averaged over all runups above still-water level, with Mase (1989) coefficients (0.88, 0.69)."),
)

# Mase (1989) power-law coefficients (R_p = H_s0 * a * xi^b), ordered high to low.
MASE_COEFFS = {
    "R_max":  (2.32, 0.77),
    "R_2":    (1.86, 0.71),
    "R_1_10": (1.70, 0.71),
    "R_1_3":  (1.38, 0.70),
    "R_mean": (0.88, 0.69),
}


@dataclass
class Result:
    L0: float; xi: float
    R_max: float; R_2: float; R_1_10: float; R_1_3: float; R_mean: float
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Estimates how high irregular waves run up a smooth, uniform, impermeable '
            'beach slope using the Mase (1989) power-law fit in the Iribarren number, '
            'returning five runup statistics: the maximum, the level exceeded by 2 percent '
            'of runups, and the averages of the highest tenth, highest third, and all '
            'runups.',
 'methods': [{'name': 'Mase (1989) power-law runup statistics',
              'when': None,
              'tag': '',
              'note': 'Coefficient pairs (a_p, b_p): R_max = (2.32, 0.77); R_2% = (1.86, '
                      '0.71); R_1/10 = (1.70, 0.71); R_1/3 = (1.38, 0.70); mean = (0.88, '
                      '0.69). Mase (1989) gentle-slope data span roughly cot 5 to 30 and '
                      'xi up to about 3.',
              'equations': [{'tex': 'L_0 = \\frac{g \\, T_p^{2}}{2 \\pi}',
                             'desc': 'Deepwater wavelength from the peak wave period.'},
                            {'tex': '\\xi = \\frac{\\tan\\theta}{\\sqrt{H_{s0} / L_0}}',
                             'desc': 'Iribarren (surf-similarity) number from foreshore '
                                     'slope and deepwater wave steepness.'},
                            {'tex': 'R_p = H_{s0} \\, a_p \\, \\xi^{b_p}',
                             'desc': 'Each runup statistic is a power law in the Iribarren '
                                     'number with statistic-specific Mase (1989) '
                                     'coefficients (a_p, b_p).'}]}],
 'symbols': [['R_p',
              'Wave runup for statistic p, vertical height above still-water level'],
             ['R_max', 'Maximum wave runup'],
             ['R_2%', 'Runup exceeded by 2% of runups'],
             ['R_1/10', 'Average of the highest 1/10 of runups'],
             ['R_1/3', 'Average of the highest 1/3 of runups (significant runup)'],
             ['R_mean', 'Mean wave runup'],
             ['H_s0', 'Deepwater significant wave height'],
             ['T_p', 'Peak energy wave period'],
             ['L_0', 'Deepwater wavelength'],
             ['xi', 'Iribarren (surf-similarity) number'],
             ['tan theta', 'Foreshore beach slope (= 1/cot theta)'],
             ['a_p, b_p', 'Statistic-specific Mase (1989) power-law coefficients'],
             ['g', 'Gravitational acceleration']],
 'references': ['Mase (1989)',
                'Hunt (1959)',
                'Walton & Ahrens (1989)',
                'Mase & Iwagaki (1984)',
                'ACES TR Chapter 5-1']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Irregular-wave runup statistics for SI inputs {Hs0, Tp, cot_theta}."""
    _validate(inp)
    Hs0 = float(inp["Hs0"]); Tp = float(inp["Tp"]); cot_theta = float(inp["cot_theta"])

    L0 = g * Tp * Tp / (2.0 * math.pi)                 # deepwater wavelength (eq 2)
    tan_theta = 1.0 / cot_theta
    xi = tan_theta / math.sqrt(Hs0 / L0)               # Iribarren number (eq 2)

    vals = {k: Hs0 * a * xi ** b for k, (a, b) in MASE_COEFFS.items()}   # eq (1)

    notes = []
    notes.append(f"Iribarren number xi = {xi:.3f}")
    # Mase (1989) gentle-slope data span roughly cot 5 to 30 and xi up to ~3
    if not (0.1 <= xi <= 3.0):
        notes.append("note: xi outside the Mase (1989) gentle-slope data range (~0.1 to 3)")
    notes.append("smooth uniform impermeable beach; see 5-2 for armored structures")

    return Result(L0=L0, xi=xi, notes="; ".join(notes), **vals)


# --- self-tests (ACES oracle + closed-form consistency) -------------------------
def _approx(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    # 1) ACES User's Guide worked example (5-1-4), exact to 0.01 ft
    r = compute({"Hs0": 4.60 * _FT, "Tp": 9.50, "cot_theta": 13.0}, g=g)
    oracle_ft = {"R_max": 8.74, "R_2": 7.11, "R_1_10": 6.50, "R_1_3": 5.29, "R_mean": 3.38}
    for k, want in oracle_ft.items():
        got = getattr(r, k) / _FT
        assert _approx(got, want, 0.01), f"{k}: got {got:.3f} ft, oracle {want} ft"
    assert _approx(r.xi, 0.771, 1e-3), r.xi

    # 2) ordering R_max > R_2% > R_1/10 > R_1/3 > R_bar
    assert r.R_max > r.R_2 > r.R_1_10 > r.R_1_3 > r.R_mean, r

    # 3) formula consistency: each statistic equals H_s0 * a * xi^b exactly
    Hs0 = 4.60 * _FT
    for k, (a, b) in MASE_COEFFS.items():
        assert _approx(getattr(r, k), Hs0 * a * r.xi ** b, 1e-12), k

    # 4) steeper slope (smaller cot) raises xi and therefore every runup statistic
    r_steep = compute({"Hs0": 4.60 * _FT, "Tp": 9.50, "cot_theta": 8.0}, g=g)
    assert r_steep.xi > r.xi and r_steep.R_max > r.R_max, (r_steep.xi, r.xi)

    print("  self-tests: PASS (ACES 5-1-4 oracle to 0.01 ft, ordering, formula, slope trend)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS:")
    for f in INPUTS:
        u = f.unit_us or f.unit_si
        v = inp[f.key] / _FT if f.unit_si == "m" else inp[f.key]
        print(f"    {f.label:34s} {f.key:10s} = {v:>10.4g} {u}")
    print("  OUTPUTS (US ft):")
    by_key = {o.key: o for o in OUTPUTS}
    for kk in ("L0", "xi", "R_max", "R_2", "R_1_10", "R_1_3", "R_mean"):
        o = by_key[kk]
        val = getattr(r, kk)
        disp = val / _FT if o.unit_si == "m" else val
        print(f"    {o.label:34s} {kk:8s} = {disp:>10.3f} {o.unit_us or ''}")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
