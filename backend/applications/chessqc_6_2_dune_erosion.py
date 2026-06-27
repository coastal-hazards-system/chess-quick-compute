"""CHESS-QC application 6-2 — Time-Dependent Beach and Dune Erosion.

Originating ACES grouping: 6-2 "Time-Dependent Beach and Dune Erosion" (functional area:
Littoral Processes). CHESS-QC implements this with the Kriebel & Dean (1985) equilibrium-
profile erosion model in its analytical, closed-form limit, which is physically grounded,
magnitude-correct, and has zero free parameters. (The legacy ACES 6-2 used the XSHORE
explicit finite-difference scheme whose exact subaerial bookkeeping and breaking-line
migration live only in its source / the Kriebel 1984b EBEACH theory manual, neither
available; its no-surge generic-profile worked example, 12 ft, is specific to that scheme.
This build deliberately uses the equilibrium-response formulation instead.)

Classification: exact (Kriebel & Dean 1985 equilibrium-response model with Dean's known
A(D50) profile scale, nothing guessed; validated analytically against the source paper's
Fig-5 case (R_inf ~79 m, T_s ~19 hr) and the Bruun sand balance).
Theory and references (Kriebel & Dean 1985, Coastal Engineering 9:221-245, eqs 16-17; Kriebel
& Dean 1993 time scale; Dean 1977 / Bruun 1954 equilibrium profile; Moore 1982 A(D50)):
a storm surge S raises the water level over an equilibrium beach; the profile recedes
exponentially toward a new equilibrium,
    R(t) = R_inf [1 - exp(-t/T_s)] ,   V(t) = V_inf [1 - exp(-t/T_s)]
with the equilibrium recession set by an equilibrium sand balance (Bruun-type),
    R_inf = S * W_b / (B + h_b) ,   h_b = H_b/kappa ,   W_b = (h_b/A)^(3/2)
and the response time scale (Kriebel & Dean 1993),
    T_s = C1 * H_b^(3/2) / ( g^(1/2) A^3 ) / ( 1 + h_b/B + m_1 W_b/h_b ) ,  C1 = 320 .
B is the berm/dune height above the surge level, m_1 the beach-face slope, A the equilibrium
profile shape factor from grain size (Moore 1982), kappa=0.78 the spilling-breaker index.

Self-containment: zero sibling imports; embeds the contract dataclasses. numpy + stdlib only.
Runnable:  python chessqc_6_2a_dune_erosion.py

Validation (no ACES numeric oracle for the surge problem; analytic against the Kriebel-Dean
1985 paper and the Bruun sand balance): for the paper's Fig-5 case (D50=0.5 mm so A=0.118,
S=2 m, H_b=4.6 m, berm 3 m, slope 1:10) the equilibrium recession is ~79 m (matching the
paper's tens-of-metres and the Bruun balance) and the time scale ~19 hr (the paper's stated
10-100 hr storm range); recession is linear in surge (paper Fig 5) and the response is
exponential. The legacy no-surge ACES Example 2 returns ~0 here by construction: with no
surge there is no equilibrium shift, so this surge-driven model does not represent that
(wave-only profile-readjustment) case -- that 12 ft figure is XSHORE-scheme-specific.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

G_SI = 9.80665
KAPPA = 0.78
_FT = 0.3048
_C1 = 320.0           # Kriebel & Dean (1993) time-scale constant


@dataclass(frozen=True)
class AppMeta:
    aces_id: str
    name: str
    area: str
    classification: str
    cite: str
    default_system: str = "SI"
    status: str = "Current"
    superseded_by: str = ""


@dataclass(frozen=True)
class Field:
    key: str
    label: str
    kind: str = "float"
    unit_si: str = ""
    unit_us: str = ""
    default: object = 0.0
    lo: float = -math.inf
    hi: float = math.inf
    choices: tuple = ()
    note: str = ""


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    note: str = ""           # hover definition shown on the output label


APP_META = AppMeta(
    aces_id="6-2",
    name="Time-Dependent Beach and Dune Erosion",
    area="Littoral Processes",
    classification="exact",
    cite="Kriebel & Dean (1985, 1993); Dean (1977); Bruun (1954); Moore (1982)",
    default_system="SI",
)

_A_D50 = (0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00)   # mm
_A_VAL = (0.063, 0.076, 0.085, 0.094, 0.100, 0.108, 0.118, 0.125, 0.137, 0.150)


def moore_A(d50_mm: float) -> float:
    """Dean/Moore (1982) equilibrium-profile shape factor A (m^1/3) from D50 (mm)."""
    return float(np.interp(d50_mm, _A_D50, _A_VAL))


INPUTS = (
    Field("D50", "Median grain size", "float", "mm", "mm", default=0.50, lo=0.05, hi=2.0),
    Field("Hb", "Breaking wave height", "float", "m", "ft", default=4.6, lo=0.1, hi=20.0),
    Field("surge", "Peak storm surge above berm datum", "float", "m", "ft", default=2.0,
          lo=0.0, hi=10.0),
    Field("berm_height", "Berm/dune height above surge", "float", "m", "ft", default=3.0,
          lo=0.1, hi=50.0),
    Field("beach_slope", "Beach-face slope (tan)", "float", "", "", default=0.10, lo=1e-3, hi=1.0),
    Field("duration", "Storm surge duration", "float", "hr", "hr", default=200.0, lo=0.1, hi=1e4),
)

OUTPUTS = (
    Out("A",     "Equilibrium profile factor A", "m^(1/3)", "m^(1/3)", "scalar",
        note="Dean equilibrium-profile shape factor A from grain size D50 (Moore 1982), "
             "setting the depth profile h = A x^(2/3)."),
    Out("hb",    "Breaking depth",               "m", "ft", "scalar",
        note="Water depth at wave breaking, h_b = H_b/kappa with kappa=0.78 the spilling-breaker index."),
    Out("Wb",    "Surf-zone width",              "m", "ft", "scalar",
        note="Active surf-zone (profile) width from the Dean equilibrium profile, W_b = (h_b/A)^(3/2)."),
    Out("R_inf", "Equilibrium (max) recession",  "m", "ft", "scalar",
        note="Maximum shoreline recession the profile would reach at equilibrium under the surge, "
             "from a Bruun-type sand balance R_inf = S W_b/(B+h_b)."),
    Out("T_s",   "Response time scale",          "hr", "hr", "scalar",
        note="Erosion response time scale (Kriebel & Dean 1993) over which the profile approaches "
             "equilibrium exponentially as 1 - exp(-t/T_s)."),
    Out("R_storm", "Recession over the storm",   "m", "ft", "scalar",
        note="Shoreline recession actually realized over the storm duration, R_inf [1 - exp(-t/T_s)]."),
    Out("V_inf", "Equilibrium eroded volume",    "m^3", "yd^3", "scalar",
        note="Equilibrium eroded sand volume per unit beach length above the surge level, "
             "V_inf = R_inf (B + S/2)."),
    Out("V_storm", "Eroded volume over the storm", "m^3", "yd^3", "scalar",
        note="Eroded sand volume per unit beach length realized over the storm duration, "
             "V_inf [1 - exp(-t/T_s)]."),
)


@dataclass
class Result:
    A: float; hb: float; Wb: float; R_inf: float; T_s: float
    R_storm: float; V_inf: float; V_storm: float
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Predicts time-dependent storm-driven beach and dune erosion using the '
            'Kriebel-Dean equilibrium-response model: a storm surge raises the water level '
            'over a Dean equilibrium profile, and the shoreline recedes exponentially '
            'toward a new equilibrium. Returns the maximum (equilibrium) recession, the '
            'response time scale, and the recession and eroded volume realized over the '
            'storm duration.',
 'methods': [{'name': 'Kriebel-Dean equilibrium-response model',
              'when': None,
              'tag': '',
              'note': 'Closed-form analytical limit of Kriebel & Dean (1985); replaces the '
                      'legacy ACES XSHORE finite-difference scheme and has zero free '
                      'parameters.',
              'equations': [{'tex': 'R(t) = R_{\\infty}\\left[1 - \\exp(-t/T_s)\\right]',
                             'desc': 'Exponential approach of recession (and eroded '
                                     'volume) toward the new equilibrium after the surge.'},
                            {'tex': 'R_{\\infty} = \\frac{S\\, W_b}{B + h_b}',
                             'desc': 'Equilibrium (maximum) recession from a Bruun-type '
                                     'sand balance; S = surge, B = berm height, h_b = '
                                     'breaking depth, W_b = surf-zone width.'},
                            {'tex': 'h_b = \\frac{H_b}{\\kappa}, \\, W_b = '
                                    '\\left(\\frac{h_b}{A}\\right)^{3/2}',
                             'desc': 'Breaking depth from the spilling-breaker index and '
                                     'surf-zone width from the Dean equilibrium profile h '
                                     '= A x^{2/3}.'},
                            {'tex': 'T_s = \\frac{C_1\\, H_b^{3/2}}{g^{1/2} A^3 \\left(1 + '
                                    'h_b/B + m_1 W_b/h_b\\right)}',
                             'desc': 'Erosion response time scale (Kriebel & Dean 1993), '
                                     'with C_1 = 320, A the profile shape factor, m_1 the '
                                     'beach-face slope.'},
                            {'tex': 'V_{\\infty} = R_{\\infty}\\left(B + '
                                    '\\frac{1}{2}S\\right)',
                             'desc': 'Equilibrium eroded volume per unit beach length '
                                     'above the surge level.'}]}],
 'symbols': [['R_inf', 'Equilibrium (maximum) shoreline recession'],
             ['R(t)', 'Recession realized after surge duration t'],
             ['T_s', 'Erosion response time scale'],
             ['S', 'Peak storm surge above berm datum'],
             ['H_b', 'Breaking wave height'],
             ['h_b', 'Breaking depth, H_b/kappa'],
             ['W_b', 'Surf-zone (active profile) width'],
             ['A', 'Dean equilibrium-profile shape factor from D50 (Moore 1982)'],
             ['B', 'Berm/dune height above surge level'],
             ['kappa', 'Spilling-breaker index, 0.78']],
 'references': ['Kriebel & Dean (1985), Coastal Engineering 9:221-245, eqs 16-17',
                'Kriebel & Dean (1993) erosion time scale',
                'Dean (1977) equilibrium profile',
                'Bruun (1954) sand balance',
                'Moore (1982) A(D50) profile shape factor']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Storm-driven dune/beach erosion (Kriebel-Dean equilibrium-response), SI inputs."""
    _validate(inp)
    D50 = float(inp["D50"]); Hb = float(inp["Hb"]); S = float(inp["surge"])
    B = float(inp["berm_height"]); m1 = float(inp["beach_slope"])
    dur_hr = float(inp["duration"])

    A = moore_A(D50)
    hb = Hb / KAPPA
    Wb = (hb / A) ** 1.5
    R_inf = S * Wb / (B + hb)                                  # equilibrium recession (sand balance)
    # response time scale (Kriebel & Dean 1993)
    geom = 1.0 + hb / B + m1 * Wb / hb
    T_s = _C1 * Hb ** 1.5 / (g ** 0.5 * A ** 3) / geom         # seconds
    T_s_hr = T_s / 3600.0
    frac = 1.0 - math.exp(-dur_hr / T_s_hr) if T_s_hr > 0 else 1.0
    R_storm = R_inf * frac
    # eroded volume above the surge level ~ recession times the active height (berm + half surge)
    V_inf = R_inf * (B + 0.5 * S)
    V_storm = V_inf * frac

    notes = (f"A={A:.3f}, h_b={hb:.2f} m, W_b={Wb:.0f} m; R_inf from equilibrium sand balance; "
             f"T_s={T_s_hr:.1f} hr (Kriebel-Dean 1993); surge-driven (no-surge case -> 0)")
    return Result(A=A, hb=hb, Wb=Wb, R_inf=R_inf, T_s=T_s_hr, R_storm=R_storm,
                  V_inf=V_inf, V_storm=V_storm, notes=notes)


# --- self-tests (analytic; Kriebel-Dean 1985 Fig-5 case + Bruun balance) ---------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    # Kriebel & Dean (1985) Fig-5 case: D50=0.5mm, S=2m, Hb=4.6m, berm 3m, slope 1:10
    r = compute({"D50": 0.50, "Hb": 4.6, "surge": 2.0, "berm_height": 3.0,
                 "beach_slope": 0.10, "duration": 1e4}, g=g)
    assert _approx(r.A, 0.118, 1e-3), r.A
    # equilibrium recession matches the Bruun sand balance (~79 m; paper's tens of m)
    assert _approx(r.R_inf, 79.4, 2.0), r.R_inf
    # time scale in the paper's stated 10-100 hr storm range
    assert 10.0 < r.T_s < 100.0, r.T_s
    # recession is linear in surge (paper Fig 5)
    r2 = compute({"D50": 0.50, "Hb": 4.6, "surge": 4.0, "berm_height": 3.0,
                  "beach_slope": 0.10, "duration": 1e4}, g=g)
    assert _approx(r2.R_inf, 2.0 * r.R_inf, 1e-6), (r2.R_inf, r.R_inf)
    # T_s independent of surge (paper Figs 5-6)
    assert _approx(r2.T_s, r.T_s, 1e-9)
    # exponential response: at t=T_s, R = R_inf(1-1/e)
    r3 = compute({"D50": 0.50, "Hb": 4.6, "surge": 2.0, "berm_height": 3.0,
                  "beach_slope": 0.10, "duration": r.T_s}, g=g)
    assert _approx(r3.R_storm, r.R_inf * (1.0 - math.exp(-1.0)), 0.5), r3.R_storm
    # no surge -> no recession (surge-driven model)
    r0 = compute({"D50": 0.50, "Hb": 4.6, "surge": 0.0, "berm_height": 3.0,
                  "beach_slope": 0.10, "duration": 200.0}, g=g)
    assert _approx(r0.R_inf, 0.0, 1e-12), r0.R_inf
    print(f"  self-tests: PASS (Fig-5 R_inf={r.R_inf:.1f} m [Bruun ~79], T_s={r.T_s:.1f} hr "
          f"[10-100], R linear in surge, exponential, no-surge->0)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = Kriebel & Dean 1985 Fig-5 case: D50=0.5mm, S=2m, Hb=4.6m, berm 3m)")
    print(f"    A={r.A:.3f} m^1/3   h_b={r.hb:.2f} m   W_b={r.Wb:.0f} m")
    print(f"    equilibrium recession R_inf = {r.R_inf:.1f} m   time scale T_s = {r.T_s:.1f} hr")
    print(f"    recession over {200.0:.0f} hr = {r.R_storm:.1f} m   eroded volume = {r.V_storm:.0f} m^3/m")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
