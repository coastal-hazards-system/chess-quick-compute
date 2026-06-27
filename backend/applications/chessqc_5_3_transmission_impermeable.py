"""CHESS-QC application 5-3 — Wave Transmission on Impermeable Structures.

Originating ACES grouping: 5-3 "Wave Transmission on Impermeable Structures" (functional
area: Wave Runup, Transmission, and Overtopping). Estimates the height of the wave
transmitted past an impermeable structure, either a sloped structure overtopped by waves
or a vertical/composite (caisson-on-berm) structure, as a transmission coefficient
K_T = H_T / H_i.

Classification: exact (Seelig 1980 / Ahrens transmission with known published coefficients,
nothing guessed; reproduces the User's Guide Examples 1-4 to 0.01 ft).
A selectable d'Angremond, van der Meer & de Jong (1996) transmission method is also provided
for sloped/low-crested structures (the CEM VI-5-2 standard, K_t bounded 0.075-0.8); Seelig
stays the default so the Examples are reproduced. The d'Angremond path is validated for its
bounds and monotonicity (no ACES oracle exists for it).
Theory and references (TR chapter 5-3, eqs 1-7 + Tables 5-3-1/2 in docs/EQUATIONS.md):
  - sloped structures (Seelig 1980, transmission by overtopping):
        H_T = K_TO * H_i,   K_TO = C * (1 - F/R),   C = 0.51 - 0.11 (B/h_s)        (1-3)
        R = runup from the 5-2 methods (rough Ahrens & McCartney; smooth Ahrens & Titus)
  - vertical / composite structures (Seelig 1976):
        K_TO = 0.5 {1 - sin[(pi/2 alpha)(F/H_i + beta)]}, clamped to [0,1] by domain   (4-7)
        alpha = 1.8 + 0.4 min(B/d_s,1);  beta = C1 beta_1 + C2 beta_2
        beta_1 = 0.1 + 0.3 min(B/d_s,1);  beta_2 = 0.1 (d1/ds<=0.3) else 0.527-0.130/(d1/ds)
        C1 = max(0, 1 - B/d_s);  C2 = min(1, B/d_s);  d1 = d_s - berm height above toe
  F = h_s - d_s is the crest freeboard (negative for a submerged structure).

Self-containment: zero sibling imports; embeds its own contract dataclasses, the Hunt
(1979) dispersion solver, and the 5-2 runup methods. stdlib + numpy only. Runnable:
    python chessqc_5_3_transmission_impermeable.py
which runs the ACES-oracle self-tests (User's Guide Examples 1-4) and prints a tabulation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

G_SI = 9.80665
_FT = 0.3048


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
    aces_id="5-3",
    name="Wave Transmission on Impermeable Structures",
    area="Wave Runup, Transmission, and Overtopping",
    classification="exact",
    cite="Seelig (1980); Seelig (1976); Ahrens & McCartney (1975); Ahrens & Titus (1985)",
    default_system="US",
)

# Default example = ACES User's Guide Example 3 (rough slope, runup + transmission).
INPUTS = (
    Field("Hi", "Incident wave height", "float", "m", "ft", default=7.50 * _FT, lo=1e-4, hi=1e3),
    Field("T", "Wave period", "float", "s", "s", default=10.0, lo=1e-2, hi=1e3),
    Field("ds", "Water depth at structure toe", "float", "m", "ft", default=10.0 * _FT,
          lo=1e-3, hi=1e4),
    Field("hs", "Structure height above toe", "float", "m", "ft", default=15.0 * _FT,
          lo=1e-3, hi=1e4),
    Field("B", "Structure crest width", "float", "m", "ft", default=7.50 * _FT, lo=0.0, hi=1e4),
    Field("structure_type", "Structure type", "choice", default="Sloped",
          choices=("Sloped", "Vertical or composite")),
    Field("transmission_method", "Transmission method (sloped)", "choice", default="Seelig",
          choices=("Seelig", "d'Angremond"),
          note="Seelig 1980 (ACES) or d'Angremond 1996 (CEM standard); sloped structures only"),
    # --- sloped-structure inputs ---
    Field("cot_theta", "Cotangent of structure slope", "float", "", "", default=3.0,
          lo=1e-3, hi=1e3, note="sloped structures"),
    Field("slope_type", "Slope type", "choice", default="Rough (riprap)",
          choices=("Rough (riprap)", "Smooth"), note="sloped structures"),
    Field("a", "Rough-slope coefficient a", "float", "", "", default=0.956, lo=0.0, hi=10.0),
    Field("b", "Rough-slope coefficient b", "float", "", "", default=0.398, lo=0.0, hi=10.0),
    Field("R_known", "Known runup (0 = compute)", "float", "m", "ft", default=0.0, lo=0.0, hi=1e4),
    # --- vertical/composite inputs ---
    Field("berm_height", "Berm height above toe", "float", "m", "ft", default=0.0, lo=0.0, hi=1e4,
          note="vertical/composite; 0 = no berm (pure vertical wall)"),
)

OUTPUTS = (
    Out("R",    "Wave runup (sloped structures)",   "m", "ft", "scalar",
        note="Maximum vertical rise of water on the sloped structure face above still-water level, from the 5-2 runup methods (zero for vertical/composite structures)."),
    Out("F",    "Crest freeboard (hs - ds)",        "m", "ft", "scalar",
        note="Height of the structure crest above the still-water level, F = h_s - d_s; negative when the structure is submerged."),
    Out("K_TO", "Transmission coefficient",         "",  "",   "scalar",
        note="Wave transmission coefficient K_TO = H_T / H_i, the fraction of incident wave height passing the structure (0 = full blocking, 1 = full transmission)."),
    Out("H_T",  "Transmitted wave height",          "m", "ft", "scalar",
        note="Height of the wave transmitted to the lee of the structure, H_T = K_TO * H_i."),
)


@dataclass
class Result:
    R: float; F: float; K_TO: float; H_T: float
    notes: str = ""


_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_length(T: float, d: float, g: float = G_SI) -> float:
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    c2 = g * d / (y + 1.0 / denom)
    return math.sqrt(c2) * T


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


def _runup(slope_type, Hi, T, ds, cot_theta, a, b, g):
    """Runup R via the 5-2 methods (rough Ahrens & McCartney; smooth Ahrens & Titus)."""
    theta = math.atan(1.0 / cot_theta)
    L0 = g * T * T / (2.0 * math.pi)
    xi = math.tan(theta) / math.sqrt(Hi / L0)
    if slope_type.startswith("Rough"):
        return Hi * a * xi / (1.0 + b * xi)
    L = wave_length(T, ds, g)
    Pi = (Hi / L) / math.tanh(2.0 * math.pi * ds / L) ** 3
    Cp = 1.002 * xi
    Cnb = 1.087 * math.sqrt(math.pi / (2.0 * theta)) + 0.775 * Pi
    if xi <= 2.0:
        C = Cp
    elif xi >= 3.5:
        C = Cnb
    else:
        C = ((3.5 - xi) / 1.5) * Cp + ((xi - 2.0) / 1.5) * Cnb
    return C * Hi


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Estimates the wave height transmitted past an impermeable coastal structure '
            'as a transmission coefficient K_TO = H_T / H_i, for either a sloped structure '
            'overtopped by waves or a vertical/composite (caisson-on-berm) structure. '
            'Returns the runup, crest freeboard, transmission coefficient, and transmitted '
            'wave height.',
 'method_key': 'transmission_method',
 'methods': [{'name': 'Seelig (ACES default)',
              'when': 'Seelig',
              'tag': 'legacy',
              'note': "ACES default; reproduces User's Guide Examples 1-4 and is the only "
                      'method available for vertical/composite structures.',
              'equations': [{'tex': 'H_T = K_{TO}\\, H_i',
                             'desc': 'Transmitted wave height from the overtopping '
                                     'transmission coefficient (eq 1).'},
                            {'tex': 'K_{TO} = C\\left(1 - \\frac{F}{R}\\right)',
                             'desc': 'Sloped structure (Seelig 1980): F = h_s - d_s is '
                                     'crest freeboard, R is the runup from the 5-2 '
                                     'methods; clamped to [0,1] (eq 2).'},
                            {'tex': 'C = 0.51 - 0.11\\,\\frac{B}{h_s}',
                             'desc': 'Crest-width coefficient for the sloped overtopping '
                                     'coefficient (eq 3).'},
                            {'tex': 'K_{TO} = 0.5\\left[1 - '
                                    '\\sin\\!\\left(\\frac{\\pi}{2\\alpha}\\left(\\frac{F}{H_i} '
                                    '+ \\beta\\right)\\right)\\right]',
                             'desc': 'Vertical/composite structure (Seelig 1976), '
                                     'evaluated within the domain -(alpha+beta) < F/H_i < '
                                     '(alpha-beta) (eqs 4 and Table 5-3-1).'},
                            {'tex': '\\alpha = 1.8 + 0.4\\,\\mathrm{min}(B/d_s,\\,1)',
                             'desc': 'Shape coefficient for the vertical/composite formula '
                                     '(Table 5-3-2).'},
                            {'tex': '\\beta = C_1\\,\\beta_1 + C_2\\,\\beta_2',
                             'desc': 'Combined offset coefficient with C_1 = max(0, 1 - '
                                     'B/d_s), C_2 = min(1, B/d_s) (eqs 6-7).'}]},
             {'name': "d'Angremond, van der Meer & de Jong",
              'when': "d'Angremond",
              'tag': 'preferred',
              'note': 'CEM VI-5-2 standard for sloped/low-crested structures; K_t bounded '
                      'to [0.075, 0.8]. Selectable for sloped structures only.',
              'equations': [{'tex': 'H_T = K_{TO}\\, H_i',
                             'desc': 'Transmitted wave height from the transmission '
                                     'coefficient.'},
                            {'tex': 'K_{TO} = -0.4\\,\\frac{F}{H_i} + '
                                    '0.64\\left(\\frac{B}{H_i}\\right)^{-0.31}\\left(1 - '
                                    '\\exp(-0.5\\,\\xi)\\right)',
                             'desc': "d'Angremond, van der Meer & de Jong (1996) "
                                     'transmission coefficient for sloped/low-crested '
                                     'structures.'},
                            {'tex': '\\xi = \\frac{\\tan\\theta}{\\sqrt{H_i / L_0}}',
                             'desc': 'Surf-similarity (Iribarren) parameter; theta is the '
                                     'structure slope angle, L_0 = g T^2 / (2 pi) the '
                                     'deepwater wavelength.'},
                            {'tex': '0.075 \\leq K_{TO} \\leq 0.8',
                             'desc': "Validity bounds applied to the d'Angremond "
                                     'coefficient.'}]}],
 'symbols': [['K_{TO}', 'Wave transmission coefficient (= H_T / H_i)'],
             ['H_T', 'Transmitted wave height behind the structure'],
             ['H_i', 'Incident wave height'],
             ['F', 'Crest freeboard, F = h_s - d_s (negative if submerged)'],
             ['R', 'Wave runup (from the 5-2 runup methods)'],
             ['C', 'Sloped-structure crest-width coefficient'],
             ['B', 'Structure crest width'],
             ['h_s', 'Structure height above the toe'],
             ['d_s', 'Water depth at the structure toe'],
             ['xi', 'Surf-similarity (Iribarren) parameter, tan(theta)/sqrt(H_i/L_0)']],
 'references': ['Seelig (1980)',
                'Seelig (1976)',
                'Ahrens & McCartney (1975)',
                'Ahrens & Titus (1985)',
                "d'Angremond, van der Meer & de Jong (1996)",
                'Cross & Sollitt (1971)',
                'CEM VI-5-2']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Transmitted wave height for SI inputs (see INPUTS)."""
    _validate(inp)
    Hi = float(inp["Hi"]); T = float(inp["T"]); ds = float(inp["ds"])
    hs = float(inp["hs"]); B = float(inp["B"])
    structure_type = str(inp["structure_type"])
    F = hs - ds
    R = 0.0
    notes = []

    if structure_type == "Sloped":
        slope_type = str(inp["slope_type"]); cot_theta = float(inp["cot_theta"])
        a = float(inp["a"]); b = float(inp["b"]); R_known = float(inp["R_known"])
        R = R_known if R_known > 0.0 else _runup(slope_type, Hi, T, ds, cot_theta, a, b, g)
        if str(inp.get("transmission_method", "Seelig")) == "d'Angremond":
            # d'Angremond, van der Meer & de Jong (1996); CEM VI-5-2. K_t in [0.075, 0.8].
            theta = math.atan(1.0 / cot_theta)
            L0 = g * T * T / (2.0 * math.pi)
            xi = math.tan(theta) / math.sqrt(Hi / L0)
            BHi = (B / Hi) if B > 0.0 else 1e-6
            K_TO = -0.4 * (F / Hi) + 0.64 * BHi ** (-0.31) * (1.0 - math.exp(-0.5 * xi))
            K_TO = min(max(K_TO, 0.075), 0.8)
            notes.append(f"sloped (d'Angremond 1996); xi={xi:.2f}, B/Hi={B / Hi:.2f}, R={R / _FT:.2f} ft")
        else:
            C = 0.51 - 0.11 * (B / hs)
            K_TO = C * (1.0 - F / R)
            K_TO = min(max(K_TO, 0.0), 1.0)
            notes.append(f"sloped (Seelig 1980); C = {C:.3f}, R = {R / _FT:.2f} ft")
    else:
        berm = float(inp["berm_height"])
        Bds = B / ds
        d1 = ds - berm                                  # depth above berm (= ds if no berm)
        d1ds = d1 / ds
        alpha = 1.8 + 0.4 * min(Bds, 1.0)
        beta1 = 0.1 + 0.3 * min(Bds, 1.0)
        beta2 = 0.1 if d1ds <= 0.3 else (0.527 - 0.130 / d1ds)
        C1 = max(0.0, 1.0 - Bds); C2 = min(1.0, Bds)
        beta = C1 * beta1 + C2 * beta2
        FHi = F / Hi
        if FHi <= -(alpha + beta):
            K_TO = 1.0
        elif FHi >= (alpha - beta):
            K_TO = 0.0
        else:
            K_TO = 0.5 * (1.0 - math.sin((math.pi / (2.0 * alpha)) * (FHi + beta)))
        if not (0.145 <= ds / Hi <= 0.5):
            notes.append("note: ds/Hi outside Seelig (1976) validity (0.145 to 0.5)")
        notes.append(f"vertical/composite (Seelig 1976); alpha={alpha:.3f}, beta={beta:.3f}")

    H_T = K_TO * Hi
    return Result(R=R, F=F, K_TO=K_TO, H_T=H_T, notes="; ".join(notes))


# --- self-tests (ACES User's Guide Examples 1-4) --------------------------------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI

    def ft(x): return x / _FT

    # Ex 1: sloped, known runup R = 15 ft -> H_T = 2.275 ft
    r1 = compute(dict(Hi=7.50 * _FT, T=10.0, ds=10.0 * _FT, hs=15.0 * _FT, B=7.50 * _FT,
                      structure_type="Sloped", slope_type="Rough (riprap)", cot_theta=3.0,
                      a=0.956, b=0.398, R_known=15.0 * _FT, berm_height=0.0), g=g)
    assert _approx(ft(r1.H_T), 2.275, 0.01), ft(r1.H_T)

    # Ex 2: vertical wall with submerged berm -> H_T = 3.798 ft
    r2 = compute(dict(Hi=7.50 * _FT, T=4.50, ds=20.0 * _FT, hs=17.50 * _FT, B=12.0 * _FT,
                      structure_type="Vertical or composite", cot_theta=3.0,
                      slope_type="Rough (riprap)", a=0.956, b=0.398, R_known=0.0,
                      berm_height=6.0 * _FT), g=g)
    assert _approx(ft(r2.H_T), 3.798, 0.01), ft(r2.H_T)

    # Ex 3: rough slope, computed runup -> R = 9.421, H_T = 1.601 ft
    r3 = compute(dict(Hi=7.50 * _FT, T=10.0, ds=10.0 * _FT, hs=15.0 * _FT, B=7.50 * _FT,
                      structure_type="Sloped", slope_type="Rough (riprap)", cot_theta=3.0,
                      a=0.956, b=0.398, R_known=0.0, berm_height=0.0), g=g)
    assert _approx(ft(r3.R), 9.421, 0.01) and _approx(ft(r3.H_T), 1.601, 0.01), (ft(r3.R), ft(r3.H_T))

    # Ex 4: smooth slope, computed runup -> R = 22.436, H_T = 2.652 ft
    r4 = compute(dict(Hi=7.50 * _FT, T=10.0, ds=10.0 * _FT, hs=15.0 * _FT, B=7.50 * _FT,
                      structure_type="Sloped", slope_type="Smooth", cot_theta=3.0,
                      a=0.956, b=0.398, R_known=0.0, berm_height=0.0), g=g)
    assert _approx(ft(r4.R), 22.436, 0.02) and _approx(ft(r4.H_T), 2.652, 0.01), (ft(r4.R), ft(r4.H_T))

    # d'Angremond (1996): selectable for sloped structures; bounds + monotonicity
    da = compute(dict(Hi=7.50 * _FT, T=10.0, ds=10.0 * _FT, hs=15.0 * _FT, B=7.50 * _FT,
                      structure_type="Sloped", slope_type="Rough (riprap)", cot_theta=3.0,
                      a=0.956, b=0.398, R_known=0.0, berm_height=0.0,
                      transmission_method="d'Angremond"), g=g)
    assert 0.075 <= da.K_TO <= 0.8 and da.H_T > 0, da.K_TO
    da_wide = compute(dict(Hi=7.50 * _FT, T=10.0, ds=10.0 * _FT, hs=15.0 * _FT, B=25.0 * _FT,
                           structure_type="Sloped", slope_type="Rough (riprap)", cot_theta=3.0,
                           a=0.956, b=0.398, R_known=0.0, berm_height=0.0,
                           transmission_method="d'Angremond"), g=g)
    assert da_wide.K_TO < da.K_TO, (da.K_TO, da_wide.K_TO)   # wider crest -> more attenuation

    print("  self-tests: PASS (ACES Examples 1-4; d'Angremond transmission consistent)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = User's Guide Example 3: rough slope, runup + transmission)")
    print(f"    runup R = {r.R/_FT:.3f} ft   freeboard F = {r.F/_FT:.2f} ft")
    print(f"    transmission coeff K_TO = {r.K_TO:.4f}   transmitted height H_T = {r.H_T/_FT:.3f} ft")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
