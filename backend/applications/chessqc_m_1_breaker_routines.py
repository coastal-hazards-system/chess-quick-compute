"""CHESS-QC application M-1 — Miscellaneous Breaker and Steepness Routines.

Originating ACES grouping: the ACES Technical Reference "Miscellaneous Routines" (a set of
shared breaker and steepness utility relations). It is renumbered M-1 here to avoid a clash
with the Harbor Design application 8-1 (Properties of Rectangular Basins), and it sorts last.
These are the Weggel (1972) breaker index
and related limiting-wave criteria that applications 3-1 (Snell) and 5-5 (wave setup) use
internally; this app exposes them directly for quick checks.

Classification: standard (the McCowan / Miche / Weggel breaker and steepness relations are
known and validated, but the transcribed TR eq-5 structure breaker is non-physical, so a
substitute relationship -- the depth-consistent eq-4 inversion -- had to be chosen, a
self-made modeling choice rather than the source's).
Theory and references (TR chapter 8-1, eqs 1-5 in docs/EQUATIONS.md):
  - maximum wave steepness (Miche 1944):     H/L = 0.142 tanh(kd)                  (1)
  - breaking height, flat slope (McCowan):   H_b = 0.78 d                          (2)
  - breaking height, finite slope            H_b = H0' [0.575 m^0.031 (H0'/L0)^-0.254]  (3)
        (Singamsetti & Wind 1980)
  - breaker depth and index (Weggel 1972):   d_b = H_b / (b - a H_b/(g T^2))       (4)
        a = 43.8 (1 - e^-19.5m),  b = 1.56 / (1 + e^-19.5m)   (gravity-explicit form,
        equivalent to the TR's US-unit a=1.36, b=1.5625; consistent across unit systems)
  - breaker height at a structure of depth d_s: the depth-consistent inversion of (4),
        H_b = b d_s / (1 + a d_s/(g T^2)).

Note on TR eq 5. The TR also gives an elaborate closed form (eq 5) for the "maximum breaker
height in the structure vicinity" with coefficients (18.5m-8), (9.25m-4) and a quadratic in
P. As transcribed in the public TR it returns non-physical values (e.g. 79 ft for a 15 ft
structure depth, growing with depth), and these utility routines have no standalone ACES
worked example to debug it against. The depth-consistent inversion of eq 4 above is used
instead for the structure breaker, which gives the physically expected result.

Self-containment: zero sibling imports; embeds the contract dataclasses and the Hunt (1979)
dispersion solver. numpy + stdlib only. Runnable:
    python chessqc_m_1_breaker_routines.py

Validation: analytic / cross-check (no standalone ACES worked example; these relations are
validated through the applications that use them, 3-1 and 5-5). McCowan is exact (0.78 d);
the Miche deepwater steepness limit is 0.142; the Weggel index has the correct slope limits
(a -> 0, b -> 0.78 at m = 0); the breaker depth/height pair are mutually consistent.
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
    aces_id="M-1",
    name="Miscellaneous Breaker and Steepness Routines",
    area="Miscellaneous Routines",
    classification="standard",
    cite="Weggel (1972); Miche (1944); McCowan (1894); Singamsetti & Wind (1980)",
    default_system="US",
)

INPUTS = (
    Field("H0", "Unrefracted deepwater wave height", "float", "m", "ft", default=6.0 * _FT,
          lo=1e-4, hi=1e3),
    Field("T", "Wave period", "float", "s", "s", default=8.0, lo=1e-2, hi=1e3),
    Field("d", "Water depth (steepness / flat breaking)", "float", "m", "ft", default=15.0 * _FT,
          lo=1e-3, hi=1e4),
    Field("m", "Nearshore slope (tan theta)", "float", "", "", default=0.05, lo=0.0, hi=1.0),
    Field("ds", "Water depth at structure", "float", "m", "ft", default=15.0 * _FT,
          lo=1e-3, hi=1e4),
)

OUTPUTS = (
    Out("L",          "Wave length at depth d",            "m", "ft", "scalar",
        note="Linear-theory wavelength at local depth d for period T, from the Hunt (1979) dispersion approximation."),
    Out("steepness_max", "Maximum steepness H/L (Miche)",  "",  "",   "scalar",
        note="Maximum stable wave steepness H/L = 0.142 tanh(kd) at depth d (Miche 1944), tending to 0.142 in deep water."),
    Out("Hmax_steep", "Maximum height at steepness limit", "m", "ft", "scalar",
        note="Largest wave height the Miche steepness limit allows at depth d, the maximum steepness times the wavelength L."),
    Out("Hb_flat",    "Breaking height, flat slope (McCowan)", "m", "ft", "scalar",
        note="Breaking wave height on a flat or unknown slope, H_b = 0.78 d (McCowan 1894)."),
    Out("Hb_sloped",  "Breaking height, finite slope (S&W)", "m", "ft", "scalar",
        note="Breaking wave height on a finite nearshore slope m, H_b = H0' [0.575 m^0.031 (H0'/L0)^-0.254] (Singamsetti & Wind 1980)."),
    Out("db_sloped",  "Breaker depth (Weggel)",            "m", "ft", "scalar",
        note="Water depth at which the sloped breaking wave breaks, d_b = H_b/(b - a H_b/(g T^2)) (Weggel 1972)."),
    Out("a_index",    "Weggel breaker index a(m)",         "",  "",   "scalar",
        note="Weggel (1972) breaker-index coefficient a = 43.8(1 - e^-19.5m), a slope-dependent factor tending to 0 as m to 0."),
    Out("b_index",    "Weggel breaker index b(m)",         "",  "",   "scalar",
        note="Weggel (1972) breaker-index coefficient b = 1.56/(1 + e^-19.5m), tending to 0.78 (the McCowan ratio) as m to 0."),
    Out("Hb_structure", "Breaker height at structure depth", "m", "ft", "scalar",
        note="Breaking wave height at a structure of depth d_s, H_b = b d_s/(1 + a d_s/(g T^2)), the depth-consistent inversion of the Weggel relation."),
)


@dataclass
class Result:
    L: float; steepness_max: float; Hmax_steep: float
    Hb_flat: float; Hb_sloped: float; db_sloped: float
    a_index: float; b_index: float; Hb_structure: float
    notes: str = ""


_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_length(T: float, d: float, g: float) -> float:
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    return math.sqrt(g * d / (y + 1.0 / denom)) * T


def weggel_index(m: float) -> tuple[float, float]:
    """Weggel (1972) breaker-index coefficients a(m), b(m), gravity-explicit form
    (a = 43.8(1-e^-19.5m), b = 1.56/(1+e^-19.5m)); used with d_b = H_b/(b - a H_b/(g T^2))."""
    e = math.exp(-19.5 * m)
    return 43.8 * (1.0 - e), 1.56 / (1.0 + e)


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes a set of limiting-wave utility relations for a quick check: the '
            'Miche maximum wave steepness, the McCowan flat-slope and Singamsetti and Wind '
            'finite-slope breaking heights, and the Weggel (1972) breaker-index '
            'coefficients with the corresponding breaker depth and the breaker height at a '
            'structure.',
 'methods': [{'name': 'Breaker and steepness utility relations',
              'when': None,
              'tag': 'standard',
              'note': 'Structure breaker uses the depth-consistent inversion of the Weggel '
                      "breaker-depth relation; the TR's closed-form eq 5 returns "
                      'non-physical values and is not used.',
              'equations': [{'tex': '\\frac{H}{L} = 0.142\\,\\tanh(k d)',
                             'desc': 'Maximum wave steepness (Miche 1944), with k = '
                                     '2\\pi/L.'},
                            {'tex': 'H_{b} = 0.78\\,d',
                             'desc': 'Breaking height on a flat or unknown slope (McCowan '
                                     '1894).'},
                            {'tex': "H_{b} = H_{0}' \\left[ 0.575\\, m^{0.031} "
                                    "\\left(\\frac{H_{0}'}{L_{0}}\\right)^{-0.254} "
                                    '\\right]',
                             'desc': 'Breaking height on a finite nearshore slope m = '
                                     'tan(theta) (Singamsetti & Wind 1980).'},
                            {'tex': 'd_{b} = \\frac{H_{b}}{b - a\\,H_{b}/(g T^{2})}, '
                                    '\\quad a = 43.8(1 - e^{-19.5 m}), \\quad b = '
                                    '\\frac{1.56}{1 + e^{-19.5 m}}',
                             'desc': 'Breaker depth and breaker-index coefficients a(m), '
                                     'b(m) (Weggel 1972), gravity-explicit form.'},
                            {'tex': 'H_{b} = \\frac{b\\,d_{s}}{1 + a\\,d_{s}/(g T^{2})}',
                             'desc': 'Breaker height at a structure of depth d_s, the '
                                     'depth-consistent inversion of the Weggel '
                                     'breaker-depth relation.'}]}],
 'symbols': [['H/L', 'wave steepness (height over wavelength)'],
             ['k', 'wave number, 2 pi / L'],
             ['d', 'local water depth (steepness / flat breaking)'],
             ['H_b', 'breaking wave height'],
             ["H_0'", 'unrefracted deepwater wave height'],
             ['L_0', 'deepwater wavelength, g T^2 / (2 pi)'],
             ['m', 'nearshore slope, tan(theta)'],
             ['d_b', 'water depth at breaking (Weggel)'],
             ['d_s', 'water depth at the structure'],
             ['a, b', 'Weggel (1972) breaker-index coefficients, functions of slope m']],
 'references': ['Weggel (1972)',
                'Miche (1944)',
                'McCowan (1894)',
                'Singamsetti & Wind (1980)',
                'ACES Technical Reference ch. 8-1 (Miscellaneous Routines)']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Breaker and steepness utility relations for SI inputs."""
    _validate(inp)
    H0 = float(inp["H0"]); T = float(inp["T"]); d = float(inp["d"])
    m = float(inp["m"]); ds = float(inp["ds"])

    L0 = g * T * T / (2.0 * math.pi)
    L = wave_length(T, d, g)
    k = 2.0 * math.pi / L

    steepness_max = 0.142 * math.tanh(k * d)                 # eq 1 (Miche)
    Hmax_steep = steepness_max * L
    Hb_flat = 0.78 * d                                       # eq 2 (McCowan)
    Hb_sloped = H0 * (0.575 * m ** 0.031 * (H0 / L0) ** (-0.254)) if m > 0 else Hb_flat  # eq 3

    a_index, b_index = weggel_index(m)
    db_sloped = Hb_sloped / (b_index - a_index * Hb_sloped / (g * T * T))   # eq 4
    # breaker height at a structure of depth ds (depth-consistent inversion of eq 4)
    Hb_structure = b_index * ds / (1.0 + a_index * ds / (g * T * T))

    notes = (f"Weggel a={a_index:.2f}, b={b_index:.3f} (m={m}); McCowan H_b=0.78d; "
             f"Miche deepwater steepness limit 0.142")
    return Result(L=L, steepness_max=steepness_max, Hmax_steep=Hmax_steep,
                  Hb_flat=Hb_flat, Hb_sloped=Hb_sloped, db_sloped=db_sloped,
                  a_index=a_index, b_index=b_index, Hb_structure=Hb_structure, notes=notes)


# --- self-tests (analytic / cross-check; no standalone ACES worked example) ------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    r = compute({"H0": 6.0 * _FT, "T": 8.0, "d": 15.0 * _FT, "m": 0.05, "ds": 15.0 * _FT}, g=g)
    ft = lambda x: x / _FT

    # McCowan flat-slope breaking is exact: H_b = 0.78 d
    assert _approx(r.Hb_flat, 0.78 * 15.0 * _FT, 1e-9), ft(r.Hb_flat)
    # Miche deepwater steepness limit -> 0.142
    deep = compute({"H0": 6.0 * _FT, "T": 8.0, "d": 1000.0 * _FT, "m": 0.05, "ds": 15.0 * _FT}, g=g)
    assert _approx(deep.steepness_max, 0.142, 1e-4), deep.steepness_max
    # Weggel index slope limits: a(0)=0, b(0)=0.78; a(large)->43.8, b->0.78125
    a0, b0 = weggel_index(0.0)
    assert _approx(a0, 0.0, 1e-9) and _approx(b0, 0.78, 1e-9), (a0, b0)
    aL, bL = weggel_index(1.0)
    assert _approx(aL, 43.8, 1e-6) and _approx(bL, 1.56, 1e-3), (aL, bL)
    # breaker depth/height consistency: a wave of height Hb_sloped breaks at db_sloped,
    # i.e. inverting eq 4 at d = db recovers Hb_sloped
    a_i, b_i = weggel_index(0.05)
    Hb_back = b_i * r.db_sloped / (1.0 + a_i * r.db_sloped / (g * 8.0 ** 2))
    assert _approx(Hb_back, r.Hb_sloped, 1e-6), (ft(Hb_back), ft(r.Hb_sloped))
    # structure breaker is physical (~0.8-1.1 of d_s here, not the eq-5 non-physical value)
    assert 0.6 * 15.0 < ft(r.Hb_structure) < 1.3 * 15.0, ft(r.Hb_structure)
    print(f"  self-tests: PASS (McCowan 0.78d exact; Miche 0.142; Weggel a/b limits; "
          f"breaker depth/height consistent; H_b@structure={ft(r.Hb_structure):.2f} ft)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    ft = lambda x: x / _FT
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"    L = {ft(r.L):.2f} ft   max steepness H/L = {r.steepness_max:.4f}   "
          f"H_max(steep) = {ft(r.Hmax_steep):.2f} ft")
    print(f"    H_b flat (McCowan) = {ft(r.Hb_flat):.2f} ft   H_b sloped (S&W) = {ft(r.Hb_sloped):.2f} ft")
    print(f"    Weggel a = {r.a_index:.2f}  b = {r.b_index:.3f}   breaker depth = {ft(r.db_sloped):.2f} ft")
    print(f"    H_b at structure (d_s) = {ft(r.Hb_structure):.2f} ft")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
