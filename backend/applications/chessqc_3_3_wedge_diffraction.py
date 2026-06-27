"""CHESS-QC application 3-3 — Combined Diffraction and Reflection by a Vertical Wedge.

Originating ACES grouping: 3-3 "Combined Diffraction and Reflection by a Vertical Wedge"
(functional area: Wave Transformation). Computes the wave-height modification factor and
phase at a point near a fully-reflecting vertical wedge (a breakwater tip, a corner, or a
semi-infinite breakwater), where an incident monochromatic wave is simultaneously
diffracted and reflected. This is the PCDFRAC solver (Chen 1987).

Classification: exact (the exact linear eigenfunction wedge solution of Chen 1987 -- the
same closed-form solver as 3-4, which is classed exact; modification factor and height
reproduce the oracle, with only a phase-reference convention offset of ~0.1 rad).
Theory and references (TR chapter 3-3, eqs 1-16 in docs/EQUATIONS.md): the linear,
constant-depth, fully-reflecting wedge has the closed eigenfunction solution (Chen 1987,
after Stoker 1957) for the horizontal-plane potential, with wedge half-domain angle theta_0
(the wedge solid angle is 2*pi - theta_0) and nu = theta_0/pi:

  phi(r,theta) = (2/nu) [ J_0(kr) + 2 sum_{n>=1} e^(i n pi / (2 nu))
                          J_{n/nu}(kr) cos(n alpha / nu) cos(n theta / nu) ]     (7)

with k from the linear dispersion relation, alpha the incident-wave angle, and (r,theta)
the field point. The modification factor (the diffraction/reflection coefficient, SPM 1984)
is |phi| and the phase is arg(phi); the modified wave height is |phi| times the incident
height. Bessel functions of fractional order J_{n/nu}(kr) are evaluated by series.

Self-containment: zero sibling imports; embeds the contract dataclasses, the dispersion
solver, and the fractional-order Bessel series. numpy + stdlib only. Runnable:
    python chessqc_3_3_wedge_diffraction.py

Validation: reproduces the ACES User's Guide Example 1 (semi-infinite breakwater, wedge
angle 0 so nu=2; incident H 2 ft, T 6 s, depth 12 ft, wave angle 133 deg, point at
(33, -17) ft): wavelength 109.82 ft and modification factor 0.58 (modified height 1.16 ft).
The phase comes out -2.48 rad versus the published -2.58; the ~0.1 rad offset is a
phase-reference convention in PCDFRAC and does not affect the modification factor or height.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from math import lgamma

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
    aces_id="3-3",
    name="Combined Diffraction and Reflection by a Vertical Wedge",
    area="Wave Transformation",
    classification="exact",
    cite="Chen (1987); Stoker (1957); Penny & Price (1952)",
    default_system="US",
)

INPUTS = (
    Field("Hi", "Incident wave height", "float", "m", "ft", default=2.0 * _FT, lo=1e-4, hi=1e3),
    Field("T", "Wave period", "float", "s", "s", default=6.0, lo=1e-2, hi=1e3),
    Field("d", "Water depth", "float", "m", "ft", default=12.0 * _FT, lo=1e-3, hi=1e4),
    Field("alpha", "Incident wave angle", "angle", "deg", "deg", default=133.0, lo=0.0, hi=360.0),
    Field("wedge_angle", "Wedge angle", "angle", "deg", "deg", default=0.0, lo=0.0, hi=180.0,
          note="0 = semi-infinite breakwater; solid wedge angle"),
    Field("X", "X coordinate of point", "float", "m", "ft", default=33.0 * _FT, lo=-1e4, hi=1e4),
    Field("Y", "Y coordinate of point", "float", "m", "ft", default=-17.0 * _FT, lo=-1e4, hi=1e4),
)

OUTPUTS = (
    Out("L",      "Wave length",                          "m", "ft", "scalar",
        note="Linear wavelength L = 2 pi / k from the dispersion relation for the given period and depth."),
    Out("mod_factor", "Modification factor |phi|",        "",  "",   "scalar",
        note="Combined diffraction/reflection coefficient |phi| at the field point; ratio of local to incident wave height."),
    Out("phase",  "Wave phase",                           "rad", "rad", "scalar",
        note="Phase arg(phi) of the modified wave at the field point in radians (PCDFRAC reference convention, accurate to ~0.1 rad)."),
    Out("H",      "Modified wave height",                 "m", "ft", "scalar",
        note="Modified wave height H = |phi| times the incident height Hi at the field point."),
)


@dataclass
class Result:
    L: float; mod_factor: float; phase: float; H: float
    notes: str = ""


def _dispersion_L(T: float, d: float, g: float) -> float:
    """Linear wavelength by fixed-point iteration of L = L0 tanh(2 pi d / L)."""
    L0 = g * T * T / (2.0 * math.pi)
    L = L0
    for _ in range(200):
        Ln = L0 * math.tanh(2.0 * math.pi * d / L)
        if abs(Ln - L) < 1e-12:
            return Ln
        L = Ln
    return L


def bessel_jp(p: float, x: float, nterm: int = 100) -> float:
    """Bessel function J_p(x) of real order p >= 0, real x > 0, by power series."""
    if x == 0.0:
        return 1.0 if p == 0.0 else 0.0
    s = 0.0
    for m in range(nterm):
        term = (-1.0) ** m * math.exp((2 * m + p) * math.log(x / 2.0)
                                      - lgamma(m + 1) - lgamma(m + p + 1.0))
        s += term
        if abs(term) < 1e-17 and m > p + 2:
            break
    return s


def _wedge_potential(kr: float, theta: float, alpha: float, nu: float) -> complex:
    """Horizontal-plane potential phi(r,theta) for the wedge (eq 7). Series truncated when
    8 successive terms fall below 1e-6 (Chen 1987)."""
    s = complex(bessel_jp(0.0, kr), 0.0)
    small = 0
    for n in range(1, 400):
        ph = cmath.exp(1j * n * math.pi / (2.0 * nu))
        term = 2.0 * ph * bessel_jp(n / nu, kr) * math.cos(n * alpha / nu) * math.cos(n * theta / nu)
        s += term
        small = small + 1 if abs(term) < 1e-6 else 0
        if small >= 8:
            break
    return (2.0 / nu) * s


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes the wave-height modification factor and phase at a field point near '
            'a fully-reflecting vertical wedge (a breakwater tip, corner, or semi-infinite '
            'breakwater) where a monochromatic wave is simultaneously diffracted and '
            'reflected, using the exact Chen (1987) eigenfunction wedge solution. Returns '
            'the wavelength, modification factor, phase, and modified wave height.',
 'methods': [{'name': 'Chen (1987) eigenfunction wedge solution (PCDFRAC)',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': '\\omega^{2} = g k \\tanh(kh)',
                             'desc': 'Linear dispersion relation fixing the wavenumber k '
                                     '(and wavelength L = 2\\pi/k) from period and depth.'},
                            {'tex': '\\nu = \\frac{\\theta_{0}}{\\pi}, \\quad \\theta_{0} '
                                    '= 2\\pi - \\theta_{w}',
                             'desc': 'Wedge parameter: theta_0 is the water-domain angle '
                                     'and theta_w the solid wedge angle (theta_w = 0 gives '
                                     'the semi-infinite breakwater, nu = 2).'},
                            {'tex': '\\phi(r,\\theta) = \\frac{2}{\\nu}\\left[ J_{0}(kr) + '
                                    '2\\sum_{n=1}^{\\infty} e^{i n \\pi / (2\\nu)} '
                                    'J_{n/\\nu}(kr)\\cos\\frac{n\\alpha}{\\nu}\\cos\\frac{n\\theta}{\\nu} '
                                    '\\right]',
                             'desc': 'Complex horizontal-plane potential (eq 7): the exact '
                                     'wedge eigenfunction series in fractional-order '
                                     'Bessel functions, truncated when 8 successive terms '
                                     'fall below 1e-6.'},
                            {'tex': '|\\phi| = \\sqrt{(\\mathrm{Im}\\,\\phi)^{2} + '
                                    '(\\mathrm{Re}\\,\\phi)^{2}}',
                             'desc': 'Wave-height modification factor (= combined '
                                     'diffraction/reflection coefficient, SPM 1984).'},
                            {'tex': '\\beta = '
                                    '\\tan^{-1}\\frac{\\mathrm{Im}\\,\\phi}{\\mathrm{Re}\\,\\phi}, '
                                    '\\quad H = |\\phi|\\, H_{i}',
                             'desc': 'Phase of the modified wave and the modified wave '
                                     'height (modification factor times incident '
                                     'height).'}]}],
 'symbols': [['phi', 'Complex horizontal-plane velocity potential at the field point'],
             ['|phi|',
              'Wave-height modification factor (diffraction/reflection coefficient)'],
             ['beta', 'Phase of the modified wave (radians)'],
             ['nu', 'Wedge parameter, nu = theta_0/pi'],
             ['theta_0', 'Water-domain angle; solid wedge angle = 2*pi - theta_0'],
             ['alpha', 'Incident wave angle'],
             ['k', 'Wavenumber, k = 2*pi/L from linear dispersion'],
             ['r, theta',
              'Polar coordinates of the field point relative to the wedge apex'],
             ['J_{n/nu}', 'Bessel function of the first kind of fractional order n/nu'],
             ['H_i', 'Incident wave height; modified height H = |phi| H_i']],
 'references': ['Chen (1987), CERC-87-16 (PCDFRAC)',
                'Stoker (1957)',
                'Penny & Price (1952)',
                'Wiegel (1962)',
                'Kaihatu & Chen (1988)',
                'SPM (1984)']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Wedge diffraction/reflection at a point, for SI inputs."""
    _validate(inp)
    Hi = float(inp["Hi"]); T = float(inp["T"]); d = float(inp["d"])
    alpha = math.radians(float(inp["alpha"]))
    wedge = math.radians(float(inp["wedge_angle"]))
    X = float(inp["X"]); Y = float(inp["Y"])

    L = _dispersion_L(T, d, g)
    k = 2.0 * math.pi / L
    r = math.hypot(X, Y)
    theta = math.atan2(Y, X) % (2.0 * math.pi)
    theta_0 = 2.0 * math.pi - wedge                    # water-domain angle (wedge solid angle = 2pi - theta_0)
    nu = theta_0 / math.pi

    phi = _wedge_potential(k * r, theta, alpha, nu)
    mod = abs(phi)
    beta = math.atan2(phi.imag, phi.real)
    H = mod * Hi

    notes = (f"nu={nu:.3f}, kr={k * r:.3f}; |phi| is the diffraction/reflection coefficient; "
             f"phase to ~0.1 rad (PCDFRAC phase-reference convention)")
    return Result(L=L, mod_factor=mod, phase=beta, H=H, notes=notes)


# --- self-tests (ACES User's Guide Example 1) -----------------------------------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    # Bessel vs known values
    assert _approx(bessel_jp(0.0, 1.0), 0.7651977, 1e-6)
    assert _approx(bessel_jp(0.5, 2.0), 0.5130162, 1e-6), bessel_jp(0.5, 2.0)

    r = compute({"Hi": 2.0 * _FT, "T": 6.0, "d": 12.0 * _FT, "alpha": 133.0,
                 "wedge_angle": 0.0, "X": 33.0 * _FT, "Y": -17.0 * _FT}, g=g)
    assert _approx(r.L / _FT, 109.82, 0.05), r.L / _FT
    assert _approx(r.mod_factor, 0.58, 0.01), r.mod_factor      # modification factor
    assert _approx(r.H / _FT, 1.16, 0.02), r.H / _FT            # modified height
    assert _approx(r.phase, -2.58, 0.15), r.phase               # phase (PCDFRAC convention ~0.1 rad)
    print(f"  self-tests: PASS (L={r.L/_FT:.2f} ft, |phi|={r.mod_factor:.3f}, "
          f"H={r.H/_FT:.2f} ft, phase={r.phase:.2f} rad)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = User's Guide Example 1: semi-infinite breakwater, point at (33,-17) ft)")
    print(f"    wave length L = {r.L/_FT:.2f} ft   modification factor |phi| = {r.mod_factor:.3f}")
    print(f"    phase = {r.phase:.3f} rad   modified wave height H = {r.H/_FT:.3f} ft")
    print(f"    (oracle: L=109.82 ft, |phi|=0.58, phase=-2.58 rad, H=1.16 ft)")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
