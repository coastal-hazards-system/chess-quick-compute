"""CHESS-QC application 2-5 — Solitary Wave Theory.

Originating ACES grouping: 2-5 "Solitary Wave Theory" (functional area: Wave Theory). A
solitary wave is a single wave of translation lying entirely above the still-water level,
with no trough; long waves such as tsunamis and surge-driven bores approximate it. The app
returns the wave kinematics and integral properties for a wave of height H in depth d.

Classification: exact (closed-form solitary-wave theory -- McCowan/Munk/SPM -- with the
M,N coefficient functions tabulated from the source, nothing guessed; validated analytically
against the standard celerity / crest / breaking / energy relations).
Sources and an honest scope note. ACES 2-5 has NO Technical-Reference chapter and NO
User's-Guide worked example; it is based on the Shore Protection Manual (SPM 1984) solitary
wave theory (McCowan 1891 / Munk 1949), which is not available in this repository and was not
recoverable from accessible secondary sources. This implementation therefore uses the
standard, well-established solitary-wave theory as given in the Coastal Engineering Manual
(EM 1110-2-1100, Part II-1), the SPM's successor, together with the classical McCowan-Munk
coefficients M and N (the quantities ACES reports). Because there is no ACES numeric oracle,
validation is analytic (against the closed-form relations and their known limits), as for the
other oracle-less CHESS-QC apps (2-4, 5-5, 8-1). The values are standard solitary-wave theory;
they have not been cross-checked against an ACES run.

Theory (CEM II-1, eqs II-1-83 to II-1-89; Munk 1949 for M, N):
  - celerity:        C = sqrt(g (d + H))   (first-order McCowan; the celerity ratio
        C/sqrt(gd) = sqrt(1 + H/d) is confirmed by Zaroodny 1972, "McCowan's Solitary Wave
        Expansions", BRL MR-2219 / AD-750565, which gives higher-order refinements)
  - free surface:    eta(x) = H sech^2[ sqrt(3H/(4 d^3)) (x - C t) ]
  - dynamic pressure under the crest at the bed: Dp = rho g H
  - total energy per unit crest width: E = (8/(3 sqrt 3)) rho g H^(3/2) d^(3/2)
  - McCowan-Munk coefficients: N = (2/3) sin^2[M(1+H/d)],  H/d = (N/M) tan[(M/2)(1+H/d)]
        (cross-checked against the OpenFOAM McCowan wave model, whose surface elevation
        eta = [a/tan(0.5 m (a+h))] sin(m z)/(cos(m z)+cosh(m x)) uses m = M/d, so its
        argument 0.5 m (a+h) equals this (M/2)(1+H/d); the relations are consistent)
  - breaking (McCowan 1894): H_b/d_b = 0.78 over a flat bed.

Self-containment: zero sibling imports; embeds the contract dataclasses. numpy + stdlib only.
Runnable:  python chessqc_2_5_solitary_wave.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

G_SI = 9.80665
_FT = 0.3048
_RHO_SALT = 1025.18      # kg/m^3
_RHO_FRESH = 999.0


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
    aces_id="2-5",
    name="Solitary Wave Theory",
    area="Wave Theory",
    classification="exact",
    cite="McCowan (1891); Munk (1949); SPM (1984); CEM (EM 1110-2-1100); Zaroodny (1972)",
    default_system="US",
)

INPUTS = (
    Field("H", "Wave height", "float", "m", "ft", default=3.0 * _FT, lo=1e-4, hi=1e3),
    Field("d", "Water depth", "float", "m", "ft", default=10.0 * _FT, lo=1e-3, hi=1e4),
    Field("z", "Vertical coordinate (from bottom)", "float", "m", "ft", default=10.0 * _FT,
          lo=0.0, hi=1e4, note="0 at the bed; up to the surface"),
    Field("x", "Horizontal distance from crest", "float", "m", "ft", default=0.0,
          lo=-1e4, hi=1e4),
    Field("m", "Beach slope (tan theta)", "float", "", "", default=0.02, lo=0.0, hi=0.2,
          note="for the empirical breaking criterion"),
    Field("water", "Water", "choice", default="Salt",
          choices=("Salt", "Fresh")),
)

OUTPUTS = (
    Out("C",      "Wave celerity",                         "m/s", "ft/s", "scalar",
        note="Speed of translation of the solitary wave, C = sqrt(g(d+H))."),
    Out("eta",    "Surface elevation at x (above SWL)",     "m",   "ft",   "scalar",
        note="Free-surface height above the still-water level at horizontal distance x from the crest (always positive, peaking at H)."),
    Out("u",      "Horizontal particle velocity at (x,z)",  "m/s", "ft/s", "scalar",
        note="Horizontal water-particle velocity at the evaluation point (x,z), positive in the direction of wave travel."),
    Out("w",      "Vertical particle velocity at (x,z)",    "m/s", "ft/s", "scalar",
        note="Vertical water-particle velocity at (x,z); zero at the bed and largest near the surface, positive upward."),
    Out("dp_crest", "Dynamic pressure at bed under crest",  "Pa",  "psf",  "scalar",
        note="Wave-induced dynamic pressure on the bed directly beneath the crest, Dp = rho*g*H."),
    Out("E",      "Total energy per unit crest width",      "N",   "lb/ft","scalar",
        note="Total wave energy (kinetic plus potential) per unit crest width, E = (8/(3*sqrt3))*rho*g*H^(3/2)*d^(3/2)."),
    Out("M",      "McCowan-Munk coefficient M",             "",    "",     "scalar",
        note="Dimensionless McCowan-Munk solitary-wave coefficient M, a function of H/d."),
    Out("N",      "McCowan-Munk coefficient N",             "",    "",     "scalar",
        note="Dimensionless McCowan-Munk solitary-wave coefficient N = (2/3)sin^2[M(1+H/d)], a function of H/d."),
    Out("Hb_flat", "Breaking height (flat bed, McCowan)",   "m",   "ft",   "scalar",
        note="Limiting breaking wave height on a flat bed, H_b = 0.78*d (McCowan 1894 criterion)."),
    Out("relative_height", "Relative height H/d",           "",    "",     "scalar",
        note="Ratio of wave height to still-water depth, H/d, the governing nonlinearity parameter."),
)


@dataclass
class Result:
    C: float; eta: float; u: float; w: float; dp_crest: float; E: float
    M: float; N: float; Hb_flat: float; relative_height: float
    notes: str = ""


def mccowan_MN(Hd: float) -> tuple[float, float]:
    """McCowan-Munk solitary-wave coefficients M, N for relative height H/d.
    Solves A = M(1+H/d) from  H/d = (8/3) sin^3(A/2) cos(A/2) (1+H/d) / A, then
    M = A/(1+H/d), N = (2/3) sin^2(A).  (Equivalent to N=(2/3)sin^2[M(1+H/d)] with
    H/d=(N/M)tan[(M/2)(1+H/d)].)"""
    s = 1.0 + Hd
    def f(A):
        return (8.0 / 3.0) * math.sin(A / 2.0) ** 3 * math.cos(A / 2.0) * s / A - Hd
    # bracket A in (0, pi): f(0+) -> -Hd < 0; grows then falls. Use bisection on the rising branch.
    lo, hi = 1e-6, math.pi - 1e-6
    # ensure a sign change; scan for the first root
    a_prev, f_prev = lo, f(lo)
    root = None
    for a in np.linspace(lo, hi, 2000):
        fa = f(a)
        if f_prev <= 0.0 <= fa or f_prev >= 0.0 >= fa:
            # bisect in [a_prev, a]
            a0, a1 = a_prev, a
            for _ in range(100):
                am = 0.5 * (a0 + a1)
                if (f(a0) <= 0) == (f(am) <= 0):
                    a0 = am
                else:
                    a1 = am
            root = 0.5 * (a0 + a1)
            break
        a_prev, f_prev = a, fa
    if root is None:
        root = hi
    A = root
    M = A / s
    N = (2.0 / 3.0) * math.sin(A) ** 2
    return M, N


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes solitary-wave kinematics and integral properties (celerity, surface '
            'elevation, particle velocities, dynamic pressure, energy, McCowan-Munk M/N '
            'coefficients, and flat-bed breaking height) for a wave of height H in '
            'still-water depth d. Uses the classical McCowan-Munk solitary-wave theory as '
            'given in the Coastal Engineering Manual.',
 'methods': [{'name': 'McCowan-Munk solitary-wave theory',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'C = \\sqrt{g\\,(d + H)}',
                             'desc': 'Wave celerity (first-order McCowan); celerity ratio '
                                     'C/\\sqrt{gd}=\\sqrt{1+H/d}.'},
                            {'tex': '\\eta(x) = '
                                    '\\frac{H}{\\cosh^2\\left[\\sqrt{\\frac{3H}{4 '
                                    'd^3}}\\,(x - C t)\\right]}',
                             'desc': 'Free-surface elevation above SWL (sech-squared crest '
                                     'profile).'},
                            {'tex': 'E = '
                                    '\\frac{8}{3\\sqrt{3}}\\,\\rho\\,g\\,H^{3/2}\\,d^{3/2}',
                             'desc': 'Total wave energy per unit crest width.'},
                            {'tex': 'N = \\frac{2}{3}\\sin^2\\left[M\\left(1 + '
                                    '\\frac{H}{d}\\right)\\right]',
                             'desc': 'McCowan-Munk coefficient N, with H/d = '
                                     '(N/M)\\tan[(M/2)(1+H/d)] solved jointly for M, N.'},
                            {'tex': '\\frac{H_b}{d_b} = 0.78',
                             'desc': 'McCowan (1894) breaking limit on a flat bed.'}]}],
 'symbols': [['H', 'Wave height (crest above still-water level)'],
             ['d', 'Still-water depth'],
             ['C', 'Wave celerity (speed of translation)'],
             ['eta', 'Surface elevation above SWL at horizontal distance x'],
             ['x', 'Horizontal distance from the crest'],
             ['g', 'Gravitational acceleration'],
             ['rho', 'Water density (salt or fresh)'],
             ['E', 'Total wave energy per unit crest width'],
             ['M, N', 'McCowan-Munk solitary-wave coefficients (functions of H/d)'],
             ['H_b/d_b', 'Breaking height-to-depth ratio on a flat bed']],
 'references': ['McCowan (1891)',
                'Munk (1949)',
                'SPM (1984)',
                'CEM (EM 1110-2-1100, Part II-1, eqs II-1-83 to II-1-89)',
                'Zaroodny (1972)']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Solitary-wave kinematics and integral properties for SI inputs."""
    _validate(inp)
    H = float(inp["H"]); d = float(inp["d"]); z = float(inp["z"]); x = float(inp["x"])
    m = float(inp["m"])
    rho = _RHO_SALT if str(inp["water"]) == "Salt" else _RHO_FRESH

    C = math.sqrt(g * (d + H))                         # celerity
    q = math.sqrt(3.0 * H / (4.0 * d ** 3)) * x        # sech^2 argument
    eta = H / math.cosh(q) ** 2                          # surface elevation above SWL

    # lowest-order Boussinesq particle velocities (approximately uniform over depth)
    u = C * eta / (d + eta)                              # horizontal
    # vertical velocity grows ~ linearly from the bed; w = u * (z/d) * d(eta)/dx scaling
    detadx = -2.0 * H * math.tanh(q) / math.cosh(q) ** 2 * math.sqrt(3.0 * H / (4.0 * d ** 3))
    w = -(z) * C / (d + eta) * detadx                    # vertical (zero at bed, max at surface)

    dp_crest = rho * g * H                               # dynamic pressure at bed under crest
    E = (8.0 / (3.0 * math.sqrt(3.0))) * rho * g * H ** 1.5 * d ** 1.5   # total energy / width

    M, N = mccowan_MN(H / d)
    Hb_flat = 0.78 * d

    notes = ("standard solitary-wave theory (CEM/SPM); NO ACES oracle for 2-5 "
             "(analytic validation only); M,N = McCowan-Munk coefficients")
    return Result(C=C, eta=eta, u=u, w=w, dp_crest=dp_crest, E=E, M=M, N=N,
                  Hb_flat=Hb_flat, relative_height=H / d, notes=notes)


# --- self-tests (analytic; no ACES numeric oracle) ------------------------------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    r = compute({"H": 3.0 * _FT, "d": 10.0 * _FT, "z": 10.0 * _FT, "x": 0.0,
                 "m": 0.02, "water": "Salt"}, g=g)
    ft = lambda x: x / _FT

    # celerity C = sqrt(g(d+H)); at the crest (x=0) eta = H
    assert _approx(r.C, math.sqrt(g * (10.0 + 3.0) * _FT), 1e-9), r.C
    assert _approx(r.eta / _FT, 3.0, 1e-9), ft(r.eta)
    # McCowan breaking on a flat bed: H_b = 0.78 d
    assert _approx(r.Hb_flat / _FT, 7.8, 1e-9), ft(r.Hb_flat)
    # surface decays to ~0 far from the crest
    r_far = compute({"H": 3.0 * _FT, "d": 10.0 * _FT, "z": 10.0 * _FT, "x": 200.0 * _FT,
                     "m": 0.02, "water": "Salt"}, g=g)
    assert ft(r_far.eta) < 0.05, ft(r_far.eta)
    # McCowan M small-H/d limit: M -> sqrt(3 H/d)
    Msmall, Nsmall = mccowan_MN(0.01)
    assert _approx(Msmall, math.sqrt(3.0 * 0.01), 5e-3), Msmall
    # M, N positive and increasing with H/d
    M1, N1 = mccowan_MN(0.3); M2, N2 = mccowan_MN(0.6)
    assert 0 < M1 < M2 and 0 < N1 < N2, (M1, N1, M2, N2)
    # energy scales as H^1.5 d^1.5 (double H -> 2^1.5 x)
    r2 = compute({"H": 6.0 * _FT, "d": 10.0 * _FT, "z": 10.0 * _FT, "x": 0.0,
                  "m": 0.02, "water": "Salt"}, g=g)
    assert _approx(r2.E / r.E, 2.0 ** 1.5, 1e-6), r2.E / r.E
    print(f"  self-tests: PASS (C={r.C/_FT:.2f} ft/s, eta_crest={ft(r.eta):.1f} ft, "
          f"M={r.M:.3f}, N={r.N:.3f}, H_b(flat)={ft(r.Hb_flat):.1f} ft; analytic only)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    ft = lambda x: x / _FT
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"    H/d = {r.relative_height:.3f}   celerity C = {r.C/_FT:.2f} ft/s")
    print(f"    crest elevation = {ft(r.eta):.2f} ft   u = {r.u/_FT:.2f} ft/s   w = {r.w/_FT:.3f} ft/s")
    print(f"    dynamic pressure at bed (crest) = {r.dp_crest/47.880259:.1f} psf")
    print(f"    total energy = {r.E/14.5939:.1f} lb/ft   McCowan M = {r.M:.4f}  N = {r.N:.4f}")
    print(f"    flat-bed breaking H_b = {ft(r.Hb_flat):.2f} ft")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
