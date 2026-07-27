"""CHESS-QC application 2-5 — Solitary Wave Theory.

Originating ACES grouping: 2-5 "Solitary Wave Theory" (functional area: Wave Theory). A
solitary wave is a single wave of translation lying entirely above the still-water level,
with no trough; long waves such as tsunamis and surge-driven bores approximate it. The app
returns the wave kinematics and integral properties for a wave of height H in depth d.

Classification: exact (CEM's McCowan/Munk solitary-wave equations, with M,N read from
and linearly interpolated along the source's Figure II-1-17; validated against the
standard celerity / crest / breaking / energy relations).
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
  - McCowan-Munk coefficients M,N: linearly interpolated from CEM Figure II-1-17,
        exactly as the source directs. They are then used in the local particle-
        velocity equations CEM II-1-92 and II-1-93.
  - breaking (Miles 1980/1981): H_b/d_b = 0.78 over a flat bed; for a
        solitary wave on a 0.01 <= m <= 0.18 slope, CEM II-1-98 gives
        H_b/d_b = 0.75 + 25m - 112m^2 + 3870m^3.

Self-containment: zero sibling imports; embeds the contract dataclasses. stdlib only.
Runnable:  python chessqc_2_5_solitary_wave.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass

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
    Out("Hb_slope", "Breaking height (sloping bed)",        "m",   "ft",   "scalar",
        note="Empirical solitary-wave breaking height H_b = d_b(0.75 + 25m - 112m^2 + 3870m^3), valid for 0.01 <= m <= 0.18 (CEM II-1-98)."),
    Out("relative_height", "Relative height H/d",           "",    "",     "scalar",
        note="Ratio of wave height to still-water depth, H/d, the governing nonlinearity parameter."),
)


@dataclass
class Result:
    C: float; eta: float; u: float; w: float; dp_crest: float; E: float
    M: float; N: float; Hb_flat: float; Hb_slope: float; relative_height: float
    notes: str = ""


_MUNK_MN = (
    # H/d, M, N — digitized from CEM Figure II-1-17 (the source supplies
    # these functions graphically rather than as an algebraic correlation).
    (0.00, 0.000, 0.000), (0.01, 0.154, 0.020), (0.02, 0.222, 0.039),
    (0.03, 0.294, 0.058), (0.04, 0.329, 0.076), (0.05, 0.370, 0.113),
    (0.10, 0.529, 0.194), (0.15, 0.629, 0.267), (0.20, 0.694, 0.330),
    (0.25, 0.743, 0.391), (0.30, 0.784, 0.443), (0.35, 0.817, 0.487),
    (0.40, 0.844, 0.525), (0.45, 0.869, 0.559), (0.50, 0.888, 0.586),
    (0.55, 0.907, 0.609), (0.60, 0.924, 0.628), (0.65, 0.940, 0.645),
    (0.70, 0.955, 0.660), (0.75, 0.968, 0.675), (0.78, 0.974, 0.682),
)


def mccowan_MN(Hd: float) -> tuple[float, float]:
    """CEM Figure II-1-17 interpolation for McCowan-Munk M and N."""
    if not 0.0 <= Hd <= 0.78:
        raise ValueError("solitary-wave theory requires 0 <= H/d <= 0.78 (CEM II-1-97)")
    for (h0, m0, n0), (h1, m1, n1) in zip(_MUNK_MN, _MUNK_MN[1:]):
        if Hd <= h1:
            f = (Hd - h0) / (h1 - h0)
            return m0 + f * (m1 - m0), n0 + f * (n1 - n0)
    return _MUNK_MN[-1][1], _MUNK_MN[-1][2]


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
                             'desc': 'McCowan-Munk coefficients M and N, interpolated '
                                     'from CEM Figure II-1-17.'},
                            {'tex': 'u = C N\\,\\frac{1+\\cos(My/d)\\cosh(Mx/d)}'
                                    '{[\\cos(My/d)+\\cosh(Mx/d)]^2}',
                             'desc': 'Horizontal particle velocity (CEM II-1-92).'},
                            {'tex': 'w = C N\\,\\frac{\\sin(My/d)\\sinh(Mx/d)}'
                                    '{[\\cos(My/d)+\\cosh(Mx/d)]^2}',
                             'desc': 'Vertical particle velocity (CEM II-1-93).'},
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
    if H / d > 0.78:
        raise ValueError("H/d exceeds the CEM solitary-wave breaking limit of 0.78")

    C = math.sqrt(g * (d + H))                         # celerity
    q = math.sqrt(3.0 * H / (4.0 * d ** 3)) * x        # sech^2 argument
    eta = H / math.cosh(q) ** 2                          # surface elevation above SWL

    M, N = mccowan_MN(H / d)
    # CEM II-1-92/93 use y measured upward from the bed. Kinematics are only
    # meaningful within the water column, so points outside it are clamped to
    # the local bed/surface before evaluating the CEM expressions.
    y = min(max(z, 0.0), d + eta)
    mx = M * x / d
    my = M * y / d
    if abs(mx) > 300.0:  # stable far-field form; both components tend to zero
        u = w = 0.0
    else:
        cmy = math.cos(my)
        chx = math.cosh(mx)
        denom = (cmy + chx) ** 2
        u = C * N * (1.0 + cmy * chx) / denom
        w = C * N * math.sin(my) * math.sinh(mx) / denom

    dp_crest = rho * g * H                               # CEM II-1-85 at q=0
    E = (8.0 / (3.0 * math.sqrt(3.0))) * rho * g * H ** 1.5 * d ** 1.5   # CEM II-1-95
    Hb_flat = 0.78 * d
    # CEM II-1-98: empirical solitary-wave limit on a sloping beach.  The
    # source reports the fit for 0.01 <= m <= 0.20 and says waves did not
    # break above about 0.18; use the physically established flat-bed limit
    # outside its experimental range rather than extrapolating the polynomial.
    if 0.01 <= m <= 0.18:
        breaker_ratio = 0.75 + 25.0 * m - 112.0 * m * m + 3870.0 * m ** 3
    else:
        breaker_ratio = 0.78
    Hb_slope = breaker_ratio * d

    notes = ("CEM II-1 solitary-wave theory; NO ACES oracle for 2-5 (analytic validation only); "
             "M,N interpolated from CEM Figure II-1-17")
    if y != z:
        notes += f"; z clamped to {y:.3g} m within the local water column"
    return Result(C=C, eta=eta, u=u, w=w, dp_crest=dp_crest, E=E, M=M, N=N,
                  Hb_flat=Hb_flat, Hb_slope=Hb_slope, relative_height=H / d, notes=notes)


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
    # CEM II-1-98 is used within its stated slope range, and is not
    # extrapolated to a flat bed or slopes where CEM reports no breaking.
    ratio_m02 = 0.75 + 25.0 * 0.02 - 112.0 * 0.02 ** 2 + 3870.0 * 0.02 ** 3
    assert _approx(r.Hb_slope / r.Hb_flat, ratio_m02 / 0.78, 1e-12), r.Hb_slope
    r_flat = compute({"H": 3.0 * _FT, "d": 10.0 * _FT, "z": 10.0 * _FT, "x": 0.0,
                      "m": 0.0, "water": "Salt"}, g=g)
    assert _approx(r_flat.Hb_slope, r_flat.Hb_flat, 1e-12), r_flat.Hb_slope
    # surface decays to ~0 far from the crest
    r_far = compute({"H": 3.0 * _FT, "d": 10.0 * _FT, "z": 10.0 * _FT, "x": 200.0 * _FT,
                     "m": 0.02, "water": "Salt"}, g=g)
    assert ft(r_far.eta) < 0.05, ft(r_far.eta)
    # CEM Figure II-1-17 reading and monotonic M,N functions.
    Msmall, Nsmall = mccowan_MN(0.30)
    assert _approx(Msmall, 0.784, 1e-12) and _approx(Nsmall, 0.443, 1e-12), (Msmall, Nsmall)
    M1, N1 = mccowan_MN(0.3); M2, N2 = mccowan_MN(0.6)
    assert 0 < M1 < M2 and 0 < N1 < N2, (M1, N1, M2, N2)
    # CEM II-1-92/93: w is zero at the crest and bed, reverses sign across
    # the crest, and the horizontal velocity decays toward zero far away.
    r_off = compute({"H": 3.0 * _FT, "d": 10.0 * _FT, "z": 5.0 * _FT, "x": 4.0 * _FT,
                     "m": 0.02, "water": "Salt"}, g=g)
    r_back = compute({"H": 3.0 * _FT, "d": 10.0 * _FT, "z": 5.0 * _FT, "x": -4.0 * _FT,
                      "m": 0.02, "water": "Salt"}, g=g)
    r_bed = compute({"H": 3.0 * _FT, "d": 10.0 * _FT, "z": 0.0, "x": 4.0 * _FT,
                     "m": 0.02, "water": "Salt"}, g=g)
    assert _approx(r_off.w, -r_back.w, 1e-12) and _approx(r_bed.w, 0.0, 1e-12), (r_off.w, r_back.w, r_bed.w)
    assert r_far.u < r.u, (r_far.u, r.u)
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
    print(f"    sloping-bed breaking H_b = {ft(r.Hb_slope):.2f} ft")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
