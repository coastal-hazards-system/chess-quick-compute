"""CHESS-QC application 5-2 — Wave Runup and Overtopping on Impermeable Structures.

Originating ACES grouping: 5-2 "Wave Runup and Overtopping on Impermeable Structures"
(functional area: Wave Runup, Transmission, and Overtopping). Computes wave runup on a
smooth or rough impermeable slope (a seawall or revetment face) and, when the runup
exceeds the crest, the resulting overtopping rate, for both monochromatic and irregular
waves, with an optional onshore-wind correction.

Classification: exact (Ahrens & McCartney 1975 / Ahrens & Titus 1985 runup + Weggel 1976
overtopping; the empirical coefficients a, b, alpha, Q*0 are supplied as known user inputs,
nothing guessed; reproduces the User's Guide Examples 1-7).
A selectable EurOtop (2018) mean-overtopping method is also provided (the modern dimensionless
standard, eqs 5.10/5.11, with a roughness factor gamma_f). Weggel stays the default so the
Examples are reproduced; the EurOtop path is validated for physical limits/monotonicity (no
ACES oracle exists for it).
Theory and references (TR chapter 5-2, eqs 1-15 in docs/EQUATIONS.md):
  - rough-slope runup:  R = H_i*a*xi/(1 + b*xi)            Ahrens & McCartney (1975)   (1)
  - smooth-slope runup: R = C*H_i, C by surf regime         Ahrens & Titus (1985)      (3-8)
        plunging   (xi<=2):   C_p  = 1.002*xi
        nonbreaking(xi>=3.5): C_nb = 1.087*sqrt(pi/(2 theta)) + 0.775*Pi     (Ahrens & Burke)
        transition (2<xi<3.5): linear blend of C_p and C_nb
        Pi = (H_i/L)/tanh^3(2 pi d_s/L)                       Goda (1983)               (8)
  - surf-similarity:    xi = tan(theta)/sqrt(H_i/L_0)                                   (2)
  - overtopping:        Q = C_w*sqrt(g*Q*0*H'0^3)*((R+F)/(R-F))^(-0.1085/alpha)  Weggel (1976) (9)
  - wind:               C_w = 1 + W_f*(F/R + 0.1)*sin(theta), W_f = U^2/1800 (U in mph)  (10-11)
  - irregular:          Q = (1/199) sum_i Q_i, R_i = sqrt(ln(1/p_i)/2)*R_s, p_i = 0.005 i  (12-14)

H_i is the incident height at the structure toe; it is deshoaled to the deepwater height
H_0 (= H_i / K_s, K_s the linear shoaling coefficient) for the overtopping term and report.
F = h_s - d_s is the crest freeboard.  The empirical coefficients a, b (rough runup) and
alpha, Q*0 (overtopping) are user inputs read from the SPM (1984) / Appendix A figures; an
optional average alpha = 0.06 - 0.01431 sin(theta) is available.

Self-containment: zero sibling imports; embeds its own Field/AppMeta/Out/Result dataclasses
and the Hunt (1979) dispersion solver. stdlib + numpy only. Runnable standalone:
    python chessqc_5_2_runup_overtopping_impermeable.py
which runs the ACES-oracle self-tests (User's Guide Examples 1-7) and prints a tabulation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

G_SI = 9.80665           # m/s^2
_FT = 0.3048
_MPH = 0.44704           # m/s per mph
_KN = 0.514444           # m/s per knot


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
    enable_if: tuple = ()    # (other_key, value): gray out (disable) unless that input == value


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    note: str = ""           # hover definition shown on the output label


APP_META = AppMeta(
    aces_id="5-2",
    name="Wave Runup and Overtopping on Impermeable Structures",
    area="Wave Runup, Transmission, and Overtopping",
    classification="exact",
    cite="Ahrens & McCartney (1975); Ahrens & Titus (1985); Weggel (1976); SPM (1984)",
    default_system="US",
)

# Default example = ACES User's Guide Example 5 (monochromatic, rough slope, runup + overtopping).
INPUTS = (
    Field("Hi", "Incident wave height at toe", "float", "m", "ft", default=7.50 * _FT,
          lo=1e-4, hi=1e3, note="significant height for irregular waves"),
    Field("T", "Wave period", "float", "s", "s", default=10.0, lo=1e-2, hi=1e3),
    Field("ds", "Water depth at structure toe", "float", "m", "ft", default=12.50 * _FT,
          lo=1e-3, hi=1e4),
    Field("cot_theta", "Cotangent of structure slope", "float", "", "", default=3.0,
          lo=1e-3, hi=1e3),
    Field("hs", "Structure height above toe", "float", "m", "ft", default=20.0 * _FT,
          lo=1e-3, hi=1e4),
    Field("wave_type", "Wave type", "choice", default="Monochromatic",
          choices=("Monochromatic", "Irregular")),
    Field("slope_type", "Slope type", "choice", default="Rough (riprap)",
          choices=("Rough (riprap)", "Smooth")),
    Field("want_overtopping", "Compute overtopping", "bool", default=True),
    Field("R_known", "Known runup (0 = compute)", "float", "m", "ft", default=0.0,
          lo=0.0, hi=1e4, note="if > 0, used directly instead of the runup formula"),
    Field("a", "Rough-slope coefficient a", "float", "", "", default=0.956, lo=0.0, hi=10.0,
          note="Ahrens & McCartney; per armor type (Appendix A)"),
    Field("b", "Rough-slope coefficient b", "float", "", "", default=0.398, lo=0.0, hi=10.0),
    Field("alpha", "Overtopping coefficient alpha", "float", "", "", default=0.076463,
          lo=1e-4, hi=1.0, note="SPM figures; or set alpha_from_slope",
          enable_if=("overtopping_method", "Weggel")),
    Field("Qstar0", "Overtopping coefficient Q*0", "float", "", "", default=0.025,
          lo=0.0, hi=10.0, enable_if=("overtopping_method", "Weggel")),
    Field("alpha_from_slope", "Use alpha = 0.06 - 0.01431 sin(theta)", "bool", default=False,
          enable_if=("overtopping_method", "Weggel")),
    Field("U", "Onshore wind velocity", "float", "km/h", "kt", default=35.0 * _KN,
          lo=0.0, hi=200.0, note="0 = no wind correction"),
    Field("KR", "Refraction coefficient", "float", "", "", default=1.0, lo=0.0, hi=1.0,
          note="H'0 = KR * H0"),
    Field("overtopping_method", "Overtopping method", "choice", default="Weggel",
          choices=("Weggel", "EurOtop"),
          note="Weggel 1976 (ACES) or EurOtop 2018 mean discharge (modern standard)"),
    Field("gamma_f", "Roughness factor (EurOtop)", "float", "", "", default=0.55, lo=0.3, hi=1.0,
          note="EurOtop only: ~0.55 rough rock, 1.0 smooth",
          enable_if=("overtopping_method", "EurOtop")),
)

OUTPUTS = (
    Out("H0",   "Deepwater wave height",            "m", "ft", "scalar",
        note="Equivalent unrefracted deepwater wave height obtained by deshoaling the incident toe height (H0 = Hi / Ks)."),
    Out("ds_H0", "Relative height ds/H0",           "",  "",   "scalar",
        note="Ratio of water depth at the structure toe to the deepwater wave height, a relative-depth indicator."),
    Out("steepness", "Wave steepness H0/gT^2",      "",  "",   "scalar",
        note="Dimensionless deepwater wave steepness H0/(g T^2) characterizing the wave."),
    Out("L0",   "Deepwater wave length",            "m", "ft", "scalar",
        note="Deepwater linear wavelength L0 = g T^2 / (2 pi)."),
    Out("xi",   "Surf-similarity (Iribarren) number", "", "",  "scalar",
        note="Surf-similarity (Iribarren) number xi = tan(theta)/sqrt(Hi/L0) classifying the breaking regime (plunging, transition, or nonbreaking)."),
    Out("R",    "Wave runup",                       "m", "ft", "scalar",
        note="Maximum vertical runup height of the wave above still-water level on the slope (significant runup R_s for irregular waves)."),
    Out("F",    "Crest freeboard (hs - ds)",        "m", "ft", "scalar",
        note="Crest freeboard F = hs - ds, the structure crest height above still-water level; positive means crest above water."),
    Out("Cw",   "Wind correction factor",           "",  "",   "scalar",
        note="Onshore-wind correction factor multiplying the overtopping rate, Cw = 1 + Wf*(F/R + 0.1)*sin(theta); 1 means no wind effect."),
    Out("Q",    "Overtopping rate per unit length", "m^3/s/m", "ft^3/s/ft", "scalar",
        note="Mean wave overtopping discharge per unit length of structure crest; 0 when runup does not reach the crest."),
)


@dataclass
class Result:
    H0: float; ds_H0: float; steepness: float; L0: float; xi: float
    R: float; F: float; Cw: float; Q: float
    notes: str = ""


# --- dispersion / shoaling (Hunt 1979 explicit Pade; TR 2-1) --------------------
_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_length(T: float, d: float, g: float = G_SI) -> float:
    """Linear wavelength L (m) at depth d via Hunt (1979) explicit celerity."""
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    c2 = g * d / (y + 1.0 / denom)
    return math.sqrt(c2) * T


def shoaling_coefficient(T: float, d: float, g: float = G_SI) -> float:
    """Linear shoaling coefficient Ks = sqrt(Cg0/Cg(d)) (deepwater -> depth d)."""
    L = wave_length(T, d, g)
    k = 2.0 * math.pi / L
    n = 0.5 * (1.0 + 2.0 * k * d / math.sinh(2.0 * k * d))
    C = L / T
    Cg = n * C
    Cg0 = 0.5 * g * T / (2.0 * math.pi)          # deepwater group velocity = C0/2
    return math.sqrt(Cg0 / Cg)


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp.get(f.key, f.default))    # optional EurOtop inputs fall back to defaults
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


def _runup(slope_type, Hi, xi, theta, Pi):
    """Runup R for the chosen slope type (eqs 1, 3-8)."""
    if slope_type.startswith("Rough"):
        a, b = _runup._a, _runup._b
        return Hi * a * xi / (1.0 + b * xi)
    # smooth slope: regime by xi
    Cp = 1.002 * xi
    Cnb = 1.087 * math.sqrt(math.pi / (2.0 * theta)) + 0.775 * Pi
    if xi <= 2.0:
        C = Cp
    elif xi >= 3.5:
        C = Cnb
    else:
        C = ((3.5 - xi) / 1.5) * Cp + ((xi - 2.0) / 1.5) * Cnb
    return C * Hi


def _overtop_rate(Cw, g, Qstar0, H0p, R, F, alpha):
    """Monochromatic overtopping rate (eq 9); 0 if runup does not reach the crest."""
    if R <= F:
        return 0.0
    base = math.sqrt(g * Qstar0 * H0p ** 3)
    return Cw * base * ((R + F) / (R - F)) ** (-0.1085 / alpha)


def _overtop_eurotop(Hm0, T, theta, Rc, gamma_f, g):
    """EurOtop (2018) mean overtopping discharge (eqs 5.10/5.11), gamma_b=gamma_beta=gamma_v=1.
    Hm0 = incident significant height at the toe; Tm-1,0 = Tp/1.1 (T treated as the peak period).
    Rc = crest freeboard (clamped at >= 0; a submerged crest uses the zero-freeboard maximum)."""
    Rc = max(Rc, 0.0)
    tan_a = math.tan(theta)
    Tm10 = T / 1.1
    xi = tan_a / math.sqrt(2.0 * math.pi * Hm0 / (g * Tm10 * Tm10))
    base = math.sqrt(g * Hm0 ** 3)
    q_break = (0.023 / math.sqrt(tan_a)) * xi * math.exp(-(2.7 * Rc / (xi * Hm0 * gamma_f)) ** 1.3)
    q_max = 0.09 * math.exp(-(1.5 * Rc / (Hm0 * gamma_f)) ** 1.3)
    return min(q_break, q_max) * base


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes wave runup on a smooth or rough impermeable structure slope '
            '(seawall/revetment face) and, when runup exceeds the crest, the resulting '
            'mean overtopping rate per unit length, for monochromatic or irregular waves '
            'with an optional onshore-wind correction. Reports the surf-similarity number, '
            'runup, freeboard, wind factor, and overtopping discharge.',
 'method_key': 'overtopping_method',
 'methods': [{'name': 'Weggel (1976) overtopping (ACES)',
              'when': 'Weggel',
              'tag': 'legacy',
              'note': "ACES default, retained to reproduce the User's Guide Examples and "
                      'back-check legacy designs.',
              'equations': [{'tex': '\\xi = \\frac{\\tan\\theta}{\\sqrt{H_i/L_0}}',
                             'desc': 'Surf-similarity (Iribarren) number; theta is the '
                                     'structure seaward-face angle, L_0 deepwater '
                                     'wavelength.'},
                            {'tex': 'R = \\frac{H_i\\,a\\,\\xi}{1 + b\\,\\xi}',
                             'desc': 'Rough-slope runup (Ahrens & McCartney 1975); a, b '
                                     'empirical per armor type.'},
                            {'tex': 'Q = '
                                    "C_w\\sqrt{g\\,Q^{*}_{0}\\,H'^{3}_{0}}\\left(\\frac{R+F}{R-F}\\right)^{-0.1085/\\alpha}",
                             'desc': 'Monochromatic overtopping rate per unit length; F = '
                                     "h_s - d_s freeboard, H'_0 unrefracted deepwater "
                                     'height.'},
                            {'tex': 'C_w = 1 + W_f\\left(\\frac{F}{R} + '
                                    '0.1\\right)\\sin\\theta',
                             'desc': 'Onshore-wind correction factor, with W_f = U^2/1800 '
                                     '(U in mph).'},
                            {'tex': 'Q = \\frac{1}{199}\\sum_{i=1}^{199} Q_i,\\, R_i = '
                                    '\\sqrt{\\frac{\\ln(1/p_i)}{2}}\\,R_s',
                             'desc': 'Irregular-wave average over a Rayleigh runup '
                                     'distribution; p_i = 0.005 i.'}]},
             {'name': 'EurOtop (2018) mean discharge (modern standard)',
              'when': 'EurOtop',
              'tag': 'preferred',
              'note': 'Modern dimensionless overtopping standard (eqs 5.10/5.11) with '
                      'roughness factor gamma_f; no ACES oracle exists, validated for '
                      'physical limits/monotonicity.',
              'equations': [{'tex': '\\xi = \\frac{\\tan\\theta}{\\sqrt{H_i/L_0}}',
                             'desc': 'Surf-similarity number used by the shared runup '
                                     'methods.'},
                            {'tex': 'R = \\frac{H_i\\,a\\,\\xi}{1 + b\\,\\xi}',
                             'desc': 'Rough-slope runup (Ahrens & McCartney 1975), '
                                     'unaffected by the overtopping-method choice.'},
                            {'tex': '\\frac{q}{\\sqrt{g\\,H_{m0}^{3}}} = '
                                    '\\frac{0.023}{\\sqrt{\\tan\\theta}}\\,\\xi\\,\\exp\\left(-\\left(\\frac{2.7\\,R_c}{\\xi\\,H_{m0}\\,\\gamma_f}\\right)^{1.3}\\right)',
                             'desc': 'EurOtop breaking-wave mean discharge (eq 5.10); R_c '
                                     'crest freeboard, gamma_f roughness factor.'},
                            {'tex': '\\frac{q}{\\sqrt{g\\,H_{m0}^{3}}} = '
                                    '0.09\\,\\exp\\left(-\\left(\\frac{1.5\\,R_c}{H_{m0}\\,\\gamma_f}\\right)^{1.3}\\right)',
                             'desc': 'EurOtop non-breaking maximum (eq 5.11); the smaller '
                                     'of the two governs.'}]}],
 'symbols': [['R', 'Wave runup height above still-water level'],
             ['H_i',
              'Incident wave height at the structure toe (significant height for irregular '
              'waves)'],
             ['xi', 'Surf-similarity (Iribarren) number'],
             ['theta', 'Structure seaward-face slope angle'],
             ['L_0', 'Deepwater wavelength'],
             ['a, b',
              'Ahrens & McCartney rough-slope runup coefficients (per armor type, Table '
              'A-3)'],
             ['Q', 'Overtopping discharge rate per unit crest length'],
             ['F', 'Crest freeboard, F = h_s - d_s (structure height minus toe depth)'],
             ['Qstar0, alpha',
              'Weggel empirical overtopping coefficients (SPM 1984 figures)'],
             ['gamma_f',
              'EurOtop roughness reduction factor (~0.55 rough rock, 1.0 smooth)']],
 'references': ['Ahrens & McCartney (1975)',
                'Ahrens & Titus (1985)',
                'Ahrens & Burke (unpublished)',
                'Weggel (1976)',
                'Goda (1983)',
                'Douglass (1986)',
                'Ahrens (1977)',
                'SPM (1984) Ch. 7',
                'EurOtop (2018)']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Runup and overtopping for SI inputs (see INPUTS)."""
    _validate(inp)
    Hi = float(inp["Hi"]); T = float(inp["T"]); ds = float(inp["ds"])
    cot_theta = float(inp["cot_theta"]); hs = float(inp["hs"])
    wave_type = str(inp["wave_type"]); slope_type = str(inp["slope_type"])
    want_ot = bool(inp["want_overtopping"]); R_known = float(inp["R_known"])
    a = float(inp["a"]); b = float(inp["b"])
    alpha = float(inp["alpha"]); Qstar0 = float(inp["Qstar0"])
    U = float(inp["U"]); KR = float(inp["KR"])

    theta = math.atan(1.0 / cot_theta)                 # structure slope angle (rad)
    L0 = g * T * T / (2.0 * math.pi)
    Ks = shoaling_coefficient(T, ds, g)
    H0 = Hi / Ks                                       # deshoal toe height to deepwater
    H0p = KR * H0                                      # unrefracted deepwater height (overtopping)
    F = hs - ds                                        # crest freeboard

    xi = math.tan(theta) / math.sqrt(Hi / L0)          # surf-similarity (eq 2)
    L = wave_length(T, ds, g)                          # incident wavelength at toe (for Goda Pi)
    Pi = (Hi / L) / math.tanh(2.0 * math.pi * ds / L) ** 3   # Goda nonlinearity (eq 8)

    if bool(inp.get("alpha_from_slope", False)):
        alpha = 0.06 - 0.01431 * math.sin(theta)

    # runup (significant runup R_s for irregular)
    _runup._a, _runup._b = a, b
    R = R_known if R_known > 0.0 else _runup(slope_type, Hi, xi, theta, Pi)

    # wind correction (eq 10-11); U stored SI (m/s) -> mph for the empirical W_f
    Cw = 1.0
    if U > 0.0 and R > 0.0:
        Wf = (U / _MPH) ** 2 / 1800.0
        Cw = 1.0 + Wf * (F / R + 0.1) * math.sin(theta)

    ot_method = str(inp.get("overtopping_method", "Weggel"))
    Q = 0.0
    if want_ot:
        if ot_method == "EurOtop":
            Q = _overtop_eurotop(Hi, T, theta, F, float(inp.get("gamma_f", 0.55)), g)
        elif wave_type == "Irregular":
            # Rayleigh runup distribution; average overtopping over 199 quantiles (eqs 12-14)
            tot = 0.0
            for i in range(1, 200):
                p = 0.005 * i
                Ri = math.sqrt(math.log(1.0 / p) / 2.0) * R
                tot += _overtop_rate(Cw, g, Qstar0, H0p, Ri, F, alpha)
            Q = tot / 199.0
        else:
            Q = _overtop_rate(Cw, g, Qstar0, H0p, R, F, alpha)

    notes = [f"xi = {xi:.3f} ({_regime(xi)})", f"freeboard F = {F / _FT:.2f} ft",
             f"overtopping: {ot_method if want_ot else 'off'}"]
    if want_ot and R <= F and wave_type != "Irregular":
        notes.append("runup below crest: no overtopping")
    return Result(H0=H0, ds_H0=ds / H0, steepness=H0 / (g * T * T), L0=L0, xi=xi,
                  R=R, F=F, Cw=Cw, Q=Q, notes="; ".join(notes))


def _regime(xi):
    return "plunging" if xi <= 2.0 else ("nonbreaking" if xi >= 3.5 else "transition")


# --- self-tests (ACES User's Guide Examples 1-7) --------------------------------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    base = dict(Hi=7.50 * _FT, T=10.0, ds=12.50 * _FT, cot_theta=3.0, hs=20.0 * _FT,
                wave_type="Monochromatic", slope_type="Rough (riprap)",
                a=0.956, b=0.398, alpha=0.076463, Qstar0=0.025, U=35.0 * _KN, KR=1.0,
                alpha_from_slope=False)

    def ft(x): return x / _FT
    def cfs(x): return x / (_FT ** 2)        # m^3/s/m (= m^2/s) -> ft^3/s/ft (= ft^2/s)

    # Ex 1: rough-slope runup only
    r1 = compute(dict(base, wave_type="Monochromatic", slope_type="Rough (riprap)",
                      want_overtopping=False, R_known=0.0), g=g)
    assert _approx(ft(r1.H0), 6.386, 0.01), ft(r1.H0)
    assert _approx(ft(r1.R), 9.421, 0.01), ft(r1.R)

    # Ex 2: smooth-slope runup only (transition regime)
    r2 = compute(dict(base, wave_type="Monochromatic", slope_type="Smooth",
                      want_overtopping=False, R_known=0.0), g=g)
    assert _approx(ft(r2.R), 21.366, 0.02), ft(r2.R)

    # Ex 3: rough-slope overtopping, known runup R = 15 ft
    r3 = compute(dict(base, slope_type="Rough (riprap)", want_overtopping=True,
                      R_known=15.0 * _FT), g=g)
    assert _approx(cfs(r3.Q), 3.565, 0.01), cfs(r3.Q)

    # Ex 4: smooth-slope overtopping, known runup R = 20 ft
    r4 = compute(dict(base, slope_type="Smooth", want_overtopping=True,
                      R_known=20.0 * _FT), g=g)
    assert _approx(cfs(r4.Q), 5.368, 0.01), cfs(r4.Q)

    # Ex 5: monochromatic rough runup + overtopping (R computed = 9.421 -> Q = 0.829)
    r5 = compute(dict(base, wave_type="Monochromatic", slope_type="Rough (riprap)",
                      want_overtopping=True, R_known=0.0), g=g)
    assert _approx(ft(r5.R), 9.421, 0.01) and _approx(cfs(r5.Q), 0.829, 0.01), (ft(r5.R), cfs(r5.Q))

    # Ex 7: irregular rough runup + overtopping (Rs = 9.421 -> Q = 0.287)
    r7 = compute(dict(base, wave_type="Irregular", slope_type="Rough (riprap)",
                      want_overtopping=True, R_known=0.0), g=g)
    assert _approx(ft(r7.R), 9.421, 0.01), ft(r7.R)
    assert _approx(cfs(r7.Q), 0.287, 0.01), cfs(r7.Q)

    # EurOtop (2018) overtopping: selectable; positive, finite, and decreasing with freeboard
    eo = compute(dict(base, want_overtopping=True, R_known=0.0, overtopping_method="EurOtop"), g=g)
    assert eo.Q > 0 and math.isfinite(eo.Q), cfs(eo.Q)
    eo_hi = compute(dict(base, want_overtopping=True, R_known=0.0,
                         overtopping_method="EurOtop", hs=30.0 * _FT), g=g)
    assert eo_hi.Q < eo.Q, (cfs(eo.Q), cfs(eo_hi.Q))     # higher crest -> less overtopping
    assert _approx(ft(eo.R), 9.421, 0.01), ft(eo.R)      # runup method unaffected

    print("  self-tests: PASS (ACES Examples 1-7; EurOtop overtopping consistent)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = User's Guide Example 5: monochromatic, rough slope, runup + overtopping)")
    print(f"    H0 = {r.H0/_FT:.3f} ft   ds/H0 = {r.ds_H0:.3f}   steepness = {r.steepness:.6f}")
    print(f"    xi = {r.xi:.3f} ({_regime(r.xi)})   runup R = {r.R/_FT:.3f} ft   freeboard F = {r.F/_FT:.2f} ft")
    print(f"    wind Cw = {r.Cw:.4f}   overtopping Q = {r.Q/_FT**2:.3f} ft^3/s/ft")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
