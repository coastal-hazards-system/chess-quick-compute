"""CHESS-QC application 5-5 — Wave Setup Across the Surf Zone.

Originating ACES grouping: 5-5 "Wave Setup" (functional area: Wave Runup,
Transmission, and Overtopping). Predicts the wave-induced change in mean water
level across the surf zone: the small set-down just seaward of breaking and the
larger set-up that raises the mean water line toward and onto the beach.

Classification: exact (closed-form radiation-stress theory).
Theory and references: radiation-stress balance of Longuet-Higgins and Stewart
(1962, 1963); breaker height by Singamsetti and Wind (1980); breaker depth and the
breaker index a(m), b(m) by Weggel (1972), as transcribed in docs/EQUATIONS.md
(TR chapter 8-1, "Miscellaneous Routines"). Explicit dispersion via Hunt (1979).

Note on the Weggel breaker index. The transcribed eq (4) db = Hb/(b - a*Hb/T^2)
with a = 1.36*(1-e^-19.5m) is written in US units (ft, s), with gravity folded into
the coefficient (43.8/32.17 = 1.36). CHESS-QC computes in SI, so this app uses the
gravity-explicit equivalent db = Hb/(b - a*Hb/(g*T^2)) with a = 43.8*(1-e^-19.5m);
the two agree exactly when g = 32.17 ft/s^2.

Self-containment: zero sibling imports; embeds its own Field/AppMeta/Out/Result
dataclasses and the Hunt dispersion solver. Runnable standalone:
    python chessqc_5_5_wave_setup.py
which runs the closed-form / consistency self-tests, then prints an ACES-style
tabulation. stdlib + numpy only.

Validation strategy (no dedicated User's Guide section for this app, and no TR
chapter of its own): the set-down at breaking reproduces the Longuet-Higgins and
Stewart closed form (and its shallow-water limit eta_b -> -gamma^2 d_b/16), the
set-up gradient reproduces the LH-S surf-zone slope, and the cross-surf-zone table
integrates from the breaker set-down up to the maximum set-up.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

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
    aces_id="5-5",
    name="Wave Setup",
    area="Wave Runup, Transmission, and Overtopping",
    classification="exact",
    cite="Longuet-Higgins & Stewart (1963); Weggel (1972); Singamsetti & Wind (1980)",
    default_system="US",
)

_FT = 0.3048
# Default example: a moderate deepwater swell on a 1:50 beach (no refraction).
INPUTS = (
    Field("H0", "Deepwater wave height", "float", "m", "ft", default=6.0 * _FT,
          lo=1e-4, hi=1e3, note="> 0 (significant or monochromatic deepwater height)"),
    Field("T",  "Wave period", "float", "s", "s", default=10.0, lo=1e-2, hi=1e3,
          note="> 0"),
    Field("m",  "Beach slope (tan beta)", "float", "", "", default=0.02, lo=1e-4, hi=1.0,
          note="bottom slope; Singamsetti & Wind valid roughly 0.02 to 0.2"),
    Field("KR", "Refraction coefficient", "float", "", "", default=1.0, lo=0.0, hi=1.0,
          note="K_R = sqrt(b0/b); 1.0 = no refraction. H'0 = K_R * H0"),
)

OUTPUTS = (
    Out("H0p",       "Unrefracted deepwater height H'0", "m",  "ft", "scalar",
        note="Unrefracted equivalent deepwater wave height H'0 = K_R * H0, used as the reference height for breaking."),
    Out("L0",        "Deepwater wave length",            "m",  "ft", "scalar",
        note="Deepwater wavelength L0 = g T^2/(2 pi) from linear dispersion."),
    Out("a_index",   "Weggel breaker coefficient a(m)",  "",   "",   "scalar",
        note="Slope-dependent Weggel (1972) breaker-index coefficient a(m) = 43.8(1 - e^-19.5m) in the gravity-explicit SI breaker-depth relation."),
    Out("b_index",   "Weggel breaker coefficient b(m)",  "",   "",   "scalar",
        note="Slope-dependent Weggel (1972) breaker-index coefficient b(m) = 1.56/(1 + e^-19.5m) in the breaker-depth relation."),
    Out("Hb",        "Breaker height",                   "m",  "ft", "scalar",
        note="Wave height at breaking on the finite beach slope (Singamsetti & Wind 1980)."),
    Out("db",        "Breaker depth",                    "m",  "ft", "scalar",
        note="Still-water depth at the breaker line (Weggel 1972)."),
    Out("gamma_b",   "Breaker index Hb/db",              "",   "",   "scalar",
        note="Breaker index gamma_b = Hb/db, the ratio of breaker height to breaker depth (usually 0.6 to 1.4)."),
    Out("setdown_b", "Set-down at breaking",             "m",  "ft", "scalar",
        note="Lowering of mean water level at the breaker line (negative); eta_b = -Hb^2 k_b/(8 sinh(2 k_b d_b))."),
    Out("setup_slope", "Set-up gradient d(eta)/dx",      "",   "",   "scalar",
        note="Cross-surf-zone slope of the mean water level d(eta)/dx = m beta/(1+beta) with beta = 3 gamma_b^2/8 (dimensionless rise per unit shoreward distance)."),
    Out("surf_width", "Surf-zone width (breaker to SWL shoreline)", "m", "ft", "scalar",
        note="Horizontal distance from the breaker line to the still-water shoreline, db/m."),
    Out("setup_swl", "Set-up at still-water shoreline",  "m",  "ft", "scalar",
        note="Mean water level rise at the still-water shoreline (set-down marched shoreward over the surf width)."),
    Out("setup_max", "Maximum set-up (at waterline)",    "m",  "ft", "scalar",
        note="Maximum mean water level rise, reached at the displaced waterline where mean depth d + eta = 0."),
    Out("dx_shore",  "Shoreline displacement",           "m",  "ft", "scalar",
        note="Horizontal distance the mean waterline moves up the beach beyond the still-water shoreline due to set-up."),
    # cross-surf-zone profile (for plots): distance shoreward of the breaker line
    Out("profile_x",   "Profile: distance shoreward of breaker", "m", "ft", "profile",
        note="Cross-shore coordinate measured shoreward from the breaker line, for the surf-zone profile plot."),
    Out("profile_bed", "Profile: bed elevation",         "m",  "ft", "profile",
        note="Bed elevation along the profile, -db + m x (negative below still-water level)."),
    Out("profile_mwl", "Profile: mean water level (set-up)", "m", "ft", "profile",
        note="Mean water level along the profile (set-down at the breaker line rising to maximum set-up at the waterline)."),
)


@dataclass
class Result:
    H0p: float; L0: float; a_index: float; b_index: float
    Hb: float; db: float; gamma_b: float
    setdown_b: float; setup_slope: float; surf_width: float
    setup_swl: float; setup_max: float; dx_shore: float
    profile_x: np.ndarray
    profile_bed: np.ndarray
    profile_mwl: np.ndarray
    notes: str = ""


# --- dispersion solver (Hunt 1979 explicit Pade; TR 2-1 eq 11) ------------------
_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_celerity(T: float, d: float, g: float = G_SI) -> float:
    """Explicit linear celerity c (m/s) via Hunt (1979), accuracy < 0.01%."""
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    c2 = g * d / (y + 1.0 / denom)
    return math.sqrt(c2)


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


def breaker_index(m: float) -> tuple[float, float]:
    """Weggel (1972) breaker-index coefficients a(m), b(m), SI / gravity-explicit.
    a = 43.8*(1 - e^-19.5m) ; b = 1.56/(1 + e^-19.5m).  (US-unit a is a/g = 1.36...)"""
    e = math.exp(-19.5 * m)
    return 43.8 * (1.0 - e), 1.56 / (1.0 + e)


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes the wave-induced change in mean water level across the surf zone — '
            'the small set-down just seaward of breaking and the larger set-up that raises '
            'the mean waterline up the beach — from deepwater wave height, period and '
            'beach slope using radiation-stress theory.',
 'methods': [{'name': 'Radiation-stress set-down and set-up (Longuet-Higgins & Stewart)',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'L_0 = \\frac{g T^2}{2 \\pi}',
                             'desc': 'Deepwater wavelength from the wave period (linear '
                                     "dispersion); H'_0 = K_R H_0 is the unrefracted "
                                     'equivalent deepwater height.'},
                            {'tex': "H_b = H_0' \\, 0.575 \\, m^{0.031} \\left( "
                                    "\\frac{H_0'}{L_0} \\right)^{-0.254}",
                             'desc': 'Breaker height on a finite slope (Singamsetti & Wind '
                                     '1980), valid roughly m = 0.02 to 0.2.'},
                            {'tex': 'd_b = \\frac{H_b}{b - a \\, H_b/(g T^2)}, \\quad a = '
                                    '43.8 (1 - e^{-19.5 m}), \\quad b = \\frac{1.56}{1 + '
                                    'e^{-19.5 m}}',
                             'desc': 'Breaker depth via the Weggel (1972) breaker index '
                                     '(gravity-explicit SI form); gives the breaker index '
                                     'gamma_b = H_b/d_b.'},
                            {'tex': '\\eta_b = - \\frac{H_b^2 k_b}{8 \\sinh(2 k_b d_b)}',
                             'desc': 'Set-down of the mean water level at the breaker line '
                                     '(Longuet-Higgins & Stewart 1963); shallow-water '
                                     'limit eta_b -> -gamma_b^2 d_b/16.'},
                            {'tex': '\\frac{d\\eta}{dx} = m \\, \\frac{\\beta}{1 + '
                                    '\\beta}, \\quad \\beta = \\frac{3 \\gamma_b^2}{8}',
                             'desc': 'Cross-surf-zone set-up gradient (LH-S); marched '
                                     'shoreward from eta_b to give set-up at the '
                                     'shoreline, maximum set-up and shoreline '
                                     'displacement.'}]}],
 'symbols': [["H_0'", "Unrefracted equivalent deepwater wave height, H_0' = K_R H_0"],
             ['K_R', 'Refraction coefficient (1.0 = no refraction)'],
             ['L_0', 'Deepwater wavelength'],
             ['m', 'Beach slope, tan(beta_slope)'],
             ['H_b', 'Breaker height'],
             ['d_b', 'Still-water depth at breaking'],
             ['gamma_b', 'Breaker index, H_b/d_b'],
             ['k_b', 'Wavenumber at the breaker line (2 pi / L_b)'],
             ['eta_b', 'Set-down of mean water level at breaking (negative)'],
             ['beta', 'Radiation-stress set-up parameter, 3 gamma_b^2 / 8']],
 'references': ['Longuet-Higgins & Stewart (1962, 1963)',
                'Weggel (1972)',
                'Singamsetti & Wind (1980)',
                'Hunt (1979)',
                'docs/EQUATIONS.md ch. 8-1 (Miscellaneous Routines)']}


def compute(inp: dict, *, g: float = G_SI, n_profile: int = 121) -> Result:
    """Wave-setup results for SI inputs {H0, T, m, KR}."""
    _validate(inp)
    H0 = float(inp["H0"]); T = float(inp["T"]); m = float(inp["m"]); KR = float(inp["KR"])

    L0 = g * T * T / (2.0 * math.pi)
    H0p = KR * H0                                   # unrefracted equivalent deepwater height

    # breaker height (Singamsetti & Wind 1980; equations.md 8-1 eq 3)
    Hb = H0p * (0.575 * m ** 0.031 * (H0p / L0) ** (-0.254))

    # breaker depth (Weggel 1972; gravity-explicit SI form of equations.md 8-1 eq 4)
    a_index, b_index = breaker_index(m)
    db = Hb / (b_index - a_index * Hb / (g * T * T))
    if db <= 0.0:
        raise ValueError("breaker depth solved non-positive; inputs out of valid range")
    gamma_b = Hb / db

    # wavelength / wavenumber at the breaker line
    Cb = wave_celerity(T, db, g)
    Lb = Cb * T
    kb = 2.0 * math.pi / Lb

    # set-down at the breaker line (Longuet-Higgins & Stewart 1963):
    #   eta_b = - H^2 k / (8 sinh(2 k d))   evaluated at breaking
    setdown_b = -(Hb * Hb * kb) / (8.0 * math.sinh(2.0 * kb * db))

    # surf-zone set-up gradient (LH-S): d(eta)/dx = m * beta/(1+beta), beta = 3 gamma^2/8
    beta = 3.0 * gamma_b * gamma_b / 8.0
    setup_slope = m * beta / (1.0 + beta)

    # geometry: x measured shoreward from the breaker line; still-water depth d(x)=db-m x.
    x_sw = db / m                                   # breaker line -> still-water shoreline
    surf_width = x_sw
    setup_swl = setdown_b + setup_slope * x_sw      # set-up at the still-water shoreline
    # actual waterline where mean depth (d + eta) = 0:
    x_wl = (db + setdown_b) / (m - setup_slope)
    setup_max = setdown_b + setup_slope * x_wl      # maximum set-up
    dx_shore = x_wl - x_sw                           # shoreline displacement up the beach

    # cross-surf-zone profile for plotting (breaker line to a little past the waterline)
    x = np.linspace(0.0, max(x_wl, x_sw) * 1.05, n_profile)
    bed = -db + m * x                               # bed elevation (negative below SWL)
    mwl = setdown_b + setup_slope * x               # mean water level
    mwl = np.minimum(mwl, np.maximum(bed, mwl[0]))  # clip MWL to the beach face past the waterline

    notes = []
    notes.append(f"breaker index gamma_b = {gamma_b:.3f}")
    if not (0.6 <= gamma_b <= 1.4):
        notes.append("WARNING: gamma_b outside the usual 0.6 to 1.4 band; check slope/steepness")
    if not (0.02 <= m <= 0.2):
        notes.append("note: beach slope outside Singamsetti & Wind range (0.02 to 0.2)")
    notes.append("set-down applies seaward of breaking; set-up across the surf zone")

    return Result(
        H0p=H0p, L0=L0, a_index=a_index, b_index=b_index, Hb=Hb, db=db, gamma_b=gamma_b,
        setdown_b=setdown_b, setup_slope=setup_slope, surf_width=surf_width,
        setup_swl=setup_swl, setup_max=setup_max, dx_shore=dx_shore,
        profile_x=x, profile_bed=bed, profile_mwl=mwl, notes="; ".join(notes))


# --- self-tests (closed forms + consistency) ------------------------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    g = G_SI
    r = compute({"H0": 6.0 * _FT, "T": 10.0, "m": 0.02, "KR": 1.0}, g=g)

    # 1) refraction: KR scales H'0 linearly
    r2 = compute({"H0": 6.0 * _FT, "T": 10.0, "m": 0.02, "KR": 0.9}, g=g)
    assert _approx(r2.H0p, 0.9 * r.H0p, 1e-9), (r2.H0p, r.H0p)

    # 2) signs / ordering: set-down negative, set-up positive and increasing shoreward
    assert r.setdown_b < 0.0, r.setdown_b
    assert r.setup_max > r.setup_swl > r.setdown_b, (r.setup_max, r.setup_swl, r.setdown_b)
    assert 0.0 < r.setup_slope < r.gamma_b, r.setup_slope  # slope is a fraction of m < 1

    # 3) breaker index in the physical band
    assert 0.6 <= r.gamma_b <= 1.4, r.gamma_b

    # 4) LH-S set-down shallow-water closed form: eta_b -> -gamma^2 d_b / 16
    #    (test on a shallow case where sinh(2kd) ~ 2kd)
    rs = compute({"H0": 4.0 * _FT, "T": 16.0, "m": 0.02, "KR": 1.0}, g=g)
    approx_shallow = -(rs.gamma_b ** 2) * rs.db / 16.0
    assert _approx(rs.setdown_b, approx_shallow, 5e-2), (rs.setdown_b, approx_shallow)

    # 5) table consistency: set-up at SWL shoreline equals gradient marched over surf width
    assert _approx(r.setup_swl, r.setdown_b + r.setup_slope * r.surf_width, 1e-9)
    # and max set-up equals gradient marched to the displaced waterline
    x_wl = r.surf_width + r.dx_shore
    assert _approx(r.setup_max, r.setdown_b + r.setup_slope * x_wl, 1e-9)

    # 6) set-up gradient closed form: m * beta/(1+beta), beta = 3 gamma^2/8
    beta = 3.0 * r.gamma_b ** 2 / 8.0
    assert _approx(r.setup_slope, 0.02 * beta / (1.0 + beta), 1e-12)

    print("  self-tests: PASS (refraction, signs/ordering, gamma band, LH-S set-down, table)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        print(f"    {f.label:30s} {f.key:4s} = {inp[f.key]:>10.4g} {f.unit_si}")
    print("  OUTPUTS:")
    by_key = {o.key: o for o in OUTPUTS}
    for kk in ("H0p", "L0", "a_index", "b_index", "Hb", "db", "gamma_b",
               "setdown_b", "setup_slope", "surf_width", "setup_swl", "setup_max", "dx_shore"):
        o = by_key[kk]
        print(f"    {o.label:34s} {kk:11s} = {getattr(r, kk):>12.5g} {o.unit_si}")
    # also show in US for the default (US) example
    print("  (US: Hb=%.2f ft, db=%.2f ft, set-down=%.3f ft, max set-up=%.3f ft, dx=%.1f ft)" % (
        r.Hb / _FT, r.db / _FT, r.setdown_b / _FT, r.setup_max / _FT, r.dx_shore / _FT))
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
