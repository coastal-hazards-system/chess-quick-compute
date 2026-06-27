"""CHESS-QC application 3-1 — Linear Wave Theory with Snell's Law.

Transforms a wave of known height/period/direction at one depth to deep water and
to a "subject" depth, using linear wave theory (Hunt 1979 dispersion), Snell's law
(O'Brien 1942) for refraction, and energy-flux conservation for shoaling. Reports
wave height / crest angle / length / celerity / group velocity / energy density /
energy flux at three locations (Known, Deep water, Subject), Ursell numbers,
deepwater steepness, and the Weggel (1972) breaker height/depth.

Classification: exact (closed-form linear wave theory + Snell's law + energy-flux shoaling;
reproduces the User's Guide Example 3-1 to the digit. The breaker height/depth, a secondary
output, uses the standard Singamsetti & Wind / Weggel index.)

Self-contained (zero sibling imports): embeds the AppMeta/Field/Out/Result
dataclasses and the Hunt dispersion solver (same kernel validated in 2-1).

Theory: TR 3-1 (eqs 1-5) for the transformation; breaker height Singamsetti & Wind
(1980) and breaker depth Weggel (1972), per TR 6-1 "Monochromatic Wave Breaking".
Fully validated against the ACES User's Guide Example 3-1, including the breaker
(Hb=12.29 ft, db=15.25 ft) -- see tests/test_manual_oracle.py.

Run:
    python chessqc_3_1_snell.py          # self-tests + tabulate the manual example
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# --- physical constants (SI) ----------------------------------------------------
G_SI = 9.80665          # m/s^2
RHO_SALT = 1025.18      # kg/m^3 (seawater)
RHO_FRESH = 999.0       # kg/m^3 (fresh water)
_FT = 0.3048            # ft -> m


# --- contract dataclasses (mirrors 2-1) -----------------------------------------
@dataclass(frozen=True)
class AppMeta:
    aces_id: str
    name: str
    area: str
    classification: str
    cite: str
    default_system: str = "SI"


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


# --- metadata --------------------------------------------------------------------
APP_META = AppMeta(
    aces_id="3-1",
    name="Linear Wave Theory with Snell's Law",
    area="Wave Transformation",
    classification="exact",
    cite="O'Brien (1942); Hunt (1979); Weggel (1972); TR 3-1",
    default_system="US",     # opens on the User's Guide Example (US units)
)

# Complete input list (ACES User's Guide 3-1). Defaults = User's Guide Example 3-1
# (H1=10 ft, T=7.5 s, d1=25 ft, alpha1=10 deg, cot phi=100, d2=20 ft).
INPUTS = (
    Field("H1",   "Wave height (known)",     "float", "m",   "ft",  default=10.0 * _FT, lo=1e-6, hi=1e4,
          note="> 0 (height at the known depth)"),
    Field("T",    "Wave period",             "float", "s",   "s",   default=7.5, lo=1e-3, hi=1e4,
          note="> 0"),
    Field("d1",   "Water depth (known)",     "float", "m",   "ft",  default=25.0 * _FT, lo=1e-6, hi=1e5,
          note="> 0"),
    Field("alpha1", "Wave crest angle (known)", "angle", "deg", "deg", default=10.0, lo=0.0, hi=90.0,
          note="angle between wave crest and depth contour, 0-90 deg"),
    Field("cot_phi", "Cotan of nearshore slope", "float", "", "", default=100.0, lo=1e-6, hi=1e6,
          note="cot(beach slope); e.g. 100 = 1:100 slope"),
    Field("d2",   "Water depth (subject)",   "float", "m",   "ft",  default=20.0 * _FT, lo=1e-6, hi=1e5,
          note="> 0 (depth to transform the wave to)"),
    Field("water", "Water type", "choice", "", "", default="Salt", choices=("Salt", "Fresh"),
          note="sets water density used for energy density/flux"),
)

# Complete output list (ACES User's Guide 3-1): per-location scalars + breaker.
OUTPUTS = (
    # Known location
    Out("H_k",  "Wave height (known)",      "m",   "ft",    "scalar",
        note="Wave height at the known (input) depth, equal to the specified H1."),
    Out("a_k",  "Wave crest angle (known)", "deg", "deg",   "scalar",
        note="Angle between the wave crest and the depth contour at the known depth (input alpha1)."),
    Out("L_k",  "Wavelength (known)",       "m",   "ft",    "scalar",
        note="Linear-theory wavelength at the known depth (Hunt 1979 dispersion)."),
    Out("C_k",  "Wave celerity (known)",    "m/s", "ft/s",  "scalar",
        note="Wave phase speed C = L/T at the known depth."),
    Out("Cg_k", "Group velocity (known)",   "m/s", "ft/s",  "scalar",
        note="Group velocity (speed of energy transport) at the known depth."),
    Out("E_k",  "Energy density (known)",   "N/m", "lb/ft", "scalar",
        note="Mean wave energy per unit surface area, E = rho g H^2/8, at the known depth."),
    Out("P_k",  "Energy flux (known)",      "N/s", "lb/s",  "scalar",
        note="Wave energy flux (power per unit crest width) P = E C_g at the known depth."),
    Out("Ur_k", "Ursell number (known)",    "",    "",      "scalar",
        note="Ursell number U_r = H L^2 / d^3 at the known depth, a check on linear-theory validity."),
    # Deep water
    Out("H_0",  "Wave height (deep)",       "m",   "ft",    "scalar",
        note="Equivalent unrefracted deepwater wave height back-transformed from the known wave."),
    Out("a_0",  "Wave crest angle (deep)",  "deg", "deg",   "scalar",
        note="Wave crest angle in deep water from Snell's law (sin alpha_0 = sin alpha_1 * C_0/C_1)."),
    Out("L_0",  "Wavelength (deep)",        "m",   "ft",    "scalar",
        note="Deepwater wavelength L_0 = g T^2 / (2 pi)."),
    Out("C_0",  "Wave celerity (deep)",     "m/s", "ft/s",  "scalar",
        note="Deepwater wave celerity C_0 = g T / (2 pi)."),
    Out("Cg_0", "Group velocity (deep)",    "m/s", "ft/s",  "scalar",
        note="Deepwater group velocity, equal to half the deepwater celerity (C_0/2)."),
    Out("E_0",  "Energy density (deep)",    "N/m", "lb/ft", "scalar",
        note="Mean wave energy per unit surface area in deep water, E = rho g H_0^2/8."),
    Out("P_0",  "Energy flux (deep)",       "N/s", "lb/s",  "scalar",
        note="Wave energy flux P = E C_g in deep water (conserved during refraction/shoaling)."),
    Out("steep_0", "Deepwater wave steepness", "", "",      "scalar",
        note="Deepwater wave steepness H_0/L_0, a dimensionless measure of wave nonlinearity."),
    # Subject location
    Out("H_2",  "Wave height (subject)",      "m",   "ft",    "scalar",
        note="Transformed wave height at the subject depth, H_0 times the refraction and shoaling coefficients."),
    Out("a_2",  "Wave crest angle (subject)", "deg", "deg",   "scalar",
        note="Wave crest angle at the subject depth from Snell's law (sin alpha_2 = sin alpha_0 * C_2/C_0)."),
    Out("L_2",  "Wavelength (subject)",       "m",   "ft",    "scalar",
        note="Linear-theory wavelength at the subject depth (Hunt 1979 dispersion)."),
    Out("C_2",  "Wave celerity (subject)",    "m/s", "ft/s",  "scalar",
        note="Wave phase speed C = L/T at the subject depth."),
    Out("Cg_2", "Group velocity (subject)",   "m/s", "ft/s",  "scalar",
        note="Group velocity (speed of energy transport) at the subject depth."),
    Out("E_2",  "Energy density (subject)",   "N/m", "lb/ft", "scalar",
        note="Mean wave energy per unit surface area, E = rho g H^2/8, at the subject depth."),
    Out("P_2",  "Energy flux (subject)",      "N/s", "lb/s",  "scalar",
        note="Wave energy flux (power per unit crest width) P = E C_g at the subject depth."),
    Out("Ur_2", "Ursell number (subject)",    "",    "",      "scalar",
        note="Ursell number U_r = H L^2 / d^3 at the subject depth, a check on linear-theory validity."),
    # Breaker (Weggel 1972)
    Out("Hb",   "Breaker height",             "m",   "ft",    "scalar",
        note="Breaking wave height for the deepwater wave on the nearshore slope (Singamsetti & Wind 1980)."),
    Out("db",   "Breaker depth",              "m",   "ft",    "scalar",
        note="Water depth at incipient breaking corresponding to Hb (Weggel 1972)."),
)


@dataclass
class Result:
    H_k: float; a_k: float; L_k: float; C_k: float; Cg_k: float; E_k: float; P_k: float; Ur_k: float
    H_0: float; a_0: float; L_0: float; C_0: float; Cg_0: float; E_0: float; P_0: float; steep_0: float
    H_2: float; a_2: float; L_2: float; C_2: float; Cg_2: float; E_2: float; P_2: float; Ur_2: float
    Hb: float; db: float
    notes: str = ""


# --- dispersion solver (Hunt 1979 explicit Pade; identical to 2-1) --------------
_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_celerity(T: float, d: float, g: float = G_SI) -> float:
    """Explicit linear celerity c (m/s) via Hunt (1979), accuracy < 0.01%."""
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    return math.sqrt(g * d / (y + 1.0 / denom))


def _wave_at(T: float, d: float, g: float):
    """Return (L, C, Cg, kd) at depth d (finite)."""
    C = wave_celerity(T, d, g)
    L = C * T
    k = 2.0 * math.pi / L
    kd = k * d
    Cg = 0.5 * C * (1.0 + 2.0 * kd / math.sinh(2.0 * kd))
    return L, C, Cg, kd


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue                                  # choice/bool fields carry no bounds
        v = inp[f.key]
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")
    if inp["d2"] <= 0 or inp["d1"] <= 0:
        raise ValueError("water depths must be > 0")


# --- breaker (TR 6-1-1/2 "Monochromatic Wave Breaking") -------------------------
def _breaker(H0: float, L0: float, T: float, m: float):
    """Breaker height + depth for a deepwater wave (H0, L0) on beach slope m=tan(phi).

    Height: Singamsetti & Wind (1980) for m>0 (TR eq 3); McCowan Hb=0.78*d for m<=0.
    Depth:  Weggel (1972) db = Hb/(b - a*Hb/T^2)  (TR eq 4),
            b = 1/(0.64*(1+e^-19.5m)),  a = 1.36*(1-e^-19m).
    The height index uses only dimensionless ratios (unit-agnostic); the Weggel depth
    coefficients are US-calibrated (g folded into 1.36), so the depth is evaluated in
    US units (ft, s) and converted back to SI. Returns (Hb, db) in SI metres."""
    if m > 0.0:
        Hb = H0 * 0.575 * m ** 0.031 * (H0 / L0) ** (-0.254)        # eq 3 (SI; ratio-based)
    else:
        Hb = 0.78 * H0            # flat/unknown slope: depth-limited (caller uses H0~d)
    b = 1.0 / (0.64 * (1.0 + math.exp(-19.5 * m)))                  # eq 4 coefficients
    a = 1.36 * (1.0 - math.exp(-19.0 * m))
    Hb_ft = Hb / _FT
    db_ft = Hb_ft / (b - a * Hb_ft / (T * T))                       # eq 4 (ft, s)
    return Hb, db_ft * _FT


# --- compute (single entry point both front-ends call) --------------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Transforms a wave of known height, period and crest angle from one depth to '
            "deep water and to a subject depth using linear wave theory, Snell's law for "
            'refraction and energy-flux conservation for shoaling. Reports height, crest '
            'angle, length, celerity, group velocity, energy density/flux and Ursell '
            'number at each location, plus the Weggel breaker height and depth.',
 'methods': [{'name': "Snell's law + energy-flux shoaling/refraction",
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': '\\frac{c}{c_{0}} = \\frac{\\sin \\alpha}{\\sin '
                                    '\\alpha_{0}}',
                             'desc': "Snell's law: ratio of celerities equals ratio of "
                                     "sines of the crest-to-contour angles (O'Brien "
                                     '1942).'},
                            {'tex': 'K_{s} = \\sqrt{\\frac{C_{g0}}{C_{g}}}',
                             'desc': 'Shoaling coefficient from energy-flux conservation '
                                     '(E = rho g H^2/8, P = E C_g).'},
                            {'tex': 'K_{r} = \\sqrt{\\frac{\\cos \\alpha_{0}}{\\cos '
                                    '\\alpha}}',
                             'desc': 'Refraction coefficient for straight, parallel depth '
                                     'contours (orthogonal spacing ratio).'},
                            {'tex': '\\frac{H}{H_{0}} = K_{r} \\cdot K_{s}',
                             'desc': 'Combined transformation: height ratio is the product '
                                     'of refraction and shoaling coefficients.'},
                            {'tex': 'H_{b} = 0.575 \\, H_{0} \\, m^{0.031} '
                                    '\\left(\\frac{H_{0}}{L_{0}}\\right)^{-0.254}',
                             'desc': 'Breaker height for finite nearshore slope m = '
                                     'tan(phi) (Singamsetti & Wind 1980).'},
                            {'tex': 'd_{b} = \\frac{H_{b}}{b - a \\, H_{b}/T^{2}}',
                             'desc': 'Breaker depth (Weggel 1972), b = '
                                     '1/(0.64(1+e^{-19.5m})), a = 1.36(1-e^{-19m}); '
                                     'US-calibrated.'}]}],
 'symbols': [['c, c_0', 'Wave celerity at the depth contour and in deep water'],
             ['alpha, alpha_0',
              'Angle between wave crest and depth contour (subject and deep water)'],
             ['C_g, C_g0', 'Group velocity at a finite depth and in deep water'],
             ['K_s', 'Shoaling coefficient, sqrt(C_g0/C_g)'],
             ['K_r', 'Refraction coefficient, sqrt(cos alpha_0 / cos alpha)'],
             ['H_0, L_0', 'Deepwater wave height and wavelength'],
             ['m', 'Nearshore beach slope, tan(phi) = 1/cot(phi)'],
             ['H_b, d_b', 'Breaker height and breaker depth'],
             ['T', 'Wave period'],
             ['U_r', 'Ursell number, H L^2 / d^3 (linearity check)']],
 'references': ["O'Brien (1942)",
                'Hunt (1979)',
                'Weggel (1972)',
                'Singamsetti & Wind (1980)',
                'McCowan (1894)',
                'SPM (1984)',
                'Dean & Dalrymple (1984)',
                'ACES TR 3-1 (eqs 1-7); breaker from TR 6-1']}


def compute(inp: dict, *, g: float = G_SI, rho: float | None = None) -> Result:
    """Snell's-law wave transformation for SI inputs
    {H1, T, d1, alpha1[deg], cot_phi, d2, water}. Water density follows the `water`
    field (Salt|Fresh) unless `rho` is given."""
    _validate(inp)
    if rho is None:
        rho = RHO_FRESH if str(inp.get("water", "Salt")) == "Fresh" else RHO_SALT
    H1 = float(inp["H1"]); T = float(inp["T"]); d1 = float(inp["d1"])
    a1 = math.radians(float(inp["alpha1"]))
    cot_phi = float(inp["cot_phi"]); d2 = float(inp["d2"])

    def energy(H, Cg):
        E = rho * g * H * H / 8.0
        return E, E * Cg

    # Known location
    L_k, C_k, Cg_k, _ = _wave_at(T, d1, g)
    # Deep water
    C_0 = g * T / (2.0 * math.pi)
    L_0 = C_0 * T
    Cg_0 = 0.5 * C_0
    # Snell to deep: sin a0 = sin a1 * C0/C1
    sin_a0 = min(1.0, math.sin(a1) * C_0 / C_k)
    a0 = math.asin(sin_a0)
    cos_a0 = math.cos(a0)
    # deepwater height (back-transform the known wave): H1 = H0 * Ks1 * Kr1
    Ks1 = math.sqrt(Cg_0 / Cg_k)
    Kr1 = math.sqrt(cos_a0 / math.cos(a1))
    H_0 = H1 / (Ks1 * Kr1)
    # Subject location
    L_2, C_2, Cg_2, _ = _wave_at(T, d2, g)
    sin_a2 = min(1.0, math.sin(a0) * C_2 / C_0)            # Snell from deep
    a2 = math.asin(sin_a2)
    Ks2 = math.sqrt(Cg_0 / Cg_2)
    Kr2 = math.sqrt(cos_a0 / math.cos(a2))
    H_2 = H_0 * Ks2 * Kr2

    E_k, P_k = energy(H1, Cg_k)
    E_0, P_0 = energy(H_0, Cg_0)
    E_2, P_2 = energy(H_2, Cg_2)
    Ur_k = H1 * L_k * L_k / (d1 ** 3)
    Ur_2 = H_2 * L_2 * L_2 / (d2 ** 3)
    steep_0 = H_0 / L_0

    m = 1.0 / cot_phi                                      # beach slope tan(phi)
    Hb, db = _breaker(H_0, L_0, T, m)

    notes = []
    if Ur_k > 26 or Ur_2 > 26:
        notes.append("Ursell > 26: linear theory questionable")
    if max(d1, d2) / L_0 > 0.5:
        notes.append("near/at deep water")
    notes.append(f"slope 1:{cot_phi:g}")

    return Result(
        H_k=H1, a_k=math.degrees(a1), L_k=L_k, C_k=C_k, Cg_k=Cg_k, E_k=E_k, P_k=P_k, Ur_k=Ur_k,
        H_0=H_0, a_0=math.degrees(a0), L_0=L_0, C_0=C_0, Cg_0=Cg_0, E_0=E_0, P_0=P_0, steep_0=steep_0,
        H_2=H_2, a_2=math.degrees(a2), L_2=L_2, C_2=C_2, Cg_2=Cg_2, E_2=E_2, P_2=P_2, Ur_2=Ur_2,
        Hb=Hb, db=db, notes="; ".join(notes),
    )


# --- self-tests + manual-example tabulation -------------------------------------
def _self_tests() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    ft = lambda x: x / _FT
    # transformation against User's Guide Example 3-1 (US)
    assert abs(ft(r.L_0) - 288.00) < 0.5, ft(r.L_0)
    assert abs(r.a_0 - 15.00) < 0.1, r.a_0
    assert abs(ft(r.H_0) - 10.68) < 0.05, ft(r.H_0)
    assert abs(ft(r.C_2) - 23.51) < 0.05 and abs(r.a_2 - 9.12) < 0.1
    assert abs(ft(r.H_2) - 10.27) < 0.05, ft(r.H_2)
    assert abs(ft(r.L_k) - 193.27) < 0.5 and abs(r.Ur_2 - 39.91) < 0.2
    assert abs(r.steep_0 - 0.04) < 0.005
    assert abs(ft(r.Hb) - 12.29) < 0.05 and abs(ft(r.db) - 15.25) < 0.05    # breaker
    print("  self-tests: PASS (matches User's Guide Example 3-1, incl. breaker)")


def _tab() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    ft = lambda x: x / _FT
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (values in US units; matches User's Guide Example 3-1)")
    print(f"  {'item':16}{'known':>10}{'deep':>10}{'subject':>10}")
    rows = [
        ("Height (ft)", ft(r.H_k), ft(r.H_0), ft(r.H_2)),
        ("Angle (deg)", r.a_k, r.a_0, r.a_2),
        ("Length (ft)", ft(r.L_k), ft(r.L_0), ft(r.L_2)),
        ("Celerity", ft(r.C_k), ft(r.C_0), ft(r.C_2)),
        ("Group vel", ft(r.Cg_k), ft(r.Cg_0), ft(r.Cg_2)),
    ]
    for name, k, o, s in rows:
        print(f"  {name:16}{k:10.2f}{o:10.2f}{s:10.2f}")
    print(f"  Ursell known={r.Ur_k:.2f}  subject={r.Ur_2:.2f}  deepwater steepness={r.steep_0:.3f}")
    print(f"  Breaker: Hb={ft(r.Hb):.2f} ft  db={ft(r.db):.2f} ft  [manual: 12.29 / 15.25]")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _tab()
