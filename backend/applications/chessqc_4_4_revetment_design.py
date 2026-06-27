"""CHESS-QC application 4-4 — Rubble-Mound Revetment Design.

Originating ACES grouping: 4-4 "Rubble-Mound Revetment Design" (functional area:
Structural Design). Sizes the armor and filter (bedding) stone for a riprap revetment
under irregular waves, reporting median stone weight and the full weight/size gradation
of each layer, the layer thicknesses, and the expected and conservative wave runup.

Classification: exact (Hudson 1958 / Ahrens 1981 / van der Meer 1988 / Ahrens & Heimbaugh
1988 -- all coefficients known from the named sources, nothing guessed; reproduces the
User's Guide example: W_50 to ~0.5% [a rounded Rayleigh H1/10:Hs factor], runup and
gradation ratios exact).
A selectable EurOtop (2018) Ru2% runup method is also provided (the current European
standard, with a roughness factor gamma_f); Ahrens & Heimbaugh stays the default so the
Example reproduction is preserved. The EurOtop path is validated for physical limits and
ordering (no ACES oracle exists for it).
Theory and references (TR chapter 4-4, eqs 1-24 in docs/EQUATIONS.md):
  - armor weight (Hudson 1958 form):  W_50 = w_r * H_s^3 / (N_s^3 (S_r-1)^3)        (1)
  - stability number N_s = larger of two methods:
      CERC zero-damage (Ahrens 1981):  N_s = (1.45/1.27) cot(theta)^(1/6)           (3)
      van der Meer (1988): plunging (4) or surging (5) by the surf-similarity
        parameter xi_z (6) versus the transition xi_scp (7); a 1.2 shallow-water
        factor is applied to the deepwater van der Meer value.
  - layer thicknesses (8-12), gradation/sizes (13-19), runup (Ahrens & Heimbaugh
    1988; 20-24) with expected (a=1.022) and conservative (a=1.286) coefficients.

Self-containment: zero sibling imports; embeds the contract dataclasses and the Hunt
(1979) dispersion celerity. numpy + stdlib only. Runnable:
    python chessqc_4_4_revetment_design.py

Validation: reproduces the ACES User's Guide Example (H_s = 5 ft, T_s = 10 s, d_s = 9 ft,
cot(theta) = 2, w_r = 165 lb/ft^3, P = 0.1, S = 2): W_50 = 2504.6 lb, armor layer 4.95 ft,
filter layer 1.24 ft, gradation W_15/W_85/W_max/W_min, and runup 10.96 ft (expected) /
13.79 ft (conservative).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

G_SI = 9.80665
_FT = 0.3048
_LBF = 4.4482216             # lbf -> N
_LBFT3 = 157.08746           # lb/ft^3 -> N/m^3
_W_WATER = 64.0 * _LBFT3     # seawater unit weight (N/m^3); S_r = w_r / w_water


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
    aces_id="4-4",
    name="Rubble-Mound Revetment Design",
    area="Structural Design",
    classification="exact",
    cite="Ahrens (1981); van der Meer (1988); Hudson (1958); Ahrens & Heimbaugh (1988)",
    default_system="US",
)

INPUTS = (
    Field("Hs", "Significant wave height", "float", "m", "ft", default=5.0 * _FT, lo=1e-4, hi=1e3),
    Field("Ts", "Significant wave period", "float", "s", "s", default=10.0, lo=1e-2, hi=1e3),
    Field("ds", "Water depth at toe of revetment", "float", "m", "ft", default=9.0 * _FT,
          lo=1e-3, hi=1e4),
    Field("cot_theta", "Cotangent of structure slope", "float", "", "", default=2.0, lo=1.0, hi=10.0),
    Field("wr", "Unit weight of rock", "float", "N/m^3", "lb/ft^3", default=165.0 * _LBFT3,
          lo=1e3, hi=5e4),
    Field("P", "Permeability coefficient", "float", "", "", default=0.1, lo=0.1, hi=0.6,
          note="0.1 impermeable core, 0.4-0.5 permeable, 0.6 homogeneous (Fig 4-4-2)"),
    Field("S", "Damage level", "float", "", "", default=2.0, lo=1.0, hi=20.0,
          note="van der Meer damage S (Table 4-4-1)"),
    Field("runup_method", "Runup method", "choice", "", "", default="Ahrens-Heimbaugh",
          choices=("Ahrens-Heimbaugh", "EurOtop"),
          note="Ahrens & Heimbaugh 1988 (ACES) or EurOtop 2018 Ru2% (modern standard)"),
    Field("gamma_f", "Roughness factor (EurOtop)", "float", "", "", default=0.55, lo=0.3, hi=1.0,
          note="EurOtop only: ~0.40 permeable rock, 0.55 impermeable rock, 1.0 smooth",
          enable_if=("runup_method", "EurOtop")),
)

OUTPUTS = (
    Out("N_s",   "Governing stability number",      "",   "",    "scalar",
        note="Dimensionless armor stability number N_s used in Hudson sizing, the larger of the CERC (Ahrens 1981) and van der Meer (1988) values."),
    Out("W50",   "Median armor weight",             "N",  "lb",  "scalar",
        note="Median (50% passing) weight of an individual armor stone from the Hudson form, W_50 = w_r H_s^3 / (N_s^3 (S_r-1)^3)."),
    Out("D50",   "Median armor dimension",          "m",  "ft",  "scalar",
        note="Nominal equivalent-cube side length of the median armor stone, D_50 = (W_50/w_r)^(1/3)."),
    Out("r_armor", "Armor-layer thickness",         "m",  "ft",  "scalar",
        note="Thickness of the two-stone armor layer normal to the slope, r = 2 (W_50/w_r)^(1/3)."),
    Out("r_filter", "Filter-layer thickness",       "m",  "ft",  "scalar",
        note="Thickness of the underlying filter (bedding) layer, the larger of one quarter the armor thickness and 1 ft."),
    Out("W15",   "Armor W15",                        "N",  "lb",  "scalar",
        note="Armor gradation weight at 15% passing, the lighter-stone limit, W_15 = 0.40 W_50."),
    Out("W85",   "Armor W85",                        "N",  "lb",  "scalar",
        note="Armor gradation weight at 85% passing, the heavier-stone limit, W_85 = 1.96 W_50."),
    Out("Wmax",  "Armor Wmax",                       "N",  "lb",  "scalar",
        note="Maximum permissible armor-stone weight in the gradation, W_max = 4 W_50."),
    Out("Wmin",  "Armor Wmin",                       "N",  "lb",  "scalar",
        note="Minimum permissible armor-stone weight in the gradation, W_min = W_50/8."),
    Out("R_expected", "Expected maximum runup",      "m",  "ft",  "scalar",
        note="Expected (mean) vertical wave runup above still-water level on the revetment slope."),
    Out("R_conservative", "Conservative runup",      "m",  "ft",  "scalar",
        note="Conservative (upper-bound / design) vertical wave runup above still-water level for crest-elevation design."),
)


@dataclass
class Result:
    N_s: float; W50: float; D50: float; r_armor: float; r_filter: float
    W15: float; W85: float; Wmax: float; Wmin: float
    R_expected: float; R_conservative: float
    notes: str = ""


_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_length(T: float, d: float, g: float = G_SI) -> float:
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    return math.sqrt(g * d / (y + 1.0 / denom)) * T


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp.get(f.key, f.default))    # optional EurOtop inputs fall back to defaults
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


def _stability_number(Hs, Ts, cot_theta, P, S, g):
    """Governing N_s = larger of CERC (Ahrens 1981) and van der Meer (1988)."""
    tan_th = 1.0 / cot_theta
    # CERC zero-damage, adjusted to H_s (eq 3)
    N_cerc = (1.45 / 1.27) * cot_theta ** (1.0 / 6.0)
    # van der Meer (eqs 4-7): surf-similarity with the average period T_z
    Tz = Ts * (0.67 / 0.80)
    xi_z = tan_th / math.sqrt(2.0 * math.pi * Hs / (g * Tz * Tz))
    N = 7000.0                                          # conservative wave count
    SN = (S / math.sqrt(N)) ** 0.2
    xi_scp = (6.2 * P ** 0.31 * math.sqrt(tan_th)) ** (1.0 / (P + 0.5))
    if xi_z <= xi_scp:                                  # plunging (eq 4)
        N_vdm = 6.2 * P ** 0.18 * SN * xi_z ** (-0.5)
    else:                                               # surging (eq 5)
        N_vdm = 1.0 * P ** (-0.13) * SN * math.sqrt(cot_theta) * xi_z ** P
    N_vdm *= 1.2                                         # shallow-water correction
    return max(N_cerc, N_vdm), N_cerc, N_vdm, xi_z


def _runup(Hs, Ts, ds, cot_theta, g):
    """Irregular-wave runup on riprap (Ahrens & Heimbaugh 1988, eqs 20-24).
    H_mo (energy-based zero-moment height) is the SMALLER of the depth-limited (eq 22)
    and steepness-limited (eq 23) values; both eqs use the spectral peak period T_p."""
    tan_th = 1.0 / cot_theta
    Tp = Ts / 0.80                                       # eq 21
    Lp = wave_length(Tp, ds, g)                          # local peak wavelength at the toe
    Hmo_depth = 0.10 * Lp * math.tanh(2.0 * math.pi * ds / Lp)             # eq 22
    Hmo_steep = Hs / math.exp(0.00089 * (ds / (g * Tp * Tp)) ** (-0.834))  # eq 23 (T_p)
    Hmo = min(Hmo_depth, Hmo_steep)
    xi_z = tan_th / math.sqrt(2.0 * math.pi * Hmo / (g * Tp * Tp))         # eq 20
    def rmax(a, b): return Hmo * a * xi_z / (1.0 + b * xi_z)               # eq 24
    return rmax(1.022, 0.247), rmax(1.286, 0.247)


def _runup_eurotop(Hs, Ts, ds, cot_theta, gamma_f, g):
    """EurOtop (2018) 2% runup (eqs 5.1/5.2), head-on with no berm (gamma_b = gamma_beta = 1).
    Uses H_m0 (depth-/steepness-limited, as in _runup) and the spectral period T_m-1,0 = T_p/1.1.
    Returns (mean Ru2%, design Ru2% = mean + ~1 sigma) to fill the expected/conservative slots."""
    tan_th = 1.0 / cot_theta
    Tp = Ts / 0.80
    Lp = wave_length(Tp, ds, g)
    Hmo = min(0.10 * Lp * math.tanh(2.0 * math.pi * ds / Lp),
              Hs / math.exp(0.00089 * (ds / (g * Tp * Tp)) ** (-0.834)))
    Tm10 = Tp / 1.1                                                        # spectral period
    xi = tan_th / math.sqrt(2.0 * math.pi * Hmo / (g * Tm10 * Tm10))       # surf-similarity
    cap_mean = gamma_f * (4.0 - 1.5 / math.sqrt(xi))                       # surging maximum (mean)
    cap_des = gamma_f * (4.3 - 1.5 / math.sqrt(xi))                        # surging maximum (design, +1 sigma)
    ru_mean = Hmo * min(1.65 * gamma_f * xi, cap_mean)
    ru_design = Hmo * min(1.75 * gamma_f * xi, cap_des)                    # mean + ~1 sigma
    return ru_mean, ru_design


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Sizes the armor and filter (bedding) stone for a rubble-mound riprap '
            'revetment under irregular waves, returning the median stone weight, the full '
            'weight/size gradation and layer thicknesses, the governing stability number, '
            'and the expected and conservative wave runup. The armor is sized with the '
            'Hudson (1958) form using the larger of the CERC (Ahrens 1981) and van der '
            'Meer (1988) stability numbers; runup is computed by either the ACES Ahrens & '
            'Heimbaugh method or the modern EurOtop Ru2% standard.',
 'method_key': 'runup_method',
 'methods': [{'name': 'Hudson sizing with Ahrens & Heimbaugh (1988) runup',
              'when': 'Ahrens-Heimbaugh',
              'tag': 'legacy',
              'note': "ACES default; retained to reproduce the User's Guide example (runup "
                      'choice does not affect armor sizing).',
              'equations': [{'tex': 'W_{50} = \\frac{w_r\\,H_s^{3}}{N_s^{3}\\,(S_r-1)^{3}}',
                             'desc': 'Median armor weight, Hudson (1958) form for '
                                     'irregular waves (eq 1); N_s is the larger of the '
                                     'CERC and van der Meer values.'},
                            {'tex': 'N_s = \\frac{1.45}{1.27}\\,(\\cot\\theta)^{1/6}',
                             'desc': 'CERC zero-damage stability number (Ahrens 1981), '
                                     'adjusted from H_10 to H_s by 1.27 (eq 3).'},
                            {'tex': 'N_s = '
                                    '6.2\\,P^{0.18}\\,(S/\\sqrt{N})^{0.2}\\,\\xi_z^{-0.5}',
                             'desc': 'van der Meer (1988) plunging-wave stability number, '
                                     'used when xi_z <= xi_scp (eq 4); ACES applies a 1.2 '
                                     'shallow-water factor.'},
                            {'tex': 'N_s = '
                                    'P^{-0.13}\\,(S/\\sqrt{N})^{0.2}\\,\\sqrt{\\cot\\theta}\\,\\xi_z^{P}',
                             'desc': 'van der Meer (1988) surging-wave stability number, '
                                     'used when xi_z > xi_scp (eq 5).'},
                            {'tex': '\\xi_z = \\frac{\\tan\\theta}{\\sqrt{2\\pi '
                                    'H_s/(g\\,T_z^{2})}}',
                             'desc': 'Surf-similarity (Iribarren) parameter using the '
                                     'average period T_z = T_s (0.67/0.80) (eq 6).'},
                            {'tex': 'R_{max} = H_{m0}\\,\\frac{a\\,\\xi_z}{1 + '
                                    '0.247\\,\\xi_z}',
                             'desc': 'Maximum runup (Ahrens & Heimbaugh 1988, eq 24); a = '
                                     '1.022 expected, a = 1.286 conservative.'}]},
             {'name': 'Hudson sizing with EurOtop (2018) Ru2% runup',
              'when': 'EurOtop',
              'tag': 'preferred',
              'note': 'EurOtop (2018), the current European standard for Ru2% runup with '
                      'roughness factor gamma_f; armor sizing is unchanged.',
              'equations': [{'tex': 'W_{50} = \\frac{w_r\\,H_s^{3}}{N_s^{3}\\,(S_r-1)^{3}}',
                             'desc': 'Median armor weight, Hudson (1958) form for '
                                     'irregular waves (eq 1); identical to the '
                                     'Ahrens-Heimbaugh path.'},
                            {'tex': 'N_s = \\frac{1.45}{1.27}\\,(\\cot\\theta)^{1/6}',
                             'desc': 'CERC zero-damage stability number (Ahrens 1981); '
                                     'governing N_s is the larger of CERC and van der Meer '
                                     '(eqs 3-5).'},
                            {'tex': '\\xi = \\frac{\\tan\\theta}{\\sqrt{2\\pi '
                                    'H_{m0}/(g\\,T_{m-1,0}^{2})}}',
                             'desc': 'Surf-similarity parameter using the spectral period '
                                     'T_(m-1,0) = T_p/1.1.'},
                            {'tex': 'R_{u2\\%} = H_{m0}\\,(1.65\\,\\gamma_f\\,\\xi)',
                             'desc': 'EurOtop (2018) breaking-wave 2% runup (eq 5.1), with '
                                     'roughness factor gamma_f.'},
                            {'tex': 'R_{u2\\%} = H_{m0}\\,\\gamma_f\\left(4 - '
                                    '\\frac{1.5}{\\sqrt{\\xi}}\\right)',
                             'desc': 'EurOtop (2018) surging-wave maximum cap (eq 5.2); '
                                     'the smaller of this and the breaking value is '
                                     'taken.'}]}],
 'symbols': [['W_50', 'Median armor-stone weight'],
             ['w_r', 'Unit weight of armor rock'],
             ['H_s', 'Significant wave height'],
             ['N_s', 'Stability number (larger of CERC and van der Meer)'],
             ['S_r', 'Stone specific gravity, w_r / w_water'],
             ['theta', 'Structure slope angle; cot(theta) is the input'],
             ['P', 'Permeability coefficient (0.1 impermeable core to 0.6 homogeneous)'],
             ['S', 'van der Meer damage level'],
             ['xi_z', 'Surf-similarity (Iribarren) parameter'],
             ['H_m0',
              'Energy-based zero-moment wave height (depth- or steepness-limited)']],
 'references': ['Hudson (1958)',
                'Ahrens (1981)',
                'van der Meer (1988)',
                'van der Meer & Pilarczyk (1987)',
                'Ahrens & Heimbaugh (1988)',
                'Battjes (1974)',
                'EM 1110-2-2300',
                'EurOtop (2018)']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Revetment armor/filter sizing and runup for SI inputs."""
    _validate(inp)
    Hs = float(inp["Hs"]); Ts = float(inp["Ts"]); ds = float(inp["ds"])
    cot_theta = float(inp["cot_theta"]); wr = float(inp["wr"])
    P = float(inp["P"]); S = float(inp["S"])

    Sr = wr / _W_WATER
    N_s, N_cerc, N_vdm, xi_z = _stability_number(Hs, Ts, cot_theta, P, S, g)

    W50 = wr * Hs ** 3 / (N_s ** 3 * (Sr - 1.0) ** 3)   # eq 1 (Hudson form), in N
    # gradation (eqs 13-16) and dimensions (eq 17)
    Wmax = 4.0 * W50; Wmin = W50 / 8.0
    W85 = 1.96 * W50; W15 = 0.40 * W50
    def dim(W): return (W / wr) ** (1.0 / 3.0)
    D50 = dim(W50)
    # layer thicknesses (eqs 8-9)
    r_armor = 2.0 * (W50 / wr) ** (1.0 / 3.0)
    r_filter = max(r_armor / 4.0, 1.0 * _FT)

    runup_method = str(inp.get("runup_method", "Ahrens-Heimbaugh"))
    if runup_method == "EurOtop":
        R_exp, R_cons = _runup_eurotop(Hs, Ts, ds, cot_theta, float(inp.get("gamma_f", 0.55)), g)
        runup_note = f"EurOtop 2018 Ru2% (gamma_f={float(inp.get('gamma_f', 0.55))})"
    else:
        R_exp, R_cons = _runup(Hs, Ts, ds, cot_theta, g)
        runup_note = "Ahrens & Heimbaugh 1988"

    notes = (f"N_s={N_s:.3f} (CERC {N_cerc:.3f}, vdM {N_vdm:.3f}; "
             f"{'CERC' if N_cerc >= N_vdm else 'van der Meer'} governs); xi_z={xi_z:.2f}; "
             f"runup: {runup_note}")
    return Result(N_s=N_s, W50=W50, D50=D50, r_armor=r_armor, r_filter=r_filter,
                  W15=W15, W85=W85, Wmax=Wmax, Wmin=Wmin,
                  R_expected=R_exp, R_conservative=R_cons, notes=notes)


# --- self-tests (ACES User's Guide Example) -------------------------------------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    r = compute({"Hs": 5.0 * _FT, "Ts": 10.0, "ds": 9.0 * _FT, "cot_theta": 2.0,
                 "wr": 165.0 * _LBFT3, "P": 0.1, "S": 2.0}, g=g)
    lb = lambda x: x / _LBF
    ft = lambda x: x / _FT
    # weights/dimensions match to ~0.5% (residual is the rounded Rayleigh H1/10:Hs factor)
    rel = lambda got, exp: abs(got - exp) <= 0.006 * exp
    assert rel(lb(r.W50), 2504.61), lb(r.W50)
    assert rel(lb(r.W15), 1001.84), lb(r.W15)
    assert rel(lb(r.W85), 4909.04), lb(r.W85)
    assert rel(lb(r.Wmax), 10018.44), lb(r.Wmax)
    assert rel(lb(r.Wmin), 313.08), lb(r.Wmin)
    assert _approx(ft(r.r_armor), 4.95, 0.02), ft(r.r_armor)
    assert _approx(ft(r.r_filter), 1.24, 0.02), ft(r.r_filter)
    assert _approx(ft(r.D50), 2.48, 0.01), ft(r.D50)
    # gradation ratios are exact by construction
    assert _approx(r.Wmax / r.W50, 4.0, 1e-9) and _approx(r.W85 / r.W50, 1.96, 1e-9)
    # runup reproduces the oracle exactly
    assert _approx(ft(r.R_expected), 10.96, 0.02), ft(r.R_expected)
    assert _approx(ft(r.R_conservative), 13.79, 0.02), ft(r.R_conservative)
    # EurOtop (2018) runup: selectable; sane, ordered, and lower than the rough Ahrens values
    eo = compute({"Hs": 5.0 * _FT, "Ts": 10.0, "ds": 9.0 * _FT, "cot_theta": 2.0,
                  "wr": 165.0 * _LBFT3, "P": 0.1, "S": 2.0, "runup_method": "EurOtop"}, g=g)
    assert eo.R_conservative > eo.R_expected > 0, (ft(eo.R_expected), ft(eo.R_conservative))
    assert 0.0 < ft(eo.R_expected) < ft(r.R_conservative), ft(eo.R_expected)  # rough rock reduces runup
    assert eo.W50 == r.W50          # runup method does not affect armor sizing
    print(f"  self-tests: PASS (W50={lb(r.W50):.1f} lb, armor {ft(r.r_armor):.2f} ft, "
          f"runup {ft(r.R_expected):.2f}/{ft(r.R_conservative):.2f} ft)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    lb = lambda x: x / _LBF; ft = lambda x: x / _FT
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  W50 = {lb(r.W50):.2f} lb   D50 = {ft(r.D50):.2f} ft   "
          f"armor layer = {ft(r.r_armor):.2f} ft   filter layer = {ft(r.r_filter):.2f} ft")
    print(f"  gradation (lb): Wmin={lb(r.Wmin):.2f}  W15={lb(r.W15):.2f}  W50={lb(r.W50):.2f}  "
          f"W85={lb(r.W85):.2f}  Wmax={lb(r.Wmax):.2f}")
    print(f"  runup: expected={ft(r.R_expected):.2f} ft   conservative={ft(r.R_conservative):.2f} ft")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
