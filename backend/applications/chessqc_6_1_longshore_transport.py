"""CHESS-QC application 6-1 — Longshore Sediment Transport (CERC method).

Originating ACES application: 6-1 "Longshore Sediment Transport" (functional area:
Littoral Processes). Estimates the potential volumetric longshore sand transport rate
from wave conditions, using the CERC energy-flux method, from breaking-wave or deepwater
inputs.

Classification: exact. Reproduces the ACES source (`Lstran.for`, subroutine SEDTRANS),
recompiled and run, to 5e-5 on both of its published examples.
Theory and references: SPM (1984) Ch. 4 (Eq. 4-49); Galvin (1979, 1980); Gravens (1988);
Wang, Kraus & Davis (1998). Equations transcribed in docs/EQUATIONS.md, TR chapter 6-1.

    Q = K * P_ls / [ (rho_s - rho) * g * a' ]              (volume per second)
    breaking:  P_ls = c_b * rho * g^1.5 * H_b^2.5 * sin(2 alpha_b)
    deepwater: P_ls = c_d * rho * g^1.5 * H_s0^2.5 * (cos alpha_0)^0.25 * sin(2 alpha_0)

with K the empirical coefficient (0.39 for field data with significant height), a' the
solids fraction (1 - porosity = 0.6), rho_s sediment density, rho water density. The
deepwater coefficient folds in the shoaling coefficient K_s = 1.3 and the breaking index
d_b = 1.28 H_b (TR 6-1).

The two coefficient vintages, and what the 25 percent really is. ACES ships its
longshore source twice, and the 2/98 revision note in `Lstran.for` states the reason:
the original code used solitary-wave approximations; Wang, Kraus & Davis measured ACES
predicting 25 percent above CERC; Galvin (1980) had already put the solitary
approximation at about that much high; so the coefficients were replaced with
linear-wave ones for consistency with SMS/Genesis. The pairs are

    deepwater  0.050389 -> 0.04031          breaking  0.088388 -> 0.07071

which are ratios of exactly 1.25. The published User's Guide examples (275,234 yd^3/yr
deepwater, 2,662,872 breaking) were computed with the solitary pair, and running the
pre-1998 source reproduces them to six figures.

This supersedes an earlier reading, recorded here and in the tests, that the examples
could only be reproduced with an effective rho_s of about 2320 kg/m^3 "below quartz and
not stated in the Technical Reference". That was wrong. ACES's sediment density is the
`5.14` in `(5.14 - RHO)` at `Lstran.for:1801`, in slug/ft^3, which is 2649 kg/m^3 -
quartz, exactly as it should be. Nothing about the density was ever unusual; the whole
discrepancy was the coefficient vintage. Both are offered, defaulting to the linear
pair, since that is what the current source computes.

Self-containment: zero sibling imports; embeds its own contract dataclasses. Runnable
standalone:
    python chessqc_6_1_longshore_transport.py
which runs the CERC-factor and both-vintage source self-tests, then prints the example.
stdlib only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

G_SI = 9.80665
_YR = 31557600.0          # seconds per Julian year (365.25 d)
_M3_TO_YD3 = 1.0 / 0.764554858


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
    aces_id="6-1",
    name="Longshore Sediment Transport",
    area="Littoral Processes",
    classification="exact",
    cite="SPM (1984) Ch.4; Galvin (1979); TR 6-1",
    default_system="US",
)

_FT = 0.3048
_METHODS = ("Deepwater wave conditions", "Breaking wave conditions")

# Energy-flux coefficients, taken from the ACES source rather than rounded. ACES
# shipped two vintages side by side and the 2/98 revision note in Lstran.for says
# exactly why: the original used solitary-wave approximations, Wang, Kraus & Davis
# measured ACES running 25 percent above CERC, Galvin (1980) had already put the
# solitary approximation at about that much high, and the coefficients were replaced
# with linear-wave ones to agree with SMS/Genesis. The pairs differ by exactly 1.25.
_LINEAR = "Linear (2/98 revision, SMS/Genesis consistent)"
_SOLITARY = "Solitary (original ACES, reproduces the User's Guide examples)"
_THEORIES = (_LINEAR, _SOLITARY)
_COEF = {                       # theory: (deepwater, breaking)
    _LINEAR: (0.04031, 0.07071),        # Lstran.for:1796, :1820
    _SOLITARY: (0.050389, 0.088388),    # lstran.forg:1771, :1795
}
INPUTS = (
    Field("method", "Wave data type", "choice", "", "", default="Deepwater wave conditions",
          choices=_METHODS, note="deepwater (H_s0, alpha_0) or breaking (H_b, alpha_b)"),
    Field("H", "Wave height", "float", "m", "ft", default=1.75 * _FT, lo=1e-4, hi=1e3,
          note="deepwater significant height H_s0, or breaker height H_b"),
    Field("angle", "Wave angle to shoreline", "angle", "deg", "deg", default=15.0, lo=0.0, hi=90.0,
          note="deepwater crest angle alpha_0, or breaker angle alpha_b"),
    Field("K", "Empirical coefficient K", "float", "", "", default=0.39, lo=0.0, hi=2.0,
          note="CERC coefficient; 0.39 for field data with significant wave height"),
    Field("rho_water", "Water density", "float", "kg/m^3", "kg/m^3", default=1025.18,
          lo=900.0, hi=1100.0, note="seawater ~1025, fresh ~1000"),
    Field("rho_sand", "Sediment density", "float", "kg/m^3", "kg/m^3", default=2650.0,
          lo=1500.0, hi=3500.0,
          note="quartz sand ~2650, which is what ACES uses: its 5.14 slug/ft^3 is "
               "2649 kg/m^3"),
    Field("wave_theory", "Energy-flux wave theory", "choice", "", "", default=_LINEAR,
          choices=_THEORIES,
          note="ACES shipped both. The solitary coefficients are 1.25x the linear ones "
               "and are what the published User's Guide examples were computed with; the "
               "2/98 revision replaced them with the linear ones"),
    Field("porosity", "Sediment porosity", "float", "", "", default=0.40, lo=0.0, hi=0.7,
          note="pore fraction; solids fraction a' = 1 - porosity (TR uses a' = 0.6)"),
)

OUTPUTS = (
    Out("Q", "Longshore transport rate", "m^3/yr", "yd^3/yr", "scalar",
        note="Potential volumetric longshore sand transport rate Q from the CERC formula, in m^3/yr (positive in the down-drift direction set by the wave angle)."),
    Out("Q_m3s", "Transport rate (volume/sec)", "m^3/s", "m^3/s", "scalar",
        note="Same potential longshore transport rate expressed as volume per second, Q in m^3/s."),
    Out("P_ls", "Longshore energy flux factor", "N/s", "lb/s", "scalar",
        note="Longshore component of the wave energy flux factor P_ls driving transport, computed from wave height and angle."),
    Out("cerc_factor", "Q / P_ls", "", "", "scalar",
        note="Ratio of transport rate to energy-flux factor, Q/P_ls in m^3/yr per W/m (~1290 for quartz sand and seawater)."),
)


@dataclass
class Result:
    Q: float; Q_m3s: float; P_ls: float; cerc_factor: float
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Estimates the potential volumetric longshore sand transport rate along a '
            'beach using the CERC energy-flux method, returning the rate in m^3/yr (and '
            'yd^3/yr) together with the longshore energy-flux factor. The user supplies '
            'either breaking-wave or deepwater wave conditions, which selects the '
            'energy-flux formula, and which of the two coefficient vintages ACES shipped '
            'to use: the linear-wave values its 2/98 revision adopted for consistency '
            'with SMS/Genesis, or the original solitary-wave values, 1.25 times larger, '
            'that the published examples were computed with.',
 'method_key': 'method',
 'methods': [{'name': 'Deepwater wave conditions',
              'when': 'Deepwater wave conditions',
              'tag': 'standard',
              'note': 'Energy flux built from deepwater significant height and crest '
                      'angle, folding in shoaling (K_s = 1.3) and the breaker index.',
              'equations': [{'tex': 'P_{ls} = 0.04 \\, \\rho \\, g^{1.5} \\, H_{s0}^{2.5} '
                                    '\\, (\\cos \\alpha_{0})^{0.25} \\sin(2 \\alpha_{0})',
                             'desc': 'Deepwater longshore energy-flux factor (TR 6-1 eq '
                                     '17; SPM 4-44 reduced).'},
                            {'tex': 'Q = \\frac{K \\, P_{ls}}{(\\rho_{s} - \\rho) \\, g '
                                    '\\, a}',
                             'desc': 'CERC volumetric longshore transport rate (TR 6-1 eq '
                                     '1), a = 1 - porosity.'}]},
             {'name': 'Breaking wave conditions',
              'when': 'Breaking wave conditions',
              'tag': 'standard',
              'note': 'Energy flux evaluated directly at the breaker line from breaker '
                      'height and breaker angle.',
              'equations': [{'tex': 'P_{ls} = 0.0707 \\, \\rho \\, g^{1.5} \\, H_{b}^{2.5} '
                                    '\\sin(2 \\alpha_{b})',
                             'desc': 'Breaking longshore energy-flux factor (TR 6-1 eq 11; '
                                     'constant = (1/16)\\sqrt{1.28}).'},
                            {'tex': 'Q = \\frac{K \\, P_{ls}}{(\\rho_{s} - \\rho) \\, g '
                                    '\\, a}',
                             'desc': 'CERC volumetric longshore transport rate (TR 6-1 eq '
                                     '1), a = 1 - porosity.'}]}],
 'symbols': [['Q', 'Potential volumetric longshore transport rate (volume per unit time)'],
             ['P_{ls}', 'Longshore energy-flux factor'],
             ['K',
              'Empirical CERC coefficient (0.39 for field data with significant wave '
              'height)'],
             ['H_b', 'Breaker significant wave height'],
             ['H_{s0}', 'Deepwater significant wave height'],
             ['c_{b}, c_{d}',
              'Energy-flux coefficients. ACES shipped two vintages: solitary-wave '
              '(0.088388 breaking, 0.050389 deepwater) and, after its 2/98 revision, '
              'linear-wave values exactly 1.25 times smaller. The published examples in '
              'the ACES manual were computed with the solitary pair'],
             ['alpha_b', 'Wave crest angle to shoreline at breaking'],
             ['alpha_0', 'Deepwater wave crest angle to shoreline'],
             ['rho_s', 'Sediment (sand) density; quartz ~2650 kg/m^3'],
             ['rho', 'Water density; seawater ~1025 kg/m^3'],
             ['a', 'Solids-to-total volume ratio, a = 1 - porosity (TR uses 0.6)']],
 'references': ['SPM (1984) Ch. 4 (Eq. 4-49)',
                'Galvin (1979), CERC TP 79-1',
                'Gravens (1988)',
                'CESM/TR 6-1']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Longshore transport rate for SI inputs. Returns Q in m^3/yr (SI internal)."""
    _validate(inp)
    method = str(inp.get("method", "Deepwater wave conditions"))
    H = float(inp["H"]); ang = math.radians(float(inp["angle"])); K = float(inp["K"])
    rho = float(inp["rho_water"]); rho_s = float(inp["rho_sand"]); por = float(inp["porosity"])
    theory = str(inp.get("wave_theory", _LINEAR))
    if theory not in _COEF:
        raise ValueError(f"unknown wave theory: {theory!r}")
    a = 1.0 - por
    if rho_s <= rho:
        raise ValueError("sediment density must exceed water density")

    c_deep, c_break = _COEF[theory]
    if method.startswith("Breaking"):
        P_ls = c_break * rho * g ** 1.5 * H ** 2.5 * math.sin(2.0 * ang)         # TR eq 11
    else:
        P_ls = (c_deep * rho * g ** 1.5 * H ** 2.5
                * math.cos(ang) ** 0.25 * math.sin(2.0 * ang))                   # TR eq 17

    Q_m3s = K * P_ls / ((rho_s - rho) * g * a)                                    # TR eq 1
    Q_m3yr = Q_m3s * _YR
    cerc = Q_m3yr / P_ls if P_ls > 0 else 0.0

    notes = [f"{method.split()[0].lower()} input; "
             f"{'linear' if theory == _LINEAR else 'solitary'} energy-flux coefficient; "
             f"a' = {a:.2f}, rho_s = {rho_s:.0f} kg/m^3"]
    if abs(ang) < 1e-9:
        notes.append("zero wave angle: no longshore transport")
    return Result(Q=Q_m3yr, Q_m3s=Q_m3s, P_ls=P_ls, cerc_factor=cerc, notes="; ".join(notes))


# --- self-tests (CERC literature factor + ACES examples) ------------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    g = G_SI
    base = dict(K=0.39, rho_water=1025.18, rho_sand=2650.0, porosity=0.40)

    # 1) physically-standard quartz reproduces the CERC literature factor ~1290 m^3/yr per W/m
    rd = compute({**base, "method": "Deepwater wave conditions", "H": 1.75 * _FT, "angle": 15.0}, g=g)
    rb = compute({**base, "method": "Breaking wave conditions", "H": 3.75 * _FT, "angle": 12.0}, g=g)
    assert _approx(rd.cerc_factor, 1287.0, 5e-3), rd.cerc_factor
    assert _approx(rb.cerc_factor, 1287.0, 5e-3), rb.cerc_factor

    # 2) angle dependence: maximum near alpha = 45 deg; zero at 0
    q0 = compute({**base, "method": "Breaking wave conditions", "H": 3.75 * _FT, "angle": 0.0}, g=g)
    q45 = compute({**base, "method": "Breaking wave conditions", "H": 3.75 * _FT, "angle": 45.0}, g=g)
    assert q0.Q == 0.0 and q45.Q > rb.Q

    # 3) against the ACES source, recompiled and run (tests/aces_oracle/fortran,
    #    `sh build.sh lst`). ACES computes in its own English-unit constants, so they
    #    are supplied here: g = 32.17 ft/s^2, rho = 1.989 and rho_s = 5.14 slug/ft^3.
    #    Quartz, standard porosity, no adjustment of any kind.
    _SLUG_FT3 = 515.378818
    aces = dict(K=0.39, rho_water=1.989 * _SLUG_FT3, rho_sand=5.14 * _SLUG_FT3,
                porosity=0.40)
    g_aces = 32.17 * _FT
    for theory, deep_yd3, break_yd3 in (
            (_SOLITARY, 275234.06, 2662871.0),      # what the User's Guide prints
            (_LINEAR, 220180.69, 2130285.0)):       # what the current source gives
        qd = compute({**aces, "wave_theory": theory,
                      "method": "Deepwater wave conditions", "H": 1.75 * _FT,
                      "angle": 15.0}, g=g_aces)
        qb = compute({**aces, "wave_theory": theory,
                      "method": "Breaking wave conditions", "H": 3.75 * _FT,
                      "angle": 12.0}, g=g_aces)
        assert _approx(qd.Q * _M3_TO_YD3, deep_yd3, 5e-5), (theory, qd.Q * _M3_TO_YD3)
        assert _approx(qb.Q * _M3_TO_YD3, break_yd3, 5e-5), (theory, qb.Q * _M3_TO_YD3)

    # 4) the two vintages differ by exactly the 25 percent the revision note describes
    a = compute({**aces, "wave_theory": _SOLITARY, "method": "Breaking wave conditions",
                 "H": 3.75 * _FT, "angle": 12.0}, g=g_aces).Q
    b = compute({**aces, "wave_theory": _LINEAR, "method": "Breaking wave conditions",
                 "H": 3.75 * _FT, "angle": 12.0}, g=g_aces).Q
    assert _approx(a / b, 1.25, 1e-4), a / b

    print("  self-tests: PASS (CERC factor 1287 with quartz; both ACES vintages "
          "reproduced to 5e-5 with quartz and ACES's own constants, no fitted density)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        vv = inp[f.key]
        sval = f"{vv:>10.4g}" if isinstance(vv, (int, float)) and f.kind != "choice" else f"{vv:>10}"
        print(f"    {f.label:28s} {f.key:10s} = {sval} {f.unit_si}")
    print("  OUTPUTS:")
    print(f"    Longshore transport rate     Q          = {r.Q:>12.0f} m^3/yr  "
          f"({r.Q * _M3_TO_YD3:.0f} yd^3/yr)")
    print(f"    Energy flux factor           P_ls       = {r.P_ls:>12.1f} W/m")
    print(f"    CERC factor (Q/P_ls)         cerc       = {r.cerc_factor:>12.0f} m^3/yr per W/m")
    print("  note: the User's Guide prints 275,234 yd^3/yr for this case. That is the "
          "pre-1998 solitary")
    print("        coefficient, 1.25x the linear one used here; select the solitary "
          "theory to match it.")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
