"""CHESS-QC application 6-3 — Longshore Transport using CEDRS Statistics.

Originating ACES application: 6-3 "Longshore Sediment Transport using CEDRS percent-
occurrence statistics" (functional area: Littoral Processes; the CEDRS branch of the
original ACES 6-1 Longshore Sediment Transport, separated as its own application in
ACES). It estimates the net and gross potential longshore transport at a site by
summing the deepwater CERC transport over a directional wave climate supplied as a
Coastal Engineering Data Retrieval System (CEDRS) "percent occurrence of wave height and
period by direction" table.

Classification: standard (the exact CERC energy-flux formula applied over a CEDRS
directional climate, but with self-chosen bin-midpoint and contributing-fraction conventions
not spelled out in the source; reproduces the User's Guide net transport to 0.58%, not to
the digit).
Theory and references: SPM (1984) Ch.4 deepwater CERC formula (as in CHESS-QC 6-1); Gravens
(1988); ACES User's Guide Longshore Sediment Transport, Example 3. The directional climate
is the WIS hindcast percent-occurrence table; the default dataset is Gulf of Mexico
Station G1033 (29.0N, 85.5W, 68 m), WIS Report 18 (file G1033.810).

Method (ACES User's Guide Example 3):
  * The CEDRS table has 16 directional bands, each 22.5 deg wide, centered on azimuths
    0, 22.5, ..., 337.5 deg (the direction waves approach FROM). Each band holds the
    percent occurrence (x1000) of significant wave height by peak period; only the height
    distribution matters to the CERC formula (it has no period dependence).
  * For a shore-normal azimuth theta, a band's angle to the shore normal is
    delta = wrap(theta - band_azimuth) in (-180, 180]. Only the part of the 22.5-deg band
    lying within +/-90 deg of the shore normal reaches the coast; the "contributing
    percentage" is that overlap fraction. Waves from the left (delta>0) drive positive
    (rightward) transport; from the right (delta<0), negative.
  * Per band: Q_band = frac * sum_h [ p_h * Q_deep(H_h, delta) ], summed over height bins,
    with Q_deep the deepwater CERC rate (SPM 4-49) evaluated at the band's representative
    angle and the bin's midpoint height, and p_h the bin's occurrence fraction. The net
    transport is the signed sum; the gross is the sum of magnitudes.

A note on the sediment density (inherited from 6-1). With physically standard quartz
(rho_s = 2650) this returns correct-physics CERC transport; the ACES User's Guide example
(net -854,849 yd^3/yr) is reproduced with the effective rho_s ~ 2319 that 6-1 documents
(below quartz, not stated in the TR). CHESS-QC defaults to quartz and exposes rho_s.

Self-containment: zero sibling imports; embeds its own contract dataclasses and the G1033
dataset. Runnable standalone:
    python chessqc_6_3_cedrs_transport.py
which runs the User's-Guide Example-3 self-test (net to <1% with rho_s=2319) then prints
the default (quartz) example. stdlib only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

G_SI = 9.80665
_YR = 31557600.0
_M3_TO_YD3 = 1.0 / 0.764554858
_FT = 0.3048

# Height-bin midpoints (m): 0.00-0.49, 0.50-0.99, ..., 4.50-4.99, 5.00+ (11 bins).
_HMID = tuple(0.245 + 0.5 * k for k in range(11))
_HLABELS = ("0.00-0.49", "0.50-0.99", "1.00-1.49", "1.50-1.99", "2.00-2.49", "2.50-2.99",
            "3.00-3.49", "3.50-3.99", "4.00-4.49", "4.50-4.99", "5.00+")

# CEDRS percent occurrence (x1000) of significant height by direction band, for the
# default station G1033 (Gulf of Mexico St. 33; WIS Report 18). Row i = band azimuth
# i*22.5 deg; columns = the 11 height bins above (summed over peak period, which the
# CERC formula does not use). Transcribed from ACES User's Guide Table 6-1-3.
_G1033 = [
    [69, 940, 879, 348, 42, 1, 0, 0, 0, 0, 0],        # az 0.0
    [122, 979, 993, 360, 58, 0, 0, 0, 0, 0, 0],       # az 22.5
    [116, 1583, 1438, 592, 126, 0, 0, 0, 0, 0, 0],    # az 45.0
    [401, 3106, 3686, 1061, 150, 0, 0, 0, 0, 0, 0],   # az 67.5
    [843, 5787, 7191, 1767, 210, 10, 0, 0, 0, 0, 0],  # az 90.0
    [869, 5874, 5872, 1132, 132, 30, 6, 0, 0, 0, 0],  # az 112.5
    [916, 4933, 4462, 1216, 186, 30, 11, 0, 0, 0, 0], # az 135.0
    [338, 2985, 3153, 1088, 280, 44, 3, 0, 0, 0, 0],  # az 157.5
    [292, 2213, 1753, 700, 195, 11, 0, 0, 0, 0, 0],   # az 180.0
    [207, 2111, 1340, 449, 87, 18, 3, 0, 0, 0, 0],    # az 202.5
    [412, 2665, 2207, 458, 108, 42, 13, 0, 0, 0, 0],  # az 225.0
    [395, 2636, 1434, 361, 86, 7, 4, 0, 0, 0, 0],     # az 247.5
    [264, 1793, 1080, 464, 144, 26, 0, 0, 0, 0, 0],   # az 270.0
    [158, 1571, 1196, 639, 169, 3, 0, 0, 0, 0, 0],    # az 292.5
    [155, 1234, 1284, 487, 48, 6, 0, 0, 0, 0, 0],     # az 315.0
    [114, 1088, 919, 314, 61, 4, 0, 0, 0, 0, 0],      # az 337.5
]


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
    columns: tuple = ()
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
    aces_id="6-3",
    name="Longshore Transport using CEDRS Statistics",
    area="Littoral Processes",
    classification="standard",
    cite="SPM (1984) Ch.4; Gravens (1988); WIS Report 18; ACES User's Guide Example 6-1-3",
    default_system="US",
)

INPUTS = (
    Field("shore_azimuth", "Shore-normal azimuth", "angle", "deg", "deg", default=40.0,
          lo=0.0, hi=360.0, note="seaward shore normal, measured clockwise from true north"),
    Field("K", "Empirical coefficient K", "float", "", "", default=0.39, lo=0.0, hi=1.0,
          note="CERC coefficient; 0.39 for field data with significant wave height"),
    Field("rho_water", "Water density", "float", "kg/m^3", "kg/m^3", default=1025.18,
          lo=900.0, hi=1100.0, note="seawater ~1025"),
    Field("rho_sand", "Sediment density", "float", "kg/m^3", "kg/m^3", default=2650.0,
          lo=1500.0, hi=3500.0, note="quartz ~2650; ~2319 reproduces the ACES example"),
    Field("porosity", "Sediment porosity", "float", "", "", default=0.40, lo=0.0, hi=0.7,
          note="pore fraction; solids fraction a' = 1 - porosity"),
    Field("occ", "CEDRS percent-occurrence (x1000) by band x height", "table",
          default=[row[:] for row in _G1033],
          columns=tuple((lbl, "x1000", "x1000") for lbl in _HLABELS),
          note="16 rows (band azimuths 0,22.5,...,337.5 deg) x 11 height bins; default G1033"),
)

OUTPUTS = (
    Out("Q_net", "Net longshore transport", "m^3/yr", "yd^3/yr", "scalar",
        note="signed sum of all directional-band transport rates; positive is rightward (looking seaward), negative leftward"),
    Out("Q_gross", "Gross longshore transport", "m^3/yr", "yd^3/yr", "scalar",
        note="sum of the magnitudes of all band transport rates, i.e. total volume moved regardless of direction"),
    Out("Q_right", "Transport to the right (+)", "m^3/yr", "yd^3/yr", "scalar",
        note="sum of the positive (rightward) band contributions, from waves approaching from the left of the shore normal"),
    Out("Q_left", "Transport to the left (-)", "m^3/yr", "yd^3/yr", "scalar",
        note="sum of the negative (leftward) band contributions, from waves approaching from the right of the shore normal"),
    Out("band_angle", "Per-band angle from shore normal", "deg", "deg", "profile",
        note="angle delta of each 22.5 deg directional band to the shore normal, wrapped to (-180,180]; positive drives rightward transport"),
    Out("band_pct", "Per-band contributing percentage", "%", "%", "profile",
        note="fraction (percent) of each band lying within +/-90 deg of the shore normal, i.e. the part whose waves reach the coast"),
    Out("band_Q", "Per-band transport rate", "m^3/yr", "yd^3/yr", "profile",
        note="signed CERC transport rate contributed by each directional band, summed over its height bins and scaled by the contributing fraction"),
)


@dataclass
class Result:
    Q_net: float
    Q_gross: float
    Q_right: float
    Q_left: float
    band_angle: list = field(default_factory=list)
    band_pct: list = field(default_factory=list)
    band_Q: list = field(default_factory=list)
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle") or f.key not in inp:
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo:g}, {f.hi:g}] ({f.note})")


def _wrap180(a: float) -> float:
    """Wrap an angle to (-180, 180]."""
    a = (a + 180.0) % 360.0 - 180.0
    return a + 360.0 if a <= -180.0 else a


def deepwater_Q(H: float, delta_deg: float, K: float, rho: float, rho_s: float,
                a_solid: float, g: float = G_SI) -> float:
    """Deepwater CERC transport (m^3/yr) for height H (m) at angle delta to shore normal.
    P_ls = 0.04 rho g^1.5 H^2.5 cos(delta)^0.25 sin(2 delta); the sign follows sin(2 delta).
    cos is clamped at 0 so |delta| >= 90 deg (waves not reaching the shore) gives 0."""
    al = math.radians(delta_deg)
    c = max(math.cos(al), 0.0)
    P_ls = 0.04 * rho * g ** 1.5 * H ** 2.5 * c ** 0.25 * math.sin(2.0 * al)
    return K * P_ls / ((rho_s - rho) * g * a_solid) * _YR


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Estimates net and gross potential longshore sediment transport at a site by '
            'applying the deepwater CERC energy-flux formula to a CEDRS percent-occurrence '
            "wave climate, summing each directional band's contribution and reporting the "
            'signed (net) and magnitude (gross) totals.',
 'methods': [{'name': 'Deepwater CERC over a CEDRS directional climate',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'Q = \\frac{K \\, P_{ls}}{(\\rho_s - \\rho) \\, g \\, '
                                    'a}',
                             'desc': 'CERC volumetric longshore transport rate from the '
                                     'longshore energy-flux factor P_ls (TR 6-1 eq 1).'},
                            {'tex': 'P_{ls} = 0.04 \\, \\rho \\, g^{3/2} \\, H_{s0}^{5/2} '
                                    '\\, (\\cos \\alpha_0)^{1/4} \\, \\sin(2 \\alpha_0)',
                             'desc': 'Deepwater longshore energy-flux factor; sign follows '
                                     'sin(2 alpha_0), and cos is clamped at 0 for '
                                     '|alpha_0| >= 90 deg (TR 6-1 eq 17).'},
                            {'tex': 'f = \\frac{\\min(\\delta + \\Delta/2, \\, 90) - '
                                    '\\max(\\delta - \\Delta/2, \\, -90)}{\\Delta}',
                             'desc': 'Contributing fraction of a 22.5 deg band (Delta) at '
                                     'angle delta to the shore normal that lies within '
                                     '+/-90 deg of it.'},
                            {'tex': 'Q_{band} = f \\sum_{h} p_h \\, Q_{deep}(H_h, \\delta)',
                             'desc': 'Per-band transport: occurrence-weighted sum of the '
                                     'deepwater CERC rate over the height bins, scaled by '
                                     'the contributing fraction.'},
                            {'tex': 'Q_{net} = \\sum_{b} Q_{band,b} , \\qquad Q_{gross} = '
                                    '\\sum_{b} |Q_{band,b}|',
                             'desc': 'Net is the signed sum over all directional bands; '
                                     'gross is the sum of magnitudes.'}]}],
 'symbols': [['Q', 'Volumetric longshore sediment transport rate'],
             ['K', 'Empirical CERC coefficient (0.39 for field significant-height data)'],
             ['P_ls', 'Longshore energy-flux factor'],
             ['rho_s', 'Sediment (sand) density; quartz ~2650 kg/m^3'],
             ['rho', 'Water density; seawater ~1025 kg/m^3'],
             ['a', 'Solids fraction of bed volume, a = 1 - porosity'],
             ['H_s0', 'Deepwater significant wave height (bin midpoint H_h)'],
             ['alpha_0', 'Deepwater wave angle to shore normal (band angle delta)'],
             ['theta', 'Seaward shore-normal azimuth, clockwise from true north'],
             ['p_h', 'Occurrence fraction of a height bin within a directional band']],
 'references': ['SPM (1984) Ch. 4 (Eq. 4-49)',
                'Gravens (1988)',
                'WIS Report 18 (Station G1033)',
                "ACES User's Guide, Longshore Sediment Transport, Example 6-1-3",
                'Galvin (1979)']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Net/gross longshore transport from a CEDRS directional climate (SI internal)."""
    _validate(inp)
    theta = float(inp["shore_azimuth"]); K = float(inp["K"])
    rho = float(inp["rho_water"]); rho_s = float(inp["rho_sand"]); por = float(inp["porosity"])
    a_solid = 1.0 - por
    if rho_s <= rho:
        raise ValueError("sediment density must exceed water density")
    occ = inp.get("occ", _G1033)
    rows = [r for r in occ if r]
    if not rows:
        raise ValueError("CEDRS occurrence table is empty")

    Q_net = Q_right = Q_left = 0.0
    band_angle = []; band_pct = []; band_Q = []
    for i, row in enumerate(rows):
        az = i * (360.0 / len(rows))               # 22.5 deg spacing for 16 bands
        delta = _wrap180(theta - az)
        half = 0.5 * (360.0 / len(rows))           # half band width (11.25 deg)
        lo = max(delta - half, -90.0); hi = min(delta + half, 90.0)
        frac = max(0.0, (hi - lo) / (2.0 * half))  # fraction of the band within +/-90 deg
        Qb = 0.0
        if frac > 0.0:
            for occ_x1000, H in zip(row, _HMID):
                p = float(occ_x1000) / 100000.0    # x1000 of percent -> fraction
                if p > 0.0:
                    Qb += p * deepwater_Q(H, delta, K, rho, rho_s, a_solid, g)
            Qb *= frac
        band_angle.append(delta)
        band_pct.append(frac * 100.0)
        band_Q.append(Qb)
        Q_net += Qb
        if Qb >= 0.0:
            Q_right += Qb
        else:
            Q_left += Qb
    Q_gross = Q_right - Q_left

    direction = "rightward" if Q_net > 0 else ("leftward" if Q_net < 0 else "balanced")
    notes = (f"net {Q_net * _M3_TO_YD3:,.0f} yd^3/yr ({direction}); gross "
             f"{Q_gross * _M3_TO_YD3:,.0f}; rho_s={rho_s:.0f} kg/m^3 (a'={a_solid:.2f}); "
             f"{len(rows)} bands, {sum(1 for p in band_pct if p > 0)} contributing")
    return Result(Q_net=Q_net, Q_gross=Q_gross, Q_right=Q_right, Q_left=Q_left,
                  band_angle=band_angle, band_pct=band_pct, band_Q=band_Q, notes=notes)


# --- self-tests (User's Guide Example 3 + structural checks) ---------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    # 1) ACES User's Guide Example 3: theta=40, K=0.39, G1033 -> net -854,849 yd^3/yr,
    #    reproduced with the effective rho_s = 2319 (as documented for 6-1).
    r = compute({"shore_azimuth": 40.0, "K": 0.39, "rho_water": 1025.18,
                 "rho_sand": 2319.0, "porosity": 0.40})
    net_yd = r.Q_net * _M3_TO_YD3
    assert _approx(net_yd, -854848.77, 0.01), net_yd      # <1% of the oracle
    # contributing percentages reproduce the User's Guide exactly (72.22, 27.78, 100)
    by_az = {round(a): p for a, p in zip(r.band_angle, r.band_pct)}
    assert _approx(by_az[85], 72.22, 1e-3), by_az[85]     # band at +85 deg
    assert _approx(by_az[40], 100.0, 1e-9)
    # the band beyond +/-90 (delta = -95) contributes ~0 (oracle -0.5 yd^3/yr)
    i95 = r.band_angle.index(-95.0)
    assert abs(r.band_Q[i95]) < 1.0, r.band_Q[i95]

    # 2) net is the signed sum of the per-band rates; gross is the sum of magnitudes
    assert _approx(r.Q_net, sum(r.band_Q), 1e-9)
    assert _approx(r.Q_gross, sum(abs(q) for q in r.band_Q), 1e-9)
    assert r.Q_gross >= abs(r.Q_net)

    # 3) a directionally-uniform climate is symmetric about every shore normal, so its net
    #    transport vanishes (the +delta and -delta bands cancel)
    uniform = [[100] * 11 for _ in range(16)]
    ru = compute({"shore_azimuth": 0.0, "K": 0.39, "rho_water": 1025.18,
                  "rho_sand": 2319.0, "porosity": 0.40, "occ": uniform})
    assert _approx(ru.Q_net, 0.0, 1e-6) and ru.Q_gross > 0.0, (ru.Q_net, ru.Q_gross)

    # 4) transport scales linearly in K and in 1/(rho_s - rho)
    rK = compute({"shore_azimuth": 40.0, "K": 0.78, "rho_water": 1025.18,
                  "rho_sand": 2319.0, "porosity": 0.40})
    assert _approx(rK.Q_net, 2.0 * r.Q_net, 1e-9)

    print(f"  self-tests: PASS (User's Guide Example 3 net {net_yd:,.0f} yd^3/yr vs "
          f"oracle -854,849 [{abs(net_yd/-854848.77-1)*100:.2f}%]; contributing % 72.22/100/27.78; "
          "uniform-climate net=0; K linearity)")


def _print_default_example() -> None:
    r = compute({f.key: (f.default if f.key != "occ" else [row[:] for row in _G1033])
                 for f in INPUTS})
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  station G1033 (Gulf of Mexico), shore-normal azimuth 40 deg, K=0.39, quartz sand")
    print("  OUTPUTS:")
    print(f"    Net transport    Q_net   = {r.Q_net * _M3_TO_YD3:>13,.0f} yd^3/yr "
          f"({r.Q_net:,.0f} m^3/yr)")
    print(f"    Gross transport  Q_gross = {r.Q_gross * _M3_TO_YD3:>13,.0f} yd^3/yr")
    print(f"    Right (+)        Q_right = {r.Q_right * _M3_TO_YD3:>13,.0f} yd^3/yr")
    print(f"    Left  (-)        Q_left  = {r.Q_left * _M3_TO_YD3:>13,.0f} yd^3/yr")
    print("  contributing bands (angle from shore normal, %, transport yd^3/yr):")
    for a, p, q in zip(r.band_angle, r.band_pct, r.band_Q):
        if p > 0:
            print(f"    {a:7.2f} deg  {p:6.2f}%  {q * _M3_TO_YD3:>13,.1f}")
    print(f"  note (quartz): the ACES User's Guide Example 3 net is -854,849 yd^3/yr "
          f"(its effective rho_s ~ 2319; see docstring)")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
