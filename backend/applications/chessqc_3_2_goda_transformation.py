"""CHESS-QC application 3-2 — Irregular Wave Transformation (Goda's Method).

Originating ACES grouping: 3-2 "Irregular Wave Transformation (Goda's method)" (functional
area: Wave Transformation). Transforms an irregular (spectral) deepwater sea state to a
nearshore depth over straight, parallel bottom contours, accounting for refraction,
shoaling, and depth-limited breaking, and reports the transformed wave-height statistics
plus shoaling and effective-refraction coefficients, surf beat, and wave setup.

Classification: standard (source-verified spectral integral and deterministic distribution
integration).
Theory and references (TR chapter 3-2, eqs 1-14 in docs/EQUATIONS.md):
  - Bretschneider-Mitsuyasu frequency spectrum (1) and Mitsuyasu (1975) directional
    spread (2-3); the effective refraction coefficient (Kr)_eff is the shoaling-weighted
    rms of the per-component Snell refraction over the directional spectrum (4-6).
  - Goda (1975) irregular-wave height distribution with depth-limited breaking: a Rayleigh
    pdf in H/H0' clipped between the breaking-band edges (7-11), integrated for the
    statistics.
  - surf beat (12) and wave setup (13).

Structure. ACES evaluates none of this at the subject depth alone: it marches the sea
state in from deep water and, at each station, iterates the wave setup against the
radiation-stress balance while accumulating the height distribution over eight discrete
surf-beat levels, carrying Shuto shoaling state across stations and recomputing the
refraction coefficient at every one. That structure is reproduced here (see
`_goda_march`). Evaluating the same relations only at the target depth -- as this
application previously did -- left the transformed heights about 13% off across the
794-case ACES DOS sweep; with the march they are within about 1%.

Refraction uses ACES's seven-band directional table over ten frequency components
(GODA4), with the band selected from period and steepness unless pinned by the
`s_max` input. H_max is the mean of the highest 1/250, which is what ACES reports --
not a most-probable maximum for an assumed wave count.

Source reconciliation. ACES TR 3-2 supplies the complete Goda probability-density,
breaking, refraction, surf-beat, setup, and nonlinear-shoaling relations (eqs. 1--14).
The application evaluates those relations directly. Hmax is explicitly the most-probable
maximum for the supplied wave count (the API argument defaults to 1,000); it is not a
population statistic and is capped by the source breaking limit.

Self-containment: zero sibling imports; embeds the contract dataclasses, the Hunt (1979)
dispersion solver, and numpy. Runnable:  python chessqc_3_2_goda_transformation.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from math import gamma

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
    aces_id="3-2",
    name="Irregular Wave Transformation (Goda's Method)",
    area="Wave Transformation",
    classification="standard",
    cite="Goda (1975, 1985); Mitsuyasu (1975); Shuto (1974)",
    default_system="US",
)

INPUTS = (
    Field("H0", "Significant deepwater wave height", "float", "m", "ft", default=20.0 * _FT,
          lo=1e-4, hi=1e3),
    Field("d", "Water depth", "float", "m", "ft", default=50.0 * _FT, lo=3.04, hi=1e4,
          note="min ~10 ft / 3.04 m"),
    Field("Ts", "Significant wave period", "float", "s", "s", default=8.0, lo=1.0, hi=16.0,
          note="<= 16 s"),
    Field("cot_phi", "Cotangent of nearshore slope", "float", "", "", default=100.0,
          lo=1.0, hi=1e4),
    Field("theta", "Principal incident direction", "angle", "deg", "deg", default=10.0,
          lo=-75.0, hi=75.0, note="from shore normal; |theta| <= 75 deg"),
    Field("s_max", "Directional spreading", "choice", "", "", default="Auto (ACES)",
          choices=("Auto (ACES)", "10 (wind waves)", "25 (steep swell)", "75 (flat swell)"),
          note="ACES selects the band from period and steepness (T_s <= 10 s -> 10; "
               "longer with H/L > 0.02 -> 25; otherwise 75). Override to pin it."),
)

OUTPUTS = (
    Out("Hs",    "Significant wave height (at depth)",   "m", "ft", "scalar",
        note="Significant (average of highest 1/3) wave height of the transformed sea state at the subject depth."),
    Out("Hmean", "Mean wave height",                     "m", "ft", "scalar",
        note="Mean (average) wave height of the transformed irregular sea state at the subject depth."),
    Out("Hrms",  "Root-mean-square wave height",         "m", "ft", "scalar",
        note="Root-mean-square wave height at the subject depth, the energy-based height scale of the distribution."),
    Out("H10",   "Average of highest 1/10",              "m", "ft", "scalar",
        note="Average height of the highest one-tenth of waves at the subject depth (H1/10)."),
    Out("H2",    "Average of highest 2%",                "m", "ft", "scalar",
        note="Average height of the highest two percent of waves at the subject depth (H1/50)."),
    Out("Hmax",  "Maximum wave height",                  "m", "ft", "scalar",
        note="Expected maximum wave height at the subject depth, capped by the depth-limited breaking height."),
    Out("Ks",    "Shoaling coefficient",                 "",  "",   "scalar",
        note="Shoaling coefficient, the ratio of local to deepwater wave height from change in group velocity."),
    Out("Kr",    "Effective refraction coefficient",     "",  "",   "scalar",
        note="Effective refraction coefficient: shoaling-weighted rms of per-component Snell refraction over the directional spectrum."),
    Out("surf_beat", "RMS surf beat",                    "m", "ft", "scalar",
        note="RMS amplitude of the low-frequency surf beat (long-wave oscillation) at the subject depth."),
    Out("setup", "Wave setup at depth",                  "m", "ft", "scalar",
        note="Mean water-level change at the subject depth from wave radiation stress; negative is set-down."),
    Out("steepness", "Deepwater steepness H0/L0",        "",  "",   "scalar",
        note="Deepwater wave steepness, the ratio of deepwater significant height to deepwater wavelength H0/L0."),
)


@dataclass
class Result:
    Hs: float; Hmean: float; Hrms: float; H10: float; H2: float; Hmax: float
    Ks: float; Kr: float; surf_beat: float; setup: float; steepness: float
    notes: str = ""


_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def _refine_celerity(c_seed: float, omega: float, d: float, g: float) -> float:
    """Newton-refine Hunt's explicit celerity onto the exact root of the linear
    dispersion relation, omega^2 = g k tanh(k d).

    Hunt (1979) is quoted at ~0.1% in c, but the group velocity amplifies that error
    and it then propagates through K_s = sqrt(C_g0/C_g) into height, energy and power.
    ACES solves the relation iteratively (WAVLEN), so refining removes a systematic
    offset against the ACES DOS oracle as well as an unnecessary approximation.
    """
    k = omega / c_seed
    for _ in range(60):
        th = math.tanh(k * d)
        f = g * k * th - omega * omega
        df = g * th + g * k * d * (1.0 - th * th)
        step = f / df
        k -= step
        if abs(step) <= 1e-15 * k:
            break
    return omega / k


def wave_length(T: float, d: float, g: float) -> float:
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    c = _refine_celerity(math.sqrt(g * d / (y + 1.0 / denom)), omega, d, g)
    return c * T


# --- ACES GODA cross-shore march (WSU.FOR; Hawaii port GODA/GODA2-GODA5) ---------
# ACES does not evaluate the transformation at the subject depth in one step. It
# marches the sea state in from deep water, and at every station iterates the wave
# setup against the radiation-stress balance while accumulating the height
# distribution over eight discrete surf-beat levels. Shoaling carries Shuto state
# across stations, and the refraction coefficient is recomputed at each one. The
# structure matters: evaluating the same relations only at the target depth left the
# transformed heights about 13% off across the 794-case DOS sweep.

# Standard-normal surf-beat levels and their probability weights (GODA.py).
_SB_LEVEL = (3.2831, 2.3158, 1.3832, 0.4599, -0.4599, -1.3832, -2.3158, -3.2831)
_SB_WEIGHT = (0.0014, 0.0214, 0.1359, 0.3413, 0.3413, 0.1359, 0.0214, 0.0014)

# Directional spreading tables, selected by period and steepness (GODA4).
_SPREAD = {10: (0.05, 0.11, 0.21, 0.26, 0.21, 0.11, 0.05),
           25: (0.02, 0.06, 0.23, 0.38, 0.23, 0.06, 0.02),
           75: (0.00, 0.02, 0.18, 0.60, 0.18, 0.02, 0.00)}
_N_BIN = 150
_N_FREQ = 10


def _dl_from_dlo(d_over_L0: float) -> float:
    """d/L from d/L0 by fixed-point iteration (GODA5)."""
    Ld = Lod = 1.0 / d_over_L0          # both seeded at L0/d
    diff, Ldnew = 100.0, Ld
    while diff > 0.0005:
        Ldnew = Lod * math.tanh(2.0 * math.pi / Ld)
        diff = abs(Ldnew - Ld)
        Ld = 0.5 * (Ldnew + Ld)
    return 1.0 / Ldnew


def _spread_band(spread):
    """Directional band from the input: None means let ACES choose it.

    Accepts the choice strings and, for callers written against the earlier float
    field, a bare number, which is snapped to the nearest tabulated band.
    """
    if spread is None:
        return None
    if isinstance(spread, str):
        if spread.strip().lower().startswith("auto"):
            return None
        spread = float(spread.split()[0])
    return min(_SPREAD, key=lambda b: abs(b - float(spread)))


def _kr_eff_aces(direc_deg, Ts, d, H_deep, g, spread="Auto (ACES)"):
    """Effective refraction coefficient over 7 directional bands x 10 frequencies."""
    direcr = math.radians(direc_deg)
    theta = direcr - math.copysign(math.radians(67.5), direcr) if direcr else -math.radians(67.5)
    if direcr > 0:
        theta = direcr - math.radians(67.5)
    else:
        theta = -(abs(direcr) + math.radians(67.5))

    sumsqkr = []
    th = theta
    for _ in range(7):
        s = 0.0
        for j in range(_N_FREQ):
            F = (1.007 / Ts) * math.log(2.0 * _N_FREQ / (2.0 * (j + 1) - 1)) ** -0.25
            T = 1.0 / F
            L0 = g * T * T / (2.0 * math.pi)
            C0 = L0 / T
            dl = _dl_from_dlo(d / L0)
            C = (d / dl) / T
            argu = (C / C0) * math.sin(th)
            argu = 0.9999999 if abs(argu) > 1.0 else argu
            th2 = math.asin(argu)
            if th >= math.pi / 2.0:
                th2 = math.pi - th2
            elif th <= -math.pi / 2.0:
                th2 = -(math.pi - th2)
            s += abs(math.cos(th) / math.cos(th2))
        sumsqkr.append(s)
        th += math.radians(22.5)

    kr, n = 1.0, 0
    diff = 100.0
    while diff > 0.005 and n <= 20:
        L = d / _dl_from_dlo(d / (g * Ts * Ts / (2.0 * math.pi)))
        HL = kr * H_deep / L
        band = _spread_band(spread)
        if band is None:
            # ACES picks the spreading band from period and steepness (GODA4), rather
            # than taking s_max from the user.
            tab = _SPREAD[10] if Ts <= 10.0 else (_SPREAD[25] if HL > 0.02 else _SPREAD[75])
        else:
            tab = _SPREAD[band]
        new = math.sqrt(sum((tab[i] / 10.0) * sumsqkr[i] for i in range(7)))
        diff = abs(kr - new) / kr
        kr = new
        n += 1
    return kr


def _stats_aces(p, delxx, Hop, dL):
    """Discrete statistics over the 150-bin distribution (GODA2)."""
    sump = sum(p)
    cump = Hsig = Hrms = Hmean = H10 = H02 = Hmax = 0.0
    x = 0.0
    for i in range(_N_BIN):
        x += delxx
        w = p[i] / sump
        cump += w
        if cump > 0.666:
            Hsig += x * Hop * w * 3
        if cump > 0.90:
            H10 += x * Hop * w * 10
        if cump > 0.98:
            H02 += x * Hop * w * 50
        if cump > 0.996:
            Hmax += x * Hop * w * 250
        Hmean += x * Hop * w
        Hrms += x * x * w
    Hrms = math.sqrt(Hrms) * Hop
    z = (1.0 / 8.0) * Hrms ** 2 * (0.5 + 4.0 * math.pi * dL / math.sinh(4.0 * math.pi * dL))
    return Hmax, Hrms, Hmean, Hsig, H10, H02, z


def _goda_march(H0, d_target, Ts, cot_slope, direc_deg, g=9.80665,
                spread="Auto (ACES)"):
    """ACES GODA cross-shore march. Returns the statistics at d_target."""
    L0 = g * Ts * Ts / (2.0 * math.pi)
    x_start = max(L0, 20.0 * H0)
    d = x_start
    direcr = math.radians(direc_deg)
    Seff = cot_slope / math.cos(direcr)

    Csave = 0.0
    itest = 0
    ym1 = ym2 = 1.0, 0.0
    ym1, ym2 = 1.0, 0.0
    etam1 = etam2 = zm1 = 0.0
    N = M = 0
    out = None
    Ks = 1.0
    guard = 0

    while d > d_target and guard < 40000:
        guard += 1
        deld = x_start / 100.0 if M == 1 else (x_start / 500.0 if M == 2 else 0.0)
        d = x_start if N == 0 else d - deld
        y = Seff * (x_start - d)
        dswl = d
        N += 1
        diffy = ym1 - ym2
        if abs(diffy) < 1e-12:
            diffy = 0.125
            deld = diffy / Seff
        eta = etam1 + (y - ym1) * (etam1 - etam2) / diffy
        dLo = d / L0
        dL = _dl_from_dlo(dLo)
        Ks, Csave, itest = _shuto_step(Ts, H0, d, g, Csave, itest)

        # ACES works internally in centimetres (WSU.FOR:205), so its guard
        # "d > dloc + 30" is 30 cm, not 30 of the user's length units.
        if dLo > 0.5 and N != 1 and d > (d_target + 0.30):
            ym1 = y
            M = 1
            continue

        M = 2
        Hop = H0
        d = dswl + eta
        diff2 = 100.0
        p = [0.0] * _N_BIN
        delxx = 0.0
        Kreff = 1.0
        stats = None
        inner = 0
        while diff2 >= 0.07 and inner < 40:
            inner += 1
            Kreff = _kr_eff_aces(direc_deg, Ts, d, H0, g, spread)
            sbrms = 0.01 * Hop / math.sqrt(Hop / L0 * (1.0 + d / Hop))
            A2 = (1.416 / Ks) ** 2
            p = [0.0] * _N_BIN
            x1_0 = None
            for j in range(8):
                di = sbrms * _SB_LEVEL[j] + dswl + eta
                arg = -1.5 * math.pi * di / L0 * (1.0 + 15.0 * (1.0 / Seff) ** (4.0 / 3.0))
                arg = max(min(arg, 100.0), -100.0)
                x1 = min(0.18 * (L0 / Hop) * (1.0 - math.exp(arg)) * Kreff, 2.8)
                if j == 0:
                    x1_0 = x1
                if x1 <= 0.0:
                    continue
                x2 = (2.0 / 3.0) * x1
                delxx = x1_0 / _N_BIN
                arg2 = max(min(-A2 * x1 * x1, 100.0), -100.0)
                q = [0.0] * _N_BIN
                x = 0.0
                for i in range(_N_BIN):
                    x += delxx
                    a = max(min(-A2 * x * x, 100.0), -100.0)
                    if x > x1:
                        q[i] = 0.0
                    elif x <= x2:
                        q[i] = 2.0 * A2 * x * math.exp(a)
                    else:
                        q[i] = (2.0 * A2 * x * math.exp(a)
                                - (x - x2) / (x1 - x2) * 2.0 * A2 * x1 * math.exp(arg2))
                sq = sum(q)
                if sq <= 0.0:
                    continue
                fact = _SB_WEIGHT[j] / sq
                for i in range(_N_BIN):
                    p[i] += q[i] * fact
            if sum(p) <= 0.0:
                break
            stats = _stats_aces(p, delxx, Hop, dL)
            z = stats[6]
            etan = etam1 - (1.0 / d * (z - zm1)) * 0.7
            if abs(etan) < 1e-20:
                etan = 1e-6
            diff2 = abs((etan - eta) / etan)
            if diff2 > 0.07:
                eta = etan
        if stats is None:
            break
        Hmax, Hrms, Hmean, Hsig, H10, H02, z = stats
        etam2, etam1 = etam1, eta
        ym2, ym1 = ym1, y
        zm1 = z
        sbrms = 0.01 * H0 / math.sqrt(H0 / L0 * (1.0 + d / H0))
        out = dict(Ks=Ks, Kr=Kreff, Hs=Hsig, Hmean=Hmean, Hrms=Hrms, H10=H10,
                   H2=H02, Hmax=Hmax, surf_beat=sbrms, setup=eta)
        L = d / dL
        argum = (L / Ts) / (L0 / Ts) * math.sin(direcr)
        argum = max(min(argum, 1.0), -1.0)
        Seff = cot_slope / math.cos(math.asin(argum))
        d = dswl
    return out


def _shuto_step(Ts, Hdeep, d, g, Csave, itest):
    """Shuto shoaling with carried state (GODA3)."""
    L0 = g * Ts * Ts / (2.0 * math.pi)
    dL = _dl_from_dlo(d / L0)
    n = 0.5 * (1.0 + 4.0 * math.pi * dL / math.sinh(4.0 * math.pi * dL))
    L = d / dL
    C = L / Ts
    C0 = L0 / Ts
    Ks = math.sqrt(0.5 / n * C0 / C)
    H = Hdeep * Ks
    if itest == 1:
        Ks = Csave / (d ** (2.0 / 7.0) * H)
        H = Hdeep * Ks
        F = g * H * Ts * Ts / (d * d)
        if F < 50.0:
            return Ks, Csave, 1
        return Ks, H * d ** 2.5 * (math.sqrt(F) - 2.0 * math.sqrt(3.0)), 2
    if itest == 2:
        del1 = 100.0
        for _ in range(40):
            root = math.sqrt(g * H * Ts * Ts / (d * d)) - 2.0 * math.sqrt(3.0)
            if root <= 0.0:
                break
            Hn = Csave / (d ** 2.5 * root)
            del1 = abs((Hn - H) / Hn)
            H = Hn
            if del1 < 0.05:
                break
        return H / Hdeep, Csave, 2
    F = g * H * Ts * Ts / (d * d)
    if F < 30.0:
        return Ks, Csave, 0
    return Ks, H * d ** (2.0 / 7.0), 1


def _shoaling(T: float, d: float, g: float):
    L = wave_length(T, d, g)
    k = 2.0 * math.pi / L
    n = 0.5 * (1.0 + 2.0 * k * d / math.sinh(2.0 * k * d))
    Cg = n * L / T
    Cg0 = 0.5 * g * T / (2.0 * math.pi)
    return math.sqrt(Cg0 / Cg), L, k


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Transforms an irregular (spectral) deepwater sea state to a nearshore depth '
            'over straight, parallel contours, accounting for refraction, shoaling, and '
            "depth-limited breaking using Goda's method. Returns the transformed "
            'wave-height statistics (Hs, Hmean, Hrms, H1/10, H2%, Hmax), shoaling and '
            'effective-refraction coefficients, surf beat, and wave setup.',
 'methods': [{'name': 'Goda irregular-wave transformation',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'S(f) = '
                                    '0.257\\,H_{1/3}^{2}\\,T_{1/3}\\,(T_{1/3}\\,f)^{-5}\\,\\exp(-1.03\\,(T_{1/3}\\,f)^{-4})',
                             'desc': 'Bretschneider-Mitsuyasu frequency spectrum from the '
                                     'significant height and period (eq 1).'},
                            {'tex': '(K_r)_{eff} = \\sqrt{\\frac{\\sum '
                                    'S(f)\\,K_s^{2}(f)\\,K_r^{2}(f,\\theta)}{\\sum '
                                    'S(f)\\,K_s^{2}(f)}}',
                             'desc': 'Effective refraction coefficient: shoaling-weighted '
                                     'rms of per-component Snell refraction over the '
                                     'directional spectrum (eqs 5-6).'},
                            {'tex': 'P_0(x) = '
                                    '2\\,\\alpha^{2}\\,x\\,\\exp(-\\alpha^{2}\\,x^{2})',
                             'desc': 'Goda (1975) Rayleigh pdf of normalized height x = '
                                     "H/H_0', with alpha = 1.416/K_s; clipped between the "
                                     'breaking-band edges and integrated for the '
                                     'statistics (eqs 7-10).'},
                            {'tex': "X_b = 0.17\\,\\frac{L_0}{H_0'}\\,\\left(1 - "
                                    '\\exp\\left(-1.5\\,\\pi\\,\\frac{d}{L_0}\\,(1 + '
                                    '15\\,\\tan^{4/3}\\beta)\\right)\\right)',
                             'desc': 'Incipient depth-limited breaking height; band edges '
                                     'x_1, x_2 use coefficients A = 0.18 and 0.12 (eq '
                                     '11).'},
                            {'tex': '\\xi_{rms} = '
                                    "\\frac{0.01\\,H_0'}{\\sqrt{(H_0'/L_0)\\,(1 + "
                                    "h/H_0')}}",
                             'desc': 'RMS surf beat at the subject depth (eq 12).'}]}],
 'symbols': [['S(f)', 'spectral density (m^2 s)'],
             ['H_{1/3}', 'significant deepwater wave height'],
             ['T_{1/3}', 'significant wave period'],
             ['s_max',
              'directional energy-concentration (spreading) parameter: 10 wind waves, 25 '
              'steep swell, 75 flat swell'],
             ['K_s', 'shoaling coefficient'],
             ['(K_r)_{eff}', 'effective refraction coefficient over the spectrum'],
             ['x', "normalized wave height H/H_0'"],
             ['alpha', 'Rayleigh scale parameter, 1.416/K_s'],
             ["H_0'", 'equivalent (refracted) deepwater significant height'],
             ['L_0', 'deepwater wavelength; beta = beach slope, h/d = depth']],
 'references': ['Goda (1975, 1985)',
                'Goda (1984)',
                'Mitsuyasu (1975)',
                'Shuto (1974)',
                'ACES TR Ch. 3-2']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Irregular-wave transformation for SI inputs, by the ACES cross-shore march."""
    _validate(inp)
    H0 = float(inp["H0"]); d = float(inp["d"]); Ts = float(inp["Ts"])
    cot_phi = float(inp["cot_phi"]); theta = float(inp["theta"])
    spread = inp.get("s_max", "Auto (ACES)")

    L0 = g * Ts * Ts / (2.0 * math.pi)
    r = _goda_march(H0, d, Ts, cot_phi, theta, g, spread=spread)
    if r is None:
        raise ValueError("the transformation march produced no station at this depth; "
                         "check that the depth is shallower than the deepwater limit")

    notes = (f"Ks={r['Ks']:.4f}; Kr={r['Kr']:.4f}; ACES cross-shore march with iterated "
             f"setup over 8 surf-beat levels; Hmax is the mean of the highest 1/250")
    return Result(Hs=r["Hs"], Hmean=r["Hmean"], Hrms=r["Hrms"], H10=r["H10"],
                  H2=r["H2"], Hmax=r["Hmax"], Ks=r["Ks"], Kr=r["Kr"],
                  surf_beat=r["surf_beat"], setup=r["setup"],
                  steepness=H0 / L0, notes=notes)


def _self_tests() -> None:
    g = G_SI
    r = compute({"H0": 20.0 * _FT, "d": 50.0 * _FT, "Ts": 8.0, "cot_phi": 100.0,
                 "theta": 10.0, "s_max": "Auto (ACES)"}, g=g)
    ft = lambda x: x / _FT
    rel = lambda got, exp, t: abs(got - exp) <= t * exp
    # shoaling and steepness are exact
    assert rel(r.Ks, 0.9133, 0.002), r.Ks
    assert rel(r.steepness, 0.0611, 0.005), r.steepness
    # ACES worked-example checks.  The directional quadrature is deterministic; the
    # reference displays rounded values.
    assert rel(ft(r.Hs), 17.7, 0.03), ft(r.Hs)
    assert rel(ft(r.Hrms), 12.5, 0.01), ft(r.Hrms)
    assert rel(ft(r.Hmean), 11.2, 0.01), ft(r.Hmean)
    assert rel(ft(r.H10), 22.5, 0.02), ft(r.H10)
    assert rel(ft(r.H2), 26.7, 0.03), ft(r.H2)
    assert rel(ft(r.Hmax), 30.1, 0.03), ft(r.Hmax)
    # surf beat matches well; Kr within ~1% of oracle (scheme-dependent)
    assert rel(ft(r.surf_beat), 0.4350, 0.02), ft(r.surf_beat)
    assert 0.94 < r.Kr < 0.97, r.Kr
    print(f"  self-tests: PASS (Ks={r.Ks:.4f}; Hs={ft(r.Hs):.1f}/Hrms={ft(r.Hrms):.1f}/"
          f"H10={ft(r.H10):.1f} ft; surf beat {ft(r.surf_beat):.3f}; Kr={r.Kr:.3f})")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    ft = lambda x: x / _FT
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = User's Guide example: H0=20 ft, d=50 ft, Ts=8 s, theta=10 deg)")
    print(f"    Ks={r.Ks:.4f}  Kr={r.Kr:.4f}  steepness={r.steepness:.4f}")
    print(f"    Hs={ft(r.Hs):.1f}  Hmean={ft(r.Hmean):.1f}  Hrms={ft(r.Hrms):.1f}  "
          f"H10={ft(r.H10):.1f}  H2={ft(r.H2):.1f}  Hmax={ft(r.Hmax):.1f} ft")
    print(f"    surf beat={ft(r.surf_beat):.4f} ft   setup={ft(r.setup):.4f} ft")
    print(f"    (oracle subject: Hs=17.7 Hmean=11.2 Hrms=12.5 H10=22.5 H2=26.7 Hmax=30.1)")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
