"""CHESS-QC application 3-2 — Irregular Wave Transformation (Goda's Method).

Originating ACES grouping: 3-2 "Irregular Wave Transformation (Goda's method)" (functional
area: Wave Transformation). Transforms an irregular (spectral) deepwater sea state to a
nearshore depth over straight, parallel bottom contours, accounting for refraction,
shoaling, and depth-limited breaking, and reports the transformed wave-height statistics
plus shoaling and effective-refraction coefficients, surf beat, and wave setup.

Classification: provisional (spectral integral + Monte-Carlo-free distribution integration).
Theory and references (TR chapter 3-2, eqs 1-14 in docs/EQUATIONS.md):
  - Bretschneider-Mitsuyasu frequency spectrum (1) and Mitsuyasu (1975) directional
    spread (2-3); the effective refraction coefficient (Kr)_eff is the shoaling-weighted
    rms of the per-component Snell refraction over the directional spectrum (4-6).
  - Goda (1975) irregular-wave height distribution with depth-limited breaking: a Rayleigh
    pdf in H/H0' clipped between the breaking-band edges (7-11), integrated for the
    statistics.
  - surf beat (12) and wave setup (13); nonlinear shoaling by Shuto (1974) (14).

Documented accuracy limitation. This app reproduces the ACES worked example's shoaling
coefficient exactly and its at-depth ("subject") significant/mean/rms heights and surf beat
to about three percent, but the high quantiles (H1/10, H1/50, Hmax) carry larger residuals
(up to ~6%) and cannot be matched to the digit from any published relation. The User's
Guide's own *plotted* deepwater distribution (Table 3-2-1) is pure Rayleigh -- at H=26.24 ft
its exceedance is 0.030 = exp(-(26.24/14.1)^2) with Hrms=14.1 ft -- so the average of the
highest tenth implied by that distribution is 1.80*Hrms = 25.4 ft, yet the *screen output*
reports H1/10 = 27.0 ft (a finite-N order-statistic average of that same Rayleigh tail
yields ~25.4 ft for any N). The screen high-quantiles therefore use an ACES-specific
finite-N largest-wave inflation (Hmax = Hrms*sqrt(ln N) with an undocumented effective wave
count N ~ 1200) that is inconsistent with the plotted distribution and is not recoverable
from the public Technical Reference. Likewise the effective refraction coefficient depends
on the exact directional-integration scheme and the default spreading parameter s_max (the
example's 0.9638 corresponds to s_max ~ 13). The transformation physics is implemented in
full; the residual is in the empirical coefficients, not the physics.

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
    classification="provisional",
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
    Field("s_max", "Directional spreading parameter", "float", "", "", default=10.0,
          lo=1.0, hi=200.0, note="10 wind waves, 25 steep swell, 75 flat swell"),
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


def wave_length(T: float, d: float, g: float) -> float:
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    return math.sqrt(g * d / (y + 1.0 / denom)) * T


def _shoaling(T: float, d: float, g: float):
    L = wave_length(T, d, g)
    k = 2.0 * math.pi / L
    n = 0.5 * (1.0 + 2.0 * k * d / math.sinh(2.0 * k * d))
    Cg = n * L / T
    Cg0 = 0.5 * g * T / (2.0 * math.pi)
    return math.sqrt(Cg0 / Cg), L, k


def _bm_spectrum(f: float, H0: float, Ts: float) -> float:
    u = Ts * f
    return 0.257 * H0 * H0 * Ts * u ** (-5) * math.exp(-1.03 * u ** (-4))


def _kr_eff(H0, Ts, d, theta_deg, s_max, g, nf=140, ndir=221):
    """Effective refraction coefficient over the directional spectrum (eqs 1-6)."""
    L0 = g * Ts * Ts / (2.0 * math.pi)
    theta_p = math.radians(theta_deg)
    fp = 0.9529 / Ts                                   # B-M spectral peak frequency
    fs = np.linspace(0.4 * fp, 3.0 * fp, nf)
    dth = np.linspace(-math.pi / 2.0, math.pi / 2.0, ndir)
    num = den = 0.0
    for f in fs:
        S = _bm_spectrum(f, H0, Ts)
        Ks, L, _ = _shoaling(1.0 / f, d, g)
        s = s_max * (f / fp) ** 5 if f <= fp else s_max * (f / fp) ** (-2.5)   # Mitsuyasu
        w0 = S * Ks * Ks
        for td in dth:
            spread = max(math.cos(td / 2.0), 0.0) ** (2.0 * s)
            a0 = theta_p + td
            if abs(a0) >= math.pi / 2.0:
                continue
            sind = math.sin(a0) * (L / L0)
            if abs(sind) >= 1.0:
                continue
            ad = math.asin(sind)
            Kr = math.sqrt(max(math.cos(a0), 1e-9) / max(math.cos(ad), 1e-9))
            num += w0 * spread * Kr * Kr
            den += w0 * spread
    return math.sqrt(num / den)


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


def _stats_at(H0e, Ks, d, L0, tanb, n_waves):
    """Integrate the breaking-clipped distribution; return (Hmean,Hrms,Hs,H10,H2,Hmax)."""
    a = 1.416 / Ks
    def Xb(A):
        return A * (L0 / H0e) * (1.0 - math.exp(
            -1.5 * math.pi * (d / L0) * (1.0 + 15.0 * tanb ** (4.0 / 3.0))))
    x1 = Xb(0.18); x2 = Xb(0.12)
    P0 = lambda x: 2.0 * a * a * x * math.exp(-a * a * x * x)
    def Pr(x):
        if x <= x2:
            return P0(x)
        if x < x1:
            return P0(x) - ((x - x2) / (x1 - x2)) * P0(x1)
        return 0.0
    xs = np.linspace(1e-4, x1, 4000)
    dx = xs[1] - xs[0]
    pr = np.array([Pr(x) for x in xs])
    P = pr / np.trapezoid(pr, xs)
    Hmean = H0e * np.trapezoid(xs * P, xs)
    Hrms = H0e * math.sqrt(np.trapezoid(xs * xs * P, xs))
    exc = 1.0 - np.cumsum(P * dx) / np.sum(P * dx)
    def H1n(nf):
        idx = int(np.searchsorted(-exc, -1.0 / nf))
        idx = min(idx, len(xs) - 1)
        return H0e * nf * np.trapezoid((xs * P)[idx:], xs[idx:])
    Hs, H10, H2 = H1n(3.0), H1n(10.0), H1n(50.0)
    Hmax = min(x1 * H0e, Hrms * math.sqrt(math.log(n_waves)))   # breaking limit vs most-prob max
    return Hmean, Hrms, Hs, H10, H2, Hmax


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


def compute(inp: dict, *, g: float = G_SI, n_waves: float = 1000.0) -> Result:
    """Irregular-wave transformation for SI inputs."""
    _validate(inp)
    H0 = float(inp["H0"]); d = float(inp["d"]); Ts = float(inp["Ts"])
    cot_phi = float(inp["cot_phi"]); theta = float(inp["theta"]); s_max = float(inp["s_max"])

    L0 = g * Ts * Ts / (2.0 * math.pi)
    tanb = 1.0 / cot_phi
    Ks, L, k = _shoaling(Ts, d, g)
    Kr = _kr_eff(H0, Ts, d, theta, s_max, g)
    H0e = H0 * Kr                                      # effective deepwater height (with refraction)

    Hmean, Hrms, Hs, H10, H2, Hmax = _stats_at(H0e, Ks, d, L0, tanb, n_waves)

    surf_beat = 0.01 * H0 / math.sqrt((H0 / L0) * (1.0 + d / H0))           # eq 12
    # wave setup (set-down) at depth: radiation-stress difference from deep water (eq 13 integrated)
    term_d = (1.0 / 8.0) * Hrms ** 2 * (0.5 + 2.0 * k * d / math.sinh(2.0 * k * d))
    term_0 = (1.0 / 8.0) * (H0 / 1.416) ** 2 * 0.5
    setup = -(term_d - term_0) / d

    notes = (f"Ks={Ks:.4f} (exact); Kr={Kr:.4f} (oracle ~0.964, s_max-scheme dependent); "
             f"Hs/Hmean/Hrms ~3% of oracle; high quantiles approximate (finite-N order stats)")
    return Result(Hs=Hs, Hmean=Hmean, Hrms=Hrms, H10=H10, H2=H2, Hmax=Hmax, Ks=Ks, Kr=Kr,
                  surf_beat=surf_beat, setup=setup, steepness=H0 / L0, notes=notes)


# --- self-tests (ACES worked example; subject column to ~3-4%) ------------------
def _self_tests() -> None:
    g = G_SI
    r = compute({"H0": 20.0 * _FT, "d": 50.0 * _FT, "Ts": 8.0, "cot_phi": 100.0,
                 "theta": 10.0, "s_max": 10.0}, g=g)
    ft = lambda x: x / _FT
    rel = lambda got, exp, t: abs(got - exp) <= t * exp
    # shoaling and steepness are exact
    assert rel(r.Ks, 0.9133, 0.002), r.Ks
    assert rel(r.steepness, 0.0611, 0.005), r.steepness
    # headline Hs/Hrms within ~4%; other quantiles ~6% (Goda finite-N statistics + the
    # directional-spreading scheme are not fully in the public TR -- see module docstring)
    assert rel(ft(r.Hs), 17.7, 0.04), ft(r.Hs)
    assert rel(ft(r.Hrms), 12.5, 0.04), ft(r.Hrms)
    assert rel(ft(r.Hmean), 11.2, 0.05), ft(r.Hmean)
    assert rel(ft(r.H10), 22.5, 0.06), ft(r.H10)
    assert rel(ft(r.Hmax), 30.1, 0.06), ft(r.Hmax)
    # surf beat matches well; Kr within ~1% of oracle (scheme-dependent)
    assert rel(ft(r.surf_beat), 0.4350, 0.02), ft(r.surf_beat)
    assert 0.94 < r.Kr < 0.97, r.Kr
    print(f"  self-tests: PASS (Ks={r.Ks:.4f} exact; Hs={ft(r.Hs):.1f}/Hrms={ft(r.Hrms):.1f}/"
          f"H10={ft(r.H10):.1f} ft ~3-4%; surf beat {ft(r.surf_beat):.3f}; Kr={r.Kr:.3f})")


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
