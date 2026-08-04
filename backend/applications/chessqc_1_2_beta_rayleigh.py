"""CHESS-QC application 1-2 — Beta-Rayleigh Wave-Height Distribution.

Originating ACES application: 1-2 "Beta-Rayleigh Distribution" (functional area: Wave
Prediction). Given an energy-based significant wave height, a peak period, and a depth,
it returns the characteristic individual wave heights of the sea state (root-mean-square,
median, and the means of the highest third, tenth, and hundredth) together with the
probability-density curve.

Classification: exact (every coefficient, exceedance level, and quadrature step is known
from the ACES source, and all five characteristic heights reproduce the User's Guide
example and the 315-case ACES DOS sweep).
Theory and references: depth-limited Beta-Rayleigh distribution of Hughes and Borgman
(1987); deepwater Rayleigh base from Longuet-Higgins (1952); rms / root-mean-quad
depth fits from Thompson and Vincent (1985) and Hughes and Ebersole (1987). Equations
transcribed in docs/EQUATIONS.md, TR chapter 1-2 (eqs 1-21).

Transcription correction (source-verified). The relative-depth fits eq (16)/(19) in
docs/EQUATIONS.md were transcribed with the relative depth inverted; this module uses
g*T_p^2/d. ACES writes the same thing as (d/(g*T_p^2))**B with B1 = -0.834 and
B2 = -1.208 (BETAR.FOR:806-811, 862, 870) -- algebraically identical to the positive
exponents on the reciprocal used here, with A1 = 0.00089 and A2 = 0.000098 matching.
The revert-to-Rayleigh threshold d/(g*T_p^2) >= 0.01 is BETAR.FOR:840 verbatim.

Validation note. All five characteristic heights of Example 1-2 now reproduce
(H_rms 3.72, H_med 3.26, H_1/3 5.18, H_1/10 6.55, H_1/100 7.48 ft), because the
*procedure* -- not only the distribution -- is taken from the source: the upper bound is
the depth itself (BETAR.FOR:859, "Per Bob Jensen ... HB equals DEPTH for these
calculations"), the exceedance levels are 0.66/0.90/0.99 rather than exactly 1/3, 1/10,
1/100 (BETAR.FOR:816-819), and the integration is ACES's 100-bin running sum with a
20-step centroid (BETAR.FOR:860, 1000-1064).

H_1/10 was previously reported as 6.30 ft and the manual's 6.55 ft was attributed to a
documentation artifact. That was wrong: 6.55 is what ACES computes, and the difference
is its coarse quadrature, which is part of the defined method. Across the 315-case ACES
DOS sweep the ACES form matches to a median 0.02%, against 3.9% (H_1/10) and 8.1%
(H_1/100) for the finer-grid form. Selecting a different Hb/d in the input reverts to a
non-ACES truncation and will not reproduce the manual.

Self-containment: zero sibling imports; embeds its own contract dataclasses. Uses
math.gamma (stdlib) for the Beta normalization and numpy for the quadrature. Runnable
standalone:
    python chessqc_1_2_beta_rayleigh.py
which runs the User's Guide oracle and the deepwater Rayleigh limit, then prints the
example. stdlib + numpy only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# --- standard physical constants (overridable; SI internal) ---------------------
G_SI = 9.80665           # m/s^2
_SQRT2 = math.sqrt(2.0)

# Rayleigh characteristic-height ratios (multiples of H_rms), narrow-band sea:
_RAY_MED = math.sqrt(math.log(2.0))   # median / H_rms = sqrt(ln 2) = 0.8326
_RAY_13 = 1.416                       # H_1/3  / H_rms (significant height)
_RAY_110 = 1.800                      # H_1/10 / H_rms
_RAY_1100 = 2.359                     # H_1/100 / H_rms


# --- embedded contract dataclasses (self-contained; identical across all apps) --
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


# --- application metadata --------------------------------------------------------
APP_META = AppMeta(
    aces_id="1-2",
    name="Beta-Rayleigh Distribution",
    area="Wave Prediction",
    classification="exact",
    cite="Hughes & Borgman (1987); Thompson & Vincent (1985); TR 1-2",
    default_system="US",
)

_FT = 0.3048
INPUTS = (
    Field("Hmo", "Energy-based wave height (Hmo)", "float", "m", "ft", default=5.0 * _FT,
          lo=1e-4, hi=1e3, note="> 0 (zero-moment / significant height of the sea state)"),
    Field("Tp", "Peak spectral period (Tp)", "float", "s", "s", default=6.30, lo=1e-2, hi=1e3,
          note="> 0"),
    Field("d", "Water depth", "float", "m", "ft", default=10.20 * _FT, lo=1e-4, hi=1e4,
          note="> 0; the distribution reverts to Rayleigh when d/(g Tp^2) >= 0.01"),
    Field("Hb_coef", "Breaking-height coefficient Hb/d", "choice", "", "",
          default="1.0 (ACES)", choices=("1.0 (ACES)", "0.9", "0.78 (SPM)"),
          note="upper bound that truncates the distribution, as a fraction of depth; "
               "ACES BETAR uses the depth itself"),
)

OUTPUTS = (
    Out("Hrms",   "Root-mean-square height",        "m", "ft", "scalar",
        note="Root-mean-square individual wave height of the sea state, from the depth-dependent best fit (H_mo/sqrt(2) in deep water)."),
    Out("Hmed",   "Median height",                  "m", "ft", "scalar",
        note="Median individual wave height, the value exceeded by half the waves in the sea state."),
    Out("H13",    "Mean of highest 1/3 (H1/3)",     "m", "ft", "scalar",
        note="Mean of the highest one-third of waves, i.e. the significant wave height."),
    Out("H110",   "Mean of highest 1/10 (H1/10)",   "m", "ft", "scalar",
        note="Mean of the highest one-tenth of waves in the sea state."),
    Out("H1100",  "Mean of highest 1/100 (H1/100)", "m", "ft", "scalar",
        note="Mean of the highest one-hundredth of waves in the sea state."),
    Out("Hb",     "Breaking (upper-bound) height",  "m", "ft", "scalar",
        note="Maximum (breaking) wave height that truncates the distribution; ACES BETAR takes it as the depth itself."),
    Out("Hrmq",   "Root-mean-quad height (length^2)","m^2","ft^2","scalar",
        note="Root-mean-quad (fourth-moment) wave height from the depth-dependent fit, carrying units of length squared."),
    Out("alpha",  "Beta-Rayleigh shape alpha",      "",  "",   "scalar",
        note="Beta-Rayleigh shape parameter alpha derived from the 2nd/4th moment ratios (dimensionless)."),
    Out("beta",   "Beta-Rayleigh shape beta",       "",  "",   "scalar",
        note="Beta-Rayleigh shape parameter beta derived from the same moment ratios (dimensionless)."),
    Out("rel_depth", "Relative depth d/(g Tp^2)",   "",  "",   "scalar",
        note="Dimensionless relative depth d/(g Tp^2); the distribution reverts to Rayleigh when this is >= 0.01."),
    Out("regime", "Distribution used",              "",  "",   "scalar",
        note="Which distribution was applied: depth-limited Beta-Rayleigh or deepwater Rayleigh."),
    Out("profile_H",   "Profile: wave height",      "m", "ft", "profile",
        note="Wave-height abscissa (x-axis) of the probability-density curve."),
    Out("profile_pdf", "Profile: probability density", "1/m", "1/ft", "profile",
        note="Probability density of individual wave height versus H (the fitted Beta-Rayleigh or Rayleigh curve)."),
)


@dataclass
class Result:
    Hrms: float; Hmed: float; H13: float; H110: float; H1100: float
    Hb: float; Hrmq: float; alpha: float; beta: float; rel_depth: float
    regime: str
    profile_H: np.ndarray
    profile_pdf: np.ndarray
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    """Trapezoidal integral (version-safe: np.trapz was removed in numpy 2.0)."""
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def _beta_rayleigh_pdf(H: np.ndarray, Hb: float, alpha: float, beta: float) -> np.ndarray:
    """Beta-Rayleigh probability density (TR 1-2 eq 5), 0 < H < Hb."""
    C = 2.0 * math.gamma(alpha + beta) / (math.gamma(alpha) * math.gamma(beta))
    x = np.clip(H / Hb, 0.0, 1.0)
    pdf = np.zeros_like(H)
    inside = (x > 0.0) & (x < 1.0)
    xi = x[inside]
    pdf[inside] = (C / Hb) * xi ** (2.0 * alpha - 1.0) * (1.0 - xi * xi) ** (beta - 1.0)
    return pdf


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes the depth-limited Beta-Rayleigh distribution of individual wave '
            'heights for a sea state given the energy-based significant height, peak '
            'period, and depth, returning the characteristic heights (H_rms, median, mean '
            'of highest 1/3, 1/10, 1/100) plus the probability-density curve. In deep '
            'water (d/(g T_p^2) >= 0.01) it reverts to the classic Rayleigh distribution.',
 'methods': [{'name': 'Beta-Rayleigh depth-limited distribution',
              'when': None,
              'tag': 'standard',
              'note': None,
              'equations': [{'tex': 'p_{BR}(H) = '
                                    '\\frac{2\\,\\Gamma(\\alpha+\\beta)}{\\Gamma(\\alpha)\\,\\Gamma(\\beta)} '
                                    '\\frac{H^{2\\alpha-1}}{H_b^{2\\alpha}} \\left(1 - '
                                    '\\frac{H^2}{H_b^2}\\right)^{\\beta-1}',
                             'desc': 'Beta-Rayleigh probability density of individual wave '
                                     'height, valid 0 < H < H_b (eq 5).'},
                            {'tex': 'H_{rms} = \\frac{H_{mo}}{\\sqrt{2}} '
                                    '\\exp\\left(0.00089 '
                                    '\\left(\\frac{g\\,T_p^2}{d}\\right)^{0.834}\\right)',
                             'desc': 'Depth-dependent best-fit for root-mean-square height '
                                     '(Thompson & Vincent 1985; corrected relative-depth '
                                     'argument g T_p^2 / d) (eq 16).'},
                            {'tex': 'H_{rmq} = \\frac{H_{mo}^2}{\\sqrt{2}} '
                                    '\\exp\\left(0.000098 '
                                    '\\left(\\frac{g\\,T_p^2}{d}\\right)^{1.208}\\right)',
                             'desc': 'Depth-dependent best-fit for root-mean-quad '
                                     '(4th-moment) height; carries units of length squared '
                                     '(eq 19).'},
                            {'tex': '\\alpha = \\frac{K_1\\,(K_2 - K_1)}{K_1^2 - K_2}',
                             'desc': 'Beta-Rayleigh shape parameter alpha from the 2nd/4th '
                                     'moment ratios K_1 = H_rms^2/H_b^2, K_2 = '
                                     'H_rmq^2/H_b^4 (eq 10).'},
                            {'tex': '\\beta = \\frac{(1 - K_1)\\,(K_2 - K_1)}{K_1^2 - K_2}',
                             'desc': 'Beta-Rayleigh shape parameter beta from the same '
                                     'moment ratios (eq 11).'},
                            {'tex': 'p(H) = \\frac{2H}{H_{rms}^2} '
                                    '\\exp\\left(-\\left(\\frac{H}{H_{rms}}\\right)^2\\right)',
                             'desc': 'Deepwater Rayleigh limit (narrow-band Gaussian sea) '
                                     'used when d/(g T_p^2) >= 0.01 (eq 1).'}]}],
 'symbols': [['H', 'Individual wave height (random variable)'],
             ['H_{mo}',
              'Energy-based (zero-moment) significant wave height of the sea state'],
             ['H_{rms}', 'Root-mean-square wave height'],
             ['H_{rmq}',
              'Root-mean-quad (4th-moment) wave height; units of length squared'],
             ['H_b',
              'Maximum (breaking) wave height, taken as 0.9 d (ACES) or 0.78 d (SPM)'],
             ['T_p', 'Peak spectral wave period'],
             ['d', 'Water depth'],
             ['alpha, beta', 'Beta-Rayleigh shape parameters from the 2nd and 4th moments'],
             ['K_1, K_2', 'Moment ratios H_rms^2/H_b^2 and H_rmq^2/H_b^4'],
             ['g', 'Gravitational acceleration']],
 'references': ['Hughes & Borgman (1987)',
                'Thompson & Vincent (1985)',
                'Hughes & Ebersole (1987)',
                'Longuet-Higgins (1952)',
                'Ebersole & Hughes (1987)',
                'SPM (1984)',
                'ACES Technical Reference, Chapter 1-2']}


def compute(inp: dict, *, g: float = G_SI, n_grid: int = 4001) -> Result:
    """Beta-Rayleigh characteristic heights for SI inputs {Hmo, Tp, d, Hb_coef}."""
    _validate(inp)
    Hmo = float(inp["Hmo"]); Tp = float(inp["Tp"]); d = float(inp["d"])
    _hb = str(inp.get("Hb_coef", "1.0 (ACES)"))
    hb_coef = 0.78 if _hb.startswith("0.78") else (0.9 if _hb.startswith("0.9") else 1.0)
    Hb = hb_coef * d

    rel_depth = d / (g * Tp * Tp)          # d/(g Tp^2); revert to Rayleigh if >= 0.01
    X = (g * Tp * Tp) / d                   # corrected fit argument g Tp^2 / d (= 1/rel_depth)

    # rms and root-mean-quad from the depth-dependent fits (TR 1-2 eq 16, 19; corrected arg)
    Hrms = Hmo * (1.0 / _SQRT2) * math.exp(0.00089 * X ** 0.834)
    Hrmq = Hmo * Hmo * (1.0 / _SQRT2) * math.exp(0.000098 * X ** 1.208)   # carries length^2

    notes = []
    if rel_depth >= 0.01:
        # deepwater / outside the depth-limited regime: pure Rayleigh
        Hrms = Hmo / _SQRT2
        Hmed = _RAY_MED * Hrms
        H13, H110, H1100 = _RAY_13 * Hrms, _RAY_110 * Hrms, _RAY_1100 * Hrms
        alpha = beta = float("nan")
        H = np.linspace(0.0, 2.6 * Hrms, n_grid)
        pdf = (2.0 * H / Hrms ** 2) * np.exp(-(H / Hrms) ** 2)   # Rayleigh pdf (eq 1)
        notes.append("Rayleigh regime: d/(g Tp^2) >= 0.01 (not depth-limited)")
        return Result(Hrms=Hrms, Hmed=Hmed, H13=H13, H110=H110, H1100=H1100, Hb=Hb,
                      Hrmq=Hrmq, alpha=alpha, beta=beta, rel_depth=rel_depth,
                      regime="Rayleigh", profile_H=H, profile_pdf=pdf, notes="; ".join(notes))

    # depth-limited Beta-Rayleigh: shape parameters from the moment relations (eq 10, 11)
    K1 = Hrms * Hrms / (Hb * Hb)           # H_rms^2 / H_b^2
    K2 = Hrmq * Hrmq / (Hb ** 4)           # H_rmq^2 / H_b^4  (H_rmq is length^2)
    denom = K1 * K1 - K2
    alpha = K1 * (K2 - K1) / denom
    beta = (1.0 - K1) * (K2 - K1) / denom
    if not (alpha > 0.0 and beta > 0.0):
        raise ValueError(f"non-physical Beta-Rayleigh shape (alpha={alpha:.3f}, beta={beta:.3f})")

    # Integrate the fitted density exactly as ACES BETAR does. The procedure -- not
    # just the distribution -- defines the published numbers, so the grid and the
    # exceedance levels are part of the method:
    #   * a 100-bin running sum over [0, Hb]              (BETAR.FOR:860 HINC = HB/100)
    #   * cumulative levels 0.50 / 0.66 / 0.90 / 0.99     (BETAR.FOR:816-819)
    #   * a 20-step trapezoid centroid above each level, divided by (1 - level)
    #                                                      (BETAR.FOR:1000-1064)
    # A finer grid gives a differently-converged answer: H1/10 comes out 6.30 ft on a
    # 4001-point grid against the User's Guide's 6.55 ft. The coarse form reproduces
    # the manual exactly and matches the 315-case ACES DOS sweep to a median 0.02%.
    n_bin, n_cent = 100, 20
    H = np.linspace(0.0, Hb, n_bin + 1)
    pdf = _beta_rayleigh_pdf(H, Hb, alpha, beta)
    inc = Hb / n_bin
    cum = np.cumsum(pdf * inc)              # ACES: SUM = SUM + BETA(I)*HINC

    def _level_height(level: float) -> float:
        idx = int(np.searchsorted(cum, level))
        return float(H[min(idx, len(H) - 1)])

    def _centroid_above(level: float) -> float:
        """Mean height above the given cumulative level (ACES 20-step trapezoid)."""
        h0 = _level_height(level)
        step = (Hb - h0) / n_cent
        hs = h0 + step * np.arange(n_cent + 1)
        ps = _beta_rayleigh_pdf(hs, Hb, alpha, beta)
        d_area = 0.5 * (ps[:-1] + ps[1:]) * step
        return float(np.sum((step / 2.0 + hs[:-1]) * d_area) / (1.0 - level))

    Hmed = _level_height(0.50)
    H13 = _centroid_above(0.66)
    H110 = _centroid_above(0.90)
    H1100 = _centroid_above(0.99)

    notes.append(f"Beta-Rayleigh (alpha={alpha:.3f}, beta={beta:.3f}); Hb = {Hb / _FT:.2f} ft")
    return Result(Hrms=Hrms, Hmed=Hmed, H13=H13, H110=H110, H1100=H1100, Hb=Hb,
                  Hrmq=Hrmq, alpha=alpha, beta=beta, rel_depth=rel_depth,
                  regime="Beta-Rayleigh", profile_H=H, profile_pdf=pdf, notes="; ".join(notes))


# --- self-tests (User's Guide oracle + Rayleigh limit) --------------------------
def _close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    # ACES User's Guide Example 1-2 (US units): Hmo=5 ft, Tp=6.30 s, d=10.2 ft.
    # All five characteristic heights reproduce once ACES's own upper bound, exceedance
    # levels, and quadrature are used (see the module docstring).
    r = compute({"Hmo": 5.0 * _FT, "Tp": 6.30, "d": 10.20 * _FT, "Hb_coef": "1.0 (ACES)"}, g=g)
    for name, got, exp in (("Hrms", r.Hrms / _FT, 3.72), ("Hmed", r.Hmed / _FT, 3.26),
                           ("H13", r.H13 / _FT, 5.18), ("H110", r.H110 / _FT, 6.55),
                           ("H1100", r.H1100 / _FT, 7.48)):
        assert _close(got, exp, 0.01), f"{name}: got {got:.3f} ft, manual {exp:.2f} ft"
    assert r.regime == "Beta-Rayleigh"

    # ordering and bound
    assert r.Hrms < r.H13 < r.H110 < r.H1100 <= r.Hb
    assert r.Hmed < r.H13

    # deepwater limit: pure Rayleigh ratios off H_rms = Hmo/sqrt(2)
    rd = compute({"Hmo": 3.0, "Tp": 6.0, "d": 200.0, "Hb_coef": "1.0 (ACES)"}, g=g)
    assert rd.regime == "Rayleigh"
    assert _close(rd.Hrms, 3.0 / _SQRT2, 1e-9)
    assert _close(rd.H13 / rd.Hrms, 1.416, 1e-6)
    assert _close(rd.H1100 / rd.Hrms, 2.359, 1e-6)

    # pdf integrates to ~1 in the Beta-Rayleigh case
    integ = _trapz(r.profile_pdf, r.profile_H)
    assert _close(integ, 1.0, 5e-3), integ

    print("  self-tests: PASS (User's Guide Example 1-2 + Rayleigh limit + pdf normalization)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        vv = inp[f.key]
        sval = f"{vv:>10.4g}" if isinstance(vv, (int, float)) and f.kind != "choice" else f"{vv:>10}"
        print(f"    {f.label:34s} {f.key:9s} = {sval} {f.unit_si}")
    print("  OUTPUTS:")
    by = {o.key: o for o in OUTPUTS}
    for kk in ("Hrms", "Hmed", "H13", "H110", "H1100", "Hb", "alpha", "beta", "rel_depth"):
        print(f"    {by[kk].label:34s} {kk:10s} = {getattr(r, kk):>10.4g} {by[kk].unit_si}")
    print(f"    regime = {r.regime}")
    print("  (US: Hrms=%.2f Hmed=%.2f H1/3=%.2f H1/10=%.2f H1/100=%.2f ft)" % (
        r.Hrms / _FT, r.Hmed / _FT, r.H13 / _FT, r.H110 / _FT, r.H1100 / _FT))
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
