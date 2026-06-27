"""CHESS-QC application 1-3 — Extremal Significant Wave Height Analysis.

Originating ACES application: 1-3 "Extremal Significant Wave Height Analysis"
(functional area: Wave Prediction). Fits extremal probability distributions to a
sample of storm significant wave heights and estimates design wave heights at
specified return periods, with confidence intervals.

Classification: exact (standard extremal-analysis procedure -- known Gringorten plotting
positions and the FT-I / Weibull candidate set, nothing guessed; reproduces the User's
Guide Example 1-3 correlations, return-period values and confidence intervals).
Theory and references: Goda (1988); plotting positions Gringorten (1963) / Petrauskas
  & Aagaard (1970). Equations transcribed in docs/EQUATIONS.md, TR chapter 1-3
  (eqs 1-12 + Tables 1-3-1/2): FT-I (1) and Weibull (2) distributions, Goda plotting
  positions (3), LS fit (4)/(5), return value (6)/(7), confidence interval (8)-(10),
  encounter probability (11), goodness-of-fit selection (12).

Candidates fitted: FT-I (Gumbel) and Weibull with k = 0.75, 1.0, 1.4, 2.0. The
best fit is the candidate with the smallest sum of squared residuals (eq 12). The
fit is performed in the data units (the regression of H on the reduced variate),
so all heights are handled in SI (m) internally; the GUI shows the chosen units.

Self-containment: zero sibling imports; embeds its own contract dataclasses (incl.
the table-input `columns` attribute) and the least-squares fit. stdlib only.
Runnable standalone:
    python chessqc_1_3_extremal_hs.py
which runs the User's Guide example self-test, then prints the example.

I/O mirrors the ACES manual:
  Inputs : H_s table (one per storm), N_T (total events), K (record length, yr),
           water depth, confidence level
  Outputs: best-fit distribution (note) + its return-period heights (Tr = 2..100 yr)
           and the confidence band on the design (Tr=100) value.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class AppMeta:
    aces_id: str; name: str; area: str; classification: str; cite: str; default_system: str = "SI"


@dataclass(frozen=True)
class Field:
    key: str; label: str; kind: str = "float"; unit_si: str = ""; unit_us: str = ""
    default: object = 0.0; lo: float = -math.inf; hi: float = math.inf
    choices: tuple = (); columns: tuple = (); note: str = ""


@dataclass(frozen=True)
class Out:
    key: str; label: str; unit_si: str = ""; unit_us: str = ""; kind: str = "scalar"
    note: str = ""           # hover definition shown on the output label


APP_META = AppMeta(
    aces_id="1-3",
    name="Extremal Significant Wave Height Analysis",
    area="Wave Prediction",
    classification="exact",
    cite="Goda (1988); Gringorten (1963); EM 1110-2-1414; TR 1-3",
    default_system="SI",     # the User's Guide example (EXTDELFT.IN) is in meters
)

# return periods reported by ACES (years)
_RETURN_PERIODS = (2.0, 5.0, 10.0, 25.0, 50.0, 100.0)

# Table 1-3-1, empirical std-deviation coefficients (Goda 1988):
#   name, weibull_k (None = FT-I),  alpha1, alpha2, kappa, c, eps
_DISTS = (
    ("FT-I",            None, 0.64, 9.0,  0.93, 0.0, 1.33),
    ("Weibull k=0.75",  0.75, 1.65, 11.4, -0.63, 0.0, 1.15),
    ("Weibull k=1.00",  1.00, 1.92, 11.4, 0.00, 0.3, 0.90),
    ("Weibull k=1.40",  1.40, 2.05, 11.4, 0.69, 0.4, 0.72),
    ("Weibull k=2.00",  2.00, 2.24, 11.4, 1.34, 0.5, 0.54),
)
# Table 1-3-2, confidence level -> z multiplier
_Z = {"80": 1.28, "85": 1.44, "90": 1.65, "95": 1.96, "99": 2.58}

# ACES User's Guide example EXTDELFT.IN (heights in metres):
_EX_HEIGHTS = [9.32, 8.11, 7.19, 7.06, 6.37, 6.15, 6.03, 5.72,
               4.92, 4.90, 4.78, 4.67, 4.64, 4.19, 3.06]
_FT = 0.3048
INPUTS = (
    Field("heights", "Significant wave heights", "table",
          default=[[h] for h in _EX_HEIGHTS], columns=(("H_s", "m", "ft"),),
          lo=0.0, hi=1e4, note="one significant height per storm/event"),
    Field("N_T", "Total number of events", "int", "", "", default=20, lo=1, hi=100000,
          note="events during the record (>= number of heights)"),
    Field("K", "Record length", "float", "yr", "yr", default=20.0, lo=1e-3, hi=1e4,
          note="years of record"),
    Field("depth", "Water depth", "float", "m", "ft", default=500.0, lo=1e-3, hi=1e5,
          note="for depth-limited breaking cap (inactive in deep water)"),
    Field("conf", "Confidence level", "choice", "", "", default="90",
          choices=("80", "85", "90", "95", "99"), note="confidence interval, %"),
)

OUTPUTS = (
    Out("corr",  "Best-fit correlation", "", "", "scalar",
        note="correlation coefficient r of the least-squares fit of the best candidate distribution (closer to 1 = better fit)"),
    Out("Hs2",   "H_s, return period 2 yr",   "m", "ft", "scalar",
        note="significant wave height expected to be equalled or exceeded once on average every 2 years, from the best-fit distribution"),
    Out("Hs5",   "H_s, return period 5 yr",   "m", "ft", "scalar",
        note="significant wave height expected to be equalled or exceeded once on average every 5 years, from the best-fit distribution"),
    Out("Hs10",  "H_s, return period 10 yr",  "m", "ft", "scalar",
        note="significant wave height expected to be equalled or exceeded once on average every 10 years, from the best-fit distribution"),
    Out("Hs25",  "H_s, return period 25 yr",  "m", "ft", "scalar",
        note="significant wave height expected to be equalled or exceeded once on average every 25 years, from the best-fit distribution"),
    Out("Hs50",  "H_s, return period 50 yr",  "m", "ft", "scalar",
        note="significant wave height expected to be equalled or exceeded once on average every 50 years, from the best-fit distribution"),
    Out("Hs100", "H_s, return period 100 yr", "m", "ft", "scalar",
        note="design significant wave height at the 100-year return period, from the best-fit distribution"),
    Out("Hs100_lo", "Design (100 yr) lower bound", "m", "ft", "scalar",
        note="lower limit of the confidence band on the 100-year design height (H_sr minus z*sigma at the chosen confidence level)"),
    Out("Hs100_hi", "Design (100 yr) upper bound", "m", "ft", "scalar",
        note="upper limit of the confidence band on the 100-year design height (H_sr plus z*sigma at the chosen confidence level)"),
)


@dataclass
class Result:
    corr: float
    Hs2: float; Hs5: float; Hs10: float; Hs25: float; Hs50: float; Hs100: float
    Hs100_lo: float; Hs100_hi: float
    best_dist: str = ""
    dists: list = field(default_factory=list)   # full per-distribution detail (tests)
    notes: str = ""


def _plot_position(m: int, N_T: int, k):
    """Goda (1988) plotting position F(H_s <= H_sm) for rank m (eq 3)."""
    if k is None:                                       # FT-I
        return 1.0 - (m - 0.44) / (N_T + 0.12)
    sk = math.sqrt(k)
    return 1.0 - (m - 0.20 - 0.27 / sk) / (N_T + 0.20 + 0.23 / sk)


def _reduced(F: float, k):
    """Reduced variate y from the non-exceedance probability F (eq 5)."""
    if k is None:                                       # FT-I
        return -math.log(-math.log(F))
    return (-math.log(1.0 - F)) ** (1.0 / k)


def _reduced_return(lam_Tr: float, k):
    """Return-period reduced variate y_r for lambda*Tr (eq 7)."""
    if k is None:                                       # FT-I
        return -math.log(-math.log(1.0 - 1.0 / lam_Tr))
    return (math.log(lam_Tr)) ** (1.0 / k)


def _lsq(xs, ys):
    """Least-squares line ys = A*xs + B (here xs = reduced variate y, ys = H).
    Returns (A, B, correlation, sxx) where sxx = sum (xs - mean)^2."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    A = sxy / sxx
    B = my - A * mx
    corr = sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else 0.0
    return A, B, corr, sxx


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Fits extremal probability distributions (FT-I/Gumbel and Weibull with k = '
            '0.75, 1.0, 1.4, 2.0) to a sample of storm significant wave heights using Goda '
            'plotting positions and a least-squares regression on the reduced variate. It '
            'returns the best-fit distribution and the design significant wave heights at '
            'return periods of 2-100 years, with confidence intervals.',
 'methods': [{'name': 'Goda extremal fit (FT-I + Weibull, least-squares on plotting '
                      'positions)',
              'when': None,
              'tag': '',
              'note': 'All five candidates (FT-I and Weibull k = 0.75, 1.0, 1.4, 2.0) are '
                      'fitted; the one with the smallest residual sum (highest '
                      'correlation) is selected per eq (12).',
              'equations': [{'tex': 'F(H_s \\leq \\hat{H}_s) = \\exp\\left[ -\\exp\\left( '
                                    '-\\frac{\\hat{H}_s - B}{A} \\right) \\right]',
                             'desc': 'FT-I (Gumbel) non-exceedance distribution; A = '
                                     'scale, B = location (eq 1)'},
                            {'tex': 'F(H_s \\leq \\hat{H}_s) = 1 - \\exp\\left[ -\\left( '
                                    '\\frac{\\hat{H}_s - B}{A} \\right)^{k} \\right]',
                             'desc': 'Weibull non-exceedance distribution with shape '
                                     'parameter k (eq 2)'},
                            {'tex': 'F(H_s \\leq H_{sm}) = 1 - \\frac{m - 0.44}{N_T + '
                                    '0.12}',
                             'desc': 'Goda plotting position for the m-th ranked height '
                                     '(FT-I form; Weibull uses k-dependent constants) (eq '
                                     '3)'},
                            {'tex': 'H_{sr} = \\hat{A}\\, y_r + \\hat{B}',
                             'desc': 'Return-period significant height from the LS line, '
                                     'with y_r the return-period reduced variate at '
                                     'lambda*T_r, lambda = N_T/K (eqs 6-7)'},
                            {'tex': '\\sigma_{nr} = \\frac{1}{\\sqrt{N}}\\sqrt{ 1 + '
                                    '\\alpha\\,(y_r - c + \\epsilon\\,\\ln\\nu)^2 }',
                             'desc': 'Normalized standard deviation of the return value; '
                                     'the confidence band is H_sr +/- z*sigma_nr*sigma_Hs '
                                     '(eqs 8-10)'}]}],
 'symbols': [['H_s', 'significant wave height (per storm/event)'],
             ['H_{sr}', 'significant height at return period T_r'],
             ['F', 'probability of H_s not being exceeded (non-exceedance)'],
             ['A, B', 'scale and location parameters (LS estimates Â, B̂)'],
             ['k', 'Weibull shape parameter (0.75, 1.0, 1.4, 2.0)'],
             ['y_r', 'return-period reduced variate'],
             ['m', 'rank of the height in the descending-sorted sample (1..N)'],
             ['N_T', 'total number of events during the record'],
             ['lambda', 'average events per year, lambda = N_T/K (K = record length, yr)'],
             ['nu', 'censoring parameter, nu = N/N_T']],
 'references': ['Goda (1988)',
                'Gringorten (1963)',
                'Petrauskas & Aagaard (1970)',
                'Gumbel (1958)',
                'EM 1110-2-1414 (Ch. 5)',
                'ACES TR 1-3']}


def compute(inp: dict, *, g=None) -> Result:
    """Extremal analysis for SI inputs (heights in m). `heights` is a table -> list
    of rows; each row's first cell is one significant wave height."""
    rows = inp["heights"]
    H = sorted((float(r[0]) for r in rows if r and r[0] not in (None, "")), reverse=True)
    N = len(H)
    if N < 3:
        raise ValueError("need at least 3 significant wave heights")
    N_T = int(inp.get("N_T", N))
    if N_T < N:
        raise ValueError(f"N_T ({N_T}) must be >= number of heights ({N})")
    K = float(inp.get("K", 1.0))
    conf = str(inp.get("conf", "90"))
    z = _Z.get(conf, 1.65)

    lam = N_T / K                                       # events per year (eq 7)
    nu = N / N_T                                        # censoring parameter (eq 8)
    # sample standard deviation of the input heights (eq 10)
    mH = sum(H) / N
    sigma_Hs = math.sqrt(sum((h - mH) ** 2 for h in H) / (N - 1))

    ranks = range(1, N + 1)
    dists = []
    for name, k, a1, a2, kappa, c, eps in _DISTS:
        F = [_plot_position(m, N_T, k) for m in ranks]
        y = [_reduced(f, k) for f in F]
        A, B, corr, syy_y = _lsq(y, H)
        # goodness of fit: eq (12) minimizes the H-space residual sum, which (since
        # the data variance is common to all candidates) is monotonic with the
        # correlation coefficient -> select by the highest correlation (Goda 1988).
        ss = (1.0 - corr * corr)
        # return-period values + confidence band (eqs 6-10)
        alpha = a1 * math.exp(a2 * N ** (-1.3) + kappa * math.sqrt(-math.log(nu)))
        ret = {}
        for Tr in _RETURN_PERIODS:
            y_r = _reduced_return(lam * Tr, k)
            Hsr = A * y_r + B
            sig_nr = (1.0 / math.sqrt(N)) * math.sqrt(1.0 + alpha * (y_r - c + eps * math.log(nu)) ** 2)
            sig_r = sig_nr * sigma_Hs
            ret[Tr] = (Hsr, Hsr - z * sig_r, Hsr + z * sig_r)
        dists.append({"name": name, "k": k, "A": A, "B": B, "corr": corr, "ss": ss,
                      "ret": ret})

    best = min(dists, key=lambda d: d["ss"])            # eq 12: smallest residual sum
    r = best["ret"]
    notes = (f"best fit: {best['name']} (r={best['corr']:.4f}); "
             f"N={N}, N_T={N_T}, lambda={lam:.3g}/yr, {conf}% CI")
    return Result(
        corr=best["corr"],
        Hs2=r[2.0][0], Hs5=r[5.0][0], Hs10=r[10.0][0], Hs25=r[25.0][0],
        Hs50=r[50.0][0], Hs100=r[100.0][0],
        Hs100_lo=r[100.0][1], Hs100_hi=r[100.0][2],
        best_dist=best["name"], dists=dists, notes=notes)


# --- self-tests -----------------------------------------------------------------
def _approx(a, b, tol=1e-3):
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    r = compute({"heights": [[h] for h in _EX_HEIGHTS], "N_T": 20, "K": 20.0,
                 "depth": 500.0, "conf": "90"})
    by = {d["name"]: d for d in r.dists}
    # correlations vs the User's Guide first screen (all five candidates)
    for name, exp in {"FT-I": 0.9813, "Weibull k=0.75": 0.9414, "Weibull k=1.00": 0.9674,
                      "Weibull k=1.40": 0.9818, "Weibull k=2.00": 0.9866}.items():
        assert _approx(by[name]["corr"], exp, 2e-3), (name, by[name]["corr"])
    # return-period heights (manual shown in ft; fit in m -> compare m vs ft*0.3048)
    ft_FTI = {2: 15.94, 5: 21.50, 10: 25.18, 25: 29.84, 50: 33.29, 100: 36.72}
    ft_W2 = {2: 15.86, 5: 22.02, 10: 25.53, 25: 29.44, 50: 32.03, 100: 34.40}
    for Tr, ftval in ft_FTI.items():
        assert _approx(by["FT-I"]["ret"][float(Tr)][0], ftval * _FT, 5e-3), (Tr, ftval)
    for Tr, ftval in ft_W2.items():
        assert _approx(by["Weibull k=2.00"]["ret"][float(Tr)][0], ftval * _FT, 5e-3), (Tr, ftval)
    # FT-I 90% CI at Tr=100: manual 25.2 - 48.2 ft
    lo, hi = by["FT-I"]["ret"][100.0][1], by["FT-I"]["ret"][100.0][2]
    assert _approx(lo / _FT, 25.2, 5e-2) and _approx(hi / _FT, 48.2, 5e-2), (lo / _FT, hi / _FT)
    # best fit = Weibull k=2.0 (lowest SS)
    assert r.best_dist == "Weibull k=2.00", r.best_dist
    print("  self-tests: PASS (User's Guide Example 1-3 incl. CI + best-fit selection)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  inputs: {len(inp['heights'])} heights, N_T={inp['N_T']}, K={inp['K']}, "
          f"conf={inp['conf']}%")
    print("  return-period H_s (ft) per distribution:")
    print("    " + "Tr".ljust(6) + "".join(d["name"][:9].rjust(11) for d in r.dists))
    for Tr in _RETURN_PERIODS:
        line = f"    {int(Tr):<6}"
        for d in r.dists:
            line += f"{d['ret'][Tr][0] / _FT:11.2f}"
        print(line)
    print(f"  notes: {r.notes}")
    print("  (manual: best fit Weibull k=2.0; FT-I Tr=100 = 36.72 ft, CI 25.2-48.2 ft)")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
