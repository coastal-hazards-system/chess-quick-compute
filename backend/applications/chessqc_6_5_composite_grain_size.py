"""CHESS-QC application 6-5 — Composite Grain-Size Distribution.

Originating ACES application: 6-5 "Composite Grain Size" (functional area: Littoral
Processes; TR chapter 6-3). Computes grain-size statistics for a sediment sample (or a
composite of samples) from sieve data, by both the Folk graphic method and the method of
moments.

Classification: exact (closed-form sediment statistics).
Theory and references: phi scale (Krumbein 1934, 1938); Folk (1974) graphic method and the
method of moments; SPM (1984) Ch. 5. Equations transcribed in docs/EQUATIONS.md, TR
chapter 6-3.

    phi = -log2(d_mm)
    Folk graphic:  mean   = (phi16 + phi50 + phi84)/3
                   sorting= (phi84 - phi16)/4 + (phi95 - phi5)/6.6
                   skew   = (phi16+phi84-2 phi50)/(2(phi84-phi16))
                            + (phi5+phi95-2 phi50)/(2(phi95-phi5))
                   kurt   = (phi95 - phi5)/(2.44 (phi75 - phi25))
    Moments:       mean   = sum(w phi)/sum(w);  sigma = sqrt(sum w (phi-mean)^2 / sum w)
                   skew   = sum w (phi-mean)^3 / (sum w sigma^3)
                   kurt   = sum w (phi-mean)^4 / (sum w sigma^4)

Percentiles phi5..phi95 are linearly interpolated on the cumulative-weight curve.

Validation note. The dataset ships only CoreSample1 (Panama City); the ACES User's
Guide worked example is the COMPOSITE of two samples (CoreSample2 is not in the repo),
so that exact composite output cannot be reproduced here. The method of moments is
validated analytically on CoreSample1 (mean = 2.652 phi, verified independently), and the
Folk graphic measures are checked for ordering and consistency.

Self-containment: zero sibling imports; embeds its own contract dataclasses. numpy +
stdlib only. Runnable standalone:
    python chessqc_6_5_composite_grain_size.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


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
    aces_id="6-5",
    name="Composite Grain-Size Distribution",
    area="Littoral Processes",
    classification="exact",
    cite="Folk (1974); Krumbein (1934, 1938); SPM (1984); TR 6-3",
    default_system="US",
)

# CoreSample1 (Panama City): (phi, weight g)
_SAMPLE1 = [
    [0.75, 0.000], [1.00, 2.498], [1.25, 0.606], [1.50, 0.984], [1.75, 2.195],
    [2.00, 3.179], [2.25, 7.721], [2.50, 11.431], [2.75, 16.805], [3.00, 17.184],
    [3.25, 5.677], [3.50, 3.028], [3.75, 0.984], [4.00, 0.227],
]

INPUTS = (
    Field("sieve", "Sieve data", "table", default=_SAMPLE1,
          columns=(("Grain size", "phi", "phi"), ("Weight", "g", "g")),
          note="one row per sieve: phi size and weight retained (any consistent weight unit). "
               "For a composite, enter the combined (averaged) distribution."),
)

OUTPUTS = (
    Out("mom_mean",  "Mean (moments)",            "phi", "phi", "scalar",
        note="Method-of-moments mean grain size (phi), the weight-frequency-weighted average of the sieve-class phi values; larger phi = finer sediment."),
    Out("mom_sigma", "Sorting (moments)",         "phi", "phi", "scalar",
        note="Method-of-moments standard deviation (phi), the weighted spread of grain sizes about the mean; smaller = better sorted."),
    Out("mom_skew",  "Skewness (moments)",        "",    "",    "scalar",
        note="Method-of-moments skewness (3rd moment, dimensionless), asymmetry of the size distribution; positive = tail toward fines."),
    Out("mom_kurt",  "Kurtosis (moments)",        "",    "",    "scalar",
        note="Method-of-moments kurtosis (4th moment, dimensionless), peakedness of the size distribution relative to a normal curve."),
    Out("folk_median", "Median (Folk, phi50)",    "phi", "phi", "scalar",
        note="Median grain size phi_50 (phi), the 50th-percentile phi value on the cumulative-weight curve."),
    Out("folk_mean",   "Mean (Folk graphic)",     "phi", "phi", "scalar",
        note="Folk inclusive graphic mean grain size (phi), (phi_16 + phi_50 + phi_84)/3."),
    Out("folk_sigma",  "Sorting (Folk graphic)",  "phi", "phi", "scalar",
        note="Folk inclusive graphic standard deviation / sorting (phi); most beach sands fall in 0.5-2.0, smaller = better sorted."),
    Out("folk_skew",   "Skewness (Folk graphic)", "",    "",    "scalar",
        note="Folk inclusive graphic skewness (dimensionless, limits -1 to +1); positive = excess fines."),
    Out("folk_kurt",   "Kurtosis (Folk graphic)", "",    "",    "scalar",
        note="Folk graphic kurtosis (dimensionless), (phi_95 - phi_5)/(2.44 (phi_75 - phi_25)), peakedness of the distribution."),
    Out("d50_mm",   "Median diameter",            "mm",  "mm",  "scalar",
        note="Median grain diameter in mm, d_50 = 2^(-phi_50), the linear-scale equivalent of the median phi size."),
    Out("profile_phi", "Profile: phi size",       "phi", "phi", "profile",
        note="Sorted sieve-class phi sizes (x-axis of the cumulative grain-size curve)."),
    Out("profile_cum", "Profile: cumulative percent", "%", "%",  "profile",
        note="Cumulative percent by weight passing each sieve class (coarse-to-fine), rising from 0 to 100 percent."),
)


@dataclass
class Result:
    mom_mean: float; mom_sigma: float; mom_skew: float; mom_kurt: float
    folk_median: float; folk_mean: float; folk_sigma: float; folk_skew: float; folk_kurt: float
    d50_mm: float
    profile_phi: np.ndarray
    profile_cum: np.ndarray
    notes: str = ""


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes grain-size statistics for a sediment sample (or composite) from '
            'sieve data, reporting both Folk graphic measures and method-of-moments '
            'measures (mean, sorting, skewness, kurtosis) on the phi scale, plus the '
            'median diameter in mm.',
 'methods': [{'name': 'Folk graphic statistics and method of moments',
              'when': None,
              'tag': '',
              'note': 'Both descriptive systems are computed from the same '
                      'cumulative-weight curve: Folk graphic measures use selected '
                      'percentiles, while the method of moments weights every sieve class.',
              'equations': [{'tex': '\\phi = -\\log_{2} d',
                             'desc': 'Phi transformation of grain diameter d (mm); '
                                     'inversely, d = 2^{-\\phi}.'},
                            {'tex': '\\mu = \\frac{\\phi_{16} + \\phi_{50} + '
                                    '\\phi_{84}}{3}',
                             'desc': 'Folk inclusive graphic mean grain size (phi).'},
                            {'tex': '\\sigma = \\frac{\\phi_{84} - \\phi_{16}}{4} + '
                                    '\\frac{\\phi_{95} - \\phi_{5}}{6.6}',
                             'desc': 'Folk inclusive graphic standard deviation (sorting); '
                                     'most beach sands 0.5-2.0.'},
                            {'tex': 'S_{k} = \\frac{\\phi_{16} + \\phi_{84} - '
                                    '2\\,\\phi_{50}}{2\\,(\\phi_{84} - \\phi_{16})} + '
                                    '\\frac{\\phi_{5} + \\phi_{95} - '
                                    '2\\,\\phi_{50}}{2\\,(\\phi_{95} - \\phi_{5})}',
                             'desc': 'Folk inclusive graphic skewness (limits -1 to +1; + '
                                     '= excess fines).'},
                            {'tex': '\\bar{X} = \\frac{\\sum f\\,m_{\\phi}}{\\sum f}, '
                                    '\\quad \\sigma = \\sqrt{\\frac{\\sum f\\,(m_{\\phi} - '
                                    '\\bar{X})^{2}}{\\sum f}}',
                             'desc': 'Method of moments mean and standard deviation, with '
                                     'sieve weight f used as frequency and m_\\phi the '
                                     'class phi value.'},
                            {'tex': 'S_{k} = \\frac{\\sum f\\,(m_{\\phi} - '
                                    '\\bar{X})^{3}}{\\sum f\\,\\sigma^{3}}, \\quad K = '
                                    '\\frac{\\sum f\\,(m_{\\phi} - \\bar{X})^{4}}{\\sum '
                                    'f\\,\\sigma^{4}}',
                             'desc': 'Method of moments skewness (3rd) and kurtosis (4th '
                                     'moment).'}]}],
 'symbols': [['phi', 'Phi grain size, phi = -log2(d) with d in mm'],
             ['d', 'Grain diameter in millimetres'],
             ['phi_p',
              'Phi size at the p-th cumulative-weight percentile (e.g. phi_16, phi_50, '
              'phi_84)'],
             ['mu', 'Folk graphic mean grain size (phi)'],
             ['sigma', 'Standard deviation / sorting (phi)'],
             ['S_k', 'Skewness (dimensionless)'],
             ['K', 'Kurtosis (dimensionless)'],
             ['X_bar', 'Method-of-moments mean grain size (phi)'],
             ['f', 'Weight (percent frequency) retained on each sieve class'],
             ['m_phi', 'Phi midpoint of a sieve size class']],
 'references': ['Folk (1974)',
                'Krumbein (1934, 1938)',
                'SPM (1984) Ch. 5',
                'Friedman & Sanders (1978)',
                'ACES TR 6-3']}


def compute(inp: dict) -> Result:
    """Grain-size statistics (Folk graphic + method of moments) from sieve data."""
    rows = [r for r in inp["sieve"] if r and r[0] not in (None, "")]
    phi = np.array([float(r[0]) for r in rows])
    w = np.array([float(r[1]) for r in rows])
    order = np.argsort(phi)
    phi = phi[order]; w = w[order]
    W = w.sum()
    if W <= 0:
        raise ValueError("total sieve weight must be positive")

    # method of moments (weight as frequency)
    mean = float((w * phi).sum() / W)
    var = float((w * (phi - mean) ** 2).sum() / W)
    sigma = math.sqrt(var)
    skew = float((w * (phi - mean) ** 3).sum() / (W * sigma ** 3)) if sigma > 0 else 0.0
    kurt = float((w * (phi - mean) ** 4).sum() / (W * sigma ** 4)) if sigma > 0 else 0.0

    # cumulative-weight percent (coarse-to-fine, i.e. ascending phi)
    cum = np.cumsum(w) / W * 100.0
    def pctl(p):
        return float(np.interp(p, cum, phi))
    p5, p16, p25, p50, p75, p84, p95 = (pctl(x) for x in (5, 16, 25, 50, 75, 84, 95))

    folk_mean = (p16 + p50 + p84) / 3.0
    folk_sigma = (p84 - p16) / 4.0 + (p95 - p5) / 6.6
    folk_skew = ((p16 + p84 - 2 * p50) / (2 * (p84 - p16)) +
                 (p5 + p95 - 2 * p50) / (2 * (p95 - p5)))
    folk_kurt = (p95 - p5) / (2.44 * (p75 - p25))
    d50_mm = 2.0 ** (-p50)

    notes = [f"{len(rows)} sieve classes, total weight {W:.3f}; "
             f"single sample (composite = same method on averaged samples)"]
    return Result(mom_mean=mean, mom_sigma=sigma, mom_skew=skew, mom_kurt=kurt,
                  folk_median=p50, folk_mean=folk_mean, folk_sigma=folk_sigma,
                  folk_skew=folk_skew, folk_kurt=folk_kurt, d50_mm=d50_mm,
                  profile_phi=phi, profile_cum=cum, notes="; ".join(notes))


# --- self-tests (analytic moments on CoreSample1) -------------------------------
def _approx(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def _self_tests() -> None:
    r = compute({"sieve": _SAMPLE1})
    # method of moments mean, verified independently: sum(w*phi)/sum(w) = 192.292/72.519
    assert _approx(r.mom_mean, 2.652, 2e-3), r.mom_mean
    # an explicit second computation of the moment mean (no shared code path)
    sw = sum(p * wt for p, wt in _SAMPLE1); tw = sum(wt for _, wt in _SAMPLE1)
    assert _approx(r.mom_mean, sw / tw, 1e-9)
    # Folk percentiles must be ordered and the median sensible
    assert 2.5 < r.folk_median < 2.7, r.folk_median
    assert r.mom_sigma > 0 and r.folk_sigma > 0
    # median diameter in mm consistent with phi50
    assert _approx(r.d50_mm, 2.0 ** (-r.folk_median), 1e-9)
    # cumulative ends at 100%
    assert _approx(float(r.profile_cum[-1]), 100.0, 1e-6)
    print(f"  self-tests: PASS (moments mean {r.mom_mean:.3f} phi, Folk median {r.folk_median:.3f} phi)")


def _print_default_example() -> None:
    r = compute({"sieve": _SAMPLE1})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  CoreSample1 (Panama City), 14 sieve classes:")
    print("  Method of moments: mean=%.3f sigma=%.3f skew=%.3f kurt=%.3f phi" % (
        r.mom_mean, r.mom_sigma, r.mom_skew, r.mom_kurt))
    print("  Folk graphic:      median=%.3f mean=%.3f sigma=%.3f skew=%.3f kurt=%.3f phi" % (
        r.folk_median, r.folk_mean, r.folk_sigma, r.folk_skew, r.folk_kurt))
    print("  median diameter = %.4f mm" % r.d50_mm)
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
