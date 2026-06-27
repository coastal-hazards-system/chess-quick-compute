"""CHESS-QC application 10-4 - Probabilistic Simulation Technique (PST).

Functional area: Coastal Hazards. Turns a Peaks-Over-Threshold sample (10-3)
into a hazard curve: response magnitude versus annual exceedance rate (AER),
with a bootstrap confidence band. The upper tail is a fitted Generalized Pareto
Distribution (GPD); the frequent range is carried empirically; the two are
spliced into one continuous curve.

Method (pure-numpy port of the PyStorm PST module):
  1. Population rate lambda_u = n_pot / record_length. record_length is given or
     auto = n_pot / events_per_year (the POT module trims to exactly
     events_per_year x effective_duration peaks, so this recovers the effective
     duration). Keep events_per_year equal to the POT target rate.
  2. Empirical AER (Weibull plotting position): sort peaks descending; rank i ->
     AER_i = i / (n_pot + 1) * lambda_u.
  3. GPD location mu (Quantile Delta Optimization): scan candidate thresholds in
     an empirical-percentile band; fit a GPD above each (method of moments, shape
     clipped to the Luceno band); score by a frequency-weighted MSE between the
     empirical AERs and the GPD-predicted magnitudes; pick mu by the WMSE-
     tolerance set (robust Tukey ceiling) with a shape-stability tie-break.
  4. Bootstrap the exceedances above mu (truncated-noise resampling) and fit a
     GPD per realization; the mean is the best estimate, the 10th/90th
     percentiles are the confidence band.
  5. Splice the GPD tail (AER < lambda_mu) onto the empirical bulk (AER >=
     lambda_mu) and report magnitudes at standard return intervals.

The GPD fit uses method of moments (xi = 1/2 (1 - m^2/v), sigma = m (1 - xi));
the PyStorm default is MLE, which needs scipy - this port stays stdlib + numpy
and is deterministic given the seed.

Input is a CSV of POT peaks (column 1 a date or index, column 2 the magnitude):
the peaks handed off from 10-3, or an uploaded file.

Classification: provisional (method-of-moments GPD and a fixed-seed bootstrap;
a faithful but simplified port of the PyStorm PST).
Theory and references: Coles (2001) GPD/POT; Nadal-Caraballo et al. PST; USACE
coastal-hazards practice.

Self-containment: zero sibling imports; embeds its own contract dataclasses.
Runnable standalone:
    python chessqc_10_4_probabilistic_simulation.py
stdlib + numpy only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime

import numpy as np

_SHAPE_CLIP_LOW = -0.5
_SHAPE_CLIP_HIGH = 0.33


# --- embedded contract dataclasses (self-contained) -----------------------------
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
    next_apps: tuple = ()


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
    show_if: tuple = ()
    data_dir: str = "water_levels"


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    group: str = ""
    x_key: str = ""
    note: str = ""           # hover definition shown on the output label


APP_META = AppMeta(
    aces_id="10-4",
    name="Probabilistic Simulation Technique",
    area="Coastal Hazards",
    classification="provisional",
    cite="Coles (2001); Nadal-Caraballo et al. PST; PyStorm PST",
    default_system="SI",
)

# A small embedded POT-like peak sample (so the app runs without a hand-off):
# 250 peaks drawn from a GPD-ish tail above 1 m, deterministic.
def _sample_csv() -> str:
    rng = np.random.default_rng(7)
    # synthetic exceedances: mu=1.0, sigma=0.35, xi=0.1
    u = rng.uniform(0, 1, 250)
    xi, sig, mu = 0.1, 0.35, 1.0
    vals = mu + sig / xi * ((1 - u) ** (-xi) - 1)
    return "value\n" + "\n".join(f"{v:.4f}" for v in vals)


_SAMPLE_CSV = _sample_csv()

# standard mean return intervals (yr) for the reporting table / plot ticks
_MRI = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                 10000], dtype=np.float64)

INPUTS = (
    Field("csv", "POT peak sample", "csv", default=_SAMPLE_CSV, choices=(),
          note="Peaks over threshold (handed off from 10-3) or an uploaded CSV "
               "with the peak magnitude in column 2."),
    Field("events_per_year", "Events per year", "float", "1/yr", "1/yr", default=10.0,
          lo=0.1, hi=365.0, note="POT target rate; sets the auto record length "
          "(n_pot / events_per_year) when record length is 0"),
    Field("record_length_years", "Record length", "float", "yr", "yr", default=0.0,
          lo=0.0, hi=1e5, note="0 = auto (n_pot / events_per_year)"),
    Field("num_simulations", "Bootstrap simulations", "int", "", "", default=300,
          lo=20, hi=2000, note="Monte-Carlo realizations for the confidence band"),
    Field("seed", "Random seed", "int", "", "", default=12345, lo=0, hi=1_000_000,
          note="fixed for reproducible confidence bands"),
    Field("threshold_min_percentile", "GPD band: min percentile", "float", "%", "%",
          default=50.0, lo=0.0, hi=99.0, note="lower bound of the GPD-location scan band"),
    Field("threshold_max_percentile", "GPD band: max percentile", "float", "%", "%",
          default=95.0, lo=1.0, hi=99.9, note="upper bound of the GPD-location scan band"),
    Field("min_exceedances", "Min exceedances", "int", "", "", default=30, lo=5, hi=100000,
          note="a candidate location must retain at least this many exceedances"),
)

OUTPUTS = (
    Out("n_pot", "POT peaks", "", "", "scalar",
        note="Number of peaks-over-threshold magnitudes in the input sample."),
    Out("gpd_threshold", "GPD location mu", "m", "ft", "scalar",
        note="Selected GPD location (threshold) above which the upper tail is fitted, chosen by quantile delta optimization."),
    Out("gpd_shape", "GPD shape xi", "", "", "scalar",
        note="Method-of-moments GPD shape parameter of the upper-tail fit (negative = bounded tail, positive = heavy tail), clipped to the stability band."),
    Out("gpd_scale", "GPD scale sigma", "m", "ft", "scalar",
        note="Method-of-moments GPD scale parameter of the upper-tail fit, setting the spread of exceedances above mu."),
    Out("lambda_u", "Base rate lambda_u", "1/yr", "1/yr", "scalar",
        note="Population (base) annual exceedance rate of the POT peaks, n_pot divided by the record length."),
    Out("lambda_mu", "Rate above mu", "1/yr", "1/yr", "scalar",
        note="Annual rate of exceedances above the GPD location mu, used to splice the GPD tail onto the empirical bulk."),
    Out("mag_10yr", "10-yr magnitude", "m", "ft", "scalar",
        note="Response magnitude at the 10-year mean return interval (annual exceedance rate 0.1/yr) read from the hazard curve."),
    Out("mag_100yr", "100-yr magnitude", "m", "ft", "scalar",
        note="Response magnitude at the 100-year mean return interval (annual exceedance rate 0.01/yr) read from the hazard curve."),
    Out("mag_500yr", "500-yr magnitude", "m", "ft", "scalar",
        note="Response magnitude at the 500-year mean return interval (annual exceedance rate 0.002/yr) read from the hazard curve."),
    # hazard curve: magnitude vs log10(mean return interval, yr)
    Out("profile_logmri", "Profile: log10 return interval", "log10 yr", "log10 yr", "profile",
        note="Hazard-curve x-axis: base-10 logarithm of the mean return interval in years (= log10 of 1/AER)."),
    Out("profile_be", "Profile: best estimate", "m", "ft", "profile", group="hc",
        note="Best-estimate hazard curve: response magnitude versus return interval (bootstrap-mean GPD tail spliced onto the empirical bulk)."),
    Out("profile_cb10", "Profile: 10% confidence", "m", "ft", "profile", group="hc",
        note="Lower confidence bound of the hazard curve, the 10th percentile of the bootstrap GPD-tail realizations."),
    Out("profile_cb90", "Profile: 90% confidence", "m", "ft", "profile", group="hc",
        note="Upper confidence bound of the hazard curve, the 90th percentile of the bootstrap GPD-tail realizations."),
    # empirical peaks (Weibull plotting positions) as markers
    Out("emp_logmri", "Empirical log10 MRI", "log10 yr", "log10 yr", "scatter_x",
        note="x-coordinate of the empirical peak markers: log10 of the return interval from the Weibull plotting position (= log10 of 1/AER_i)."),
    Out("emp_mag", "Empirical peaks", "m", "ft", "scatter", group="hc", x_key="emp_logmri",
        note="Observed POT peak magnitudes plotted at their Weibull empirical return intervals."),
)


@dataclass
class Result:
    n_pot: float
    gpd_threshold: float
    gpd_shape: float
    gpd_scale: float
    lambda_u: float
    lambda_mu: float
    mag_10yr: float
    mag_100yr: float
    mag_500yr: float
    profile_logmri: np.ndarray
    profile_be: np.ndarray
    profile_cb10: np.ndarray
    profile_cb90: np.ndarray
    emp_logmri: np.ndarray
    emp_mag: np.ndarray
    notes: str = ""


# --- parsing --------------------------------------------------------------------
def _parse_values(text: str) -> np.ndarray:
    """Read the magnitude column (column 2 if present, else column 1) from CSV
    text; skip the header and blank/non-numeric rows. Returns the values."""
    vals: list[float] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(",")
        cell = parts[1].strip() if len(parts) >= 2 else parts[0].strip()
        try:
            v = float(cell)
        except ValueError:
            continue
        if math.isfinite(v):
            vals.append(v)
    return np.asarray(vals, dtype=np.float64)


# --- GPD primitives (method of moments; closed-form ICDF) -----------------------
def _fit_gpd_mom(data: np.ndarray, loc: float):
    """Method-of-moments GPD above `loc`; shape clipped, scale recomputed.
    Returns (xi, scale) or None for a degenerate sample."""
    ex = np.asarray(data, dtype=np.float64) - loc
    if ex.size < 2:
        return None
    m = float(np.mean(ex))
    v = float(np.var(ex, ddof=1))
    if not (v > 0.0 and m > 0.0):
        return None
    xi = 0.5 * (1.0 - m * m / v)
    xi = min(max(xi, _SHAPE_CLIP_LOW), _SHAPE_CLIP_HIGH)
    scale = max(m * (1.0 - xi), 1e-12)
    return xi, scale


def _gpd_ppf(p, xi: float, loc: float, scale: float):
    """GPD inverse CDF (quantile) at non-exceedance probability p."""
    p = np.clip(p, 0.0, 1.0 - 1e-12)
    if abs(xi) < 1e-8:
        return loc - scale * np.log(1.0 - p)
    return loc + scale / xi * ((1.0 - p) ** (-xi) - 1.0)


def _local_dispersion(series: np.ndarray, window: int) -> np.ndarray:
    """Scaled MAD of `series` over a +/-window neighbourhood (shape-stability)."""
    n = series.size
    disp = np.full(n, np.inf)
    for i in range(n):
        a, b = max(0, i - window), min(n, i + window + 1)
        seg = series[a:b]
        seg = seg[np.isfinite(seg)]
        if seg.size:
            disp[i] = 1.4826 * float(np.median(np.abs(seg - np.median(seg))))
    return disp


def _validate(inp: dict) -> None:
    for fld in INPUTS:
        if fld.kind not in ("float", "int"):
            continue
        v = float(inp.get(fld.key, fld.default))
        if not (fld.lo <= v <= fld.hi):
            raise ValueError(f"{fld.label} ({fld.key}) = {v} outside [{fld.lo}, {fld.hi}]")


# --- compute --------------------------------------------------------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Turns a Peaks-Over-Threshold peak sample into a coastal hazard curve '
            '(response magnitude versus annual exceedance rate) by carrying the frequent '
            'range empirically and splicing on a Generalized Pareto Distribution upper '
            'tail, with a bootstrap confidence band. It reports magnitudes at standard '
            'mean return intervals (e.g. 10-, 100-, 500-yr) plus the fitted GPD '
            'parameters.',
 'methods': [{'name': 'Probabilistic Simulation Technique (empirical bulk + GPD tail)',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': '\\lambda_u = \\frac{n_{pot}}{T_R}',
                             'desc': 'Population (base) exceedance rate from the POT count '
                                     'over the record length T_R (auto = n_pot / '
                                     'events-per-year).'},
                            {'tex': '\\mathrm{AER}_i = \\frac{i}{n_{pot} + 1}\\, '
                                    '\\lambda_u',
                             'desc': 'Empirical annual exceedance rate of the rank-i '
                                     '(descending) peak via the Weibull plotting '
                                     'position.'},
                            {'tex': '\\xi = \\frac{1}{2}\\left(1 - \\frac{m^2}{v}\\right)',
                             'desc': 'Method-of-moments GPD shape from the mean m and '
                                     'variance v of the exceedances above the location mu '
                                     '(clipped to the Luceno band).'},
                            {'tex': '\\sigma = m\\,(1 - \\xi)',
                             'desc': 'Method-of-moments GPD scale from the exceedance mean '
                                     'and the fitted shape.'},
                            {'tex': 'x(p) = \\mu + \\frac{\\sigma}{\\xi}\\left[(1 - '
                                    'p)^{-\\xi} - 1\\right]',
                             'desc': 'GPD quantile (inverse CDF) giving the tail magnitude '
                                     'at non-exceedance probability p = 1 - AER/lambda_mu; '
                                     'reduces to mu - sigma ln(1-p) as xi -> 0.'},
                            {'tex': '\\lambda_{\\mu} = \\frac{n_{\\mu}}{T_R}',
                             'desc': 'Rate of exceedances above the GPD location mu, used '
                                     'to map AER to the GPD probability scale and to '
                                     'splice the tail onto the empirical bulk at AER = '
                                     'lambda_mu.'}]}],
 'symbols': [['lambda_u', 'Population (base) annual exceedance rate of POT peaks (1/yr)'],
             ['lambda_mu', 'Annual rate of exceedances above the GPD location mu (1/yr)'],
             ['n_pot', 'Number of peaks-over-threshold in the sample'],
             ['T_R',
              'Effective record length in years (given or auto = n_pot / events-per-year)'],
             ['AER', 'Annual exceedance rate of a peak (1/yr); return interval = 1/AER'],
             ['mu', 'GPD location (threshold) selected by quantile delta optimization'],
             ['xi', 'GPD shape parameter (clipped to the stability band)'],
             ['sigma', 'GPD scale parameter'],
             ['m, v', 'Sample mean and variance of the exceedances above mu'],
             ['p', 'Non-exceedance probability on the GPD, p = 1 - AER/lambda_mu']],
 'references': ['Coles (2001), An Introduction to Statistical Modeling of Extreme Values '
                '(GPD/POT)',
                'Nadal-Caraballo et al., Probabilistic Simulation Technique (PST)',
                'PyStorm PST module',
                'USACE coastal-hazards practice']}


def compute(inp: dict) -> Result:
    _validate(inp)
    values = _parse_values(str(inp.get("csv", _SAMPLE_CSV)))
    n_pot = values.size
    if n_pot < 10:
        raise ValueError(f"need at least 10 POT peaks for PST; got {n_pot}")

    epy = float(inp["events_per_year"])
    rec = float(inp["record_length_years"])
    record_length = rec if rec > 0 else n_pot / epy
    lambda_u = n_pot / record_length

    v_desc = np.sort(values)[::-1]
    ranks = np.arange(1, n_pot + 1)
    weibull_aer = (ranks / (n_pot + 1)) * lambda_u

    # --- Step 3: QDO scan for the GPD location mu ---
    pmin = float(inp["threshold_min_percentile"])
    pmax = float(inp["threshold_max_percentile"])
    min_exc = int(inp["min_exceedances"])
    band_lo = float(np.percentile(v_desc, pmin))
    band_hi = float(np.percentile(v_desc, pmax))
    cand = np.linspace(float(v_desc.min()), band_hi, 50)

    wmse = np.full(cand.size, np.nan)
    shp = np.full(cand.size, np.nan)
    n_exc = np.zeros(cand.size, dtype=np.int64)
    for i, th in enumerate(cand):
        mask = v_desc > th
        pot = v_desc[mask]; aer = weibull_aer[mask]
        n_exc[i] = pot.size
        if np.unique(pot).size <= 1:
            continue
        fit = _fit_gpd_mom(pot, th)
        if fit is None:
            continue
        xi, sc = fit
        lam_mu = pot.size / record_length
        prob = np.clip(1.0 - aer / lam_mu, 1e-12, 1.0 - 1e-12)
        pred = _gpd_ppf(prob, xi, th, sc)
        wm = (aer < 1.0) & np.isfinite(pred)
        if not np.any(wm):
            continue
        w = 1.0 / aer[wm]
        wmse[i] = float(np.sum(w * (pot[wm] - pred[wm]) ** 2) / np.sum(w))
        shp[i] = xi

    in_band = (cand >= band_lo - 1e-9) & (cand <= band_hi + 1e-9)
    band_mask = in_band & (n_exc >= min_exc)
    wmse_band = np.where(band_mask, wmse, np.nan)
    if not np.any(np.isfinite(wmse_band)):
        raise ValueError("GPD-location search failed: no in-band candidate kept "
                         f">= {min_exc} exceedances with a finite fit. Lower the "
                         "minimum exceedances or raise the max percentile.")
    best = float(np.nanmin(wmse_band))
    finite = wmse_band[np.isfinite(wmse_band)]
    q1, q3 = np.percentile(finite, [25.0, 75.0])
    fence = q3 + 1.5 * (q3 - q1)
    inlier = finite[finite <= fence]
    upper = float(np.max(inlier)) if inlier.size else float(np.max(finite))
    ceiling = best + 0.05 * (upper - best)
    sel = np.where(wmse_band <= ceiling)[0]
    if sel.size == 0:
        sel = np.array([int(np.nanargmin(wmse_band))])
    xi_disp = _local_dispersion(shp, 3)
    best_idx = int(sel[int(np.argmin(xi_disp[sel]))])     # tie-break: min xi-dispersion
    threshold = float(cand[best_idx])

    above = v_desc > threshold
    pot_above = v_desc[above]
    aer_above = weibull_aer[above]
    pot_below = v_desc[~above]
    aer_below = weibull_aer[~above]
    if pot_above.size < 2:
        raise ValueError("fewer than 2 peaks exceed the GPD location; widen the band")
    lambda_mu = pot_above.size / record_length

    # --- Step 4: truncated-noise bootstrap of the exceedances ---
    n_sims = int(inp["num_simulations"])
    rng = np.random.default_rng(int(inp["seed"]))
    psort = np.sort(pot_above)[::-1]
    delta = np.append(np.diff(psort), 0.0)            # spacing to the next-smaller value
    n_ab = psort.size

    # dense AER grid for the GPD tail: from just below lambda_mu down to 1e-4
    aer_gpd = np.logspace(math.log10(min(lambda_mu * 0.999, lambda_u)), -4.0, 240)
    aer_gpd = aer_gpd[aer_gpd < lambda_mu]
    q_gpd = 1.0 - aer_gpd / lambda_mu                  # GPD non-exceedance prob

    ens = np.full((n_sims, aer_gpd.size), np.nan)
    for j in range(n_sims):
        idx = rng.integers(0, n_ab, n_ab)
        noise = _trunc_norm(rng, -1.0, 1.0, n_ab)
        sample = np.sort(psort[idx] + delta[idx] * noise)[::-1]
        fit = _fit_gpd_mom(sample, threshold)
        if fit is None:
            continue
        xi, sc = fit
        ens[j, :] = _gpd_ppf(q_gpd, xi, threshold, sc)

    be = np.nanmean(ens, axis=0)
    cb10 = np.nanpercentile(ens, 10, axis=0)
    cb90 = np.nanpercentile(ens, 90, axis=0)

    # best-estimate GPD parameters (fit on the actual exceedances)
    fit0 = _fit_gpd_mom(pot_above, threshold)
    xi0, sc0 = fit0 if fit0 else (float("nan"), float("nan"))

    # --- Step 5: splice GPD tail + empirical bulk; report key return periods ---
    # full curve AER (descending so magnitude increases): GPD tail then bulk
    curve_aer = np.concatenate([aer_gpd, aer_below])
    curve_be = np.concatenate([be, pot_below])
    curve_c10 = np.concatenate([cb10, pot_below])
    curve_c90 = np.concatenate([cb90, pot_below])
    order = np.argsort(curve_aer)[::-1]               # AER descending
    curve_aer = curve_aer[order]; curve_be = curve_be[order]
    curve_c10 = curve_c10[order]; curve_c90 = curve_c90[order]

    logmri = np.log10(1.0 / curve_aer)               # x-axis: log10 return interval (yr)

    def _at(mri):
        aer = 1.0 / mri
        if aer >= curve_aer[0] or aer <= curve_aer[-1]:
            return float("nan")
        return float(np.interp(-math.log10(aer), -np.log10(curve_aer), curve_be))
    mag10, mag100, mag500 = _at(10.0), _at(100.0), _at(500.0)

    emp_logmri = np.log10((n_pot + 1) / ranks / lambda_u)   # = log10(1/weibull_aer)

    notes = (f"n_pot={n_pot}, record {record_length:.1f} yr; mu={threshold:.3f} "
             f"(xi={xi0:+.3f}, sigma={sc0:.3f}); lambda_u={lambda_u:.2f}, "
             f"lambda_mu={lambda_mu:.2f}/yr; {n_sims} bootstraps")

    return Result(
        n_pot=float(n_pot), gpd_threshold=threshold, gpd_shape=xi0, gpd_scale=sc0,
        lambda_u=lambda_u, lambda_mu=lambda_mu,
        mag_10yr=mag10, mag_100yr=mag100, mag_500yr=mag500,
        profile_logmri=logmri, profile_be=curve_be,
        profile_cb10=curve_c10, profile_cb90=curve_c90,
        emp_logmri=emp_logmri, emp_mag=v_desc, notes=notes,
    )


def _trunc_norm(rng, lo: float, hi: float, n: int) -> np.ndarray:
    """Standard normal truncated to [lo, hi] by rejection."""
    out = rng.standard_normal(n)
    bad = (out < lo) | (out > hi)
    while bad.any():
        out[bad] = rng.standard_normal(int(bad.sum()))
        bad = (out < lo) | (out > hi)
    return out


# --- self-tests -----------------------------------------------------------------
def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}
    r = compute(base)
    assert r.n_pot == 250
    assert np.isfinite(r.gpd_threshold) and np.isfinite(r.gpd_shape)
    assert _SHAPE_CLIP_LOW - 1e-9 <= r.gpd_shape <= _SHAPE_CLIP_HIGH + 1e-9
    # confidence band is ordered and the curve increases toward rarer events
    fin = np.isfinite(r.profile_be)
    assert np.all(r.profile_cb10[fin] <= r.profile_be[fin] + 1e-6)
    assert np.all(r.profile_be[fin] <= r.profile_cb90[fin] + 1e-6)
    # rarer (higher MRI) -> larger magnitude (monotone best estimate)
    bex = r.profile_be[fin]
    assert bex[-1] >= bex[0] - 1e-6, (bex[0], bex[-1])
    # key return-period magnitudes increase with return period
    assert r.mag_10yr <= r.mag_100yr <= r.mag_500yr, (r.mag_10yr, r.mag_100yr, r.mag_500yr)
    # determinism (fixed seed)
    r2 = compute(base)
    assert abs(r2.mag_100yr - r.mag_100yr) < 1e-9
    print("  self-tests: PASS (GPD fit, band order, monotone curve, return periods, determinism)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"    n_pot         = {int(r.n_pot)}")
    print(f"    GPD mu        = {r.gpd_threshold:.3f} m  (xi={r.gpd_shape:+.3f}, sigma={r.gpd_scale:.3f})")
    print(f"    lambda_u/mu   = {r.lambda_u:.2f} / {r.lambda_mu:.2f} /yr")
    print(f"    10/100/500 yr = {r.mag_10yr:.2f} / {r.mag_100yr:.2f} / {r.mag_500yr:.2f} m")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
