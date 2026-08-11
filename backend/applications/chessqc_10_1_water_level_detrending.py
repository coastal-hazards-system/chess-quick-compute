"""CHESS-QC application 10-1 - Water Level Detrending.

Functional area: Coastal Hazards. Estimates the long-term sea-level trend of a
water-level record and removes it, leaving a series that varies about the tidal
datum rather than about a moving mean.

The preferred estimator is the NOAA CO-OPS 053 (Zervas 2009) procedure, Steps
1-8, as implemented by the PyStorm water-level-detrending engine:

  1-2. Average the record within each calendar month and keep only months that
       clear an explicit completeness floor (months are never interpolated
       across). Each retained month is stamped at its midpoint in decimal years.
  3-4. Regress the monthly means on centered time plus twelve calendar-month
       indicators - and, where the operator declares one, a level-shift dummy
       per confirmed step discontinuity:

           y_i = b*tau_i + sum_j m_j*D_ij [+ sum_k d_k*f_ik] + eps_i

       There is no separate intercept: the monthly constants absorb both the
       level and the mean seasonal cycle, which makes the slope invariant to the
       datum the series is expressed in.
  5-6. Fit by ordinary least squares, then estimate the lag-1 autocorrelation of
       the residuals over CONSECUTIVE month pairs only, so gaps neither pair
       across nor bias rho toward zero.
    7. Cochrane-Orcutt: quasi-difference the consecutive pairs by rho, refit,
       re-estimate rho on the original-unit residuals, and iterate to
       convergence. Serially correlated monthly means otherwise make the trend
       look far better determined than it is.
    8. Take the standard error from the converged GLS normal equations with
       df = N* - p, and report the 95% confidence interval t(0.975, df)*se.

The fit needs no pre-specified origin: it centers on the mean month internally,
and the slope is unaffected by that choice. The National Tidal Datum Epoch
enters only when the trend is subtracted,

    detrended(t) = y(t) - b*(t - t_0),   t_0 = (epoch start + epoch end) / 2

so the detrended series sits on the epoch datum (t_0 = 1992.0 for the 1983-2001
NTDE - the epoch's center year, not the midpoint of its calendar span). Ordinary
least squares on the raw samples, and a directly supplied slope (a published
regional or global rate), are offered as alternatives.

Input is a CSV record with column 1 a date and column 2 the water level (meters):
either a bundled NOAA station or a user-supplied file. Blank water-level cells
(gaps) and a header row are ignored. Hourly records and ready-made monthly-mean
series are both accepted; the completeness floor scales with the record's own
cadence. The full record resolution is used for the fit; the returned profile
series are decimated only for display when the record is large.

Classification: standard (the published CO-OPS 053 estimator, reproducing NOAA's
own trends for the bundled stations to well within their confidence intervals).
Theory and references: Zervas (2009), "Sea Level Variations of the United States
1854-2006", NOAA Technical Report NOS CO-OPS 053; NOAA CO-OPS tidal datum
(NTDE) conventions; Cochrane & Orcutt (1949).

Self-containment: zero sibling imports; embeds its own contract dataclasses.
Runnable standalone:
    python chessqc_10_1_water_level_detrending.py
stdlib + numpy only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime

import numpy as np


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
    next_apps: tuple = ()      # workflow "Next" targets: ((id, label), ...) carrying the series


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
    enable_if: tuple = ()    # (other_key, value): gray out (disable) unless that input == value
    show_if: tuple = ()        # (other_key, value): show this field only when inp[other_key] == value


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    group: str = ""        # profile panel id; profiles sharing a group plot together
    note: str = ""           # hover definition shown on the output label


# --- application metadata --------------------------------------------------------
APP_META = AppMeta(
    aces_id="10-1",
    name="Water Level Detrending",
    area="Coastal Hazards",
    classification="standard",
    cite="Zervas (2009) NOAA CO-OPS 053; Cochrane & Orcutt (1949); NTDE datum convention",
    default_system="SI",
    next_apps=(("10-2", "Non-Tidal Residual"), ("10-3", "Peaks Over Threshold")),
)

# Bundled NOAA CO-OPS stations (id|label); the front-ends fetch
# data/water_levels/<id>.csv. Records are full-resolution hourly water level (m).
STATIONS = (
    "8518750|The Battery, NY",
    "8651370|Duck, NC",
    "8724580|Key West, FL",
    "8761724|Grand Isle, LA",
    "8771450|Galveston Pier 21, TX",
)

_FIT_COOPS = "CO-OPS 053 (monthly GLS)"
_FIT_OLS = "Ordinary least squares"
_FIT_GIVEN = "Specified slope"

# Embedded monthly-mean sample (the default, so the app computes without a file):
# 30 years rising at 3 mm/yr, carrying the mean seasonal cycle and AR(1)-correlated
# monthly noise that the CO-OPS 053 estimator exists to handle. Deterministic (fixed
# generator seed), so the shipped example never moves.
_SEASONAL = (-0.08, -0.09, -0.05, 0.01, 0.03, 0.05,
             0.06, 0.08, 0.07, 0.02, -0.03, -0.07)     # m, mean ~ 0


def _sample_csv(slope: float = 0.003, rho: float = 0.4, sigma: float = 0.03,
                years: int = 30, start: int = 1991, base: float = 0.10) -> str:
    rng = np.random.default_rng(0)
    n = years * 12
    innov = rng.normal(0.0, sigma, n)
    rows = ["date,water_level_m"]
    eps = 0.0
    for i in range(n):
        yr, mo = start + i // 12, 1 + i % 12
        eps = rho * eps + innov[i]
        t = yr + (mo - 0.5) / 12.0
        v = base + slope * (t - (start + years / 2.0)) + _SEASONAL[mo - 1] + eps
        rows.append(f"{yr}-{mo:02d}-15,{v:.4f}")
    return "\n".join(rows)


_SAMPLE_CSV = _sample_csv()

# Cap on the number of points returned in the display profiles (the fit always
# uses every valid sample; only the plotted series is strided when larger).
_PLOT_MAX = 4000

_DAYS_TO_EPOCH = 719163          # datetime(1970,1,1).toordinal()
_SECONDS_PER_DAY = 86400.0

INPUTS = (
    Field("csv", "Water-level record", "csv", default=_SAMPLE_CSV, choices=STATIONS,
          note="Select a bundled NOAA station or upload your own CSV "
               "(column 1 = date, column 2 = water level in m). Header and blank "
               "water-level rows are ignored. Hourly records and ready-made "
               "monthly means are both accepted."),
    Field("fit_mode", "Trend estimator", "choice", default=_FIT_COOPS,
          choices=(_FIT_COOPS, _FIT_OLS, _FIT_GIVEN),
          note="CO-OPS 053 fits monthly means with a seasonal cycle and AR(1) "
               "errors (NOAA's published procedure); ordinary least squares "
               "regresses the raw samples on time; or supply a slope directly."),
    Field("min_valid_days", "Month completeness floor", "float", "d", "d", default=28.0,
          lo=1.0, hi=31.0,
          note="Days' worth of valid samples a calendar month must hold to enter "
               "the fit; incomplete months are excluded, never interpolated.",
          enable_if=("fit_mode", _FIT_COOPS)),
    Field("level_shifts", "Level shifts", "list", "yr", "yr", default=(),
          note="Decimal years of confirmed step discontinuities (gauge relocation, "
               "datum shift), as a JSON list, e.g. [1947.5]. Each adds a dummy that "
               "absorbs the step instead of letting it bias the slope.",
          enable_if=("fit_mode", _FIT_COOPS)),
    Field("ntde_start", "NTDE start year", "int", "yr", "yr", default=1983,
          lo=1800, hi=2100,
          note="National Tidal Datum Epoch start year; with the end year it sets "
               "the epoch center the detrended series is referenced to. It plays "
               "no part in the trend fit."),
    Field("ntde_end", "NTDE end year", "int", "yr", "yr", default=2001,
          lo=1800, hi=2100, note="National Tidal Datum Epoch end year (inclusive)."),
    Field("slope_value", "Specified slope", "float", "mm/yr", "in/yr", default=0.003,
          lo=-1.0, hi=1.0, note="used only when the estimator is a supplied slope",
          show_if=("fit_mode", _FIT_GIVEN)),
)

OUTPUTS = (
    Out("slope_per_year", "Sea-level trend (slope)", "mm/yr", "in/yr", "scalar",
        note="Fitted (or specified) linear sea-level rate; positive means the water "
             "level is rising over time."),
    Out("ci_halfwidth", "Trend 95% CI half-width", "mm/yr", "in/yr", "scalar",
        note="Half-width of the 95 percent confidence interval on the trend, "
             "t(0.975, df) times the GLS standard error; zero when the slope is "
             "supplied rather than fitted."),
    Out("rho_ar1", "Residual AR(1) rho", "", "", "scalar",
        note="Converged lag-1 autocorrelation of the monthly residuals, estimated "
             "over consecutive month pairs only; the Cochrane-Orcutt correction "
             "widens the confidence interval by roughly sqrt((1+rho)/(1-rho))."),
    Out("seasonal_range", "Fitted seasonal range", "m", "ft", "scalar",
        note="Peak-to-trough range of the twelve fitted calendar-month constants, "
             "i.e. the mean annual cycle removed from the trend estimate."),
    Out("n_months", "Complete months in fit", "", "", "scalar",
        note="Calendar months clearing the completeness floor and entering the "
             "regression (N); incomplete months are excluded."),
    Out("dof", "Degrees of freedom", "", "", "scalar",
        note="Degrees of freedom of the GLS fit, N* - p, with N* the consecutive "
             "quasi-differenced month pairs and p the regression columns."),
    Out("epoch_year", "Epoch center (datum) year", "yr", "yr", "scalar",
        note="Decimal year the detrended series is referenced to: the center of "
             "the National Tidal Datum Epoch, where detrended equals observed."),
    Out("total_trend", "Total trend over record", "m", "ft", "scalar",
        note="Net change in level attributable to the trend across the record, equal "
             "to slope times the record length."),
    Out("record_years", "Record length", "yr", "yr", "scalar",
        note="Time span of the record in years, from the first to the last valid "
             "sample."),
    Out("n_samples", "Samples in record", "", "", "scalar",
        note="Number of finite (date, water-level) samples in the record; "
             "blank-value gaps are excluded."),
    Out("rms_residual", "RMS residual about trend", "m", "ft", "scalar",
        note="Root-mean-square scatter of the detrended series about its mean; the "
             "typical residual variability remaining after the trend is removed."),
    Out("converged", "Fit converged", "", "", "scalar",
        note="yes when the Cochrane-Orcutt iteration met its tolerance on rho "
             "within the iteration cap."),
    # vertical marker on the plots at the epoch center; plotted on the x (year) axis
    Out("pivot_line", "Epoch center", "yr", "yr", "vline",
        note="Vertical marker at the tidal-datum epoch center, where the detrended "
             "series coincides with the observed record."),
    Out("profile_year", "Profile: year", "yr", "yr", "profile",
        note="Decimal calendar year (x-axis) for the plotted series, decimated for "
             "display when the record is large."),
    # Panel 1 (group "obs"): observed level with the fitted linear trend on top.
    Out("profile_original", "Profile: observed", "m", "ft", "profile", group="obs",
        note="Observed water level over time; NaN at gaps so the line breaks across "
             "missing data."),
    Out("profile_trend", "Profile: linear trend", "m", "ft", "profile", group="obs",
        note="Fitted linear-trend line overlaid on the observed level, passing through "
             "the datum at the epoch center."),
    # Panel 2 (group "detr"): detrended level with the horizontal datum.
    Out("profile_detrended", "Profile: detrended", "m", "ft", "profile", group="detr",
        note="Water level with the linear trend removed, referenced to the epoch "
             "center; the residual variability about the trend."),
    Out("profile_datum", "Profile: datum", "m", "ft", "profile", group="detr",
        note="Horizontal reference level (mean of the detrended series) drawn across "
             "the detrended panel."),
    # full-resolution detrended series for the workflow hand-off (emitted only when
    # the `handoff` input is set); carried into 10-2 / 10-3. Not shown or plotted.
    Out("handoff_csv", "handoff", "", "", "data",
        note="Full-resolution detrended series (year, level) passed in-memory to the "
             "10-2 / 10-3 workflow apps; not displayed."),
)


@dataclass
class Result:
    slope_per_year: float
    ci_halfwidth: float
    rho_ar1: float
    seasonal_range: float
    n_months: float
    dof: float
    epoch_year: float
    total_trend: float
    record_years: float
    n_samples: float
    rms_residual: float
    converged: str
    pivot_line: float
    profile_year: np.ndarray
    profile_original: np.ndarray
    profile_trend: np.ndarray
    profile_detrended: np.ndarray
    profile_datum: np.ndarray
    handoff_csv: str = ""
    notes: str = ""


# --- Student-t quantile (keeps the module stdlib + numpy; no scipy) --------------
def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function (modified Lentz)."""
    tiny = 1e-30
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-16:
            break
    return h


def _betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(lbeta)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _t_sf(t: float, df: float) -> float:
    """Upper-tail probability P(T > t) of Student's t with df degrees of freedom."""
    p_two = _betainc(0.5 * df, 0.5, df / (df + t * t))
    return 0.5 * p_two if t >= 0.0 else 1.0 - 0.5 * p_two


def _t_ppf(p: float, df: float) -> float:
    """Quantile of Student's t: the t with P(T <= t) = p (bisection on the CDF)."""
    if not (0.0 < p < 1.0) or df <= 0:
        return float("nan")
    lo, hi = 0.0, 1.0
    while _t_sf(hi, df) > 1.0 - p and hi < 1e6:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _t_sf(mid, df) > 1.0 - p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --- date / CSV parsing ---------------------------------------------------------
_YEAR_BOUNDS: dict[int, tuple] = {}


def _decimal_year(s: str) -> float:
    """Convert a 'YYYY-MM-DD HH:MM' / 'YYYY-MM-DD' / 'YYYY-MM' string to a
    fractional calendar year (leap-year aware). Raises ValueError on a bad date."""
    return _parse_stamp(s)[0]


def _parse_stamp(s: str) -> tuple:
    """'YYYY-MM-DD HH:MM' -> (decimal year, absolute seconds, year, month).

    A bare number is read as a decimal year (the hand-off format), from which the
    calendar year and month are recovered. Absolute seconds come from the proleptic
    Gregorian day number, so pre-1970 records (some gauges start in the 1850s) are
    handled without platform timestamp limits."""
    s = s.strip()
    sep = " " if " " in s else ("T" if "T" in s else "")
    datepart, timepart = (s.split(sep, 1) if sep else (s, ""))
    if "-" not in datepart:          # bare (decimal) calendar year
        t = float(datepart)
        y = int(math.floor(t))
        frac = t - y
        st, span = _year_bounds(y)
        secs = st * _SECONDS_PER_DAY + frac * span
        mo = min(12, 1 + int(frac * 12.0))
        return t, secs, y, mo
    p = datepart.split("-")
    y = int(p[0])
    mo = int(p[1]) if len(p) > 1 and p[1] else 1
    d = int(p[2]) if len(p) > 2 and p[2] else 1
    hh = mm = 0
    if timepart:
        tp = timepart.split(":")
        hh = int(tp[0]) if tp[0] else 0
        mm = int(tp[1]) if len(tp) > 1 and tp[1] else 0
    ordinal = datetime(y, mo, d).toordinal() - _DAYS_TO_EPOCH
    secs = ordinal * _SECONDS_PER_DAY + hh * 3600.0 + mm * 60.0
    st, span = _year_bounds(y)
    return y + (secs - st * _SECONDS_PER_DAY) / span, secs, y, mo


def _year_bounds(y: int) -> tuple:
    """(day number of Jan 1 since the epoch, seconds in the year) for year y."""
    if y not in _YEAR_BOUNDS:
        st = datetime(y, 1, 1).toordinal() - _DAYS_TO_EPOCH
        en = datetime(y + 1, 1, 1).toordinal() - _DAYS_TO_EPOCH
        _YEAR_BOUNDS[y] = (st, (en - st) * _SECONDS_PER_DAY)
    return _YEAR_BOUNDS[y]


def _parse_csv(text: str) -> tuple:
    """Parse CSV text -> (decimal_year, seconds, year, month, value) arrays, sorted
    by time. Column 1 is the date, column 2 the value. Rows with a valid date but a
    blank/non-numeric value are kept with value = NaN, so gaps in the record stay as
    gaps (the fit ignores them and the plot breaks the line across them). The header
    row and rows with an unparseable date are dropped."""
    years: list = []
    secs: list = []
    yy: list = []
    mm: list = []
    vals: list = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(",")
        if len(parts) < 2:
            continue
        ds = parts[0].strip()
        if not ds:
            continue
        try:
            t, sec, y, mo = _parse_stamp(ds)
        except (ValueError, IndexError):
            continue  # header or unparseable date
        try:
            v = float(parts[1].strip())
            if not math.isfinite(v):
                v = math.nan
        except ValueError:
            v = math.nan          # blank / non-numeric -> gap
        years.append(t); secs.append(sec); yy.append(y); mm.append(mo); vals.append(v)
    t = np.asarray(years, dtype=np.float64)
    sec = np.asarray(secs, dtype=np.float64)
    y_arr = np.asarray(yy, dtype=np.int64)
    m_arr = np.asarray(mm, dtype=np.int64)
    v = np.asarray(vals, dtype=np.float64)
    if int(np.isfinite(v).sum()) < 2:
        raise ValueError(
            "need at least 2 valid (date, water level) rows to fit a trend; "
            f"parsed {len(years)} rows, {int(np.isfinite(v).sum())} with values"
        )
    order = np.argsort(sec, kind="stable")
    return t[order], sec[order], y_arr[order], m_arr[order], v[order]


# --- CO-OPS 053 Steps 1-2: monthly means with a completeness floor ---------------
def _monthly_means(sec: np.ndarray, yr: np.ndarray, mo: np.ndarray, v: np.ndarray,
                   min_valid_days: float) -> tuple:
    """Calendar-month means of the record -> (t, mean, serial, month) for the months
    that clear the floor. `t` is the decimal-year midpoint of each month. The floor
    is expressed in days of valid data and converted with the record's own modal
    sampling step, so hourly records and ready-made monthly means both work."""
    fin = np.isfinite(v)
    if int(fin.sum()) < 2:
        raise ValueError("no valid samples to average into monthly means")
    d = np.diff(np.unique(sec))
    d = d[d > 0]
    if d.size == 0:
        raise ValueError("need at least 2 distinct sample times to infer the cadence")
    step = float(np.median(d))
    n_required = max(1, int(round(min_valid_days * _SECONDS_PER_DAY / step)))

    serial_all = yr * 12 + (mo - 1)
    serial, inv = np.unique(serial_all[fin], return_inverse=True)
    counts = np.bincount(inv)
    sums = np.bincount(inv, weights=v[fin])
    keep = counts >= n_required
    if int(keep.sum()) == 0:
        raise ValueError(
            f"no calendar month holds {min_valid_days:g} days of valid data "
            f"({n_required} samples at the record's {step / 3600.0:.3g} h cadence); "
            "lower the completeness floor")
    serial = serial[keep]
    means = (sums[keep] / counts[keep])
    m_of = (serial % 12) + 1
    y_of = serial // 12
    t = np.empty(serial.size, dtype=np.float64)
    for i in range(serial.size):
        y_i, m_i = int(y_of[i]), int(m_of[i])
        start = datetime(y_i, m_i, 1).toordinal() - _DAYS_TO_EPOCH
        end = (datetime(y_i + 1, 1, 1) if m_i == 12
               else datetime(y_i, m_i + 1, 1)).toordinal() - _DAYS_TO_EPOCH
        mid = (start + end) / 2.0 * _SECONDS_PER_DAY
        y0, span = _year_bounds(y_i)
        t[i] = y_i + (mid - y0 * _SECONDS_PER_DAY) / span
    return t, means, serial, m_of


# --- CO-OPS 053 Steps 3-8: seasonal + trend regression with AR(1) GLS ------------
def _fit_seasonal_trend(t: np.ndarray, y: np.ndarray, serial: np.ndarray,
                        month: np.ndarray, shifts: tuple = (),
                        ar_tol: float = 1e-4, ar_max_iter: int = 25) -> dict:
    """Cochrane-Orcutt AR(1) GLS fit of trend + 12 monthly constants (+ level-shift
    dummies). Returns the slope, its GLS standard error and CI half-width, the
    converged rho, and the fit's shape diagnostics."""
    n = t.size
    if n < 24:
        raise ValueError(
            f"need at least 24 complete months to fit trend + seasonal cycle; got {n}")
    if np.unique(month).size < 12:
        raise ValueError("record does not cover all 12 calendar months; the "
                         "seasonal constants are unidentifiable")

    # design: centered time, 12 calendar-month indicators (no separate intercept,
    # which keeps the slope invariant to the datum), then one dummy per level shift
    masks = []
    for s in shifts:
        mask = t < float(s)
        if mask.all() or not mask.any():
            raise ValueError(f"level shift {float(s):g} lies outside the fitted record")
        masks.append(mask)
    p = 13 + len(masks)
    X = np.zeros((n, p))
    t_center = float(np.mean(t))
    X[:, 0] = t - t_center
    for j in range(12):
        X[:, 1 + j] = (month == j + 1).astype(float)
    for k, mask in enumerate(masks):
        X[:, 13 + k] = mask.astype(float)

    consecutive = np.zeros(n, dtype=bool)
    consecutive[1:] = np.diff(serial) == 1
    idx = np.flatnonzero(consecutive)
    n_star = int(idx.size)
    if n_star <= p:
        raise ValueError(
            f"only {n_star} consecutive month pairs for {p} regression columns; "
            "the record is too fragmented for the AR(1) GLS fit")

    def lstsq(A, b):
        return np.linalg.lstsq(A, b, rcond=None)[0]

    def lag1(resid):
        den = float(np.dot(resid[idx - 1], resid[idx - 1]))
        if den == 0.0:
            return 0.0
        return float(np.clip(float(np.dot(resid[idx], resid[idx - 1])) / den, -0.99, 0.99))

    beta_ols = lstsq(X, y)                                   # Step 5
    resid_ols = y - X @ beta_ols
    se_ols = float(np.sqrt(float(resid_ols @ resid_ols) / (n - p)
                           * np.linalg.pinv(X.T @ X)[0, 0]))
    rho = lag1(resid_ols)                                    # Step 6

    converged, iterations = False, 0
    for iterations in range(1, ar_max_iter + 1):             # Step 7
        beta = lstsq(X[idx] - rho * X[idx - 1], y[idx] - rho * y[idx - 1])
        rho_new = lag1(y - X @ beta)
        if abs(rho_new - rho) < ar_tol:
            rho = rho_new
            converged = True
            break
        rho = rho_new

    ys = y[idx] - rho * y[idx - 1]                           # Step 8
    Xs = X[idx] - rho * X[idx - 1]
    beta = lstsq(Xs, ys)
    resid_star = ys - Xs @ beta
    dof = n_star - p
    sigma2 = float(resid_star @ resid_star) / dof
    se = float(np.sqrt(sigma2 * np.linalg.pinv(Xs.T @ Xs)[0, 0]))
    monthly = beta[1:13]
    return {
        "slope": float(beta[0]), "slope_ols": float(beta_ols[0]), "se": se,
        "se_ols": se_ols, "ci": float(_t_ppf(0.975, dof) * se), "dof": int(dof),
        "rho": float(rho), "iterations": int(iterations), "converged": converged,
        "n_months": int(n), "n_star": n_star,
        "seasonal_range": float(np.ptp(monthly)), "monthly": monthly,
        "shift_coeffs": beta[13:],
    }


def _validate(inp: dict) -> None:
    for fld in INPUTS:
        if fld.kind not in ("float", "int"):
            continue
        v = float(inp.get(fld.key, fld.default))
        if not (fld.lo <= v <= fld.hi):
            raise ValueError(f"{fld.label} ({fld.key}) = {v} outside [{fld.lo}, {fld.hi}]")
    if int(inp.get("ntde_end", 2001)) < int(inp.get("ntde_start", 1983)):
        raise ValueError("NTDE end year must be >= start year")


def _decimate(*arrays: np.ndarray, nmax: int = _PLOT_MAX) -> tuple:
    """Stride large arrays down to <= nmax points for display (keeps the last
    point). Returns the arrays unchanged when already small enough."""
    n = len(arrays[0])
    if n <= nmax:
        return arrays
    stride = int(math.ceil(n / nmax))
    idx = np.arange(0, n, stride)
    if idx[-1] != n - 1:
        idx = np.append(idx, n - 1)
    return tuple(a[idx] for a in arrays)


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Estimates the long-term sea-level trend of a water-level record and removes '
            'it. The preferred estimator is the NOAA CO-OPS 053 procedure: calendar-month '
            'means that clear a completeness floor are regressed on centered time plus '
            'twelve monthly constants (and any declared level shifts), then '
            'Cochrane-Orcutt iteration corrects the AR(1) serial correlation of the '
            'residuals and sets the confidence interval. The fit needs no pre-specified '
            'origin; the tidal-datum epoch enters only when the trend is subtracted.',
 'method_key': 'fit_mode',
 'methods': [{'name': 'CO-OPS 053 monthly GLS (seasonal + AR(1))',
              'when': 'CO-OPS 053 (monthly GLS)',
              'tag': 'preferred',
              'note': "NOAA's published sea-level trend procedure (Zervas 2009, Steps "
                      '1-8): monthly means, seasonal constants, Cochrane-Orcutt AR(1) '
                      'GLS, and a t-based confidence interval.',
              'equations': [{'tex': 'y_i = b\\,\\tau_i + \\sum_{j=1}^{12} m_j D_{ij} + '
                                    '\\sum_k d_k f_{ik} + \\varepsilon_i',
                             'desc': 'Seasonal-plus-trend model on the complete monthly '
                                     'means: centered time, twelve calendar-month '
                                     'indicators (which absorb the level, so no separate '
                                     'intercept), and one dummy per declared level shift'},
                            {'tex': '\\tau_i = t_i - \\bar{t}',
                             'desc': 'Time centered on the mean month of the fit; the '
                                     'slope does not depend on this origin, so none has '
                                     'to be specified'},
                            {'tex': '\\rho = \\frac{\\sum_{i \\in C} e_i\\,e_{i-1}}{\\sum_{i '
                                    '\\in C} e_{i-1}^{2}}',
                             'desc': 'Lag-1 residual autocorrelation over the set C of '
                                     'consecutive month pairs only, so gaps neither pair '
                                     'across nor bias rho toward zero'},
                            {'tex': 'y_i^{*} = y_i - \\rho\\,y_{i-1}, \\quad X_i^{*} = X_i '
                                    '- \\rho\\,X_{i-1}',
                             'desc': 'Cochrane-Orcutt quasi-differencing of the '
                                     'consecutive pairs; refit, re-estimate rho, and '
                                     'iterate to convergence'},
                            {'tex': 's_b = \\sqrt{\\sigma^{2}\\left[(X^{*\\top}X^{*})^{-1}'
                                    '\\right]_{11}}, \\quad \\nu = N^{*} - p',
                             'desc': 'Standard error of the trend from the converged GLS '
                                     'normal equations, with N* quasi-differenced rows and '
                                     'p regression columns'},
                            {'tex': 'b \\pm t_{0.975,\\nu}\\, s_b',
                             'desc': '95 percent confidence interval on the trend'},
                            {'tex': 'y_d(t) = y(t) - b\\,(t - t_0), \\quad t_0 = '
                                    '\\frac{y_s + y_e}{2}',
                             'desc': 'Detrended series, pivoted at the epoch center so it '
                                     'stays on the tidal datum (1992.0 for the 1983-2001 '
                                     'NTDE)'}]},
             {'name': 'Ordinary least squares',
              'when': 'Ordinary least squares',
              'tag': 'standard',
              'note': 'Straight regression of every sample on time: no seasonal cycle and '
                      'no serial-correlation correction, so it carries no confidence '
                      'interval. Useful for short records that cannot support the monthly '
                      'fit.',
              'equations': [{'tex': 'b = \\frac{\\sum (t_i - \\bar{t})(y_i - '
                                    '\\bar{y})}{\\sum (t_i - \\bar{t})^{2}}',
                             'desc': 'Least-squares slope about the data centroid, '
                                     'computed on the finite samples only'},
                            {'tex': 'y_d(t) = y(t) - b\\,(t - t_0)',
                             'desc': 'Detrended series, pivoted at the epoch center'}]},
             {'name': 'Specified slope',
              'when': 'Specified slope',
              'tag': '',
              'note': 'Applies a slope supplied by the operator (a published regional or '
                      'global rate) instead of fitting one.',
              'equations': [{'tex': 'y_d(t) = y(t) - b\\,(t - t_0)',
                             'desc': 'Detrended series using the supplied rate b'}]}],
 'symbols': [['b', 'Sea-level trend (slope), m/yr'],
             ['t_i', 'Decimal calendar year of monthly mean i (its midpoint)'],
             ['y_i', 'Monthly mean water level i, m'],
             ['bar t', 'Mean time of the fitted months (the internal centering)'],
             ['m_j', 'Fitted constant for calendar month j, m'],
             ['D_ij', 'Indicator: 1 when month i is calendar month j'],
             ['d_k', 'Level-shift coefficient k (pre-shift level relative to post), m'],
             ['f_ik', 'Level-shift indicator: 1 before shift k, 0 after'],
             ['e_i', 'Regression residual of month i, m'],
             ['rho', 'Lag-1 autocorrelation of the residuals'],
             ['N', 'Complete months entering the fit'],
             ['N^*', 'Consecutive quasi-differenced month pairs'],
             ['p', 'Regression columns (13 + level shifts)'],
             ['nu', 'Degrees of freedom, N* - p'],
             ['s_b', 'GLS standard error of the trend, m/yr'],
             ['t_0', 'Tidal-datum epoch center the detrended series is referenced to'],
             ['y_s, y_e', 'NTDE start and end years'],
             ['y_d', 'Detrended water level, m']],
 'references': ['Zervas, C. (2009), Sea Level Variations of the United States 1854-2006, '
                'NOAA Technical Report NOS CO-OPS 053 (Steps 1-8 of the trend procedure)',
                'Cochrane, D. & Orcutt, G. H. (1949), Application of least squares '
                'regression to relationships containing autocorrelated error terms',
                'NOAA CO-OPS tidal datum (NTDE) conventions']}


def compute(inp: dict) -> Result:
    """Detrend a water-level record (SI inputs). Returns the trend with its
    confidence interval and the observed, trend and detrended series for plotting."""
    _validate(inp)
    t, sec, yr, mo, y = _parse_csv(str(inp.get("csv", _SAMPLE_CSV)))
    fin = np.isfinite(y)                    # gaps (blank cells) are NaN -> excluded
    tf, yf = t[fin], y[fin]

    # The epoch center is where the trend is removed, not where it is fitted:
    # CO-OPS 053 takes the mean of the two epoch year labels (1992.0 for 1983-2001).
    epoch = (float(inp["ntde_start"]) + float(inp["ntde_end"])) / 2.0

    mode = str(inp.get("fit_mode", _FIT_COOPS))
    fit = None
    if mode == _FIT_GIVEN:
        slope = float(inp["slope_value"])
    elif mode == _FIT_OLS:
        tc = tf - float(np.mean(tf))
        denom = float(np.dot(tc, tc))
        if denom == 0.0:
            raise ValueError("record times are constant - cannot fit a trend")
        slope = float(np.dot(tc, yf) / denom)
    else:
        shifts = tuple(float(s) for s in (inp.get("level_shifts") or ()))
        mt, my, serial, month = _monthly_means(
            sec, yr, mo, y, float(inp.get("min_valid_days", 28.0)))
        fit = _fit_seasonal_trend(mt, my, serial, month, shifts)
        slope = fit["slope"]

    trend_comp = slope * (t - epoch)       # the linear component removed (zero at t0)
    detrended = y - trend_comp             # NaN propagates at gaps
    datum = float(np.mean(yf - slope * (tf - epoch)))   # mean of the finite detrended
    trend_line = datum + trend_comp        # absolute line overlaying the observations

    resid = (yf - slope * (tf - epoch)) - datum
    rms = float(np.sqrt(np.mean(resid * resid)))
    record_years = float(tf[-1] - tf[0])
    total_trend = slope * record_years
    n = int(fin.sum())

    pt, po, ptr, pdt = _decimate(t, y, trend_line, detrended)
    pdatum = np.full(len(pt), datum)
    gaps = len(t) - n
    if fit is not None:
        notes = [
            f"CO-OPS 053 monthly GLS: {slope * 1000:.2f} +/- {fit['ci'] * 1000:.2f} mm/yr "
            f"(95% CI); rho={fit['rho']:.3f} after {fit['iterations']} Cochrane-Orcutt "
            f"pass{'es' if fit['iterations'] != 1 else ''}, {fit['n_months']} complete "
            f"months, {fit['dof']} dof; seasonal range {fit['seasonal_range'] * 100:.1f} cm; "
            f"referenced to the {epoch:.1f} epoch center"
        ]
        if not fit["converged"]:
            notes.append("rho did not meet the iteration tolerance")
    else:
        notes = [
            f"{'supplied' if mode == _FIT_GIVEN else 'ordinary least-squares'} slope "
            f"{slope * 1000:.2f} mm/yr, no confidence interval; referenced to the "
            f"{epoch:.1f} epoch center"
        ]
    notes.append(f"n={n} samples over {record_years:.1f} yr"
                 + (f" ({gaps} gaps skipped)" if gaps else ""))
    if len(t) != len(pt):
        notes.append(f"plot strided to {len(pt)} of {len(t)} points (fit uses all)")

    handoff = ""
    if inp.get("handoff"):           # full-resolution detrended series for the next app
        handoff = "year,detrended\n" + "\n".join(
            f"{tt:.6f},{'' if not math.isfinite(vv) else format(vv, '.6f')}"
            for tt, vv in zip(t.tolist(), detrended.tolist()))

    return Result(
        slope_per_year=slope,
        ci_halfwidth=(fit["ci"] if fit else 0.0),
        rho_ar1=(fit["rho"] if fit else 0.0),
        seasonal_range=(fit["seasonal_range"] if fit else 0.0),
        n_months=float(fit["n_months"] if fit else 0),
        dof=float(fit["dof"] if fit else 0),
        epoch_year=epoch, total_trend=total_trend,
        record_years=record_years, n_samples=float(n), rms_residual=rms,
        converged=("yes" if (fit is None or fit["converged"]) else "no"),
        pivot_line=epoch,
        profile_year=pt, profile_original=po, profile_trend=ptr,
        profile_detrended=pdt, profile_datum=pdatum, handoff_csv=handoff,
        notes="; ".join(notes),
    )


# --- self-tests -----------------------------------------------------------------
def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}

    # Student-t quantiles against published table values (the CI depends on them).
    for df, want in ((1, 12.7062), (10, 2.228139), (30, 2.042272), (1000, 1.962339)):
        got = _t_ppf(0.975, df)
        assert abs(got - want) < 1e-4, (df, got, want)

    # The embedded sample rises at 3 mm/yr under a seasonal cycle and AR(1) noise;
    # the CO-OPS estimator recovers the rate and covers it with its own interval.
    r = compute(base)
    assert abs(r.slope_per_year - 0.003) < 5e-4, r.slope_per_year
    assert abs(r.slope_per_year - 0.003) < r.ci_halfwidth, (r.slope_per_year, r.ci_halfwidth)
    assert r.ci_halfwidth > 0.0 and r.converged == "yes"
    assert 0.2 < r.rho_ar1 < 0.6, r.rho_ar1          # injected rho = 0.4
    assert r.n_months == 360 and r.dof == 360 - 1 - 13, (r.n_months, r.dof)
    # the fitted monthly constants reproduce the injected seasonal cycle
    assert abs(r.seasonal_range - (max(_SEASONAL) - min(_SEASONAL))) < 0.03, r.seasonal_range

    # Datum invariance: a constant offset moves the level, never the slope. This is
    # what the intercept-free monthly-constant design buys.
    shifted = "date,v\n" + "\n".join(
        f"{ln.split(',')[0]},{float(ln.split(',')[1]) + 10.0:.4f}"
        for ln in _SAMPLE_CSV.splitlines()[1:])
    r_off = compute({**base, "csv": shifted})
    assert abs(r_off.slope_per_year - r.slope_per_year) < 1e-12
    assert abs(r_off.ci_halfwidth - r.ci_halfwidth) < 1e-9

    # The epoch center is the mean of the two year labels (1992.0 for 1983-2001),
    # it is where the detrended series meets the observed record, and the fit does
    # not depend on it: moving the epoch shifts the datum, never the slope.
    assert abs(r.epoch_year - 1992.0) < 1e-12, r.epoch_year
    r_ep = compute({**base, "ntde_start": 2002, "ntde_end": 2020})
    assert abs(r_ep.slope_per_year - r.slope_per_year) < 1e-12
    assert abs(r_ep.epoch_year - 2011.0) < 1e-12

    # rho is estimated over consecutive month pairs only: punching out every 7th
    # month (coprime with 12, so all calendar months survive) leaves it intact.
    rows = _SAMPLE_CSV.splitlines()
    gapped = "\n".join([rows[0]] + [ln for i, ln in enumerate(rows[1:]) if i % 7 != 0])
    r_gap = compute({**base, "csv": gapped})
    assert abs(r_gap.rho_ar1 - r.rho_ar1) < 0.08, (r_gap.rho_ar1, r.rho_ar1)

    # A level-shift dummy absorbs a step that otherwise biases the slope.
    step_at = 2006.0
    stepped = "date,v\n" + "\n".join(
        f"{ln.split(',')[0]},{float(ln.split(',')[1]) + (0.15 if _decimal_year(ln.split(',')[0]) >= step_at else 0.0):.4f}"
        for ln in rows[1:])
    r_bias = compute({**base, "csv": stepped})
    r_fix = compute({**base, "csv": stepped, "level_shifts": [step_at]})
    assert abs(r_bias.slope_per_year - 0.003) > 5 * abs(r_fix.slope_per_year - 0.003)
    assert abs(r_fix.slope_per_year - 0.003) < 5e-4, r_fix.slope_per_year

    # Alternative estimators: OLS runs without the monthly machinery, and a supplied
    # slope is used verbatim. Neither carries a confidence interval.
    r_ols = compute({**base, "fit_mode": _FIT_OLS})
    assert abs(r_ols.slope_per_year - 0.003) < 1e-3 and r_ols.ci_halfwidth == 0.0
    r_sp = compute({**base, "fit_mode": _FIT_GIVEN, "slope_value": 0.005})
    assert abs(r_sp.slope_per_year - 0.005) < 1e-12

    # Gaps (blank values) stay NaN in the plotted series so the line breaks.
    gap_csv = "date,v\n" + "\n".join(
        (f"{ln.split(',')[0]}," if i in (5, 6) else ln) for i, ln in enumerate(rows[1:]))
    rg = compute({**base, "csv": gap_csv})
    assert rg.n_samples == 358, rg.n_samples
    assert np.isnan(rg.profile_original).sum() == 2
    assert np.isfinite(rg.profile_trend).all()            # trend line spans the gap

    # Decimation caps the plotted series but not the fit.
    many = "date,v\n" + "\n".join(
        f"{1900 + h // 8766}-{1 + (h // 730) % 12:02d}-15 {h % 24:02d}:00,{0.10 + 1e-6 * h}"
        for h in range(_PLOT_MAX + 500))
    rm = compute({**base, "csv": many, "fit_mode": _FIT_OLS})
    assert rm.n_samples == _PLOT_MAX + 500
    assert len(rm.profile_year) <= _PLOT_MAX + 1
    print("  self-tests: PASS (t quantiles, CO-OPS 053 slope + CI, datum invariance, "
          "epoch independence, gap-safe rho, level shift, OLS/supplied, gaps, decimation)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  embedded sample: {int(r.n_samples)} monthly means over {r.record_years:.0f} yr")
    print(f"    trend        = {r.slope_per_year * 1000:.2f} +/- "
          f"{r.ci_halfwidth * 1000:.2f} mm/yr (95% CI)")
    print(f"    AR(1) rho    = {r.rho_ar1:.3f}   dof = {int(r.dof)}   "
          f"months = {int(r.n_months)}")
    print(f"    seasonal     = {r.seasonal_range * 100:.1f} cm peak to trough")
    print(f"    epoch center = {r.epoch_year:.1f}")
    print(f"    RMS residual = {r.rms_residual:.4f} m")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
