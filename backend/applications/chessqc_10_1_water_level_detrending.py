"""CHESS-QC application 10-1 - Water Level Detrending.

Functional area: Coastal Hazards. Removes the long-term linear sea-level trend
from a water-level record by least-squares regression of level on time, so that
the residual series reflects variability about the trend rather than the trend
itself.

Two trend references (centerings) are offered:

  - NTDE midpoint (pivot): time is centered on the midpoint of the National
    Tidal Datum Epoch (NTDE), so the fitted trend is zero at that epoch and the
    detrended level equals the observed level there. This matches NOAA tidal
    datum conventions. The epoch spans [start, end + 1) in years, e.g. the
    official 1983-2001 NTDE covers [1983, 2002) with midpoint 1992.5.
  - Record mean (no pivot): time is centered on the mean of the record, so the
    detrended series is referenced to the record-average level.

The least-squares slope is identical for both references; only the constant
offset of the detrended series (the reference level) differs. The slope may also
be supplied directly (for example a published regional or global rate) instead
of being fitted from the record.

Input is a CSV record with column 1 a date and column 2 the water level (meters):
either a bundled NOAA station or a user-supplied file. Blank water-level cells
(gaps) and a header row are ignored. The full record resolution is used for the
fit; the returned profile series are decimated only for display when the record
is large.

Classification: standard (ordinary least-squares linear trend; a standard,
well-defined statistical procedure under the linear-trend assumption).
Theory and references: Zervas (2009), "Sea Level Variations of the United States
1854-2006", NOAA Technical Report NOS CO-OPS 053; NOAA CO-OPS tidal datum
(NTDE) conventions.

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
    group: str = ""        # profile panel id; profiles sharing a group plot together


# --- application metadata --------------------------------------------------------
APP_META = AppMeta(
    aces_id="10-1",
    name="Water Level Detrending",
    area="Coastal Hazards",
    classification="standard",
    cite="Zervas (2009) NOAA CO-OPS 053; NTDE datum convention",
    default_system="SI",
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

# A small embedded annual-mean sample (default, so the app computes without a
# file): level rises at exactly 3 mm/yr about year 2005, base 0.10 m.
_SAMPLE_CSV = "date,water_level_m\n" + "\n".join(
    f"{yr}-07-01,{0.10 + 0.003 * (yr - 2005):.4f}" for yr in range(1990, 2021)
)

# Cap on the number of points returned in the display profiles (the fit always
# uses every valid sample; only the plotted series is strided when larger).
_PLOT_MAX = 4000

INPUTS = (
    Field("csv", "Water-level record", "csv", default=_SAMPLE_CSV, choices=STATIONS,
          note="Select a bundled NOAA station or upload your own CSV "
               "(column 1 = date, column 2 = water level in m). Header and blank "
               "water-level rows are ignored."),
    Field("method", "Trend reference", "choice", default="NTDE midpoint (pivot)",
          choices=("NTDE midpoint (pivot)", "Record mean (no pivot)"),
          note="Pivot the trend at the NTDE midpoint (NOAA datum convention) or "
               "center it on the record mean."),
    Field("ntde_start", "NTDE start year", "int", "yr", "yr", default=1983,
          lo=1800, hi=2100, note="National Tidal Datum Epoch start year (inclusive); "
          "used only with the NTDE midpoint pivot"),
    Field("ntde_end", "NTDE end year", "int", "yr", "yr", default=2001,
          lo=1800, hi=2100, note="National Tidal Datum Epoch end year (inclusive)"),
    Field("fit_mode", "Slope source", "choice", default="Fit (least squares)",
          choices=("Fit (least squares)", "Specified slope"),
          note="Fit the slope from the record, or apply a supplied slope."),
    Field("slope_value", "Specified slope", "float", "m/yr", "ft/yr", default=0.003,
          lo=-1.0, hi=1.0, note="used only when Slope source = Specified slope"),
)

OUTPUTS = (
    Out("slope_per_year", "Linear trend (slope)", "m/yr", "ft/yr", "scalar"),
    Out("pivot_year", "Pivot (reference) year", "yr", "yr", "scalar"),
    Out("total_trend", "Total trend over record", "m", "ft", "scalar"),
    Out("record_years", "Record length", "yr", "yr", "scalar"),
    Out("n_samples", "Samples used in fit", "", "", "scalar"),
    Out("rms_residual", "RMS residual about trend", "m", "ft", "scalar"),
    Out("profile_year", "Profile: year", "yr", "yr", "profile"),
    # Panel 1 (group "obs"): observed level with the fitted linear trend on top.
    Out("profile_original", "Profile: observed", "m", "ft", "profile", group="obs"),
    Out("profile_trend", "Profile: linear trend", "m", "ft", "profile", group="obs"),
    # Panel 2 (group "detr"): detrended level with the horizontal datum.
    Out("profile_detrended", "Profile: detrended", "m", "ft", "profile", group="detr"),
    Out("profile_datum", "Profile: datum", "m", "ft", "profile", group="detr"),
)


@dataclass
class Result:
    slope_per_year: float
    pivot_year: float
    total_trend: float
    record_years: float
    n_samples: float
    rms_residual: float
    profile_year: np.ndarray
    profile_original: np.ndarray
    profile_trend: np.ndarray
    profile_detrended: np.ndarray
    profile_datum: np.ndarray
    notes: str = ""


# --- date / CSV parsing ---------------------------------------------------------
_YEAR_BOUNDS: dict[int, tuple] = {}


def _decimal_year(s: str) -> float:
    """Convert a 'YYYY-MM-DD HH:MM' / 'YYYY-MM-DD' / 'YYYY-MM' string to a
    fractional calendar year (leap-year aware). Raises ValueError on a bad date."""
    s = s.strip()
    sep = " " if " " in s else ("T" if "T" in s else "")
    datepart, timepart = (s.split(sep, 1) if sep else (s, ""))
    p = datepart.split("-")
    y = int(p[0])
    mo = int(p[1]) if len(p) > 1 and p[1] else 1
    d = int(p[2]) if len(p) > 2 and p[2] else 1
    hh = mm = 0
    if timepart:
        tp = timepart.split(":")
        hh = int(tp[0]) if tp[0] else 0
        mm = int(tp[1]) if len(tp) > 1 and tp[1] else 0
    dt = datetime(y, mo, d, hh, mm)
    if y not in _YEAR_BOUNDS:
        st = datetime(y, 1, 1)
        _YEAR_BOUNDS[y] = (st, (datetime(y + 1, 1, 1) - st).total_seconds())
    st, span = _YEAR_BOUNDS[y]
    return y + (dt - st).total_seconds() / span


def _parse_csv(text: str) -> tuple[np.ndarray, np.ndarray]:
    """Parse CSV text -> (decimal_year, value) arrays, sorted by time. Column 1 is
    the date, column 2 the value. Rows with a valid date but a blank/non-numeric
    value are kept with value = NaN, so gaps in the record stay as gaps (the trend
    fit ignores them and the plot breaks the line across them). The header row and
    rows with an unparseable date are dropped."""
    years: list[float] = []
    vals: list[float] = []
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
            t = _decimal_year(ds)
        except (ValueError, IndexError):
            continue  # header or unparseable date
        try:
            v = float(parts[1].strip())
            if not math.isfinite(v):
                v = math.nan
        except ValueError:
            v = math.nan          # blank / non-numeric -> gap
        years.append(t)
        vals.append(v)
    t = np.asarray(years, dtype=np.float64)
    y = np.asarray(vals, dtype=np.float64)
    if int(np.isfinite(y).sum()) < 2:
        raise ValueError(
            "need at least 2 valid (date, water level) rows to fit a trend; "
            f"parsed {len(years)} rows, {int(np.isfinite(y).sum())} with values"
        )
    order = np.argsort(t, kind="stable")
    return t[order], y[order]


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
def compute(inp: dict) -> Result:
    """Detrend a water-level record (SI inputs). Returns slope and the original,
    trend, and detrended series for plotting."""
    _validate(inp)
    t, y = _parse_csv(str(inp.get("csv", _SAMPLE_CSV)))
    fin = np.isfinite(y)                    # gaps (blank cells) are NaN -> excluded from the fit
    tf, yf = t[fin], y[fin]
    tbar = float(np.mean(tf))

    method = str(inp.get("method", "NTDE midpoint (pivot)")).lower()
    use_pivot = "midpoint" in method or "ntde" in method
    if use_pivot:
        pivot = (float(inp["ntde_start"]) + float(inp["ntde_end"]) + 1.0) / 2.0
    else:
        pivot = tbar

    # Slope by ordinary least squares about the data centroid (intercept absorbed),
    # so the fitted rate is independent of the chosen reference. The pivot affects
    # only the offset of the detrended series, not the slope. Computed on finite
    # samples only (gaps excluded).
    if "specified" in str(inp.get("fit_mode", "")).lower():
        slope = float(inp["slope_value"])
    else:
        tc = tf - tbar
        denom = float(np.dot(tc, tc))
        if denom == 0.0:
            raise ValueError("record times are constant - cannot fit a trend")
        slope = float(np.dot(tc, yf) / denom)

    trend_comp = slope * (t - pivot)       # the linear component removed (zero at pivot)
    detrended = y - trend_comp             # NaN propagates at gaps
    datum = float(np.mean(yf - slope * (tf - pivot)))   # mean of the finite detrended series
    # absolute trend line that overlays the observed data (slope through the centroid)
    trend_line = datum + trend_comp

    resid = (yf - slope * (tf - pivot)) - datum
    rms = float(np.sqrt(np.mean(resid * resid)))
    record_years = float(tf[-1] - tf[0])
    total_trend = slope * record_years
    n = int(fin.sum())

    pt, po, ptr, pdt = _decimate(t, y, trend_line, detrended)
    pdatum = np.full(len(pt), datum)
    gaps = len(t) - n
    notes = [
        f"{'NTDE midpoint pivot' if use_pivot else 'record-mean centered'} "
        f"at {pivot:.2f}; slope {slope:.5f} m/yr ({slope * 1000:.2f} mm/yr); "
        f"n={n} samples over {record_years:.1f} yr"
        + (f" ({gaps} gaps skipped)" if gaps else "")
    ]
    if len(t) != len(pt):
        notes.append(f"plot strided to {len(pt)} of {len(t)} points (fit uses all)")

    return Result(
        slope_per_year=slope, pivot_year=pivot, total_trend=total_trend,
        record_years=record_years, n_samples=float(n), rms_residual=rms,
        profile_year=pt, profile_original=po, profile_trend=ptr,
        profile_detrended=pdt, profile_datum=pdatum, notes="; ".join(notes),
    )


# --- self-tests -----------------------------------------------------------------
def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}

    # The embedded sample rises at ~3 mm/yr (exactly linear in calendar year;
    # decimal-year leap spacing leaves a ~1e-8 wobble). Both references recover it.
    r_mid = compute({**base, "method": "NTDE midpoint (pivot)"})
    r_ord = compute({**base, "method": "Record mean (no pivot)"})
    assert abs(r_mid.slope_per_year - 0.003) < 1e-6, r_mid.slope_per_year
    assert abs(r_ord.slope_per_year - 0.003) < 1e-6, r_ord.slope_per_year

    # Slope is reference-independent; the two references give the same rate.
    assert abs(r_mid.slope_per_year - r_ord.slope_per_year) < 1e-12
    assert abs(r_mid.pivot_year - 1992.5) < 1e-9, r_mid.pivot_year
    mean_year = float(np.mean(r_ord.profile_year))
    assert abs(r_ord.pivot_year - mean_year) < 1e-9

    # Only the detrended offset differs, by slope*(pivot_mid - mean_year).
    expected = r_mid.slope_per_year * (r_mid.pivot_year - mean_year)
    diff = r_mid.profile_detrended - r_ord.profile_detrended
    assert np.allclose(diff, expected, atol=1e-9), (diff[0], expected)

    # The detrended series has negligible residual for a near-linear input.
    assert r_ord.rms_residual < 1e-4, r_ord.rms_residual

    # Specified-slope override is used verbatim.
    r_sp = compute({**base, "fit_mode": "Specified slope", "slope_value": 0.005})
    assert abs(r_sp.slope_per_year - 0.005) < 1e-12

    # Decimation caps the plotted series but not the fit.
    many = "date,v\n" + "\n".join(
        f"2000-01-01 {h:02d}:00,{0.0 + 1e-5 * h}" for h in range(24)
    )  # small, no decimation
    rm = compute({**base, "csv": "date,v\n" + "\n".join(
        f"{1900 + i // 12}-{i % 12 + 1:02d}-15,{0.001 * i}" for i in range(_PLOT_MAX + 500))})
    assert rm.n_samples == _PLOT_MAX + 500
    assert len(rm.profile_year) <= _PLOT_MAX + 1
    _ = many

    # Gaps (blank values) are excluded from the fit but kept as NaN in the plotted
    # series so the line breaks across them.
    gap_csv = ("date,v\n2000-01-01,0.10\n2001-01-01,\n2002-01-01,0.106\n"
               "2003-01-01,0.109\n2004-01-01,0.112")
    rg = compute({**base, "csv": gap_csv, "method": "Record mean (no pivot)"})
    assert rg.n_samples == 4, rg.n_samples                 # 5 rows, 1 blank -> 4 in the fit
    assert len(rg.profile_original) == 5                   # all timestamps kept for the plot
    assert np.isnan(rg.profile_original).sum() == 1        # the gap stays NaN
    assert np.isnan(rg.profile_detrended).sum() == 1
    assert np.isfinite(rg.profile_trend).all()            # trend line spans the gap
    print(f"  self-tests: PASS (slope recovery, pivot offset, override, decimation, gaps)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  embedded sample: {int(r.n_samples)} annual values over {r.record_years:.0f} yr")
    print(f"    slope        = {r.slope_per_year:.5f} m/yr ({r.slope_per_year * 1000:.2f} mm/yr)")
    print(f"    pivot year   = {r.pivot_year:.2f}")
    print(f"    total trend  = {r.total_trend:.4f} m over the record")
    print(f"    RMS residual = {r.rms_residual:.5f} m")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
