"""CHESS-QC application 10-3 - Peaks Over Threshold (POT).

Functional area: Coastal Hazards. Extracts independent storm peaks from a
continuous water-level / non-tidal-residual series: an automatically chosen
percentile threshold is raised until the declustered exceedance rate matches a
target average number of events per year, then the exceedances are declustered
into one peak per independent storm and rank-trimmed to exactly
round(target x effective-duration) peaks. These peaks are the input to the
Probabilistic Simulation Technique (10-4).

Method (after the PyStorm peaks_over_threshold module):
  - Effective duration = (non-NaN samples) x (median sampling interval) in years;
    gaps do not count toward the rate. For the hourly records this application is
    normally driven with, that is the PyStorm form (non-NaN samples) / (365.25 x
    24), and it stays correct for records held at any other resolution.
  - Iterative threshold search: the threshold is raised from a start percentile
    in small percentile steps; at each level the exceedances are declustered and
    the event rate is measured against the effective duration. The highest
    threshold whose rate is still >= target is kept (converged when that rate
    lands in [target, target + tolerance]).
  - Declustering: "hydrograph" groups consecutive exceedances separated by more
    than the inter-event window and keeps each group's maximum; "peak_gap" drops
    a sample that lies within the window of, and is not larger than, the
    preceding exceedance.
  - Rank-trim the converged peaks to the largest round(target x eff-duration),
    so the retained count is deterministic.

Input is a CSV record (column 1 a date, column 2 the value): either a bundled
NOAA station or a user-supplied file (the detrended water level or non-tidal
residual from 10-1 / 10-2 in the typical workflow). Blank values are gaps.

Classification: standard (a standard declustering / threshold-selection
procedure; deterministic given the parameters).
Theory and references: Coles (2001), extreme value POT; USACE coastal-hazards
practice. Pure-numpy port of the PyStorm POT kernel.

Self-containment: zero sibling imports; embeds its own contract dataclasses.
Runnable standalone:
    python chessqc_10_3_peaks_over_threshold.py
stdlib + numpy only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np

_HOURS_PER_YEAR = 365.25 * 24.0
_SECONDS_PER_YEAR = _HOURS_PER_YEAR * 3600.0


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
    data_dir: str = ""       # bundled-record folder under data/ for a `csv` field
    show_if: tuple = ()


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    group: str = ""
    x_key: str = ""        # for kind "scatter": the output key giving the marker x-values
    note: str = ""           # hover definition shown on the output label


APP_META = AppMeta(
    aces_id="10-3",
    name="Peaks Over Threshold",
    area="Coastal Hazards",
    classification="standard",
    cite="Coles (2001); USACE coastal-hazards practice; PyStorm POT",
    default_system="SI",
    next_apps=(("10-4", "Probabilistic Simulation Technique"),),
)

STATIONS = (
    "8518750|The Battery, NY",
    "8651370|Duck, NC",
    "8735180|Dauphin Island, AL",
    "8761724|Grand Isle, LA",
    "8771450|Galveston Pier 21, TX",
    "8774770|Rockport, TX",
    "9755371|San Juan, PR",
)

# Small embedded sample: 6 years of daily residual-like values built from an annual
# cycle plus four incommensurate short-period components, which produces a dense,
# quasi-random population of independent storm peaks - enough for the threshold scan
# to converge on the default 10 events/yr without a data file.
_SAMPLE_START = datetime(1990, 1, 1)


def _sample_value(d: int) -> float:
    return (0.20
            + 0.12 * math.sin(2.0 * math.pi * d / 365.25 - 1.2)    # seasonal
            + 0.10 * math.sin(2.0 * math.pi * d / 29.5 + 0.7)
            + 0.08 * math.sin(2.0 * math.pi * d / 11.7 + 1.9)
            + 0.06 * math.sin(2.0 * math.pi * d / 7.3 + 2.8)
            + 0.05 * math.sin(2.0 * math.pi * d / 3.1 + 0.4))


_SAMPLE_CSV = "date,value\n" + "\n".join(
    f"{(_SAMPLE_START + timedelta(days=d)).strftime('%Y-%m-%d')},{_sample_value(d):.4f}"
    for d in range(6 * 365)
)

_PLOT_MAX = 4000

INPUTS = (
    Field("csv", "Water-level / NTR record", "csv", default=_SAMPLE_CSV, choices=STATIONS,
          data_dir="detrended_levels",
          note="Select a bundled NOAA station (its detrended hourly water level, the "
               "series 10-1 produces) or upload your own record: Parquet written by "
               "the CHESS/PyStorm engines, or a CSV with column 1 = date and column "
               "2 = value. Typically the detrended level (10-1) or the non-tidal "
               "residual (10-2)."),
    Field("target_events_per_year", "Target events per year", "float", "1/yr", "1/yr",
          default=10.0, lo=0.1, hi=365.0,
          note="average number of independent peaks per year to retain (matches PST)"),
    # held in SI seconds like every other duration in the suite (the "hr" unit is the
    # display unit; the front-ends divide by 3600 to show it)
    Field("interevent_hours", "Inter-event window", "float", "hr", "hr",
          default=48.0 * 3600.0, lo=1.0 * 3600.0, hi=2160.0 * 3600.0,
          note="minimum separation between independent events"),
    Field("method", "Declustering method", "choice", default="hydrograph",
          choices=("hydrograph", "peak_gap"),
          note="hydrograph: group + per-group max; peak_gap: sequential gap filter"),
    Field("start_percentile", "Start percentile", "float", "%", "%", default=75.0,
          lo=0.0, hi=99.9, note="series percentile where the upward threshold scan begins"),
    Field("step_size", "Percentile step", "float", "%", "%", default=0.05,
          lo=0.001, hi=5.0, note="percentile increment per scan step (smaller = finer, slower)"),
    Field("tolerance", "Convergence tolerance", "float", "1/yr", "1/yr", default=0.25,
          lo=0.001, hi=10.0, note="accept a rate in [target, target + tolerance]"),
)

OUTPUTS = (
    Out("threshold", "Threshold", "m", "ft", "scalar",
        note="auto-selected exceedance threshold u; the series level above which samples are counted as exceedances."),
    Out("n_peaks", "Peaks retained", "", "", "scalar",
        note="number of independent declustered storm peaks kept after rank-trimming to round(target x effective duration)."),
    Out("events_per_year", "Effective rate", "1/yr", "1/yr", "scalar",
        note="achieved declustered exceedance rate = retained peaks divided by the effective record duration (events per year)."),
    Out("final_percentile", "Threshold percentile", "%", "%", "scalar",
        note="series percentile p at which the converged threshold sits; higher percentile means a more extreme threshold."),
    Out("eff_duration", "Effective duration", "yr", "yr", "scalar",
        note="effective record length in years from the count of valid (non-NaN) samples; gaps do not count toward the rate."),
    Out("converged", "Converged", "", "", "scalar",
        note="yes if the achieved rate landed within [target, target + tolerance], otherwise no."),
    Out("profile_year", "Profile: year", "yr", "yr", "profile",
        note="decimal-year time axis for the plotted water-level / NTR series (decimated for display only)."),
    Out("profile_series", "Profile: series", "m", "ft", "profile", group="ts",
        note="the input water-level / non-tidal-residual series plotted against time."),
    Out("profile_threshold", "Profile: threshold", "m", "ft", "profile", group="ts",
        note="the selected threshold u drawn as a constant horizontal level across the series panel."),
    # peaks as scatter markers (own x), drawn on the series panel
    Out("peaks_t", "Peak times", "yr", "yr", "scatter_x",
        note="decimal-year times of the extracted independent storm peaks (marker x-values)."),
    Out("peaks_v", "Peaks", "m", "ft", "scatter", group="ts", x_key="peaks_t",
        note="values of the declustered storm peaks (one maximum per independent event), shown as markers on the series."),
    # the declustered peaks as a year,value series for the hand-off into 10-4 PST
    Out("handoff_csv", "handoff", "", "", "data",
        note="declustered peaks as a year,value CSV passed in memory to the 10-4 PST application."),
)


@dataclass
class Result:
    threshold: float
    n_peaks: float
    events_per_year: float
    final_percentile: float
    eff_duration: float
    converged: str
    profile_year: np.ndarray
    profile_series: np.ndarray
    profile_threshold: np.ndarray
    peaks_t: np.ndarray
    peaks_v: np.ndarray
    handoff_csv: str = ""
    notes: str = ""


# --- date / CSV parsing (shared style with 10-1) --------------------------------
_YEAR_BOUNDS: dict[int, tuple] = {}


def _decimal_year(s: str) -> float:
    s = s.strip()
    sep = " " if " " in s else ("T" if "T" in s else "")
    datepart, timepart = (s.split(sep, 1) if sep else (s, ""))
    if "-" not in datepart:          # bare (decimal) calendar year, e.g. a hand-off "year"
        return float(datepart)
    p = datepart.split("-")
    y = int(p[0]); mo = int(p[1]) if len(p) > 1 and p[1] else 1
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
    """CSV text -> (decimal_year, value) sorted by time; blank/non-numeric values
    are kept as NaN gaps. Header / unparseable-date rows are dropped."""
    years: list[float] = []
    vals: list[float] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(",")
        if len(parts) < 2 or not parts[0].strip():
            continue
        try:
            t = _decimal_year(parts[0].strip())
        except (ValueError, IndexError):
            continue
        try:
            v = float(parts[1].strip())
            if not math.isfinite(v):
                v = math.nan
        except ValueError:
            v = math.nan
        years.append(t); vals.append(v)
    t = np.asarray(years, dtype=np.float64)
    y = np.asarray(vals, dtype=np.float64)
    if int(np.isfinite(y).sum()) < 2:
        raise ValueError(f"need at least 2 valid (date, value) rows; got {int(np.isfinite(y).sum())}")
    order = np.argsort(t, kind="stable")
    return t[order], y[order]


# --- declustering segmenters (pure-numpy; ported from PyStorm) -------------------
def _segment_hydrograph(values, times, exceed_idx, interevent):
    if exceed_idx.size == 0:
        return np.empty(0, dtype=np.int64)
    t_exc = times[exceed_idx]
    breaks = np.where(np.diff(t_exc) > interevent)[0] + 1
    starts = np.concatenate(([0], breaks))
    ends = np.concatenate((breaks, [exceed_idx.size]))
    out = np.empty(starts.size, dtype=np.int64)
    for g, (s, e) in enumerate(zip(starts, ends)):
        block = exceed_idx[s:e]
        out[g] = int(block[np.argmax(values[block])])
    return out


def _segment_peak_gap(values, times, exceed_idx, interevent):
    if exceed_idx.size == 0:
        return np.empty(0, dtype=np.int64)
    keep = np.ones(exceed_idx.size, dtype=bool)
    for k in range(1, exceed_idx.size):
        dt = times[exceed_idx[k]] - times[exceed_idx[k - 1]]
        if dt < interevent and values[exceed_idx[k]] <= values[exceed_idx[k - 1]]:
            keep[k] = False
    return exceed_idx[keep]


def _validate(inp: dict) -> None:
    for fld in INPUTS:
        if fld.kind not in ("float", "int"):
            continue
        v = float(inp.get(fld.key, fld.default))
        if not (fld.lo <= v <= fld.hi):
            raise ValueError(f"{fld.label} ({fld.key}) = {v} outside [{fld.lo}, {fld.hi}]")


def _decimate(*arrays, nmax=_PLOT_MAX):
    n = len(arrays[0])
    if n <= nmax:
        return arrays
    stride = int(math.ceil(n / nmax))
    idx = np.arange(0, n, stride)
    if idx[-1] != n - 1:
        idx = np.append(idx, n - 1)
    return tuple(a[idx] for a in arrays)


# --- compute --------------------------------------------------------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Extracts independent storm peaks from a continuous water-level or '
            'non-tidal-residual series by raising a percentile threshold until the '
            'declustered exceedance rate matches a target number of events per year, then '
            'declustering and rank-trimming to a deterministic peak count for input to the '
            '10-4 PST.',
 'method_key': 'method',
 'methods': [{'name': 'Hydrograph declustering',
              'when': 'hydrograph',
              'tag': 'preferred',
              'note': 'Default: groups consecutive exceedances into one storm hydrograph '
                      "and keeps that group's single maximum.",
              'equations': [{'tex': 'D_{eff} = N \\, \\Delta t',
                             'desc': 'Effective duration in years: the count of non-NaN '
                                     'samples times the median sampling interval (gaps '
                                     'do not count toward the rate).'},
                            {'tex': 'u = y_{(k)}, \\quad k = \\lfloor (1 - p/100)(N - 1) '
                                    '\\rfloor',
                             'desc': 'Threshold u is the value at percentile p of the '
                                     'series (descending order index k).'},
                            {'tex': '\\lambda = \\frac{N_p}{D_{eff}}',
                             'desc': 'Declustered exceedance rate: declustered peak count '
                                     'over effective duration.'},
                            {'tex': '\\mu \\leq \\lambda \\leq \\mu + \\epsilon',
                             'desc': 'Convergence: the highest threshold whose rate is '
                                     'still at least the target mu, within tolerance '
                                     'epsilon.'},
                            {'tex': 't_{i+1} - t_i > \\tau',
                             'desc': 'A new event starts where consecutive exceedances are '
                                     'separated by more than the inter-event window tau; '
                                     'the group maximum is the peak.'},
                            {'tex': 'N_{keep} = \\mathrm{round}(\\mu \\cdot D_{eff})',
                             'desc': 'Rank-trim the converged peaks to the largest '
                                     'round(target x effective-duration) for a '
                                     'deterministic count.'}]},
             {'name': 'Peak-gap declustering',
              'when': 'peak_gap',
              'tag': '',
              'note': 'Sequential gap filter: drops a later sample that lies within the '
                      'window of, and is not larger than, the preceding exceedance.',
              'equations': [{'tex': 'D_{eff} = N \\, \\Delta t',
                             'desc': 'Effective duration in years: the count of non-NaN '
                                     'samples times the median sampling interval (gaps '
                                     'do not count toward the rate).'},
                            {'tex': 'u = y_{(k)}, \\quad k = \\lfloor (1 - p/100)(N - 1) '
                                    '\\rfloor',
                             'desc': 'Threshold u is the value at percentile p of the '
                                     'series (descending order index k).'},
                            {'tex': '\\lambda = \\frac{N_p}{D_{eff}}',
                             'desc': 'Declustered exceedance rate: declustered peak count '
                                     'over effective duration.'},
                            {'tex': '\\mu \\leq \\lambda \\leq \\mu + \\epsilon',
                             'desc': 'Convergence: the highest threshold whose rate is '
                                     'still at least the target mu, within tolerance '
                                     'epsilon.'},
                            {'tex': 't_k - t_{k-1} < \\tau, \\, y_k \\leq y_{k-1}',
                             'desc': 'Drop exceedance k when it falls within the window '
                                     'tau of, and is no larger than, the preceding '
                                     'exceedance.'},
                            {'tex': 'N_{keep} = \\mathrm{round}(\\mu \\cdot D_{eff})',
                             'desc': 'Rank-trim the converged peaks to the largest '
                                     'round(target x effective-duration) for a '
                                     'deterministic count.'}]}],
 'symbols': [['u', 'selected exceedance threshold (m)'],
             ['p', 'series percentile defining the current threshold (%)'],
             ['N', 'number of valid (non-NaN) samples in the series'],
             ['Delta t', 'median sampling interval of the record, in years'],
             ['D_{eff}', 'effective record duration in years (valid samples x sampling interval)'],
             ['N_p', 'number of declustered peaks retained'],
             ['lambda', 'declustered exceedance rate (events per year)'],
             ['mu', 'target average events per year'],
             ['epsilon', 'convergence tolerance on the rate (1/yr)'],
             ['tau', 'inter-event window in years = window seconds / (365.25 x 24 x 3600)'],
             ['N_{keep}',
              'deterministic retained peak count = round(target x effective duration)']],
 'references': ['Coles (2001), An Introduction to Statistical Modeling of Extreme Values '
                '(POT)',
                'USACE coastal-hazards practice',
                'PyStorm peaks_over_threshold module']}


def compute(inp: dict) -> Result:
    """Extract declustered peaks over an auto-selected threshold (SI inputs)."""
    _validate(inp)
    t_all, y_all = _parse_csv(str(inp.get("csv", _SAMPLE_CSV)))
    fin = np.isfinite(y_all)
    t = t_all[fin]; y = y_all[fin]
    n = t.size

    target = float(inp["target_events_per_year"])
    interevent = float(inp["interevent_hours"]) / _SECONDS_PER_YEAR  # seconds -> years
    method = str(inp.get("method", "hydrograph"))
    start_pct = float(inp["start_percentile"])
    step = float(inp["step_size"])
    tol = float(inp["tolerance"])
    segmenter = _segment_hydrograph if method == "hydrograph" else _segment_peak_gap

    # Effective duration: the valid-sample count times the record's own sampling
    # interval, so gaps do not count toward the rate. The interval is the median
    # timestamp spacing, which for the hourly records this application is normally
    # driven with reduces to the count-based form n / (365.25 x 24).
    dt_year = float(np.median(np.diff(t))) if n > 1 else 0.0
    if not (dt_year > 0.0):
        dt_year = 1.0 / _HOURS_PER_YEAR
    eff_dur = n * dt_year                          # effective duration (yr)
    sorted_desc = np.sort(y)[::-1]
    max_iter = int((100.0 - start_pct) / step) + 1

    # one-sided scan: keep the highest-threshold state whose rate is still >= target
    ge = None; first = None
    pct = start_pct
    for _ in range(max_iter):
        if pct >= 100.0:
            break
        frac = 1.0 - pct / 100.0
        k = max(0, min(int(np.floor(frac * (n - 1))), n - 1))
        threshold = float(sorted_desc[k])
        exceed = np.flatnonzero(y > threshold)
        if exceed.size:
            peak_idx = segmenter(y, t, exceed, interevent)
            if peak_idx.size:
                rate = peak_idx.size / eff_dur
                if first is None:
                    first = (threshold, peak_idx, rate, pct)
                if rate >= target:
                    ge = (threshold, peak_idx, rate, pct)
        pct += step

    why = ""
    if ge is not None:
        threshold, peak_idx, rate, fpct = ge
        converged = rate <= target + tol
        if not converged:
            why = (f"the rate at the highest qualifying threshold ({rate:.2f}/yr) "
                   f"overshoots the tolerance: lower the percentile step or widen it")
    elif first is not None:
        # the target rate is out of reach for this record: fall back to the lowest
        # scanned threshold, which yields the most events (the closest achievable
        # rate), rather than the most extreme threshold, which yields almost none
        threshold, peak_idx, rate, fpct = first
        converged = False
        why = (f"{target:.3g}/yr is out of reach for this record: even at the p"
               f"{start_pct:g} start threshold it holds only {rate:.2f} independent "
               f"events/yr over {eff_dur:.2f} yr")
    else:
        # no exceedances at any scanned percentile (e.g. a flat / near-constant
        # series); report zero peaks rather than failing
        threshold = float(sorted_desc[0])
        peak_idx = np.empty(0, dtype=np.int64)
        rate, fpct, converged = 0.0, start_pct, False
        why = "no sample exceeds any scanned threshold (a flat or near-constant series)"

    # rank-trim to exactly round(target * eff_dur) largest peaks (deterministic)
    n_keep = int(round(target * eff_dur))
    if n_keep >= 1 and peak_idx.size > n_keep:
        top = np.argsort(y[peak_idx], kind="stable")[::-1][:n_keep]
        peak_idx = np.sort(peak_idx[top])

    peaks_t = t[peak_idx]
    peaks_v = y[peak_idx]
    eff_rate = peak_idx.size / eff_dur

    pyear, pser = _decimate(t_all, y_all)
    pthr = np.full(len(pyear), threshold)
    notes = (f"{method} declustering, inter-event "
             f"{float(inp['interevent_hours']) / 3600.0:.0f} h; "
             f"threshold {threshold:.3f} at p{fpct:.2f}; {peak_idx.size} peaks over "
             f"{eff_dur:.1f} yr ({eff_rate:.2f}/yr); "
             f"{'converged' if converged else 'did not converge'}"
             + ("" if converged or not why else f" - {why}"))

    handoff = ""
    if inp.get("handoff"):           # the peaks (year,value) for 10-4 PST
        handoff = "year,peak\n" + "\n".join(
            f"{tt:.6f},{vv:.6f}" for tt, vv in zip(peaks_t.tolist(), peaks_v.tolist()))

    return Result(
        threshold=threshold, n_peaks=float(peak_idx.size), events_per_year=eff_rate,
        final_percentile=fpct, eff_duration=eff_dur,
        converged=("yes" if converged else "no"),
        profile_year=pyear, profile_series=pser, profile_threshold=pthr,
        peaks_t=peaks_t, peaks_v=peaks_v, handoff_csv=handoff, notes=notes,
    )


# --- self-tests -----------------------------------------------------------------
def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}
    r = compute(base)
    assert r.n_peaks >= 1, r.n_peaks
    assert np.isfinite(r.threshold)
    # peaks must exceed the threshold and be declustered (spacing >= inter-event)
    assert np.all(r.peaks_v > r.threshold - 1e-9), "peaks below threshold"
    assert r.peaks_t.size == int(r.n_peaks)
    assert np.allclose(r.profile_threshold, r.threshold)
    # rank-trim never exceeds round(target * eff_dur)
    assert r.n_peaks <= round(base["target_events_per_year"] * r.eff_duration) + 0
    # a higher target retains at least as many peaks as a lower one
    r_lo = compute({**base, "target_events_per_year": 2.0})
    r_hi = compute({**base, "target_events_per_year": 20.0})
    assert r_hi.n_peaks >= r_lo.n_peaks, (r_lo.n_peaks, r_hi.n_peaks)
    # peak_gap method also runs and yields declustered peaks
    r_pg = compute({**base, "method": "peak_gap"})
    assert r_pg.n_peaks >= 1

    # the default (daily) sample spans its full 6 years and converges on the target
    assert abs(r.eff_duration - 6.0) < 0.02, r.eff_duration
    assert abs(r.events_per_year - base["target_events_per_year"]) < 0.02
    assert r.converged == "yes", r.notes

    # effective duration follows the record's own sampling interval: the same
    # signal sampled hourly and daily reports the same duration in years
    hourly = "date,v\n" + "\n".join(
        f"{(datetime(2000, 1, 1) + timedelta(hours=h)).strftime('%Y-%m-%d %H:%M')},"
        f"{math.sin(h / 7.0):.4f}" for h in range(24 * 200))
    daily = "date,v\n" + "\n".join(
        f"{(datetime(2000, 1, 1) + timedelta(days=d)).strftime('%Y-%m-%d')},"
        f"{math.sin(d / 3.0):.4f}" for d in range(200))
    r_h = compute({**base, "csv": hourly})
    r_d = compute({**base, "csv": daily})
    assert abs(r_h.eff_duration - 200 / 365.25) < 0.01, r_h.eff_duration
    assert abs(r_d.eff_duration - 200 / 365.25) < 0.01, r_d.eff_duration

    # an unreachable target falls back to the lowest scanned threshold (the most
    # events available), not to the most extreme one (which would yield almost none)
    r_un = compute({**base, "target_events_per_year": 300.0})
    assert r_un.converged == "no"
    assert r_un.final_percentile == base["start_percentile"], r_un.final_percentile
    assert r_un.n_peaks > compute({**base, "target_events_per_year": 2.0}).n_peaks
    print(f"  self-tests: PASS (threshold search, declustering, rank-trim, methods, "
          f"sampling interval, unreachable-target fallback)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"    threshold        = {r.threshold:.3f} m at p{r.final_percentile:.2f}")
    print(f"    peaks retained   = {int(r.n_peaks)}")
    print(f"    effective rate   = {r.events_per_year:.2f} /yr over {r.eff_duration:.2f} yr")
    print(f"    converged        = {r.converged}")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
