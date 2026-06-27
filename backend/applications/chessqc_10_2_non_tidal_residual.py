"""CHESS-QC application 10-2 - Non-Tidal Residual (NTR).

Functional area: Coastal Hazards. The non-tidal residual is the (detrended)
observed water level minus the predicted astronomical tide on the same time
grid:

    NTR(t) = WL(t) - tide(t)

It isolates the storm/meteorological surge from the deterministic tide and is the
usual input to Peaks Over Threshold (10-3). The tide is sampled onto the
water-level timestamps by linear interpolation in time (an exact no-op on
matched hourly stamps; it only bridges a missing tide hour). Water-level gaps
(NaN) are preserved so the NTR series carries the same coverage as the record.

Typical workflow: detrend the water level (10-1), then subtract the tide here.
Standalone, this app subtracts the tide from whatever water-level series is
supplied (bundled NOAA station or uploaded CSV).

Classification: standard (a direct, well-defined subtraction with time
alignment).
Theory and references: NOAA CO-OPS tide/residual practice; PyStorm NTR engine.

Self-containment: zero sibling imports; embeds its own contract dataclasses.
Runnable standalone:
    python chessqc_10_2_non_tidal_residual.py
stdlib + numpy only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np


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
    show_if: tuple = ()
    data_dir: str = "water_levels"   # bundled-station folder under data/ for a "csv" field


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
    aces_id="10-2",
    name="Non-Tidal Residual",
    area="Coastal Hazards",
    classification="standard",
    cite="NOAA CO-OPS tide/residual practice; PyStorm NTR",
    default_system="SI",
    next_apps=(("10-3", "Peaks Over Threshold"),),
)

STATIONS = (
    "8518750|The Battery, NY",
    "8651370|Duck, NC",
    "8724580|Key West, FL",
    "8761724|Grand Isle, LA",
    "8771450|Galveston Pier 21, TX",
)


def _sample(kind: str) -> str:
    """Small hourly sample: tide = 0.6 m semidiurnal; WL = tide + a slow surge."""
    base = datetime(2000, 1, 1)
    rows = ["date,value"]
    for h in range(24 * 30):
        ts = (base + timedelta(hours=h)).strftime("%Y-%m-%d %H:%M")
        tide = 0.6 * math.sin(2 * math.pi * h / 12.42)
        if kind == "tide":
            val = tide
        else:  # water level = tide + slow surge residual
            val = tide + 0.10 + 0.08 * math.sin(2 * math.pi * h / 240.0)
        rows.append(f"{ts},{val:.4f}")
    return "\n".join(rows)


_SAMPLE_WL = _sample("wl")
_SAMPLE_TIDE = _sample("tide")
_PLOT_MAX = 4000

INPUTS = (
    Field("csv_wl", "Water-level record", "csv", default=_SAMPLE_WL, choices=STATIONS,
          data_dir="water_levels",
          note="Observed (ideally detrended, via 10-1) water level: a bundled NOAA "
               "station or an uploaded CSV (column 1 = date, column 2 = level in m)."),
    Field("csv_tide", "Tide prediction", "csv", default=_SAMPLE_TIDE, choices=STATIONS,
          data_dir="tide_predictions",
          note="Predicted astronomical tide on the same grid: a bundled NOAA station "
               "or an uploaded CSV (column 1 = date, column 2 = tide in m)."),
)

OUTPUTS = (
    Out("mean_ntr", "Mean NTR", "m", "ft", "scalar",
        note="Time-mean of the non-tidal residual (WL - tide) over the overlapping finite samples, in m."),
    Out("rms_ntr", "RMS NTR", "m", "ft", "scalar",
        note="Root-mean-square of the non-tidal residual, a measure of typical surge magnitude, in m."),
    Out("max_ntr", "Maximum NTR", "m", "ft", "scalar",
        note="Peak (largest) non-tidal residual over the overlapping record, in m."),
    Out("n_samples", "Samples", "", "", "scalar",
        note="Count of overlapping finite (non-gap) NTR samples used in the statistics."),
    Out("profile_year", "Profile: year", "yr", "yr", "profile",
        note="Decimal-year time axis of the plotted series, shared by the water-level, tide, and NTR curves."),
    Out("profile_wl", "Profile: water level", "m", "ft", "profile", group="obs",
        note="Observed (ideally detrended) water-level series at each timestamp, in m."),
    Out("profile_tide", "Profile: tide", "m", "ft", "profile", group="obs",
        note="Predicted astronomical tide linearly interpolated onto the water-level timestamps, in m."),
    Out("profile_ntr", "Profile: NTR", "m", "ft", "profile", group="ntr",
        note="Non-tidal residual series, NTR(t) = WL(t) - tide(t), the isolated storm/meteorological surge, in m."),
    # full-resolution NTR series for the workflow hand-off into 10-3 (emitted only
    # when the `handoff` input is set). Not shown or plotted.
    Out("handoff_csv", "handoff", "", "", "data",
        note="Full-resolution year,NTR CSV passed in-memory to 10-3 (Peaks Over Threshold); not displayed."),
)


@dataclass
class Result:
    mean_ntr: float
    rms_ntr: float
    max_ntr: float
    n_samples: float
    profile_year: np.ndarray
    profile_wl: np.ndarray
    profile_tide: np.ndarray
    profile_ntr: np.ndarray
    handoff_csv: str = ""
    notes: str = ""


# --- date / CSV parsing (shared style with 10-1 / 10-3) -------------------------
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


def _parse_csv(text: str, keep_nan: bool) -> tuple[np.ndarray, np.ndarray]:
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
        if not keep_nan and math.isnan(v):
            continue
        years.append(t); vals.append(v)
    t = np.asarray(years, dtype=np.float64)
    y = np.asarray(vals, dtype=np.float64)
    order = np.argsort(t, kind="stable")
    return t[order], y[order]


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
ABOUT = {'summary': 'Computes the non-tidal residual (storm/meteorological surge) by subtracting '
            'the predicted astronomical tide from the observed water level on a common '
            'time grid, and reports its mean, RMS, and maximum. The tide is linearly '
            'interpolated onto the water-level timestamps and water-level gaps are '
            'preserved.',
 'methods': [{'name': 'Non-tidal residual by tide subtraction',
              'when': None,
              'tag': 'standard',
              'note': None,
              'equations': [{'tex': '\\mathrm{NTR}(t) = \\mathrm{WL}(t) - '
                                    '\\mathrm{tide}(t)',
                             'desc': 'Non-tidal residual: observed (detrended) water level '
                                     'minus the predicted astronomical tide on the same '
                                     'time grid.'},
                            {'tex': '\\mathrm{tide}(t_i) = \\mathrm{tide}(t_a) + '
                                    '\\frac{t_i - t_a}{t_b - t_a}\\,[\\mathrm{tide}(t_b) - '
                                    '\\mathrm{tide}(t_a)]',
                             'desc': 'Linear interpolation of the tide onto each '
                                     'water-level timestamp t_i between bracketing tide '
                                     'stamps t_a, t_b (a no-op on matched hourly stamps).'},
                            {'tex': '\\bar{\\eta} = \\frac{1}{N}\\sum_{i=1}^{N} \\eta_i',
                             'desc': 'Mean NTR over the N overlapping (finite) samples '
                                     'eta_i = NTR(t_i).'},
                            {'tex': '\\eta_{\\mathrm{rms}} = '
                                    '\\sqrt{\\frac{1}{N}\\sum_{i=1}^{N} \\eta_i^{2}}',
                             'desc': 'Root-mean-square NTR, a measure of typical surge '
                                     'magnitude.'},
                            {'tex': '\\eta_{\\max} = \\max_{i} \\, \\eta_i',
                             'desc': 'Maximum NTR (peak surge) over the overlapping '
                                     'record.'}]}],
 'symbols': [['NTR', 'Non-tidal residual (storm/meteorological surge), in metres'],
             ['WL(t)', 'Observed, ideally detrended water level at time t'],
             ['tide(t)', 'Predicted astronomical tide at time t'],
             ['t_i', 'Water-level timestamp onto which the tide is interpolated'],
             ['eta_i', 'NTR sample value at timestamp t_i'],
             ['N', 'Number of overlapping finite (non-gap) samples'],
             ['eta_rms', 'Root-mean-square NTR'],
             ['eta_max', 'Maximum (peak) NTR']],
 'references': ['NOAA CO-OPS tide/residual practice', 'PyStorm NTR engine']}


def compute(inp: dict) -> Result:
    """NTR = water level - tide, with the tide interpolated onto the WL grid."""
    t_wl, v_wl = _parse_csv(str(inp.get("csv_wl", _SAMPLE_WL)), keep_nan=True)
    t_td, v_td = _parse_csv(str(inp.get("csv_tide", _SAMPLE_TIDE)), keep_nan=False)
    if t_wl.size < 2:
        raise ValueError("water-level record needs at least 2 valid rows")
    if t_td.size < 2:
        raise ValueError("tide record needs at least 2 valid rows")

    # linear interpolation of tide onto the WL timestamps; outside the tide span -> NaN
    tide_at_wl = np.interp(t_wl, t_td, v_td, left=np.nan, right=np.nan)
    ntr = v_wl - tide_at_wl

    fin = np.isfinite(ntr)
    if int(fin.sum()) < 1:
        raise ValueError("no overlapping water-level and tide samples")
    nf = ntr[fin]
    mean_ntr = float(np.mean(nf))
    rms_ntr = float(np.sqrt(np.mean(nf * nf)))
    max_ntr = float(np.max(nf))

    pyear, pwl, ptide, pntr = _decimate(t_wl, v_wl, tide_at_wl, ntr)
    notes = (f"NTR = WL - tide; {int(fin.sum())} overlapping samples; "
             f"mean {mean_ntr:.3f} m, RMS {rms_ntr:.3f} m, max {max_ntr:.3f} m")
    handoff = ""
    if inp.get("handoff"):           # full-resolution NTR series for the next app
        handoff = "year,ntr\n" + "\n".join(
            f"{tt:.6f},{'' if not math.isfinite(vv) else format(vv, '.6f')}"
            for tt, vv in zip(t_wl.tolist(), ntr.tolist()))
    return Result(
        mean_ntr=mean_ntr, rms_ntr=rms_ntr, max_ntr=max_ntr, n_samples=float(fin.sum()),
        profile_year=pyear, profile_wl=pwl, profile_tide=ptide, profile_ntr=pntr,
        handoff_csv=handoff, notes=notes,
    )


# --- self-tests -----------------------------------------------------------------
def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}
    r = compute(base)
    # the embedded surge residual is 0.10 + 0.08*sin(...), so mean NTR ~ 0.10 m
    assert abs(r.mean_ntr - 0.10) < 0.02, r.mean_ntr
    assert r.n_samples == 24 * 30, r.n_samples
    assert r.profile_ntr.shape == r.profile_year.shape
    # NTR = WL - tide pointwise (on the sample, no gaps)
    assert np.allclose(r.profile_ntr, r.profile_wl - r.profile_tide, equal_nan=True)
    # a constant tide offset shifts NTR by the opposite amount
    flat = "date,value\n" + "\n".join(
        f"2000-01-{d + 1:02d},1.0" for d in range(5))
    wl = "date,value\n" + "\n".join(
        f"2000-01-{d + 1:02d},{1.5 + 0.1 * d}" for d in range(5))
    r2 = compute({"csv_wl": wl, "csv_tide": flat})
    assert np.allclose(r2.profile_ntr, np.array([0.5, 0.6, 0.7, 0.8, 0.9]), atol=1e-9), r2.profile_ntr
    print("  self-tests: PASS (NTR subtraction, alignment, mean surge)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"    samples   = {int(r.n_samples)}")
    print(f"    mean NTR  = {r.mean_ntr:.3f} m")
    print(f"    RMS NTR   = {r.rms_ntr:.3f} m")
    print(f"    max NTR   = {r.max_ntr:.3f} m")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
