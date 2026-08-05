"""CHESS-QC application 6-2 — Time-Dependent Beach and Dune Erosion (XSHORE).

Originating ACES application: 6-2 "Numerical Simulation of Time-Dependent Beach and Dune
Erosion" (functional area: Littoral Processes). A transcription of the XSHORE explicit
finite-difference scheme ACES runs (`KDMAIN.FOR` subroutine `KD1`, with `KD2`-`KD5B`),
after Kriebel (1984b, EBEACH) and Kriebel & Dean (1985).

Classification: exact. Reproduces the ACES source, recompiled and run on the decks it
ships with (tests/aces_oracle/fortran, target `kd`).

A note on what this replaced. Until this rework CHESS-QC 6-2 carried the Kriebel & Dean
analytical equilibrium-response model instead, on the stated grounds that the XSHORE
scheme "live[s] only in its source / the Kriebel 1984b EBEACH theory manual, neither
available". The source is available and compiles, so 6-2 now follows it. The analytical
model is retained, unregistered, as `_kriebel_dean_1985.py`: it is sound and answers a
different question (equilibrium recession for a design surge), and nothing was gained by
deleting it.

Theory. The profile is carried as bed elevation z(x) on a fixed 4 ft grid, x seaward
from a landward baseline, z relative to the still-water datum. Within the surf zone the
scheme drives the profile toward the equilibrium energy dissipation per unit volume,

    D(x)   = K_diss * ( |z_i|^2.5 - |z_i+1|^2.5 ) / ( z_bar * dx )      dissipation
    D_eq   = (5 rho gamma^2 g^1.5 / 24) * A^1.5                         equilibrium
    q(x)   = K_trans * ( D(x) - D_eq )                                  transport
    dz/dt  = -(1/dx) * dq/dx                                            continuity

with gamma = 0.78 the breaker index, K_trans = 0.0011445 ft^4/lb the transport
coefficient, K_diss = rho g^1.5 gamma^2 / 8, and A the equilibrium profile shape factor
from the mean grain size (Moore 1982, in ACES's own two-branch fit).

Each hour the breaking depth is found from the offshore wave by Snell's law with
straight and parallel contours (`KD3`, using ACES's own linear wave theory), which sets
the seaward boundary; the landward boundary is where the profile is tangent to the
equilibrium form, z_tan = -A (2A / 3 m_beach)^2. Transport is smoothed with a
five-point stencil, extended linearly across the beach face to close continuity, and
integrated with a sub-hourly step

    dt = 0.5 dx^2 / ( alpha sqrt(-z_break) ) ,  alpha = (5/16) K_trans rho g^1.5 gamma^2

halved on the fly wherever the local bed step exceeds 0.5 ft, as the source does.

Units: the application contract is SI like every other CHESS-QC application - inputs
arrive in SI and results are returned in SI. XSHORE is stated and validated in US
customary units and its grid interval is 4 ft exactly, so compute() converts to feet on
entry and back on exit. ACES's own constants are used inside: g = 32.17 ft/s^2 and
rho = 1.989 slug/ft^3.

The 0, 5, 10 and 15 ft contour changes are reported at those elevations whatever the
display system, because that is what ACES does: `KD1` converts the horizontal distances
by the unit factor but leaves the contour levels themselves fixed.

Self-containment: zero sibling imports; embeds the contract dataclasses, the generic
profile builder, the natural cubic spline, ACES's linear wave theory and the marcher.
numpy + stdlib only. Runnable: python chessqc_6_2_dune_erosion.py

Validation, against the recompiled source on the decks ACES ships:

  - `XSHORE2.IN`, generic profile, one wave record, no water-level forcing;
  - `XSHORE3.IN`, generic profile with a six-point surge on a 4 h interval and four
    wave records on a 5 h one, which ACES steps onto the shorter interval first;
  - `XSHORE4.IN`, a 53-point surveyed profile, resampled by cubic spline;
  - `XSHORE6.IN`, a 200-point surveyed profile.

The tabulated volume change and all four contour changes agree on every one of them to
better than half of the source's printed digit (0.005). The first three are asserted in
the self-tests, which carry the source's tables; the fourth is checked by the harness in
tests/aces_oracle.

`XSHORE1.IN` is the one deck not reproduced. It drives its water level from a
37-constituent tidal synthesis, which this application does not carry - the water level
is a supplied series here - so that deck is out of scope rather than wrong. It is
recorded in FINDINGS.md.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# ACES's own constants, English units (KDMAIN.FOR:420-424 and ACEDRLIB.FOR:910).
G_US = 32.17          # ft/s^2
RHO_US = 1.989        # slug/ft^3, salt water
GAMMA = 0.78          # breaker index, wave height over depth
KTRANS = 0.0011445    # ft^4/lb, transport coefficient
DEGRAD = 57.29578     # degrees per radian, as ACES defines it
_DX = 4               # ft, the fixed grid interval XSHORE works on
_MAXDEP = 45          # ft, how deep the generic profile is carried (KD2.FOR:17)

_FT = 0.3048          # m per ft
_FT2 = _FT * _FT
_FT3 = _FT2 * _FT


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
    enable_if: tuple = ()    # (other_key, value): gray out unless that input == value


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    note: str = ""


APP_META = AppMeta(
    aces_id="6-2",
    name="Time-Dependent Beach and Dune Erosion",
    area="Littoral Processes",
    classification="exact",
    cite="Kriebel (1984b, EBEACH); Kriebel & Dean (1985); Moore (1982); "
         "Dean (1977); ACES Technical Reference Ch. 6-2",
    default_system="US",
)

_GENERIC = "Generic (dune, berm and slopes)"
_SURVEYED = "Surveyed profile"
_PROFILE_TYPES = (_GENERIC, _SURVEYED)

# ACES deck XSHORE2.IN: a generic profile, one wave record, no water-level forcing.
_EX_WAVES = ((8.0, 8.0, 10.0),)           # height ft, period s, angle deg
# ACES deck XSHORE4.IN, a surveyed profile: (x ft, bed elevation ft)
_EX_SURVEY = (
    (0, 14.1), (4, 13.4), (11.2, 13.1), (25.1, 10.6), (45, 15),
    (54.4, 14.1), (75.7, 12.5), (105.2, 12.9), (139.6, 13.5), (163.9, 12.5),
    (189.4, 10.5), (205.5, 8.6), (242.5, 4.3), (281.6, 2.3), (320.6, 1.1),
    (374.7, 0.4), (393.7, 0.2), (421.3, -0.5), (453.6, -3.1), (497.3, -6.9),
    (539.2, -7), (577.4, -6.6), (626.7, -7.6), (638.7, -8.7), (672.6, -9.8),
    (721.9, -9.7), (735.7, -8.8), (772.9, -9), (821.8, -7.6), (883.9, -7.5),
    (957, -10), (975.6, -11.3), (998, -12), (1028, -13.4), (1076, -16.1),
    (1120, -18.1), (1153, -19.2), (1190, -20.5), (1226, -21.5), (1285, -22.9),
    (1316, -23), (1372, -24.4), (1421, -25.5), (1485, -26.5), (1532, -27.2),
    (1585, -28.3), (1625, -29.2), (1682, -30.1), (1723, -30.4), (1777, -31.1),
    (1821, -31.5), (1870, -32.2), (1916, -32.4),
)

INPUTS = (
    Field("profile_type", "Profile definition", "choice", "", "", default=_GENERIC,
          choices=_PROFILE_TYPES,
          note="build a generic profile from dune, berm and slope descriptors, or "
               "supply a surveyed one"),
    Field("dune_elev", "Dune crest elevation", "float", "m", "ft",
          default=20.0 * _FT, lo=0.0, hi=200.0 * _FT, enable_if=("profile_type", _GENERIC),
          note="elevation of the dune top above the datum"),
    Field("dune_base_elev", "Dune base elevation", "float", "m", "ft",
          default=6.0 * _FT, lo=0.0, hi=200.0 * _FT, enable_if=("profile_type", _GENERIC),
          note="elevation at the foot of the dune face; must be at least 5 ft below the crest"),
    Field("dune_width", "Dune crest width", "float", "m", "ft",
          default=50.0 * _FT, lo=0.0, hi=1000.0 * _FT, enable_if=("profile_type", _GENERIC)),
    Field("berm_elev", "Berm elevation", "float", "m", "ft",
          default=6.0 * _FT, lo=0.0, hi=200.0 * _FT, enable_if=("profile_type", _GENERIC)),
    Field("berm_width", "Berm width", "float", "m", "ft",
          default=100.0 * _FT, lo=0.0, hi=2000.0 * _FT, enable_if=("profile_type", _GENERIC)),
    Field("slope_dune", "Dune face slope (run per rise)", "float", "", "",
          default=2.0, lo=0.1, hi=1000.0, enable_if=("profile_type", _GENERIC),
          note="1 on this many; 2 means 1:2"),
    Field("slope_beach", "Beach face slope (run per rise)", "float", "", "",
          default=10.0, lo=0.1, hi=1000.0, enable_if=("profile_type", _GENERIC)),
    Field("slope_nearshore", "Nearshore slope (run per rise)", "float", "", "",
          default=20.0, lo=0.1, hi=1000.0, enable_if=("profile_type", _GENERIC)),
    Field("survey", "Surveyed profile", "matrix", "m", "m",
          default=tuple((round(x * _FT, 6), round(z * _FT, 6)) for x, z in _EX_SURVEY),
          enable_if=("profile_type", _SURVEYED),
          note="rows of (distance seaward from the baseline, bed elevation above datum); "
               "resampled onto the 4 ft grid by natural cubic spline"),
    Field("grain_size", "Mean grain size", "float", "mm", "mm", default=0.22,
          lo=0.01, hi=10.0,
          note="sets the equilibrium profile shape factor A (Moore 1982)"),
    Field("gauge_depth", "Depth of the wave records", "float", "m", "ft",
          default=60.0 * _FT, lo=5.0 * _FT, hi=9999.0 * _FT,
          note="still-water depth at which the wave conditions below were measured"),
    Field("waves", "Wave records", "matrix", "m", "m",
          default=tuple((round(h * _FT, 6), t, a) for h, t, a in _EX_WAVES),
          note="rows of (height, period s, crest angle to the shoreline deg), one per "
               "wave interval; heights are SI on the contract"),
    Field("wave_dt_hr", "Wave record interval", "float", "h", "h", default=20.0,
          lo=1.0, hi=1000.0),
    Field("length_hr", "Length of simulation", "float", "h", "h", default=20.0,
          lo=1.0, hi=1000.0),
    Field("out_interval_hr", "Tabular output interval", "float", "h", "h", default=2.0,
          lo=1.0, hi=1000.0,
          note="ignored when explicit output times are given below"),
    Field("out_times", "Explicit output times", "list", "h", "h", default=(),
          note="report at these hours instead of at a fixed interval; empty to use the "
               "interval. ACES takes one or the other, never both"),
    Field("water_level", "Still-water level series", "list", "m", "ft", default=(),
          note="surge or tide elevation above the datum, at the interval below; empty "
               "for none. ACES applies it as a uniform shift of the whole profile"),
    Field("wl_dt_hr", "Water-level interval", "float", "h", "h", default=1.0,
          lo=1.0, hi=1000.0),
)

OUTPUTS = (
    # The plotted family is the profile itself: the front-ends give every "profile"
    # output one shared horizontal axis, and the hourly table below is on a different
    # one. The table's endpoints are therefore reported as scalars, and the full
    # series stay on the Result for anyone driving compute() directly.
    Out("profile_X", "Distance seaward from baseline", "m", "ft", "profile",
        note="Horizontal coordinate of the 4 ft computational grid, measured seaward from the landward baseline."),
    Out("z_initial", "Initial bed elevation", "m", "ft", "profile",
        note="Bed elevation above the still-water datum at the start of the simulation."),
    Out("z_final", "Final bed elevation", "m", "ft", "profile",
        note="Bed elevation above the datum at the end of the simulation, after the profile has responded to the wave record."),
    Out("vol_total", "Change in subaerial volume", "m^3/m", "yd^3/ft", "scalar",
        note="Change in the volume of sand above the still-water level over the whole simulation, per unit length of shoreline; negative is erosion."),
    Out("recession", "Shoreline recession", "m", "ft", "scalar",
        note="Landward migration of the 0 ft contour over the whole simulation, reported positive for erosion."),
    Out("contour_5_final", "Contour change, +5 ft", "m", "ft", "scalar",
        note="Net migration of the +5 ft contour. ACES fixes these levels at 0, 5, 10 and 15 ft whatever the display units."),
    Out("contour_10_final", "Contour change, +10 ft", "m", "ft", "scalar",
        note="Net migration of the +10 ft contour."),
    Out("contour_15_final", "Contour change, +15 ft", "m", "ft", "scalar",
        note="Net migration of the +15 ft contour."),
    Out("A_shape", "Equilibrium profile shape factor", "m^(1/3)", "ft^(1/3)", "scalar",
        note="Dean's A(D50) from the mean grain size, the scale of the equilibrium profile h = A x^(2/3)."),
    Out("D_eq", "Equilibrium energy dissipation", "", "", "scalar",
        note="Equilibrium wave energy dissipation per unit water volume in the surf zone, the target the scheme drives toward."),
)


@dataclass
class Result:
    profile_X: np.ndarray; z_initial: np.ndarray; z_final: np.ndarray
    # the hourly table ACES prints; on its own time axis, so not declared as plotted
    # outputs (see the note on OUTPUTS)
    hours: np.ndarray; vol_change: np.ndarray
    contour_0: np.ndarray; contour_5: np.ndarray
    contour_10: np.ndarray; contour_15: np.ndarray
    vol_total: float; recession: float
    contour_5_final: float; contour_10_final: float; contour_15_final: float
    A_shape: float; D_eq: float
    notes: str = ""


# --- ACES linear wave theory (ACECLIB.FOR) --------------------------------------
_WAVLEN_COEFF = (.66667, .35550, .16084, .06320, .02174,
                 .00654, .00171, .00039, .00011)


def wavlen(t: float, d: float, g: float = G_US) -> float:
    """Wavelength by ACES's nine-term series (ACECLIB.FOR:492)."""
    twopi = 2.0 * math.pi
    y = (twopi / t) ** 2 * d / g
    s = 1.0
    for n, c in enumerate(_WAVLEN_COEFF, 1):
        s += c * y ** n
    s = 1.0 / s
    s = 1.0 / (y + s)
    return t * math.sqrt(g * d * s)


def lwtgen(d: float, t: float, g: float = G_US):
    """Celerity, group celerity and wavelength at depth d (ACECLIB.FOR:186)."""
    twopi = 2.0 * math.pi
    c0 = g * t / twopi
    cg0 = 0.5 * c0
    l0 = c0 * t
    l = wavlen(t, d, g)
    reldep = d / l
    if reldep > 0.5:
        return c0, c0, cg0, cg0, l0, l0, reldep
    c = l / t
    arg = 4.0 * math.pi * reldep
    arg = arg / math.sinh(arg)
    cg = 0.5 * c * (1.0 + arg)
    return c, c0, cg, cg0, l, l0, reldep


def lwttws(alpha0: float, c: float, cg: float, c0: float, h0: float):
    """Snell refraction and shoaling to a transitional depth (ACECLIB.FOR:293).

    Returns (alpha deg, H, K_r, K_s)."""
    arg = (c / c0) * math.sin(alpha0 / DEGRAD)
    arg = max(-1.0, min(1.0, arg))
    alpha = math.asin(arg) * DEGRAD
    ks = math.sqrt(0.5 * c0 / cg)
    kr = math.sqrt(math.cos(alpha0 / DEGRAD) / math.cos(alpha / DEGRAD))
    return alpha, h0 * ks * kr, kr, ks


def breaking_point(hd: float, period: float, dang: float, d: float, g: float = G_US):
    """Breaking wave height and the depth it occurs in (KD3.FOR).

    Bisects on depth until the shoaled and refracted height matches the depth-limited
    height 0.78 d. Note that the crest angle passed to LWTTWS is the deepwater angle
    and the source never updates it: LWTTWS returns the local angle in a separate
    variable, which KD3 discards. Returns (H_b, d_b)."""
    c, c0, cg, cg0, _l, _l0, _r = lwtgen(d, period, g)
    ddeep = d
    dshal = 1.0
    h = hd
    ang = dang
    cg = cg0                                   # KD3.FOR:17, before the iteration
    ibit = 0
    hb = 0.78 * d
    for _ in range(200):                       # the source loops on GOTO; bound it here
        hb = 0.78 * d
        if abs((hb - h) / h) <= 0.002:
            return hb, d
        if h < hb:                             # label 110: step shoreward
            ddeep = d
            d = 0.5 * (dshal + d)
            c, c0, cg, cg0, _l, _l0, _r = lwtgen(d, period, g)
            _alpha, h, _kr, _ks = lwttws(ang, c, cg, c0, hd)
            ibit = 1
            continue
        dshal = d
        if ibit == 0:
            # first pass and the wave is already broken at the gauge: take the
            # depth-limited height here and stop (KD3.FOR:26-34)
            ang = dang
            c, c0, cg, cg0, _l, _l0, _r = lwtgen(d, period, g)
            _alpha, h, _kr, _ks = lwttws(ang, c, cg, c0, hd)
            return hb, d
        d = 0.5 * (ddeep + d)
        c, c0, cg, cg0, _l, _l0, _r = lwtgen(d, period, g)
        _alpha, h, _kr, _ks = lwttws(ang, c, cg, c0, hd)
    return hb, d


# --- profile construction -------------------------------------------------------
def generic_profile(ddune, wdune, sdune, dbase, dberm, wberm, sbeach, snear, dx=_DX):
    """Build the generic dune/berm/nearshore profile (KD2.FOR).

    Slopes are rise per run here, the reciprocal of the run-per-rise the inputs carry.
    Returns (x, z) as 1-based lists with index 0 unused, and the point count."""
    n = 1001
    x = [0.0] * (n + 6)
    z = [0.0] * (n + 6)
    for j in range(2, n + 1):
        x[j] = x[j - 1] + dx

    xend = 0
    if wdune > 0.0 and ddune > 0.0:
        xend = int(wdune / dx) + 1
        for j in range(1, xend + 1):
            z[j] = ddune

    if ddune > 0.0:
        xbegin = xend + 1
        h = ddune - dbase
        w = h / sdune
        xend = xend + int(w / dx)
        if math.fmod(w + 0.1, float(dx)) >= 2:
            xend += 1
        step = h / ((xend - xbegin) + 1)
        sdune = h / (((xend - xbegin) + 1) * dx)
        if xbegin > 1:
            for j in range(xbegin, xend + 1):
                z[j] = z[j - 1] - step
        else:
            z[1] = ddune
            for j in range(xbegin + 1, xend + 1):
                z[j] = z[j - 1] - step

    xbegin = xend + 1
    if wberm > 0.0:
        xend = xend + int(wberm / dx)
        if math.fmod(wberm, float(dx)) >= 2.0:
            xend += 1
    h = dbase - dberm
    step = 0.0
    if (xend - xbegin) + 1 > 0:
        step = h / ((xend - xbegin) + 1)
    if xbegin == 1:
        z[1] = dbase
        for j in range(xbegin + 1, xend + 1):
            z[j] = z[j - 1] - step
    else:
        for j in range(xbegin, xend + 1):
            z[j] = z[j - 1] - step

    xbegin = xend + 1
    if sbeach > 0.0 and dberm > 0.0:
        w = dberm / sbeach
        xend = xend + int(w / dx) + 1
        if math.fmod(w + 0.1, float(dx)) >= 2:
            xend += 1
        step = dberm / ((xend - xbegin) + 1)
        sbeach = dberm / (((xend - xbegin) + 1) * dx)
        if xbegin == 1:
            z[1] = dberm
            for j in range(xbegin + 1, xend + 1):
                z[j] = z[j - 1] - step
        else:
            for j in range(xbegin, xend + 1):
                z[j] = z[j - 1] - step

    xbegin = xend + 1
    w = _MAXDEP / snear
    bottom = int(w / dx)
    xend = xbegin + bottom
    step = dx * snear
    for j in range(xbegin, xend + 1):
        z[j] = z[j - 1] - step
    return x, z, xend, sdune, sbeach


def _spline_natural(xs, ys):
    """Natural cubic spline coefficients, KD5A with IEND = 1 (which is what KD5 uses).

    Returns per-interval (a3, a2, a1, a0) as in the source's A(i, 1:4)."""
    n = len(xs)
    s = [0.0] * n
    if n >= 3:
        nm2 = n - 2
        a1 = [0.0] * (nm2 + 1)
        a2 = [0.0] * (nm2 + 1)
        a3 = [0.0] * (nm2 + 1)
        a4 = [0.0] * (nm2 + 1)
        dx1 = xs[1] - xs[0]
        dy1 = (ys[1] - ys[0]) / dx1 * 6.0
        for i in range(1, nm2 + 1):
            dx2 = xs[i + 1] - xs[i]
            dy2 = (ys[i + 1] - ys[i]) / dx2 * 6.0
            a1[i] = dx1
            a2[i] = 2.0 * (dx1 + dx2)
            a3[i] = dx2
            a4[i] = dy2 - dy1
            dx1, dy1 = dx2, dy2
        for i in range(2, nm2 + 1):
            a2[i] -= a1[i] / a2[i - 1] * a3[i - 1]
            a4[i] -= a1[i] / a2[i - 1] * a4[i - 1]
        a4[nm2] /= a2[nm2]
        for i in range(2, nm2 + 1):
            j = (n - 1) - i
            a4[j] = (a4[j] - a3[j] * a4[j + 1]) / a2[j]
        for i in range(1, nm2 + 1):
            s[i] = a4[i]
        s[0] = 0.0
        s[n - 1] = 0.0
    coef = []
    for i in range(n - 1):
        h = xs[i + 1] - xs[i]
        coef.append(((s[i + 1] - s[i]) / 6.0 / h,
                     s[i] / 2.0,
                     (ys[i + 1] - ys[i]) / h - h * (2.0 * s[i] + s[i + 1]) / 6.0,
                     ys[i]))
    return coef


def resample_profile(xs, zs, dx=_DX):
    """Resample a surveyed profile onto the dx grid by natural cubic spline (KD5)."""
    coef = _spline_natural(list(xs), list(zs))
    pts = int(xs[-1] / dx) + 1
    n = 1001
    x = [0.0] * (n + 6)
    z = [0.0] * (n + 6)
    for i in range(pts):
        xx = float(i) * dx
        x[i + 1] = xx
        val = 0.0
        for k in range(len(xs) - 1):
            if xs[k] <= xx <= xs[k + 1]:
                t = xx - xs[k]
                a3, a2, a1, a0 = coef[k]
                val = a3 * t ** 3 + a2 * t ** 2 + a1 * t + a0
                break
        z[i + 1] = val
    return x, z, pts


def synchronise(stint, wtint, simtim, surgel, waves):
    """Put the water-level series and the wave records on one interval (KD4.FOR).

    ACES calls this twice: once to spline the water level onto the shorter interval,
    and once to step the wave records onto it. Returns
    (surgel, waves, numrec, tint)."""
    nsr = int(simtim / stint) + 1
    nwr = int(simtim / wtint)
    surgel = list(surgel)
    waves = list(waves)
    tint = wtint
    if stint > wtint:
        # resample the water level by the same cubic spline the profile uses
        xint = wtint
        nnr = int(simtim / xint) + 1
        t = [(i) * stint for i in range(nsr)]
        y = [surgel[i] if i < len(surgel) else 0.0 for i in range(nsr)]
        coef = _spline_natural(t, y)
        out = []
        for i in range(nnr):
            xx = float(i) * xint
            val = 0.0
            for kk in range(len(t) - 1):
                if t[kk] <= xx <= t[kk + 1]:
                    d = xx - t[kk]
                    a3, a2, a1, a0 = coef[kk]
                    val = a3 * d ** 3 + a2 * d ** 2 + a1 * d + a0
                    break
            out.append(val)
        surgel = out
        tint = wtint
        nsr -= 1
    elif stint != wtint:
        # step the wave records onto the shorter interval, holding each until the
        # next one is closer (KD4.FOR:53-72)
        hhour = wtint
        tint = stint
        j = 1
        held = []
        for i in range(1, nsr + 1):
            shour = i * stint
            held.append(waves[min(j, len(waves)) - 1])
            if ((hhour - shour) < stint
                    and (shour + stint) * 10 >= hhour * 10 + wtint * 5):
                j += 1
                if j > nwr:
                    j = nwr
                hhour += wtint
        waves = held
    numrec = nwr
    if nsr > nwr:
        numrec = nsr
    return surgel, waves, numrec, tint


def shape_factor(gs_mm: float) -> float:
    """Dean's A(D50), in ft^(1/3), by ACES's two-branch fit (KDMAIN.FOR:465)."""
    if gs_mm < 0.28:
        return (10.0 ** (1.082 * math.log10(gs_mm) - 0.256)) / 0.3048 ** 0.3333
    return (10.0 ** (0.303 * math.log10(gs_mm) - 0.686)) / 0.3048 ** 0.3333


# --- 'Method & equations' panel content ----------------------------------------
ABOUT = {
    'summary':
        'Marches a beach and dune profile through a storm by the XSHORE explicit '
        'finite-difference scheme, the one ACES 6-2 runs. Sand is moved across the '
        'surf zone in proportion to how far the local wave energy dissipation per unit '
        'water volume departs from its equilibrium value, and the profile is advanced '
        'by continuity on a fixed 4 ft grid with a sub-hourly step that halves itself '
        'wherever the bed steepens. Reports the profile at each output time, the '
        'change in subaerial sand volume, and the migration of the 0, 5, 10 and 15 ft '
        'contours. The scheme is stated in US customary units and the grid interval is '
        '4 ft exactly; inputs and results are carried in SI and converted at that '
        'boundary.',
    'method_key': 'profile_type',
    'methods': [{
        'name': 'XSHORE explicit finite-difference scheme',
        'when': None,
        'tag': '',
        'note': None,
        'equations': [
            {'tex': r'D_i = \frac{K_D \left( |z_i|^{5/2} - |z_{i+1}|^{5/2} \right)}'
                    r'{\bar{z}\,\Delta x}, \quad K_D = \tfrac{1}{8}\rho\,\gamma^{2} g^{3/2}',
             'desc': 'Wave energy dissipation per unit volume of water across the cell, '
                     'from the depth-limited height H = gamma z.'},
            {'tex': r'D_{eq} = \frac{5\,\rho\,\gamma^{2} g^{3/2}}{24}\,A^{3/2}',
             'desc': 'The equilibrium dissipation the profile is driven toward, set by '
                     'the grain size alone through the shape factor A.'},
            {'tex': r'q_i = K_t \left( D_i - D_{eq} \right)',
             'desc': 'Cross-shore transport, offshore where dissipation exceeds '
                     'equilibrium and onshore where it falls short. K_t = 0.0011445 '
                     'ft^4/lb.'},
            {'tex': r'\frac{\partial z}{\partial t} = -\frac{1}{\Delta x}'
                    r'\frac{\partial q}{\partial x}',
             'desc': 'Continuity: the bed rises or falls with the gradient of transport.'},
            {'tex': r'z_{tan} = -A \left( \frac{2A}{3 m_b} \right)^{2}',
             'desc': 'Landward limit of the surf-zone calculation, the depth at which '
                     'the equilibrium profile is tangent to the linear beach face of '
                     'slope m_b.'},
            {'tex': r'\Delta t = \frac{\Delta x^{2}}{2\,\alpha \sqrt{-z_b}}, \quad '
                    r'\alpha = \tfrac{5}{16} K_t\,\rho\,g^{3/2}\gamma^{2}',
             'desc': 'Stable sub-hourly step for the explicit scheme, from the breaking '
                     'depth; halved locally wherever the bed steps down by more than '
                     '0.5 ft.'},
            {'tex': r'A = 10^{\,1.082 \log_{10} D_{50} - 0.256} \;\; (D_{50} < 0.28\ \mathrm{mm})',
             'desc': 'Equilibrium profile shape factor from mean grain size (Moore 1982), '
                     'in ACES\'s two-branch fit; the other branch is '
                     '10^(0.303 log10 D50 - 0.686).'},
        ]}],
    'symbols': [
        ['z(x)', 'Bed elevation above the still-water datum (ft internally)'],
        ['x', 'Distance seaward from the landward baseline'],
        ['q', 'Cross-shore sediment transport rate per unit width'],
        ['D', 'Wave energy dissipation per unit volume of water in the surf zone'],
        ['D_{eq}', 'Equilibrium value of that dissipation'],
        ['A', 'Equilibrium profile shape factor from grain size'],
        ['\\gamma', 'Breaker index, wave height over water depth, 0.78'],
        ['K_t', 'Transport coefficient, 0.0011445 ft^4/lb'],
        ['m_b', 'Beach face slope, rise over run'],
        ['z_b', 'Bed elevation at the breaking point, the seaward calculation limit'],
        ['\\Delta x', 'Grid interval, fixed at 4 ft'],
    ],
    'references': ['Kriebel (1984b), EBEACH users manual vol II: theory and background',
                   'Kriebel & Dean (1985), Coastal Engineering 9:221-245',
                   'Moore (1982)',
                   'Dean (1977)',
                   'ACES Technical Reference Ch. 6-2'],
}


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} out of range: {v} not in [{f.lo}, {f.hi}]")
    if not inp.get("waves"):
        raise ValueError("at least one wave record is required")


def compute(inp: dict, *, g: float = G_US) -> Result:
    """Run the XSHORE march. Inputs are SI; results are returned in SI."""
    _validate(inp)
    ptype = str(inp.get("profile_type", _GENERIC))
    gs = float(inp["grain_size"])
    dd = float(inp["gauge_depth"]) / _FT
    wtint = float(inp["wave_dt_hr"])
    simtim = float(inp["length_hr"])
    deltat = float(inp["out_interval_hr"])
    out_times = sorted(float(v) for v in (inp.get("out_times") or ()))
    waves = [(float(h) / _FT, float(t), float(a)) for h, t, a in inp["waves"]]
    wl = [float(v) / _FT for v in (inp.get("water_level") or ())]
    wl_dt = float(inp.get("wl_dt_hr", 1.0))

    dx = _DX
    rho = RHO_US
    kdiss = 0.125 * rho * GAMMA ** 2 * g ** 1.5
    alpha_s = (5.0 / 16.0) * KTRANS * rho * g ** 1.5 * GAMMA ** 2
    a_shape = shape_factor(gs)
    deq = ((5.0 * rho * GAMMA ** 2 * g ** 1.5) / 24.0) * a_shape ** 1.5

    ddune = float(inp["dune_elev"]) / _FT
    dbase = float(inp["dune_base_elev"]) / _FT
    dberm = float(inp["berm_elev"]) / _FT
    sbeach = 1.0 / float(inp["slope_beach"]) if ptype == _GENERIC else 0.0
    sdune = 0.0
    if ptype == _GENERIC:
        x, z0, dpts, sdune, sbeach = generic_profile(
            ddune, float(inp["dune_width"]) / _FT, 1.0 / float(inp["slope_dune"]),
            dbase, dberm, float(inp["berm_width"]) / _FT, sbeach,
            1.0 / float(inp["slope_nearshore"]), dx)
    else:
        rows = [(float(a) / _FT, float(b) / _FT) for a, b in inp["survey"]]
        x, z0, dpts = resample_profile([r[0] for r in rows], [r[1] for r in rows], dx)

    n = len(z0)
    # depth[0] is the original profile, depth[1] and depth[2] the working pair, exactly
    # as DEPTH(0:2, NPTS) is used in the source
    depth = [list(z0), [0.0] * n, [0.0] * n]
    qsed = [0.0] * n
    qtemp = [0.0] * n
    diss = [0.0] * n

    # landward reference point and the initial contour positions (KDMAIN.FOR:500-548)
    cnt = 0
    nds = 0
    # KD1 accumulates the local face slopes ONTO whatever KD2 left in these, then
    # divides the total by the count, so the builder's value is the first term of the
    # average rather than being discarded (KDMAIN.FOR:507, :544)
    sdune_acc = sdune if ptype == _GENERIC else 0.0
    sbeach_acc = sbeach
    ft0 = ft5 = ft10 = ft15 = 0.0
    for i in range(dpts, 0, -1):
        if dberm != 0.0:
            if depth[0][i] + 0.01 >= dberm and depth[0][i + 1] < dberm and nds == 0:
                nds = i
            if depth[0][i] < dberm and depth[0][i + 1] > 0 and nds == 0:
                sbeach_acc += (depth[0][i] - depth[0][i + 1]) / dx
                cnt += 1
        elif dbase != 0.0:
            if depth[0][i] + 0.01 >= dbase and depth[0][i + 1] < dbase and nds == 0:
                nds = i
            if depth[0][i] < dbase and depth[0][i + 1] > 0 and nds == 0:
                sbeach_acc += (depth[0][i] - depth[0][i + 1]) / dx
                cnt += 1
        else:
            if depth[0][i] + 0.01 >= ddune and depth[0][i + 1] < ddune and nds == 0:
                nds = i
            if depth[0][i] < ddune and depth[0][i + 1] > dbase and nds == 0:
                sdune_acc += (depth[0][i] - depth[0][i + 1]) / dx
                cnt += 1
        for lev, name in ((0.0, "0"), (5.0, "5"), (10.0, "10"), (15.0, "15")):
            if depth[0][i] >= lev > depth[0][i + 1]:
                v = x[i + 1] - ((x[i + 1] - x[i]) * (lev - depth[0][i + 1])
                                / (depth[0][i] - depth[0][i + 1]))
                if name == "0":
                    ft0 = v
                elif name == "5":
                    ft5 = v
                elif name == "10":
                    ft10 = v
                else:
                    ft15 = v
        depth[1][i] = depth[0][i]
        depth[2][i] = depth[0][i]
    if cnt:
        if dberm != 0:
            sbeach = sbeach_acc / cnt
        else:
            sdune_acc /= cnt
            if ptype == _SURVEYED and sbeach == 0.0:
                sbeach = sdune_acc

    # put the water level and the wave records on one interval, as KD1 does before
    # it begins (KDMAIN.FOR:566-582)
    isurge = bool(wl)
    mintint = min(wl_dt, wtint) if isurge else wtint
    tint = mintint
    nrec = len(waves)
    if isurge and mintint != wl_dt:
        wl, _w, _n, tint = synchronise(wl_dt, mintint, simtim, wl, waves)
    if wtint != mintint:
        _s, waves, nrec, tint = synchronise(mintint, wtint, simtim, wl, waves)
    nwr = max(nrec, len(waves)) if wtint != mintint else len(waves)
    kflg = 2 if dberm == 0 else 1
    oldflg = 0
    ckslp = 0.0
    proht = 0.0
    ktime = 0
    vol0 = 0.0
    nstop = 0
    ftz = ftv = ftx = ftxv = 0.0
    hours = []; volc = []; c0 = []; c5 = []; c10 = []; c15 = []
    next_out = out_times[0] if out_times else deltat
    out_k = 0
    msg = ""

    for i1 in range(1, nwr + 1):
        hd1, period, dang = waves[i1 - 1]
        d = dd
        dtan = -a_shape * ((2.0 * a_shape) / (3.0 * sbeach)) ** 2
        hb1, db_depth = breaking_point(hd1, period, dang, d, g)
        db = -(hb1 / 0.78)

        # water level: ACES shifts the whole profile once per wave record, not per
        # hour, and the record index is what selects the level (KDMAIN.FOR:625-646)
        if isurge:
            eta = wl[i1 - 1] if i1 - 1 < len(wl) else (wl[-1] if wl else 0.0)
            surge = proht - eta
            proht = eta
            posval = 0
            for i2 in range(1, dpts + 1):
                depth[1][i2] += surge
                depth[2][i2] = depth[1][i2]
                if depth[1][i2] > 0.0:
                    posval += 1
            if posval < 1:
                msg = "entire profile submerged; check the water levels and beach data"
                break
            ddune += surge
            dberm += surge
            dbase += surge

        for _i2 in range(int(round(tint))):
            ktime += 1
            if ktime > simtim:
                break

            # boundaries: landward tangent point, seaward breaking point
            ndtan = 1
            i3 = 1
            while True:
                if depth[1][i3] >= dtan:
                    ndtan = i3
                if depth[1][i3] < db and depth[1][i3 - 1] >= db:
                    break
                i3 += 1
                if i3 > dpts:
                    break
            nbreak = i3 - 1

            if depth[1][nds] + 0.1 >= dbase and dberm > 0:
                kflg = 1
                top = depth[1][nds]
            else:
                kflg = 2
                top = ddune

            while True:
                if oldflg != kflg:
                    nstop = 0                      # DO 65 falling through leaves I3 = 0
                    for i3 in range(ndtan, 0, -1):
                        if depth[1][i3] + 0.1 >= top:
                            nstop = i3
                            break
                    ckslp = (depth[1][nstop] - depth[1][nstop + 1]) * 0.90
                oldflg = kflg
                nds = nstop
                if ndtan - nstop < 2:
                    if kflg == 1:
                        kflg = 2
                        top = ddune
                        continue
                    msg = "entire profile submerged; check the water levels and beach data"
                break
            if msg:
                break

            if abs(depth[1][nstop] - depth[1][nstop + 1]) > ckslp:
                nstop -= 1
            if nstop < 1:
                nstop = 1

            dt = (0.5 * dx ** 2) / (alpha_s * (-db) ** 0.5)
            stbtime = int(round(3600.0 / (dt - 0.5)))
            io, ix = 1, 2

            i3 = 1
            while True:
                restart = False
                for i4 in range(ndtan, nbreak + 1):
                    disnum = abs(depth[io][i4]) ** 2.5 - abs(depth[io][i4 + 1]) ** 2.5
                    disavg = (depth[io][i4] + depth[io][i4 + 1]) / 2.0
                    if depth[io][i4 - 1] - depth[io][i4] < -0.5:
                        dt = dt / 2.0
                        if dt >= 1.5:
                            stbtime = int(round(3600.0 / (dt - 0.5)))
                            io, ix = (2, 1) if io == 1 else (1, 2)
                            if i3 > 1:
                                i3 -= 1
                            restart = True
                            break
                    diss[i4] = (kdiss * disnum) / (disavg * dx)
                    qtemp[i4] = KTRANS * (diss[i4] - deq)
                if restart:
                    continue

                for i4 in range(ndtan, nbreak + 4):
                    qsed[i4] = (0.1 * qtemp[i4 - 2] + 0.15 * qtemp[i4 - 1]
                                + 0.5 * qtemp[i4] + 0.15 * qtemp[i4 + 1]
                                + 0.1 * qtemp[i4 + 2])
                    diss[i4] = deq + qsed[i4] / KTRANS

                if nstop < 1 or ndtan == nstop:
                    msg = "landward profile boundary could not be established"
                    break
                qfshor = qsed[ndtan] / (ndtan - nstop)
                qsed[nstop] = 0.0
                for i4 in range(nstop + 1, ndtan):
                    qsed[i4] = qsed[i4 - 1] + qfshor

                for i4 in range(nstop + 1, nbreak + 4):
                    depth[ix][i4] = depth[io][i4] - (dt / dx) * (qsed[i4] - qsed[i4 - 1])

                io, ix = (2, 1) if io == 1 else (1, 2)
                i3 += 1
                if i3 > stbtime:
                    break
            if msg:
                break

            if io == 2:
                for i3 in range(1, dpts + 1):
                    depth[1][i3], depth[2][i3] = depth[2][i3], depth[1][i3]
            else:
                for i3 in range(1, dpts + 1):
                    depth[2][i3] = depth[1][i3]

            # subaerial volume and the four contour positions
            vol = 0.0
            for m in range(2, nbreak + 5):
                qtemp[m] = 0.0
                if ktime == 1 and depth[0][m] > 0:
                    vol0 += (depth[0][m] + depth[0][m - 1]) * dx / 2
                if depth[1][m] >= -proht:
                    vol += (depth[1][m] + depth[1][m - 1] + 2 * proht) * dx / 2
                if depth[1][m] < -proht and depth[1][m - 1] > -proht:
                    vol += (((depth[1][m - 1] + proht) * dx)
                            / (depth[1][m - 1] - depth[1][m])
                            * ((depth[1][m - 1] + proht) / 2))
                for lev in (0.0, 5.0, 10.0, 15.0):
                    if depth[1][m] + proht >= lev > depth[1][m + 1] + proht:
                        v = x[m + 1] - ((x[m + 1] - x[m])
                                        * (lev - (depth[1][m + 1] + proht))
                                        / (depth[1][m] - depth[1][m + 1]))
                        if lev == 0.0:
                            ftz = v
                        elif lev == 5.0:
                            ftv = v
                        elif lev == 10.0:
                            ftx = v
                        else:
                            ftxv = v
            vchange = vol - vol0

            if ktime + 1e-9 >= next_out:
                hours.append(float(ktime))
                volc.append(vchange / 27.0)          # ft^3/ft -> yd^3/ft, as ACES reports
                c0.append(ftz - ft0)
                c5.append(ftv - ft5)
                c10.append(ftx - ft10)
                c15.append(ftxv - ft15)
                if out_times:
                    out_k += 1
                    next_out = (out_times[out_k] if out_k < len(out_times)
                                else float("inf"))
                else:
                    next_out += deltat
        if msg or ktime >= simtim:
            break

    z_final = [depth[1][i] + proht for i in range(1, dpts + 1)]
    z_init = [depth[0][i] for i in range(1, dpts + 1)]
    xs = [x[i] for i in range(1, dpts + 1)]

    _YD3_FT = 0.764554858 / _FT          # yd^3/ft -> m^3/m
    notes = (f"{dpts} points at {dx} ft; A={a_shape:.4f} ft^(1/3), "
             f"D_eq={deq:.4f}; {nwr} wave record(s) over {ktime} h")
    if msg:
        notes += f"; stopped: {msg}"
    return Result(
        profile_X=np.array(xs) * _FT, z_initial=np.array(z_init) * _FT,
        z_final=np.array(z_final) * _FT,
        hours=np.array(hours), vol_change=np.array(volc) * _YD3_FT,
        contour_0=np.array(c0) * _FT, contour_5=np.array(c5) * _FT,
        contour_10=np.array(c10) * _FT, contour_15=np.array(c15) * _FT,
        vol_total=(volc[-1] if volc else 0.0) * _YD3_FT,
        recession=-(c0[-1] if c0 else 0.0) * _FT,
        contour_5_final=(c5[-1] if c5 else 0.0) * _FT,
        contour_10_final=(c10[-1] if c10 else 0.0) * _FT,
        contour_15_final=(c15[-1] if c15 else 0.0) * _FT,
        A_shape=a_shape * _FT ** (1.0 / 3.0), D_eq=deq, notes=notes)


# --- self-tests (against the recompiled ACES source) ----------------------------
def _approx(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


# The ACES source, recompiled and run on the decks it ships with
# (tests/aces_oracle/fortran, `sh build.sh kd`). Each row is
# (hour, volume change yd^3/ft, then the 0 / +5 / +10 / +15 ft contour changes ft),
# read from the source's own kd_profile.txt.
_SRC_XSHORE2 = (
    (2, -0.01, -0.18, -0.18, 0.00, 0.00),
    (4, -0.41, -1.96, -1.96, 0.00, 0.00),
    (6, -0.89, -4.14, -4.14, 0.00, 0.00),
    (8, -1.51, -5.93, -5.93, -0.09, -0.09),
    (10, -2.13, -7.32, -6.88, -0.27, -0.27),
    (12, -2.84, -7.99, -7.98, -0.48, -0.48),
    (14, -3.58, -9.13, -9.13, -0.69, -0.69),
    (16, -4.33, -10.30, -10.30, -0.91, -0.91),
    (18, -4.96, -11.57, -11.27, -1.10, -1.10),
    (20, -5.67, -12.41, -12.35, -1.31, -1.31),
)
# XSHORE3: the same generic profile with a six-point surge on a 4 h interval and four
# wave records on a 5 h one, which ACES steps onto the shorter interval before it runs
_SRC_XSHORE3 = (
    (1, 0.08, 0.21, 0.22, 0.00, 0.00),
    (2, -0.01, -0.18, -0.18, 0.00, 0.00),
    (4, -0.41, -1.96, -1.96, 0.00, 0.00),
    (6, -1.94, -1.97, -9.64, 0.00, 0.00),
    (8, -2.81, 0.52, -10.92, -0.26, -0.26),
    (10, -3.50, 20.63, -12.80, -0.74, -0.74),
    (12, -3.92, 31.47, -14.47, -1.15, -1.15),
    (14, -6.22, 31.47, -44.57, -1.78, -1.78),
    (16, -5.61, 31.47, -51.36, -1.47, -1.47),
    (18, -5.64, 7.50, -47.42, -1.28, -1.28),
    (20, -7.23, -1.34, -47.56, -1.29, -1.29),
)
_XSHORE3_WAVES = ((8.0, 8.0, 10.0), (5.0, 5.0, 45.0), (3.0, 4.0, 30.0), (12.0, 10.0, 0.0))
_XSHORE3_WL = (0.0, 3.0, 5.0, 7.0, 5.0, 2.0)
_XSHORE3_TIMES = (1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0)
# XSHORE4: the 53-point surveyed profile that ships as this application's default
_SRC_XSHORE4 = (
    (2, -1.68, -17.54, -2.41, 0.00, 0.00),
    (4, -2.76, -34.67, -3.98, 0.00, 0.00),
    (6, -3.50, -48.52, -5.09, 0.00, 0.00),
    (8, -4.07, -58.17, -5.96, -0.18, 0.00),
    (10, -4.54, -64.89, -6.69, -0.84, 0.00),
    (12, -4.94, -69.85, -7.29, -1.39, 0.00),
    (14, -5.28, -73.61, -7.79, -1.86, 0.00),
    (16, -5.59, -76.76, -8.25, -2.33, 0.00),
    (18, -5.85, -79.29, -8.65, -2.73, 0.00),
    (20, -6.07, -81.33, -8.99, -3.09, 0.00),
)
_YD3_PER_M3 = 0.764554858 / _FT       # yd^3/ft -> m^3/m


def _check_deck(name, inp, table):
    """Compare one run against the source's table; returns the worst residuals."""
    r = compute(inp)
    got = {round(float(h), 6): k for k, h in enumerate(r.hours)}
    worst_v = worst_c = 0.0
    for hr, vol, ca, cb, cc, cd in table:
        k = got.get(float(hr))
        assert k is not None, (name, hr, list(got))
        worst_v = max(worst_v, abs(r.vol_change[k] / _YD3_PER_M3 - vol))
        for series, exp in ((r.contour_0, ca), (r.contour_5, cb),
                            (r.contour_10, cc), (r.contour_15, cd)):
            worst_c = max(worst_c, abs(series[k] / _FT - exp))
    # the source prints two decimals, so half a digit is 0.005 in each column
    assert worst_v <= 0.006, (name, "volume", worst_v)
    assert worst_c <= 0.006, (name, "contour", worst_c)
    return worst_v, worst_c


def _self_tests() -> None:
    # 1) the generic profile builder reproduces the source's own listing
    x, z, dpts, _sd, _sb = generic_profile(20.0, 50.0, 1 / 2.0, 6.0, 6.0, 100.0,
                                           1 / 10.0, 1 / 20.0)
    assert dpts == 287, dpts
    assert _approx(z[1], 20.0, 1e-9) and _approx(z[13], 20.0, 1e-9), (z[1], z[13])
    assert _approx(z[20], 6.0, 1e-6), z[20]        # dune face lands on the berm
    assert _approx(x[dpts], 1144.0, 1e-9), x[dpts]
    assert _approx(z[dpts], -45.2, 1e-6), z[dpts]

    # 2) the shape factor. ACES's lower branch agrees with Dean's tabulated
    #    A ~ 0.100 m^(1/3) at D50 = 0.20 mm to about 3 percent; its upper branch runs
    #    well above the table (0.156 against 0.125 m^(1/3) at 0.40 mm), which is ACES's
    #    own fit and is reproduced here rather than corrected.
    assert _approx(shape_factor(0.20) * _FT ** (1.0 / 3.0), 0.100, 0.005),         shape_factor(0.20) * _FT ** (1.0 / 3.0)
    assert shape_factor(0.27) > shape_factor(0.22) > shape_factor(0.20)

    # 3) the breaking point, against the source's own iteration on the default case
    hb, db = breaking_point(8.0, 8.0, 10.0, 60.0)
    assert _approx(hb, 8.73463, 5e-5), hb
    assert _approx(db, 11.19824, 5e-5), db

    # 4) the march, against the source on three of the decks it ships with
    base = {f.key: f.default for f in INPUTS}
    generic = dict(base, profile_type=_GENERIC)
    w2, c2 = _check_deck("XSHORE2", generic, _SRC_XSHORE2)

    x3 = dict(generic,
              waves=tuple((h * _FT, t, a) for h, t, a in _XSHORE3_WAVES),
              wave_dt_hr=5.0, water_level=tuple(v * _FT for v in _XSHORE3_WL),
              wl_dt_hr=4.0, out_times=_XSHORE3_TIMES)
    w3, c3 = _check_deck("XSHORE3", x3, _SRC_XSHORE3)

    x4 = dict(base, profile_type=_SURVEYED, dune_elev=14.1 * _FT)
    w4, c4 = _check_deck("XSHORE4", x4, _SRC_XSHORE4)

    wv = max(w2, w3, w4)
    wc = max(c2, c3, c4)
    print(f"  self-tests: PASS (generic profile {dpts} points to 1144 ft; breaking "
          f"H_b={hb:.4f} ft in {db:.4f} ft; XSHORE2/3/4 volume within {wv:.4f} yd^3/ft "
          f"and contours within {wc:.4f} ft of the source, print resolution 0.005)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = the ACES deck XSHORE2.IN: generic 20 ft dune on a 6 ft berm,")
    print("   0.22 mm sand, one 8 ft / 8 s wave record at 10 deg in 60 ft, 20 hours)")
    _YD3_FT = 0.764554858 / _FT
    print(f"    grid: {len(r.profile_X)} points to {r.profile_X[-1] / _FT:.0f} ft;  "
          f"A = {r.A_shape / _FT ** (1/3):.4f} ft^(1/3)")
    print(f"    {'hour':>6} {'volume change':>16} {'0 ft':>9} {'+5 ft':>9} "
          f"{'+10 ft':>9} {'+15 ft':>9}")
    print(f"    {'':>6} {'(yd^3/ft)':>16} {'(ft)':>9} {'(ft)':>9} {'(ft)':>9} {'(ft)':>9}")
    for k in range(len(r.hours)):
        print(f"    {r.hours[k]:6.0f} {r.vol_change[k] / _YD3_FT:16.2f} "
              f"{r.contour_0[k] / _FT:9.2f} {r.contour_5[k] / _FT:9.2f} "
              f"{r.contour_10[k] / _FT:9.2f} {r.contour_15[k] / _FT:9.2f}")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
