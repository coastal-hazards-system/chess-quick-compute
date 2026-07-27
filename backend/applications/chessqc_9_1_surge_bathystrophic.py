"""CHESS-QC application — Bathystrophic Storm Surge (Bodine 1971).

First Quick Compute tool beyond the original 34 ACES applications (functional area:
Storm Surge). Estimates open-coast hurricane surge along a single cross-shelf traverse
by the quasi-1D bathystrophic method of Bodine (1971), CERC TM-35 (verification:
Pararas-Carayannis 1975, TM-50).

Classification: standard for the published Bodine TM-35 numerical method. The original
appendix supplies the Chesapeake 17-station bathymetry and the radius, wind-speed, and wind-
direction curve ordinates used by its FORTRAN calculation; they are encoded below in the
source-curve option. Holland and Myers remain optional parametric sensitivity paths. The
bathystrophic approximation itself remains explicitly screening-only, as Bodine explains
its restricted geometry/quadrant validity and need for local calibration.

Method (Bodine 1971): the surge is built by integrating, along a traverse from the
shelf edge to the shore as the storm passes, the reduced vertically-integrated
equations of motion:
  setup      dSx/dx = k W^2 cos(theta) / (g D)            (onshore wind setup, eq 15/23)
             dSy/dx = f V / (g D)                          (bathystrophic/Coriolis, eq 24)
  transport  dV/dt  = k W^2 sin(theta) - K V|V| / D^2      (alongshore flux, eq 16)
solved by the finite-difference analogs (eqs 25, 26, 33) with the flux limiter (36).
The total still-water rise at the shore is the composite (eq 34):
  S = Sx + Sy + Se(initial rise) + SA(astronomical tide) + S_dp(pressure) + Sw(wave setup).
Stresses: wind tau_s = rho k W^2 (Van Dorn 1953, k = K1 + K2(1-Wc/W)^2, W>=Wc=14 kt,
x WKCOR=1.1); bottom tau_b/rho = K V|V|/D^2 (K ~ 0.0025).

Wind field, three selectable paths:
  * Bodine TM-35 source curves: the original appendix curve ordinates and its English-unit
    finite-difference form; this reproduces the documented Chesapeake reference calculation.
  * Holland (1980): p(r) = Pc + dP exp(-(R/r)^B);  default, modern standard.
  * Myers (1954) / Bodine: p(r) = Pc + dP exp(-R/r)  == Holland with B = 1.
  Gradient wind  V_gr(r) = sqrt( (B dP/rho_a)(R/r)^B exp(-(R/r)^B) + (r f/2)^2 ) - r f/2,
  reduced to the surface (~0.9), given a forward-speed asymmetry and an inflow angle,
  then split into onshore/alongshore components that drive the integrator. The Holland
  shape factor B (default 1.5) may be overridden by an explicit Vmax via
  B = rho_a e Vmax^2 / dP.

Validation: Bodine TM-35 Chesapeake Bay Entrance example (Pc=27.57, Pn=29.92 inHg,
R=35 nm, Vf=22 kt, lat 37 deg, K=0.0025) -> reported peak 13.41 ft at 17 h. The source-
curve path is checked against that published result (within the precision of its rounded curve
ordinates); the parametric paths
are checked separately by their analytic wind/pressure identities. self-contained, numpy +
stdlib only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

G_SI = 9.80665             # m/s^2
RHO_W = 1025.0            # kg/m^3 sea water
RHO_AIR = 1.20            # kg/m^3
OMEGA = 7.2921159e-5      # rad/s (earth rotation)
_E = math.e


@dataclass(frozen=True)
class AppMeta:
    aces_id: str; name: str; area: str; classification: str; cite: str; default_system: str = "SI"
    status: str = "Current"          # Current | Screening only | Superseded
    superseded_by: str = ""          # newer method, if any (surfaced in the docs)


@dataclass(frozen=True)
class Field:
    key: str; label: str; kind: str = "float"; unit_si: str = ""; unit_us: str = ""
    default: object = 0.0; lo: float = -math.inf; hi: float = math.inf
    choices: tuple = (); columns: tuple = (); note: str = ""
    enable_if: tuple = ()    # (other_key, value): gray out (disable) unless that input == value


@dataclass(frozen=True)
class Out:
    key: str; label: str; unit_si: str = ""; unit_us: str = ""; kind: str = "scalar"
    note: str = ""           # hover definition shown on the output label


APP_META = AppMeta(
    aces_id="9-1",
    name="Bathystrophic Storm Surge",
    area="Storm Surge",
    classification="standard",
    cite="Bodine (1971) TM-35; Holland (1980); Myers (1954); TR/CERC",
    default_system="US",
    status="Screening only",
    superseded_by="ADCIRC (risk assessment)",
)

_WIND_MODELS = ("Bodine TM-35 source curves", "Holland (1980)", "Myers / Bodine (1954)")
_NM, _FT, _KT, _MPH, _MI = 1852.0, 0.3048, 0.514444, 0.44704, 1609.344
_INHG = 3386.389

# Default = Bodine TM-35 Chesapeake Bay Entrance example (US units; stored SI).
# Bathymetry: distance from shore (nm) vs depth below SWL (ft), shelf edge -> shore.
_BATHY = [(62.0, 600.0), (60.0, 350.0), (55.0, 240.0), (50.0, 178.0),
          (45.0, 154.0), (40.0, 130.0), (35.0, 115.0), (30.0, 102.0),
          (25.0, 84.0), (20.0, 71.0), (15.0, 64.0), (10.0, 56.0),
          (5.0, 44.0), (3.5, 40.0), (2.0, 36.0), (1.0, 33.0), (0.0, 30.0)]

# Bodine TM-35 appendix, p. 52: curve abscissae and ordinates. The original program
# linearly interpolates them at 1-nmi spacing before its time/space march.
_TM35_XR = (0, 252, 290, 320, 333, 345, 355, 365, 370, 375, 378, 380, 385,
            390, 400, 410, 420, 433, 450, 510, 765)
_TM35_R = (377, 129, 92.5, 65.5, 54, 45, 39, 34, 32.5, 32, 31.5, 32, 32.5,
           35, 40, 46, 53, 64, 80, 137, 389)
_TM35_XW = (0, 90, 140, 180, 210, 252, 270, 290, 300, 310, 320, 340, 345,
            360, 369, 380, 390, 402, 410, 420, 433, 450, 472, 490, 540,
            598, 669, 765)
_TM35_W = (0, 10, 17, 24, 30, 40, 46, 53, 57.5, 63, 70, 86.5, 90, 97, 100,
           102, 102, 100, 97, 90, 80, 70, 60, 54, 41, 30, 20, 10)
_TM35_XT = (0, 190, 230, 260, 280, 300, 310, 320, 330, 340, 350, 370, 380,
            390, 395, 400, 405, 410, 430, 450, 470, 490, 510, 560, 640, 765)
_TM35_THETA = (125, 123, 121, 118, 115.5, 111, 107, 102, 95.5, 87, 77, 50,
               35, 15, 0, 338, 330.5, 326.5, 316.5, 308, 301, 296, 293,
               288.8, 285, 280)
# Continuous equivalent used by the source program's special 360-degree crossing logic.
_TM35_THETA_CONT = _TM35_THETA[:15] + tuple(v - 360.0 for v in _TM35_THETA[15:])

INPUTS = (
    Field("bathy", "Traverse bathymetry", "table",
          default=[[d * _NM, z * _FT] for d, z in _BATHY],
          columns=(("Distance from shore", "km", "nm"), ("Depth below SWL", "m", "ft")),
          note="shelf edge -> shore; one row per station"),
    Field("Pc", "Central pressure", "float", "hPa", "inHg", default=27.57 * _INHG,
          lo=800 * 100.0, hi=1020 * 100.0, note="storm central pressure"),
    Field("Pn", "Peripheral pressure", "float", "hPa", "inHg", default=29.92 * _INHG,
          lo=980 * 100.0, hi=1030 * 100.0, note="ambient/peripheral pressure"),
    Field("R", "Radius of maximum winds", "float", "km", "nm", default=35.0 * _NM,
          lo=2e3, hi=2e5, note="> 0"),
    Field("Vf", "Forward speed", "float", "km/h", "kt", default=22.0 * _KT,
          lo=0.0, hi=40.0, note="storm translation speed"),
    Field("track_offset", "Track offset from traverse", "float", "km", "nm",
          default=35.0 * _NM, lo=0.0, hi=5e5,
          note="alongshore distance from the traverse to the storm landfall/track"),
    Field("lat", "Latitude", "angle", "deg", "deg", default=37.0, lo=0.0, hi=80.0,
          note="for the Coriolis parameter"),
    Field("wind_model", "Wind model", "choice", "", "", default="Bodine TM-35 source curves",
          choices=_WIND_MODELS, note="Bodine's published source curves, or Holland/Myers parametric sensitivity paths"),
    Field("B_holland", "Holland B (peakedness)", "float", "", "", default=1.5,
          lo=0.5, hi=2.5, note="Holland shape factor; locked to 1.0 for Myers/Bodine",
          enable_if=("wind_model", "Holland (1980)")),
    Field("Vmax", "Max wind (optional)", "float", "km/h", "mph", default=0.0,
          lo=0.0, hi=120.0, note="if > 0, overrides B via B = rho_a e Vmax^2 / dP",
          enable_if=("wind_model", "Holland (1980)")),
    Field("rho_air", "Air density", "float", "kg/m^3", "kg/m^3", default=RHO_AIR,
          lo=1.0, hi=1.3, note="ambient air density"),
    Field("K_bottom", "Bottom friction coefficient", "float", "", "", default=0.0025,
          lo=1e-4, hi=2e-2, note="bed friction K (~0.002-0.005)"),
    Field("Se", "Initial water-level rise", "float", "m", "ft", default=0.5 * _FT,
          lo=-3.0, hi=10.0, note="initial setup at start of computation"),
    Field("SA", "Astronomical tide", "float", "m", "ft", default=2.5 * _FT,
          lo=-5.0, hi=10.0, note="astronomical tide above MSL datum"),
    Field("dt", "Time step", "float", "hr", "hr", default=0.5 * 3600.0,
          lo=60.0, hi=7200.0, note="integration time step (stored in seconds)"),
    Field("n_steps", "Number of time steps", "int", "", "", default=62, lo=10, hi=400,
          note="storm is swept past the traverse over these steps"),
)

OUTPUTS = (
    Out("peak_surge", "Peak surge at shore", "m", "ft", "scalar",
        note="Peak total still-water rise at the shore, summing the wind, bathystrophic, pressure, initial-rise and tide components (eq 34)."),
    Out("S_wind", "  wind setup (Sx)", "m", "ft", "scalar",
        note="Onshore wind-stress setup component Sx at the shore at the time of peak surge (eq 25)."),
    Out("S_bathy", "  bathystrophic setup (Sy)", "m", "ft", "scalar",
        note="Bathystrophic (Coriolis) setup component Sy at the shore, from the alongshore wind-driven transport (eq 26)."),
    Out("S_press", "  pressure setup", "m", "ft", "scalar",
        note="Inverse-barometer pressure setup at the shore, from the local pressure deficit relative to the peripheral pressure."),
    Out("Vmax_out", "Max wind speed (30 ft)", "km/h", "mph", "scalar",
        note="Maximum surface (30 ft) wind speed at the radius of maximum winds, including half the storm forward speed."),
    Out("B_used", "Holland B used", "", "", "scalar",
        note="Holland shape (peakedness) factor actually used in the wind field; 1.0 for Myers/Bodine or when set by an explicit Vmax."),
    Out("t_peak", "Time of peak", "s", "hr", "scalar",
        note="Elapsed time from the start of the storm sweep at which the peak shore surge occurs."),
    Out("profile_X", "Profile: distance from shore", "m", "nm", "profile",
        note="Cross-shelf distance from the shore for each traverse station, ordered shore to shelf edge."),
    Out("profile_eta", "Profile: setup at peak", "m", "ft", "profile",
        note="Total still-water setup along the traverse at the time of peak surge (wind plus bathystrophic plus pressure plus initial rise and tide)."),
    Out("profile_u", "Profile: still-water depth", "m", "ft", "profile",
        note="Undisturbed still-water depth below SWL at each traverse station."),
    Out("profile_w", "Profile: total depth at peak", "m", "ft", "profile",
        note="Total water depth at peak surge, the still-water depth plus the setup at each traverse station."),
)


@dataclass
class Result:
    peak_surge: float; S_wind: float; S_bathy: float; S_press: float
    Vmax_out: float; B_used: float; t_peak: float
    profile_X: np.ndarray; profile_eta: np.ndarray
    profile_u: np.ndarray; profile_w: np.ndarray
    notes: str = ""


# --- parametric wind field (Holland 1980; Myers 1954 == Holland B=1) ------------
def _gradient_wind(r, R, dP, B, f, rho_a):
    """Gradient wind speed (m/s) at radius r (m). Holland (1980); Myers = B=1."""
    r = max(r, 1.0)
    x = (R / r) ** B
    inside = (B * dP / rho_a) * x * math.exp(-x) + (r * f / 2.0) ** 2
    return math.sqrt(max(inside, 0.0)) - r * f / 2.0


def _pressure_at(r, R, Pc, dP, B):
    """Surface pressure (Pa) at radius r, Holland/Myers profile."""
    return Pc + dP * math.exp(-((R / max(r, 1.0)) ** B))


def _validate(inp):
    for fdef in INPUTS:
        if fdef.kind not in ("float", "int", "angle") or fdef.key not in inp:
            continue
        v = float(inp[fdef.key])
        if not (fdef.lo <= v <= fdef.hi):
            raise ValueError(f"{fdef.label} ({fdef.key}) = {v} outside [{fdef.lo:g}, {fdef.hi:g}]")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Estimates open-coast hurricane storm surge along a single cross-shelf '
            "traverse using Bodine's (1971) quasi-1D bathystrophic method, integrating the "
            'vertically-averaged equations of motion as a parametric storm sweeps past. '
            'Returns the peak still-water rise at the shore and its wind, bathystrophic, '
            'and pressure components.',
 'method_key': 'wind_model',
 'methods': [{'name': 'Bodine TM-35 source-curve wind field',
              'when': 'Bodine TM-35 source curves',
              'tag': 'source reference',
              'note': 'Replays the original TM-35 appendix curve ordinates for radius, '
                      '30-ft wind speed, and wind angle, with Bodine\'s published '
                      'finite-difference equations and English-unit constants. Use this '
                      'for the documented Chesapeake reference case; use a parametric path '
                      'or externally supplied forcing for another storm.',
              'equations': [{'tex': 'S_p = 1.14(P_n-P_c)\\left(1-e^{-R/r}\\right)',
                             'desc': 'Myers/Bodine pressure setup in feet when pressure is '
                                     'in inches Hg (TM-35 eq. 33).'},
                            {'tex': 'V^{n+1}=\\frac{V^n+\\frac{\\Delta t}{2}(B^n+B^{n+1})}'
                                    '{1+(C_3/D)^2K|V^n|\\Delta t}',
                             'desc': 'Bodine\'s finite-difference alongshore-transport '
                                     'update, with the source flux limiter.'},
                            {'tex': 'S=S_x+S_y+S_e+S_A+S_p',
                             'desc': 'TM-35 composite still-water level (no wave setup at '
                                     'the Chesapeake Bay mouth reference).'}]},
             {'name': 'Holland (1980) wind field',
              'when': 'Holland (1980)',
              'tag': 'preferred',
              'note': 'Modern standard; shape factor B (default 1.5) adjustable, or set by '
                      'an explicit Vmax via B = rho_a e Vmax^2 / dP.',
              'equations': [{'tex': 'p(r) = P_c + \\Delta P \\, '
                                    '\\exp\\left(-\\left(\\frac{R}{r}\\right)^{B}\\right)',
                             'desc': 'Holland radial surface-pressure profile (deficit dP '
                                     '= Pn - Pc).'},
                            {'tex': 'V_{gr}(r) = \\sqrt{\\frac{B \\, \\Delta '
                                    'P}{\\rho_a}\\left(\\frac{R}{r}\\right)^{B}\\exp\\left(-\\left(\\frac{R}{r}\\right)^{B}\\right) '
                                    '+ \\left(\\frac{r f}{2}\\right)^{2}} - \\frac{r f}{2}',
                             'desc': 'Gradient-wind speed; surface-reduced and split into '
                                     'onshore/alongshore components.'},
                            {'tex': '\\frac{dS_x}{dx} = \\frac{k \\, W^{2} \\cos\\theta}{g '
                                    'D}',
                             'desc': 'Onshore wind setup gradient (eq 15/23).'},
                            {'tex': '\\frac{dS_y}{dx} = \\frac{f V}{g D}',
                             'desc': 'Bathystrophic (Coriolis) setup gradient from '
                                     'alongshore transport (eq 24).'},
                            {'tex': '\\frac{dV}{dt} = k \\, W^{2} \\sin\\theta - \\frac{K '
                                    'V |V|}{D^{2}}',
                             'desc': 'Alongshore transport balance: wind stress minus '
                                     'bottom friction (eq 16).'},
                            {'tex': 'S = S_x + S_y + S_e + S_A + S_{dp} + S_w',
                             'desc': 'Composite still-water rise at the shore (eq 34).'}]},
             {'name': 'Myers (1954) / Bodine wind field',
              'when': 'Myers / Bodine (1954)',
              'tag': 'legacy',
              'note': 'Holland profile with B locked to 1; retained to reproduce the '
                      'original Bodine TM-35 graphical examples.',
              'equations': [{'tex': 'p(r) = P_c + \\Delta P \\, '
                                    '\\exp\\left(-\\frac{R}{r}\\right)',
                             'desc': 'Myers exponential pressure profile (Holland with B = '
                                     '1).'},
                            {'tex': 'V_{gr}(r) = \\sqrt{\\frac{\\Delta '
                                    'P}{\\rho_a}\\left(\\frac{R}{r}\\right)\\exp\\left(-\\frac{R}{r}\\right) '
                                    '+ \\left(\\frac{r f}{2}\\right)^{2}} - \\frac{r f}{2}',
                             'desc': 'Gradient-wind speed for B = 1.'},
                            {'tex': '\\frac{dS_x}{dx} = \\frac{k \\, W^{2} \\cos\\theta}{g '
                                    'D}',
                             'desc': 'Onshore wind setup gradient (eq 15/23).'},
                            {'tex': '\\frac{dS_y}{dx} = \\frac{f V}{g D}',
                             'desc': 'Bathystrophic (Coriolis) setup gradient (eq 24).'},
                            {'tex': '\\frac{dV}{dt} = k \\, W^{2} \\sin\\theta - \\frac{K '
                                    'V |V|}{D^{2}}',
                             'desc': 'Alongshore transport balance (eq 16).'},
                            {'tex': 'S = S_x + S_y + S_e + S_A + S_{dp} + S_w',
                             'desc': 'Composite still-water rise at the shore (eq 34).'}]}],
 'symbols': [['S', 'Total still-water surge at the shore'],
             ['S_x', 'Onshore wind setup component'],
             ['S_y', 'Bathystrophic (Coriolis) setup component'],
             ['S_{dp}', 'Pressure (inverse-barometer) setup component'],
             ['W', 'Surface wind speed driving the stress'],
             ['theta', 'Wind direction relative to the traverse (onshore vs alongshore)'],
             ['V', 'Vertically-integrated alongshore transport (flux)'],
             ['D', 'Total water depth (still-water depth plus setup)'],
             ['k', 'Surface wind-stress coefficient (Van Dorn 1953)'],
             ['K', 'Bottom friction coefficient (~0.0025)'],
             ['B', 'Holland shape (peakedness) factor; B = 1 for Myers/Bodine'],
             ['dP', 'Central pressure deficit, Pn - Pc'],
             ['f', 'Coriolis parameter, 2 omega sin(lat)'],
             ['R', 'Radius of maximum winds']],
 'references': ['Bodine (1971) TM-35',
                'Pararas-Carayannis (1975) TM-50',
                'Holland (1980)',
                'Myers (1954)',
                'Van Dorn (1953)']}


def _interp_source(x: float, xp: tuple, fp: tuple) -> float:
    """TM-35 integer-table lookup after its linear interpolation of source ordinates.

    FORTRAN stores the x=0 ordinate at array index 1, then assigns LD=DI by integer
    truncation. Thus lookup index LD corresponds to curve coordinate LD-1.
    """
    x = max(float(math.trunc(x) - 1), 0.0)
    if x <= xp[0]:
        return float(fp[0])
    if x >= xp[-1]:
        # The published curves end in the remote, low-forcing portion of the storm.
        # Do not extend a graphical curve beyond its documented support.
        return 0.0 if fp is _TM35_W else float(fp[-1])
    return float(np.interp(x, xp, fp))


def _compute_tm35_source(inp: dict) -> Result:
    """Bodine TM-35 Appendix source-curve calculation (original English-unit scheme)."""
    rows = sorted(([float(c[0]) / _NM, float(c[1]) / _FT] for c in inp["bathy"] if c), reverse=True)
    if len(rows) < 3:
        raise ValueError("need at least 3 bathymetry stations")
    x_nm = np.array([r[0] for r in rows]); depth_ft = np.array([max(r[1], 1.0) for r in rows])
    N = len(x_nm)
    dt_hr = float(inp["dt"]) / 3600.0; n_steps = int(inp["n_steps"])
    vf_kt = float(inp["Vf"]) / _KT
    Pc = float(inp["Pc"]) / _INHG; Pn = float(inp["Pn"]) / _INHG
    R_nm = float(inp["R"]) / _NM; K = float(inp.get("K_bottom", 0.0025))
    Se = float(inp.get("Se", 0.0)) / _FT; SA = float(inp.get("SA", 0.0)) / _FT
    sinphi = math.sin(math.radians(float(inp["lat"])))
    C1, C2, C3, WKCOR = 203.0, 106.0, 5280.0, 1.10  # TM-35 card 7

    prev_v = np.zeros(N - 1); prev_bn = np.zeros(N - 1)
    # Source arrays SXP/SYP/SPP are reach values (after each reach), not station values.
    prev_sx = np.zeros(N - 1); prev_sy = np.zeros(N - 1); prev_sp = np.zeros(N - 1)
    peak = -math.inf; best = None
    for n in range(n_steps):
        hour = n * dt_hr
        r_nm = np.empty(N); wind_mph = np.empty(N); wind_raw = np.empty(N); theta = np.empty(N)
        for i, dist in enumerate(x_nm):
            # TM-35 program listing, p. 48: DISTC=HOUR*VF; DI=DIST+DISTC+1.0 nmi.
            coordinate = dist + hour * vf_kt + 1.0
            r_nm[i] = _interp_source(coordinate, _TM35_XR, _TM35_R)
            wind_raw[i] = _interp_source(coordinate, _TM35_XW, _TM35_W)
            wind_mph[i] = wind_raw[i]
            theta[i] = math.radians(_interp_source(coordinate, _TM35_XT, _TM35_THETA_CONT))
            if dist < 1.0:
                wind_mph[i] *= 0.89
            elif abs(dist - 1.0) < 1e-12:
                wind_mph[i] *= 0.945
        sp = 1.14 * (Pn - Pc) * (1.0 - np.exp(-R_nm / np.maximum(r_nm, 1e-9)))
        wk = np.full(N, 1.1e-6)
        # The listing deliberately evaluates Van Dorn's coefficient with unreduced W(LD),
        # while WWX/WWY use the land-interference-reduced WIND value.
        mask = wind_raw > 16.0
        wk[mask] += 2.5e-6 * (1.0 - 16.0 / wind_raw[mask]) ** 2
        wwx = wind_mph * wind_mph * np.cos(theta)
        wwy = wind_mph * wind_mph * np.sin(theta)
        sx_reach = np.zeros(N - 1); sy_reach = np.zeros(N - 1)
        sp_reach = np.zeros(N - 1); new_v = np.zeros(N - 1); new_bn = np.zeros(N - 1)
        sum_sx = 0.0; sum_sy = 0.0
        for i in range(N - 1):
            dx = x_nm[i] - x_nm[i + 1]
            an = wk[i] * (wwx[i] + wwx[i + 1]) * WKCOR
            bn = wk[i] * (wwy[i] + wwy[i + 1]) * 0.5 * WKCOR
            spn = 0.5 * (sp[i] + sp[i + 1])
            if n == 0:  # TM-35 lines 445-450 initialization
                sxp = syp = 0.0
                spp, vp, bp, sap = spn, 0.0, bn, SA
            else:
                sxp, syp = prev_sx[i], prev_sy[i]
                spp, vp, bp, sap = prev_sp[i], prev_v[i], prev_bn[i], SA
            dts = 0.5 * (depth_ft[i] + depth_ft[i + 1]) + Se + sxp + syp
            dtn = dts + SA + spn
            dth = dts + 0.5 * (SA + sap) + 0.5 * (spp + spn)
            dtn = max(dtn, 1.0); dth = max(dth, 1.0)
            da = C3 / dth
            vn = (vp + 0.5 * (bn + bp) * dt_hr) / (1.0 + da * da * K * abs(vp) * dt_hr)
            limit = math.sqrt(abs(bn) / (K * da * da)) if bn else 0.0
            vn = max(-limit, min(limit, vn))
            new_v[i], new_bn[i] = vn, bn
            sum_sx += C1 * dx * an / dtn
            sum_sy += C2 * dx * (2.0 * sinphi) * vn / dtn
            sx_reach[i], sy_reach[i], sp_reach[i] = sum_sx, sum_sy, spn
        total = sum_sx + sum_sy + Se + SA + sp_reach[-1]
        if total > peak:
            peak = total
            best = (n + 1, sx_reach.copy(), sy_reach.copy(), sp_reach.copy())
        prev_v, prev_bn = new_v, new_bn
        prev_sx, prev_sy, prev_sp = sx_reach, sy_reach, sp_reach

    step, sx_reach, sy_reach, sp_reach = best
    # Convert reach-end cumulative quantities to a station profile; the seaward
    # boundary has zero wind/bathystrophic setup and its local pressure setup.
    sx = np.concatenate(([0.0], sx_reach)); sy = np.concatenate(([0.0], sy_reach))
    sp = np.concatenate(([sp_reach[0]], sp_reach))
    eta = sx + sy + Se + SA + sp
    return Result(peak_surge=peak * _FT, S_wind=sx_reach[-1] * _FT,
                  S_bathy=sy_reach[-1] * _FT, S_press=sp_reach[-1] * _FT,
                  Vmax_out=102.0 * _MPH, B_used=1.0,
                  t_peak=step * dt_hr * 3600.0, profile_X=x_nm[::-1] * _NM,
                  profile_eta=eta[::-1] * _FT, profile_u=depth_ft[::-1] * _FT,
                  profile_w=(depth_ft + eta)[::-1] * _FT,
                  notes=("Bodine TM-35 Appendix source curves and original finite-difference "
                         f"scheme; peak at {step * dt_hr:.1f} h"))


def compute(inp: dict, *, g: float = G_SI, rho_w: float = RHO_W) -> Result:
    """Bathystrophic surge along a traverse. SI inputs; the GUI converts at the edge."""
    _validate(inp)
    if str(inp.get("wind_model", "Bodine TM-35 source curves")) == "Bodine TM-35 source curves":
        return _compute_tm35_source(inp)
    rows = sorted(([float(c[0]), float(c[1])] for c in inp["bathy"] if c), reverse=True)
    if len(rows) < 3:
        raise ValueError("need at least 3 bathymetry stations")
    X = np.array([r[0] for r in rows])           # distance from shore (m), seaward->shore
    d = np.array([max(r[1], 1.0) for r in rows])  # depth below SWL (m), >0
    N = len(X)

    Pc, Pn = float(inp["Pc"]), float(inp["Pn"])
    dP = max(Pn - Pc, 1.0)                        # pressure deficit (Pa)
    R = float(inp["R"]); Vf = float(inp["Vf"]); Y0 = float(inp["track_offset"])
    lat = math.radians(float(inp["lat"]))
    f = 2.0 * OMEGA * math.sin(lat)
    rho_a = float(inp.get("rho_air", RHO_AIR))
    K = float(inp.get("K_bottom", 0.0025))
    Se = float(inp.get("Se", 0.0)); SA = float(inp.get("SA", 0.0))
    dt = float(inp["dt"]); n_steps = int(inp["n_steps"])
    model = str(inp.get("wind_model", "Holland (1980)"))

    # Holland B: Myers locks B=1; an explicit Vmax overrides B
    if model.startswith("Myers"):
        B = 1.0
    else:
        B = float(inp.get("B_holland", 1.5))
        Vmax_in = float(inp.get("Vmax", 0.0))
        if Vmax_in > 0.0:
            B = rho_a * _E * Vmax_in * Vmax_in / dP
    B = min(max(B, 0.3), 3.0)

    WKCOR, BETA, RED = 1.1, math.radians(22.0), 0.865  # stress corr, inflow angle, surface reduction (Bodine eq 30)
    Wc = 14.0 * _KT

    def van_dorn_k(W):
        if W <= Wc:
            return 1.1e-6 * WKCOR
        return (1.1e-6 + 2.5e-6 * (1.0 - Wc / W) ** 2) * WKCOR

    # storm sweeps shoreward along a track offset Y0 alongshore; center cross-shore
    # position X_c(t) goes from well seaward to past landfall.
    travel = Vf * dt * n_steps
    Xc0 = 0.65 * travel + X[0]                    # start seaward of the shelf edge
    Wg_max = RED * _gradient_wind(R, R, dP, B, f, rho_a)   # surface wind at RMW

    def wind_at(Xi, Xc):
        """Onshore & alongshore wind components (m/s) at traverse station Xi."""
        dx, dy = (Xi - Xc), -Y0                   # vector storm-center -> station (seaward+, along)
        r = math.hypot(dx, dy)
        Wg = RED * _gradient_wind(r, R, dP, B, f, rho_a)
        if r < 1.0:
            return 0.0, 0.0, 0.0
        # cyclonic (CCW) tangential + inflow toward center
        tx, ty = dy / r, dx / r                   # tangential CCW (rotate radial +90)
        ix, iy = -dx / r, -dy / r                 # inward radial
        wx = Wg * (math.cos(BETA) * tx + math.sin(BETA) * ix)
        wy = Wg * (math.cos(BETA) * ty + math.sin(BETA) * iy)
        # forward-speed asymmetry: storm motion (shoreward, -X) added in proportion to the
        # local rotational strength so it vanishes away from the storm (Bodine adds 0.5 Vf
        # to the peak wind).
        wx += -0.5 * Vf * (Wg / Wg_max)
        Wmag = math.hypot(wx, wy)
        return -wx, wy, Wmag                      # onshore(+ toward shore = -X), alongshore, |W|

    V = np.zeros(N)                               # alongshore transport at each reach (m^2/s)
    peak = -1e9
    best = None
    for n in range(n_steps + 1):
        Xc = Xc0 - Vf * dt * n
        Won = np.zeros(N); Wal = np.zeros(N); Sdp = np.zeros(N)
        for i in range(N):
            won, wal, wm = wind_at(X[i], Xc)
            kk = van_dorn_k(wm)
            Won[i] = kk * wm * won               # A = k W^2 cos(theta)  (signed onshore)
            Wal[i] = kk * wm * wal               # B = k W^2 sin(theta)  (signed alongshore)
            r_i = math.hypot(X[i] - Xc, Y0)
            Sdp[i] = (Pn - _pressure_at(r_i, R, Pc, dP, B)) / (rho_w * g)   # inverse barometer
        # march shoreward, accumulating setup; update transport (semi-implicit)
        Sx = np.zeros(N); Sy = np.zeros(N)
        for i in range(N - 1):
            dxr = X[i] - X[i + 1]                 # reach length (m) > 0
            D = 0.5 * (d[i] + d[i + 1]) + 0.5 * (Sx[i] + Sy[i] + Sdp[i] + Se + SA)
            D = max(D, 1.0)
            Bavg = 0.5 * (Wal[i] + Wal[i + 1])
            Vn = (V[i] + Bavg * dt) / (1.0 + K * abs(V[i]) * dt / (D * D))
            lim = math.sqrt(D * D * abs(Bavg) / K) if Bavg != 0 else 1e30   # eq 36 flux limit
            Vn = max(min(Vn, lim), -lim)
            V[i] = Vn
            Aavg = 0.5 * (Won[i] + Won[i + 1])
            Sx[i + 1] = Sx[i] + dxr * Aavg / (g * D)            # eq 25
            Sy[i + 1] = Sy[i] + dxr * f * Vn / (g * D)          # eq 26
        S_total = Sx[-1] + Sy[-1] + Sdp[-1] + Se + SA
        if S_total > peak:
            peak = S_total
            best = (n * dt, Sx.copy(), Sy.copy(), Sdp.copy(), Sx[-1], Sy[-1], Sdp[-1])

    t_peak, Sxp, Syp, Sdpp, Sxs, Sys, Sdps = best
    Vmax_out = RED * _gradient_wind(R, R, dP, B, f, rho_a) + 0.5 * Vf
    eta = Sxp + Syp + Sdpp + Se + SA
    Dtot = d + eta
    notes = [f"{model}, B={B:.2f}; dP={dP/_INHG:.2f} inHg; peak at t={t_peak/3600:.1f} h",
             f"components at shore (ft): wind {Sxs/_FT:.2f}, bathystrophic {Sys/_FT:.2f}, "
             f"pressure {Sdps/_FT:.2f}"]
    return Result(peak_surge=peak, S_wind=Sxs, S_bathy=Sys, S_press=Sdps,
                  Vmax_out=Vmax_out, B_used=B, t_peak=t_peak,
                  profile_X=X[::-1], profile_eta=eta[::-1],
                  profile_u=d[::-1], profile_w=Dtot[::-1],
                  notes="; ".join(notes))


# --- self-tests -----------------------------------------------------------------
def _approx(a, b, tol=1e-3):
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    # 1) Myers (B=1) pressure == Holland(B=1); Holland B>1 is more peaked
    dP = 2.35 * _INHG
    assert _approx(_pressure_at(35000.0, 35000.0, 96000.0, dP, 1.0),
                   96000.0 + dP * math.exp(-1.0), 1e-9)
    # 2) cyclostrophic Vmax relation: at f=0, V_gr(R) = sqrt(B dP/(rho_a e))
    for B in (1.0, 1.5, 2.0):
        vg = _gradient_wind(35000.0, 35000.0, dP, B, 0.0, 1.20)
        assert _approx(vg, math.sqrt(B * dP / (1.20 * _E)), 1e-6), (B, vg)
    # 3) pressure setup ~ 1.1-1.14 ft per inHg of deficit (inverse barometer)
    s_per_inhg = (_INHG) / (RHO_W * G_SI) / _FT
    assert 1.0 < s_per_inhg < 1.2, s_per_inhg
    # 4) Chesapeake example, Myers/Bodine (B=1): max wind ~ Bodine's Vx=102 mph, and
    #    peak surge in the ballpark of Bodine's 13.4 ft (parametric vs graphical isovel).
    inp = {fd.key: fd.default for fd in INPUTS}
    rm = compute({**inp, "wind_model": "Myers / Bodine (1954)"})
    assert _approx(rm.Vmax_out / _MPH, 102.0, 0.05), rm.Vmax_out / _MPH      # within 5%
    assert 10.0 < rm.peak_surge / _FT < 20.0, rm.peak_surge / _FT           # Bodine ~13.4
    assert rm.B_used == 1.0 and rm.S_wind > 0 and rm.S_press > 0
    # 5) The recovered TM-35 source curve produces the documented Chesapeake
    #    reference magnitude; Holland remains available and explicit Vmax overrides B.
    rs = compute(inp)
    assert _approx(rs.peak_surge / _FT, 13.41, 0.06), rs.peak_surge / _FT
    assert _approx(rs.t_peak / 3600.0, 17.0, 0.06), rs.t_peak / 3600.0
    assert rs.B_used == 1.0
    rh = compute({**inp, "wind_model": "Holland (1980)"})
    assert rh.B_used == 1.5
    rv = compute({**inp, "wind_model": "Holland (1980)", "Vmax": 50.0})
    assert _approx(rv.B_used, 1.20 * _E * 50.0 ** 2 / (rv_dP := (float(inp["Pn"]) - float(inp["Pc"]))), 1e-6)
    print(f"  self-tests: PASS (TM-35 source curves, wind models, pressure setup, Vmax->B; "
          f"source peak {rs.peak_surge/_FT:.1f} ft; Myers Vmax {rm.Vmax_out/_MPH:.0f} mph)")


def _print_default_example() -> None:
    inp = {fd.key: fd.default for fd in INPUTS}
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  Chesapeake Bay Entrance traverse (Bodine TM-35 example):")
    for wm in ("Bodine TM-35 source curves", "Myers / Bodine (1954)", "Holland (1980)"):
        r = compute({**inp, "wind_model": wm})
        print(f"    [{wm:22s}] peak {r.peak_surge/_FT:6.2f} ft  "
              f"(Sx {r.S_wind/_FT:.1f}, Sy {r.S_bathy/_FT:.1f}, Sp {r.S_press/_FT:.1f}); "
              f"Vmax {r.Vmax_out/_MPH:5.1f} mph; B {r.B_used:.2f}")
    print(f"  Bodine TM-35 reported reference: peak surge 13.41 ft at 17 h; the source-curve")
    print(f"  path preserves its appendix forcing and finite-difference scheme. Parametric paths are sensitivity tools.")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
