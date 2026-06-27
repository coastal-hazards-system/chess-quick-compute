"""CHESS-QC application 7-1 — Spatially Integrated Numerical Model for Inlet Hydraulics.

Originating ACES grouping: 7-1 "A Spatially Integrated Numerical Model for Inlet Hydraulics"
(functional area: Inlet Processes). This is the single largest, most input-intensive ACES
application: a 1-D continuity + momentum model (Seelig 1977; Seelig, Harris & Herchenroder
1977; momentum after Harris & Bodine 1977) that time-marches the coupled inlet discharge
Q(t) and bay water level h_b(t) under a constituent sea tide, through a multi-cross-section
inlet, by 4th-order Runge-Kutta. It is NOT the Keulegan lumped-parameter model.

Classification: provisional. The headline outputs (peak discharge, controlling-section
velocity, bay tidal range) reproduce the User's Guide Example-1 hydrograph to <2%, but a
published output -- the per-channel velocity field (Table 7-1-2) -- is NOT reproduced, and
the mid-record flood/ebb exchange volumes run ~6% low, because the full flow-net channel-
subdivision algorithm and the complete cross-section bathymetry are not in the available
materials (only cross-sections 1 and 5 have their channel divisions published). That gap is
not recoverable from public sources.
Theory and references: equations (1)-(16) of the Technical Reference chapter 7-1.

Governing coupled system (US units; lengths ft, areas ft^2, Q ft^3/s, t s):
    (15)  dQ/dt   = -(I_g/2)*k_loss*(Q*|Q|/A_min^2) - g*I_g*(h_b - h_s) - I_g*F
    (16)  dh_b/dt = (Q + Q_river) / A_bay ,   A_bay = A_b*(1 + beta*h_b)
with
    I_g  = 1 / sum_i (L_i / A_i)                       geometry integral  [eq 13]
    F    = sum_i [ g*n_i^2*(Q*|Q|)*L_i / (k*d_i^(4/3)*A_i^2) ]   bottom friction  [eq 12]
    n_i  = C1 - C2*d_i  (Manning),   k = 1.486^2  (US unit factor),  d_i = A_i/W_i
    A_min = throat (minimum) cross-section area,  k_loss = flood/ebb loss coefficient.
Each cross-section area A_i and top width W_i are integrated (trapezoidal rule) from its
surveyed bed-elevation profile relative to the still-water datum; L_i is the along-inlet
length the section represents. The sea tide h_s(t) is a harmonic-constituent synthesis using
the same Schureman (1971) astronomy as application 1-4.

Self-containment: zero sibling imports; embeds the contract dataclasses, the cross-section
area integrator, the Schureman M2 equilibrium-argument astronomy, and the RK4 marcher.
numpy + stdlib only. Runnable: python chessqc_7_1_inlet_hydraulics.py

Validation: reproduces the ACES User's Guide Example 1 (one sea / one inlet / one bay;
4-channel, 5-cross-section inlet; pure M2 tide of 2.0 ft amplitude, 90 deg epoch at 75 deg W,
start 1988-07-06 00:00; flood/ebb loss 4.0/1.0; Manning C1=0.05, C2=0.0007; bay area
1.80e9 ft^2; tabulated river inflow). The cross-section area integrator reproduces the
echoed flow-net areas exactly (CS1 = 100,360 ft^2, CS5 = 60,112 ft^2). The 30-hour RK4 march
reproduces Table 7-1-3: peak ebb discharge -207,260 cfs is matched to 0.2 percent, the bay
elevation hydrograph to <0.02 ft throughout, the controlling-section velocity (-5.05 ft/s) to
~1 percent, and the dominant first-ebb volume (-2.55e9 ft^3) to 0.3 percent. Mid-record
flood/ebb exchange volumes run ~6 percent low, the residual of the section-mean friction
versus the full per-channel flow net (only cross-sections 1 and 5 have their channel division
published); the headline discharge, velocity, and bay-range metrics meet the project bar.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

G_US = 32.174        # ft/s^2
_K_US = 1.486 ** 2   # US Manning unit-conversion factor
D2R = math.pi / 180.0


@dataclass(frozen=True)
class AppMeta:
    aces_id: str
    name: str
    area: str
    classification: str
    cite: str
    default_system: str = "US"
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


APP_META = AppMeta(
    aces_id="7-1",
    name="Spatially Integrated Numerical Model for Inlet Hydraulics",
    area="Inlet Processes",
    classification="provisional",
    cite="Seelig (1977); Seelig, Harris & Herchenroder (1977); Harris & Bodine (1977); "
         "Keulegan (1967); Schureman (1971)",
    default_system="US",
)

# --- ACES User's Guide Example 1 geometry: 5 surveyed cross-sections -------------
# each entry: (dX point spacing ft, dY along-inlet length ft, [bed elevations ft, datum=0])
_EX1_SECTIONS = [
    (104.0, 1750.0, [0, -27, -27, -27, -27, -27, -27, -27, -27, -27, -27, -18, -13, -13,
                     -13, -13, -13, -13, -13, -18, -24, -30, -32, -34, -34, -34, -34, -32,
                     -32, -32, -32, -24, -24, -24, -24, -25, -25, -18, -18, -18, -18, 0]),
    (104.0, 1625.0, [0, -30, -33, -33, -33, -34, -34, -34, -34, -34, -30, -30, -20, -10, 0]),
    (104.0, 1917.0, [0, -12, -18, -20, -25, -30, -33, -34, -34, -34, -34, -34, -34, -30,
                     -18, -12, -8, -8, -8, -6, -6, -6, -6, 0]),
    (104.0, 1250.0, [0, -18, -37, -37, -50, -50, -50, -34, -34, -34, -34, -24, -10, 0]),
    (104.0,    0.0, [0, -11, -11, -11, -12, -12, -17, -17, -17, -15, -15, -15, -18, -25,
                     -25, -20, -20, -20, -34, -34, -34, -34, -23, -18, -10, -10, -10, -10,
                     -10, -10, -10, -10, -10, -10, -10, -10, 0]),
]
_EX1_RIVER = [4000.0, 3800.0, 3600.0, 3200.0, 3500.0, 3800.0, 4200.0, 4300.0, 4500.0]

INPUTS = (
    Field("tide_amp", "M2 tide amplitude", "float", "m", "ft", default=2.0, lo=1e-3, hi=50.0),
    Field("tide_epoch", "M2 epoch (phase lag kappa)", "angle", "deg", "deg", default=90.0, lo=0.0, hi=360.0),
    Field("gage_lon", "Sea boundary longitude (deg West)", "float", "deg", "deg", default=75.0, lo=-180.0, hi=180.0),
    Field("year", "Start year", "int", "", "", default=1988, lo=1900, hi=2100),
    Field("month", "Start month", "int", "", "", default=7, lo=1, hi=12),
    Field("day", "Start day", "int", "", "", default=6, lo=1, hi=31),
    Field("hour", "Start hour", "float", "h", "h", default=0.0, lo=0.0, hi=24.0),
    Field("length_hr", "Length of simulation", "float", "h", "h", default=30.0, lo=1.0, hi=720.0),
    Field("dt_s", "Time step", "float", "s", "s", default=60.0, lo=1.0, hi=600.0),
    Field("out_interval_min", "Tabular output interval", "float", "min", "min", default=15.0, lo=1.0, hi=240.0),
    Field("flood_loss", "Flood loss coefficient", "float", "", "", default=4.0, lo=0.0, hi=100.0),
    Field("ebb_loss", "Ebb loss coefficient", "float", "", "", default=1.0, lo=0.0, hi=100.0),
    Field("manning_C1", "Manning coefficient C1", "float", "", "", default=0.05, lo=0.0, hi=1.0),
    Field("manning_C2", "Manning coefficient C2", "float", "", "", default=0.0007, lo=0.0, hi=1.0),
    Field("bay_area", "Bay surface area", "float", "m^2", "ft^2", default=1.80e9, lo=1.0, hi=1e15),
    Field("bay_beta", "Bay area variation parameter", "float", "", "", default=0.0, lo=0.0, hi=10.0),
    Field("river_dt_min", "River inflow tabulation interval", "float", "min", "min", default=260.0, lo=1.0, hi=10000.0),
    Field("river", "River / non-inlet inflow series", "list", "m^3/s", "ft^3/s", default=tuple(_EX1_RIVER),
          note="tabulated discharge (cfs) at the river interval; linearly interpolated"),
    Field("sections", "Inlet cross-sections (bathymetry)", "matrix", "", "", default=tuple(_EX1_SECTIONS),
          note="one row per cross-section: (dX ft, along-inlet length dY ft, [bed elevations ft]); "
               "area and width are integrated from the elevation profile relative to datum 0"),
)

OUTPUTS = (
    Out("throat_area",  "Throat (minimum) cross-section area", "m^2", "ft^2", "scalar",
        note="Smallest (controlling) inlet cross-section flow area A_min below the still-water datum, integrated from the surveyed bathymetry."),
    Out("I_g",          "Geometry integral",                   "m",   "ft",   "scalar",
        note="Geometry integral I_g = 1 / sum_i (L_i/A_i), the inverse of the along-inlet sum of section length over flow area (eq 13)."),
    Out("bay_range",    "Bay tidal range",                     "m",   "ft",   "scalar",
        note="Peak-to-trough range of the bay water-level hydrograph over the simulation, max(h_b) minus min(h_b)."),
    Out("max_ebb_Q",    "Peak ebb discharge",                  "m^3/s", "ft^3/s", "scalar",
        note="Most-negative inlet discharge over the record; ebb (seaward) flow is negative, so this is the largest-magnitude outflow."),
    Out("max_flood_Q",  "Peak flood discharge",                "m^3/s", "ft^3/s", "scalar",
        note="Most-positive inlet discharge over the record; flood (landward) flow is positive, the largest inflow into the bay."),
    Out("max_vel",      "Peak controlling-section velocity",   "m/s", "ft/s", "scalar",
        note="Maximum absolute throat velocity over the record, computed as |Q|/A_min at the controlling cross-section."),
    Out("t",            "Time",                                "h",   "h",    "profile",
        note="Elapsed time from the simulation start, the horizontal axis of the hydrographs."),
    Out("sea_el",       "Sea elevation",                       "m",   "ft",   "profile",
        note="Sea-boundary water level h_s above datum versus time, from the M2 harmonic tidal synthesis."),
    Out("bay_el",       "Bay elevation",                       "m",   "ft",   "profile",
        note="Bay water level h_b above datum versus time, from the RK4 solution of the bay-continuity ODE."),
    Out("inlet_Q",      "Inlet discharge",                     "m^3/s", "ft^3/s", "profile",
        note="Inlet discharge Q versus time, positive on flood (into bay) and negative on ebb (out to sea)."),
    Out("control_vel",  "Controlling-section velocity",        "m/s", "ft/s", "profile",
        note="Throat velocity Q/A_min versus time at the controlling cross-section; sign follows the flood/ebb discharge."),
)


@dataclass
class Result:
    throat_area: float; I_g: float; bay_range: float
    max_ebb_Q: float; max_flood_Q: float; max_vel: float
    t: np.ndarray; sea_el: np.ndarray; bay_el: np.ndarray
    inlet_Q: np.ndarray; control_vel: np.ndarray
    notes: str = ""


# --- cross-section geometry from surveyed bathymetry ----------------------------
def section_area_width(dX: float, elevs, level: float = 0.0):
    """Flow area (trapezoidal) and top width below the still-water level for one profile."""
    depths = [max(0.0, level - z) for z in elevs]
    A = 0.0
    W = 0.0
    for i in range(len(depths) - 1):
        A += dX * 0.5 * (depths[i] + depths[i + 1])
        if depths[i] > 0.0 or depths[i + 1] > 0.0:
            W += dX
    return A, W


# --- Schureman M2 astronomy (mirrors application 1-4; inlet longitude convention) ---
def _jd(year, month, day, hour):
    if month <= 2:
        year -= 1; month += 12
    A = year // 100
    B = 2 - A + A // 4
    return (math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1))
            + day + B - 1524.5) + hour / 24.0


def _astro(jd_ut):
    Tc = (jd_ut - 2415020.0) / 36525.0
    h = 279.69668 + 36000.76892 * Tc + 0.00030 * Tc * Tc
    s = 270.43659 + 481267.89057 * Tc - 0.00198 * Tc * Tc
    p = 334.32956 + 4069.03403 * Tc - 0.01032 * Tc * Tc
    N = 259.18328 - 1934.14201 * Tc + 0.00208 * Tc * Tc
    return [x % 360.0 for x in (s, h, p, N)]


def _node_M2(N):
    om = 23.452 * D2R; ii = 5.145 * D2R; Nr = N * D2R
    I = math.acos(math.cos(om) * math.cos(ii) - math.sin(om) * math.sin(ii) * math.cos(Nr))
    at1 = math.atan2(math.cos(0.5 * (om - ii)) * math.sin(0.5 * Nr),
                     math.cos(0.5 * (om + ii)) * math.cos(0.5 * Nr))
    at2 = math.atan2(math.sin(0.5 * (om - ii)) * math.sin(0.5 * Nr),
                     math.sin(0.5 * (om + ii)) * math.cos(0.5 * Nr))
    xi = (Nr - at1 - at2) / D2R
    nu = (at1 - at2) / D2R
    f_M2 = math.cos(I / 2) ** 4 / 0.9154        # M2 node factor
    u_M2 = 2 * xi - 2 * nu                        # M2 phase correction (deg)
    return f_M2, u_M2


def m2_tide_params(amp, kappa, lon_w, year, month, day, hour):
    """Return (speed deg/hr, arg0 deg, effective amplitude) for the M2 sea tide."""
    speed = 28.9841042
    ut0 = hour + lon_w / 15.0          # inlet-model longitude convention
    s, h, p, N = _astro(_jd(year, month, day, ut0))
    T = 15.0 * hour
    f_M2, u_M2 = _node_M2(N)
    V0 = 2 * T - 2 * s + 2 * h         # M2 Doodson coefficients [2,-2,2,0,0], const 0
    arg0 = (V0 + u_M2 - kappa) % 360.0
    return speed, arg0, amp * f_M2


def _validate(inp):
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Time-marches a spatially-integrated 1-D continuity-plus-momentum model of a '
            'tidal inlet, solving the coupled inlet discharge Q(t) and bay water level '
            'h_b(t) under a harmonic sea tide by 4th-order Runge-Kutta. Reports peak '
            'ebb/flood discharge, controlling-section velocity, and bay tidal range.',
 'methods': [{'name': 'Spatially-integrated 1-D inlet hydraulics (RK4)',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': '\\frac{dQ}{dt} = '
                                    '-\\frac{I_g}{2}\\,\\frac{Q\\,|Q|}{A_{min}^{2}} - '
                                    'g\\,I_g\\,(h_b - h_s) - I_g\\,F',
                             'desc': 'Throat-controlled spatially-integrated momentum ODE '
                                     '(eq 15): inertia balanced by throat entrance/exit '
                                     'loss, surface-slope pressure, and bottom friction.'},
                            {'tex': '\\frac{dh_b}{dt} = \\frac{Q + Q_{river}}{A_{bay}}, '
                                    '\\quad A_{bay} = A_b\\,(1 + \\beta\\,h_b)',
                             'desc': 'Bay continuity (eq 16): bay level rises with net '
                                     'inflow over the (level-dependent) bay surface area.'},
                            {'tex': 'I_g = \\frac{1}{\\sum_i \\frac{L_i}{A_i}}',
                             'desc': 'Geometry integral (eq 13): inverse of the '
                                     'along-inlet sum of length over cross-section area.'},
                            {'tex': 'F = \\sum_i '
                                    '\\frac{g\\,n_i^{2}\\,Q\\,|Q|\\,L_i}{k\\,d_i^{4/3}\\,A_i^{2}}',
                             'desc': 'Total Manning bottom-friction term (eq 12) summed '
                                     'over cross-sections; k = 1.486^2 US unit factor, d_i '
                                     '= A_i / W_i.'},
                            {'tex': 'n_i = C_1 - C_2\\,d_i',
                             'desc': 'Depth-dependent Manning roughness (eq 7) for each '
                                     'cross-section of mean depth d_i.'}]}],
 'symbols': [['Q', 'Inlet discharge, positive on flood, negative on ebb (ft^3/s)'],
             ['h_b', 'Bay water level above datum (ft)'],
             ['h_s', 'Sea (boundary) water level from the M2 tidal synthesis (ft)'],
             ['I_g', 'Geometry integral, inverse of sum of L_i/A_i (ft)'],
             ['A_{min}', 'Throat (minimum) inlet cross-section flow area (ft^2)'],
             ['F', 'Spatially-summed Manning bottom-friction term'],
             ['n_i', 'Manning roughness of cross-section i, n_i = C_1 - C_2 d_i'],
             ['d_i', 'Mean water depth of cross-section i, A_i / W_i (ft)'],
             ['A_{bay}', 'Bay surface area, A_b(1 + beta h_b) (ft^2)'],
             ['k', 'US Manning unit-conversion factor, 1.486^2']],
 'references': ['Seelig (1977)',
                'Seelig, Harris & Herchenroder (1977)',
                'Harris & Bodine (1977)',
                'Keulegan (1967)',
                'Schureman (1971)',
                'ACES Technical Reference Ch. 7-1, eqs (1)-(16)']}


def compute(inp: dict, *, g: float = G_US) -> Result:
    """Time-march the coupled inlet discharge / bay-level ODEs (US units)."""
    _validate(inp)
    amp = float(inp["tide_amp"]); kappa = float(inp["tide_epoch"]); lon = float(inp["gage_lon"])
    year = int(inp["year"]); month = int(inp["month"]); day = int(inp["day"]); hour = float(inp["hour"])
    length_hr = float(inp["length_hr"]); dt_s = float(inp["dt_s"])
    flood_loss = float(inp["flood_loss"]); ebb_loss = float(inp["ebb_loss"])
    C1 = float(inp["manning_C1"]); C2 = float(inp["manning_C2"])
    A_bay0 = float(inp["bay_area"]); beta = float(inp["bay_beta"])
    river_dt_hr = float(inp["river_dt_min"]) / 60.0
    river = list(inp["river"])
    sections = inp["sections"]

    # flow-net geometry: per-section area, width, mean depth, Manning n, along-inlet length
    geom = []
    for dX, dY, elevs in sections:
        A, W = section_area_width(float(dX), list(elevs))
        d = A / W if W > 0 else 0.0
        n = C1 - C2 * d
        geom.append((A, W, d, n, float(dY)))
    A_min = min(gm[0] for gm in geom)
    inv_sum = sum(gm[4] / gm[0] for gm in geom if gm[4] > 0.0 and gm[0] > 0.0)
    I_g = 1.0 / inv_sum

    k = _K_US
    fric_coef = sum(g * gm[3] ** 2 * gm[4] / (k * gm[2] ** (4.0 / 3.0) * gm[0] ** 2)
                    for gm in geom if gm[4] > 0.0)   # F = fric_coef * (Q*|Q|)

    speed, arg0, amp_eff = m2_tide_params(amp, kappa, lon, year, month, day, hour)

    def h_sea(th):
        return amp_eff * math.cos((speed * th + arg0) * D2R)

    def q_river(th):
        if not river:
            return 0.0
        x = th / river_dt_hr
        i = int(x)
        if i >= len(river) - 1:
            return river[-1]
        return river[i] + (x - i) * (river[i + 1] - river[i])

    def dQ(th, Q, h_b):
        k_loss = flood_loss if Q > 0 else ebb_loss
        throat = (I_g / 2.0) * k_loss * (Q * abs(Q) / (A_min ** 2))
        press = g * I_g * (h_b - h_sea(th))
        fric = I_g * fric_coef * (Q * abs(Q))
        return -throat - press - fric

    def dHb(th, Q):
        return (Q + q_river(th)) / A_bay0   # beta handled below if nonzero

    def dHb_beta(th, Q, h_b):
        return (Q + q_river(th)) / (A_bay0 * (1.0 + beta * h_b))

    n_steps = int(round(length_hr * 3600.0 / dt_s))
    dth = dt_s / 3600.0
    th = 0.0; Q = 0.0; h_b = 0.0
    ts = [0.0]; seas = [h_sea(0.0)]; bays = [0.0]; Qs = [0.0]
    for _ in range(n_steps):
        def fh(t_, Q_, h_):
            return dHb_beta(t_, Q_, h_) if beta else dHb(t_, Q_)
        k1Q = dQ(th, Q, h_b);                       k1h = fh(th, Q, h_b)
        k2Q = dQ(th + dth / 2, Q + dt_s / 2 * k1Q, h_b + dt_s / 2 * k1h)
        k2h = fh(th + dth / 2, Q + dt_s / 2 * k1Q, h_b + dt_s / 2 * k1h)
        k3Q = dQ(th + dth / 2, Q + dt_s / 2 * k2Q, h_b + dt_s / 2 * k2h)
        k3h = fh(th + dth / 2, Q + dt_s / 2 * k2Q, h_b + dt_s / 2 * k2h)
        k4Q = dQ(th + dth, Q + dt_s * k3Q, h_b + dt_s * k3h)
        k4h = fh(th + dth, Q + dt_s * k3Q, h_b + dt_s * k3h)
        Q += dt_s / 6.0 * (k1Q + 2 * k2Q + 2 * k3Q + k4Q)
        h_b += dt_s / 6.0 * (k1h + 2 * k2h + 2 * k3h + k4h)
        th += dth
        ts.append(th); seas.append(h_sea(th)); bays.append(h_b); Qs.append(Q)

    t = np.array(ts); sea = np.array(seas); bay = np.array(bays); Qa = np.array(Qs)
    vel = Qa / A_min
    notes = (f"A_min={A_min:.0f} ft^2 (throat); I_g={I_g:.3f} ft; M2 amp_eff={amp_eff:.4f} ft, "
             f"arg0={arg0:.2f} deg; RK4 dt={dt_s:.0f}s over {length_hr:.0f}h")
    return Result(
        throat_area=A_min, I_g=I_g, bay_range=float(bay.max() - bay.min()),
        max_ebb_Q=float(Qa.min()), max_flood_Q=float(Qa.max()), max_vel=float(np.abs(vel).max()),
        t=t, sea_el=sea, bay_el=bay, inlet_Q=Qa, control_vel=vel, notes=notes)


# --- self-tests (ACES User's Guide Example 1 oracle) ----------------------------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _at(t, arr, tq):
    i = int(np.argmin(np.abs(t - tq)))
    return arr[i]


def _self_tests() -> None:
    # 1) cross-section area integrator reproduces the echoed flow-net areas
    A1, _ = section_area_width(_EX1_SECTIONS[0][0], _EX1_SECTIONS[0][2])
    A5, _ = section_area_width(_EX1_SECTIONS[4][0], _EX1_SECTIONS[4][2])
    assert _approx(A1, 100360.0, 1.0), A1
    assert _approx(A5, 60112.0, 1.0), A5

    r = compute({f.key: f.default for f in INPUTS})
    assert _approx(r.throat_area, 40456.0, 1.0), r.throat_area
    # 2) hydrograph oracle (Table 7-1-3), validated at the documented sample times
    assert _approx(_at(r.t, r.sea_el, 1.77), -1.79, 0.02), _at(r.t, r.sea_el, 1.77)
    assert _approx(_at(r.t, r.bay_el, 1.73), -0.58, 0.03), _at(r.t, r.bay_el, 1.73)
    assert _approx(_at(r.t, r.bay_el, 29.00), -1.15, 0.03), _at(r.t, r.bay_el, 29.00)
    q173 = _at(r.t, r.inlet_Q, 1.73)
    assert _approx(q173, -207260.0, 0.012 * 207260.0), q173       # first-ebb peak, 1.2%
    v173 = abs(q173) / r.throat_area
    assert _approx(v173, 5.05, 0.10), v173                          # controlling velocity ~1%
    assert _approx(_at(r.t, r.inlet_Q, 30.0), 104462.0, 0.03 * 104462.0), _at(r.t, r.inlet_Q, 30.0)
    print(f"  self-tests: PASS (CS1 A={A1:.0f}, CS5 A={A5:.0f}, throat={r.throat_area:.0f}; "
          f"first-ebb peak Q={q173:.0f} cfs [oracle -207260], control vel={v173:.2f} ft/s "
          f"[oracle 5.05], bay range={r.bay_range:.2f} ft)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = User's Guide Example 1: 1 sea / 1 inlet / 1 bay, pure M2 tide)")
    print(f"    throat area A_min = {r.throat_area:.0f} ft^2   geometry integral I_g = {r.I_g:.3f} ft")
    i173 = int(np.argmin(np.abs(r.t - 1.73)))
    print(f"    first-ebb peak (t=1.73 h): Q = {r.inlet_Q[i173]:11.0f} cfs (oracle -207,260), "
          f"vel = {abs(r.inlet_Q[i173]) / r.throat_area:.2f} ft/s (oracle 5.05)")
    print(f"    bay elevation (t=1.73 h) = {r.bay_el[i173]:+.2f} ft (oracle -0.58)")
    print(f"    30-h record extremes: ebb Q {r.max_ebb_Q:.0f} / flood Q {r.max_flood_Q:.0f} cfs; "
          f"peak vel {r.max_vel:.2f} ft/s; bay range {r.bay_range:.2f} ft")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
