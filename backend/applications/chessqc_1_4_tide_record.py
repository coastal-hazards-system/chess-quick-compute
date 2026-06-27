"""CHESS-QC application 1-4 — Constituent Tide Record Generation.

Originating ACES application: 1-4 "Constituent Tide Record Generation" (functional area:
Wave Prediction). Predicts a water-level time series at a gage from the harmonic tidal
constituents (amplitude and epoch per constituent) by the classical harmonic method.

Classification: exact (closed-form harmonic synthesis with exact Schureman astronomical
node factors and equilibrium arguments -- no empirical fitting; reproduces the User's Guide
Table 1-4-1 to the stated tolerance).
Theory and references: Schureman (1971) [reprint of C&GS Special Pub. 98, 1940]; harmonic
method after Lord Kelvin. Equations transcribed in docs/EQUATIONS.md, TR chapter 1-4
(eqs 1-2); constituent speeds from ACES Appendix Table A-5.

    h(t) = H0 + sum_n  f_n * A_n * cos[ a_n*t + (V0+u)_n - kappa_n ]

A_n (amplitude) and kappa_n (epoch / phase lag) are user inputs; a_n (speed) is from
Table A-5; f_n (node factor) and (V0+u)_n (equilibrium argument) are computed
astronomically from the start date/time and the gage longitude.

Convention (pinned by reproducing ACES User's Guide Example 1-4 to 0.04 ft over a 120 h
record, zero free parameters): the equilibrium argument uses T = 15*H (mean solar time
measured from local midnight, H = start hour of day), and the slow astronomical
longitudes s,h,p,N,p1 are evaluated at UT = H - longitude_west/15. (This is the ACES
longitude handling; it is what reproduces the published Buzzards Bay record.)

Self-containment: zero sibling imports; embeds its own contract dataclasses and the
Schureman astronomy. Runnable standalone:
    python chessqc_1_4_tide_record.py
which reproduces the User's Guide Example 1-4 oracle, then prints a summary.
stdlib + numpy only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

D2R = math.pi / 180.0


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
    note: str = ""           # hover definition shown on the output label


# --- application metadata --------------------------------------------------------
APP_META = AppMeta(
    aces_id="1-4",
    name="Constituent Tide Record Generation",
    area="Wave Prediction",
    classification="exact",
    cite="Schureman (1971); Table A-5; TR 1-4",
    default_system="US",
)

# Canonical constituent table (ACES Appendix Table A-5 order). Each entry:
#   name: (speed deg/hr, [cT,cs,ch,cp,cp1], const deg, f-key, [(u-term, sign), ...])
# V0 = cT*T + cs*s + ch*h + cp*p + cp1*p1 + const ;  u = sum(sign * u_term)
_C = {
 "M2":(28.9841042,[2,-2,2,0,0],0,   "M2",  [("M2",1)]),
 "S2":(30.0000000,[2,0,0,0,0],0,    "1",   []),
 "N2":(28.4397295,[2,-3,2,1,0],0,   "M2",  [("M2",1)]),
 "K1":(15.0410686,[1,0,1,0,0],90,   "K1",  [("K1",1)]),
 "M4":(57.9682084,[4,-4,4,0,0],0,   "M2^2",[("M2",2)]),
 "O1":(13.9430356,[1,-2,1,0,0],-90, "O1",  [("O1",1)]),
 "M6":(86.9523127,[6,-6,6,0,0],0,   "M2^3",[("M2",3)]),
 "MK3":(44.0251729,[3,-2,3,0,0],90, "M2K1",[("M2",1),("K1",1)]),
 "S4":(60.0000000,[4,0,0,0,0],0,    "1",   []),
 "MN4":(57.4238337,[4,-5,4,1,0],0,  "M2^2",[("M2",2)]),
 "NU2":(28.5125831,[2,-3,4,-1,0],0, "M2",  [("M2",1)]),
 "S6":(90.0000000,[6,0,0,0,0],0,    "1",   []),
 "MU2":(27.9682084,[2,-4,4,0,0],0,  "M2",  [("M2",1)]),
 "2N2":(27.8953548,[2,-4,2,2,0],0,  "M2",  [("M2",1)]),
 "OO1":(16.1391017,[1,2,1,0,0],90,  "OO1", [("OO1",1)]),
 "LAMBDA2":(29.4556253,[2,-1,2,-1,0],180,"M2",[("M2",1)]),
 "S1":(15.0000000,[1,0,0,0,0],0,    "1",   []),
 "M1":(14.4966939,[1,-1,1,0,0],90,  "O1",  [("O1",1)]),
 "J1":(15.5854433,[1,1,1,-1,0],90,  "J1",  [("J1",1)]),
 "MM":(0.5443747,[0,1,0,-1,0],0,    "MM",  []),
 "SSA":(0.0821373,[0,0,2,0,0],0,    "1",   []),
 "SA":(0.0410686,[0,0,1,0,0],0,     "1",   []),
 "MSF":(1.0158958,[0,2,-2,0,0],0,   "M2",  [("M2",1)]),
 "MF":(1.0980331,[0,2,0,0,0],0,     "MF",  [("MF",1)]),
 "RHO1":(13.4715145,[1,-3,3,-1,0],-90,"O1",[("O1",1)]),
 "Q1":(13.3986609,[1,-3,1,1,0],-90, "O1",  [("O1",1)]),
 "T2":(29.9589333,[2,0,-1,0,1],0,   "1",   []),
 "R2":(30.0410667,[2,0,1,0,-1],180, "1",   []),
 "2Q1":(12.8542862,[1,-4,1,2,0],-90,"O1",  [("O1",1)]),
 "P1":(14.9589314,[1,0,-1,0,0],-90, "1",   []),
 "2SM2":(31.0158958,[2,2,-2,0,0],0, "M2",  [("M2",1)]),
 "M3":(43.4761563,[3,-3,3,0,0],0,   "M3",  [("M3",1)]),
 "L2":(29.5284789,[2,-1,2,-1,0],180,"M2",  [("M2",1)]),
 "2MK3":(42.9271398,[3,-4,3,0,0],-90,"M2^2K1",[("M2",2),("K1",-1)]),
 "K2":(30.0821373,[2,0,2,0,0],0,    "K2",  [("K2",1)]),
 "M8":(115.9364169,[8,-8,8,0,0],0,  "M2^4",[("M2",4)]),
 "MS4":(58.9841042,[4,-2,2,0,0],0,  "M2",  [("M2",1)]),
}
CANON = list(_C.keys())   # 37 constituents in Table A-5 order; row index = CANON.index + 1

_FT = 0.3048
# ACES User's Guide Example 1-4 (Buzzards Bay Entrance, MA): 25 constituents.
# row = [Table A-5 index (1-based), amplitude (SI m), epoch (deg)]
_EX = {
 "M2":(1.621,269.90),"S2":(0.303,283.60),"N2":(0.447,245.10),"K1":(0.262,114.00),
 "M4":(0.266,136.70),"O1":(0.221,123.90),"M6":(0.070,241.90),"MK3":(0.045,138.00),
 "MN4":(0.113,82.20),"NU2":(0.077,262.20),"MU2":(0.070,225.00),"2N2":(0.071,225.70),
 "LAMBDA2":(0.011,276.30),"S1":(0.038,55.30),"M1":(0.016,119.00),"J1":(0.017,109.00),
 "SSA":(0.037,44.60),"SA":(0.112,151.60),"Q1":(0.045,112.60),"T2":(0.018,283.60),
 "P1":(0.091,123.80),"L2":(0.045,294.70),"2MK3":(0.039,159.00),"K2":(0.091,274.20),
 "MS4":(0.076,231.00),
}
_DEFAULT_ROWS = [[float(CANON.index(nm) + 1), amp * _FT, ep] for nm, (amp, ep) in _EX.items()]

INPUTS = (
    Field("year", "Start year", "int", "", "", default=1989, lo=1800, hi=2200),
    Field("month", "Start month", "int", "", "", default=1, lo=1, hi=12),
    Field("day", "Start day", "int", "", "", default=10, lo=1, hi=31),
    Field("hour", "Start hour (local)", "float", "hr", "hr", default=10.0, lo=0.0, hi=24.0,
          note="hour of day at the gage (0 to 24)"),
    Field("length_hr", "Record length", "float", "hr", "hr", default=120.0, lo=0.1, hi=1e5),
    Field("interval_min", "Output interval", "float", "min", "min", default=15.0, lo=0.1, hi=1440.0),
    Field("H0", "Mean water level above datum", "float", "m", "ft", default=1.79 * _FT,
          lo=-1e3, hi=1e3, note="datum offset (e.g. MLLW)"),
    Field("gage_lon", "Gage longitude", "float", "deg", "deg", default=70.62, lo=-180.0, hi=360.0,
          note="degrees West (positive)"),
    Field("constituents", "Constituents", "table", default=_DEFAULT_ROWS,
          columns=(("Table A-5 #", "", ""), ("Amplitude", "m", "ft"), ("Epoch", "deg", "deg")),
          note="one row per constituent: Table A-5 index (1=M2,2=S2,3=N2,4=K1,5=M4,6=O1,...), "
               "amplitude, epoch in degrees. See CANON / module docstring for the full index list."),
)

OUTPUTS = (
    Out("h_max", "Maximum elevation", "m", "ft", "scalar",
        note="highest predicted water level above the prediction datum over the record (high water)."),
    Out("h_min", "Minimum elevation", "m", "ft", "scalar",
        note="lowest predicted water level above the prediction datum over the record (low water)."),
    Out("range", "Tidal range (max - min)", "m", "ft", "scalar",
        note="peak-to-trough range of the predicted water level, h_max minus h_min."),
    Out("h_start", "Elevation at start", "m", "ft", "point",
        note="predicted water level above datum at the first sample (start date/time)."),
    Out("n_constituents", "Number of active constituents", "", "", "scalar",
        note="count of harmonic constituents with nonzero amplitude included in the synthesis."),
    Out("profile_t", "Profile: time", "hr", "hr", "profile",
        note="elapsed time from the start epoch at each sample, in hours."),
    Out("profile_h", "Profile: tide elevation", "m", "ft", "profile",
        note="predicted water level above datum versus time, the harmonic sum h(t) = H0 + sum f_n A_n cos(...)."),
)


@dataclass
class Result:
    h_max: float; h_min: float; range: float; h_start: float; n_constituents: float
    profile_t: np.ndarray
    profile_h: np.ndarray
    notes: str = ""


# --- Schureman astronomy --------------------------------------------------------
def _jd(year: int, month: int, day: int, hour: float) -> float:
    """Julian date (UT) for a Gregorian calendar date and fractional hour."""
    if month <= 2:
        year -= 1; month += 12
    A = year // 100
    B = 2 - A + A // 4
    return (math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1))
            + day + B - 1524.5) + hour / 24.0


def _astro(jd_ut: float):
    """Mean longitudes s,h,p,N,p1 (deg) via Schureman 1900-epoch polynomials."""
    Tc = (jd_ut - 2415020.0) / 36525.0
    h = 279.69668 + 36000.76892 * Tc + 0.00030 * Tc * Tc
    s = 270.43659 + 481267.89057 * Tc - 0.00198 * Tc * Tc
    p = 334.32956 + 4069.03403 * Tc - 0.01032 * Tc * Tc
    N = 259.18328 - 1934.14201 * Tc + 0.00208 * Tc * Tc
    p1 = 281.22083 + 1.71902 * Tc + 0.00045 * Tc * Tc
    return [x % 360.0 for x in (s, h, p, N, p1)]


def _node(N: float):
    """Schureman node-factor auxiliaries (I, xi, nu, nu', nu'') in degrees, from node N."""
    om = 23.452 * D2R; ii = 5.145 * D2R; Nr = N * D2R
    I = math.acos(math.cos(om) * math.cos(ii) - math.sin(om) * math.sin(ii) * math.cos(Nr))
    at1 = math.atan2(math.cos(0.5 * (om - ii)) * math.sin(0.5 * Nr),
                     math.cos(0.5 * (om + ii)) * math.cos(0.5 * Nr))
    at2 = math.atan2(math.sin(0.5 * (om - ii)) * math.sin(0.5 * Nr),
                     math.sin(0.5 * (om + ii)) * math.cos(0.5 * Nr))
    xi = Nr - at1 - at2; nu = at1 - at2
    nup = math.atan2(math.sin(2 * I) * math.sin(nu), math.sin(2 * I) * math.cos(nu) + 0.3347)
    nupp = 0.5 * math.atan2(math.sin(I) ** 2 * math.sin(2 * nu),
                            math.sin(I) ** 2 * math.cos(2 * nu) + 0.0727)
    deg = 1.0 / D2R
    return I, xi * deg, nu * deg, nup * deg, nupp * deg


def _factors(N: float):
    """Return (f_dict, u_dict) of node factors and phase corrections (deg) for node N."""
    I, xi, nu, nup, nupp = _node(N)
    Ir = I; sinI = math.sin(Ir); cosh2 = math.cos(Ir / 2) ** 2
    f = {
        "1": 1.0,
        "M2": math.cos(Ir / 2) ** 4 / 0.9154,
        "O1": sinI * cosh2 / 0.3800,
        "K1": math.sqrt(0.8965 * math.sin(2 * Ir) ** 2
                        + 0.6001 * math.sin(2 * Ir) * math.cos(nu * D2R) + 0.1006),
        "K2": math.sqrt(19.0444 * sinI ** 4 + 2.7702 * sinI ** 2 * math.cos(2 * nu * D2R) + 0.0981),
        "J1": math.sin(2 * Ir) / 0.7214,
        "OO1": sinI * math.sin(Ir / 2) ** 2 / 0.0164,
        "MF": sinI ** 2 / 0.1578,
        "MM": (2.0 / 3.0 - sinI ** 2) / 0.5021,
        "M3": math.cos(Ir / 2) ** 6 / 0.8758,
    }
    u = {"M2": 2 * xi - 2 * nu, "O1": 2 * xi - nu, "K1": -nup, "K2": -2 * nupp,
         "J1": -nu, "OO1": -2 * xi - nu, "MF": -2 * xi, "MM": 0.0, "M3": 3 * xi - 3 * nu,
         "1": 0.0}
    return f, u


def _fval(fkey: str, f: dict) -> float:
    if fkey in f:
        return f[fkey]
    if fkey == "M2^2": return f["M2"] ** 2
    if fkey == "M2^3": return f["M2"] ** 3
    if fkey == "M2^4": return f["M2"] ** 4
    if fkey == "M2K1": return f["M2"] * f["K1"]
    if fkey == "M2^2K1": return f["M2"] ** 2 * f["K1"]
    return 1.0


def _validate(inp: dict) -> None:
    for fld in INPUTS:
        if fld.kind not in ("float", "int"):
            continue
        v = float(inp[fld.key])
        if not (fld.lo <= v <= fld.hi):
            raise ValueError(f"{fld.label} ({fld.key}) = {v} outside [{fld.lo}, {fld.hi}]")


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Synthesizes a water-level time series at a tide gage from harmonic tidal '
            'constituents (per-constituent amplitude and epoch), summing each constituent '
            'with astronomically computed node factors and equilibrium arguments. Returns '
            'the elevation profile plus its max, min, range, and starting value.',
 'methods': [{'name': 'Classical harmonic tide synthesis (Schureman)',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'h = H_0 + \\sum_{n=1}^{N} f_n \\, A_n \\cos\\left[ '
                                    'a_n t + (V_0 + u)_n - \\kappa_n \\right]',
                             'desc': 'Full harmonic synthesis (TR 1-4-2 eq 2): tide '
                                     'elevation at time t as the sum of N constituents '
                                     'above datum H_0.'},
                            {'tex': 'h_n = A_n \\cos\\left( a_n t + \\alpha_n \\right)',
                             'desc': 'Contribution of a single harmonic constituent n (TR '
                                     '1-4-1 eq 1); alpha_n is its phase at the initial '
                                     'epoch.'},
                            {'tex': '(V_0)_n = c_T T + c_s s + c_h h + c_p p + c_{p_1} p_1 '
                                    '+ C_n',
                             'desc': 'Equilibrium argument from the astronomical '
                                     "longitudes, with the constituent's integer "
                                     'multipliers from Table A-5.'},
                            {'tex': 'T = 15 \\, H',
                             'desc': 'Hour-angle term in mean solar time measured from '
                                     'local midnight (H = start hour of day), the ACES '
                                     'longitude convention.'}]}],
 'symbols': [['h', 'tide elevation at time t (above prediction datum)'],
             ['H_0', 'mean water level above the prediction datum'],
             ['A_n', 'amplitude of constituent n (user input)'],
             ['kappa_n', 'local phase lag (epoch) of constituent n, degrees (user input)'],
             ['a_n', 'angular speed of constituent n, deg/hr (Table A-5)'],
             ['f_n', 'node factor of constituent n (Schureman astronomy)'],
             ['(V_0+u)_n',
              'local equilibrium argument of constituent n at the initial epoch'],
             ['t', 'time measured from the initial epoch'],
             ['s, h, p, p_1',
              'mean longitudes of moon, sun, lunar perigee, and solar perigee'],
             ['N', 'number of active constituents']],
 'references': ['Schureman (1971), C&GS Special Pub. 98 (reprint of 1940)',
                'Harris (1981), WES SR-7',
                'Headquarters DA (1989), EM 1110-2-1414 Ch. 2',
                'ACES Appendix Table A-5 (constituent speeds)']}


def compute(inp: dict) -> Result:
    """Constituent tide record for SI inputs. Returns the elevation time series."""
    _validate(inp)
    year = int(inp["year"]); month = int(inp["month"]); day = int(inp["day"])
    hour = float(inp["hour"]); length = float(inp["length_hr"])
    interval = float(inp["interval_min"]) / 60.0
    H0 = float(inp["H0"]); lon = float(inp["gage_lon"])

    # equilibrium argument: T from local midnight; slow longitudes at UT = hour - lon_west/15
    ut0 = hour - lon / 15.0
    s, h, p, N, p1 = _astro(_jd(year, month, day, ut0))
    T = 15.0 * hour
    f, u = _factors(N)

    terms = []
    for row in inp["constituents"]:
        if not row or len(row) < 3:
            continue
        idx = int(round(float(row[0]))); amp = float(row[1]); kappa = float(row[2])
        if amp == 0.0 or not (1 <= idx <= len(CANON)):
            continue
        name = CANON[idx - 1]
        speed, c, const, fkey, uc = _C[name]
        V0 = c[0] * T + c[1] * s + c[2] * h + c[3] * p + c[4] * p1 + const
        un = sum(sign * u[t] for t, sign in uc)
        arg0 = (V0 + un - kappa) % 360.0
        terms.append((speed, _fval(fkey, f) * amp, arg0))

    n = max(2, int(round(length / interval)) + 1)
    t = np.linspace(0.0, length, n)
    hsum = np.full(n, H0)
    for speed, amp, arg0 in terms:
        hsum = hsum + amp * np.cos((speed * t + arg0) * D2R)

    notes = [f"{len(terms)} active constituents; start {year:04d}-{month:02d}-{day:02d} "
             f"{hour:05.2f}h, {length:.0f} h record"]
    return Result(h_max=float(hsum.max()), h_min=float(hsum.min()),
                  range=float(hsum.max() - hsum.min()), h_start=float(hsum[0]),
                  n_constituents=float(len(terms)), profile_t=t, profile_h=hsum,
                  notes="; ".join(notes))


# --- self-tests (User's Guide Example 1-4 oracle) -------------------------------
def _self_tests() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    # ACES User's Guide Table 1-4-1 (Buzzards Bay), elevations in ft at given hours
    oracle = {0.0: 4.26, 0.25: 4.35, 0.5: 4.39, 0.75: 4.38, 1.0: 4.32,
              118.5: 0.65, 119.0: 0.38, 119.75: 0.01, 120.0: -0.08}
    t = r.profile_t; hft = r.profile_h / _FT
    worst = 0.0
    for tt, exp in oracle.items():
        i = int(round(tt / 0.25))
        got = hft[i]
        worst = max(worst, abs(got - exp))
        assert abs(got - exp) <= 0.06, f"t={tt}: got {got:.2f} ft, manual {exp:.2f} ft"
    assert r.n_constituents == 25
    print(f"  self-tests: PASS (User's Guide Example 1-4, max dev {worst:.3f} ft over the record)")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  Buzzards Bay Entrance, MA: {int(r.n_constituents)} constituents, "
          f"{len(r.profile_t)} samples")
    print("  elevation (ft): start=%.2f  max=%.2f  min=%.2f  range=%.2f" % (
        r.h_start / _FT, r.h_max / _FT, r.h_min / _FT, r.range / _FT))
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
