"""CHESS-QC application 1-6 — Holland Hurricane Wind Model.

Originating ACES application: 1-6 "Holland Hurricane Wind Model" (functional area:
Wave Prediction; a later ACES addition). Given the two Holland (1980) profile
parameters and the storm's central / peripheral pressures, it reconstructs the radial
pressure profile and the gradient- and cyclostrophic-wind profiles of a tropical
cyclone, and reports the maximum wind speed.

Classification: exact (the Holland 1980 model is fully specified -- its A, B, p_c, p_n
parameters are known inputs, nothing guessed; gradient-wind and pressure profiles validated
analytically, including the closed-form cyclostrophic V_max).
Theory and references: Holland (1980) "An analytic model of the wind and pressure
profiles in hurricanes", Mon. Wea. Rev. 108, 1212-1218. ACES help manual
(ACESManual.rtf, "Holland hurricane wind model"). Shares the gradient-wind engine
used by CHESS-QC 9-1 (Bathystrophic Storm Surge).

Holland's two-parameter model:
    radial pressure         p(r) = p_c + dP * exp[ -(R_max / r)^B ],   dP = p_n - p_c
    radius of maximum wind  R_max = A^(1/B)               (so (R_max/r)^B = A / r^B)
    gradient wind           V_gr(r) = sqrt[ (B/rho_a)(R_max/r)^B dP e^{-(R_max/r)^B}
                                            + (r f / 2)^2 ] - r f / 2
    cyclostrophic wind      V_c(r)  = sqrt[ (B/rho_a)(R_max/r)^B dP e^{-(R_max/r)^B} ]
                            (the f -> 0 limit; the maximum, at r = R_max, is the
                             classic V_max = sqrt( B dP / (rho_a e) ).)
A is a length-scaling parameter; B (~0.5-2.5, "peakedness") controls how sharply the
wind peaks at R_max. The user fixes any two of {A, B, R_max} and the third follows
from R_max = A^(1/B); CHESS-QC exposes a `solve_for` selector for this.

Self-containment: zero sibling imports; embeds its own contract dataclasses. Runnable
standalone:
    python chessqc_1_6_holland_hurricane.py
which runs the analytic self-tests (cyclostrophic V_max, pressure limits, the
R_max <-> A,B identity) then prints the default example. numpy + stdlib only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

G_SI = 9.80665
RHO_AIR = 1.15            # kg/m^3 (warm, humid tropical air)
OMEGA = 7.2921159e-5      # rad/s (earth rotation)
_E = math.e
_HPA = 100.0             # Pa per hPa (= mb)
_KM = 1000.0
_MPS_TO_KT = 1.0 / 0.514444
_MPS_TO_MPH = 1.0 / 0.44704


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


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    note: str = ""           # hover definition shown on the output label


APP_META = AppMeta(
    aces_id="1-6",
    name="Holland Hurricane Wind Model",
    area="Wave Prediction",
    classification="exact",
    cite="Holland (1980) Mon. Wea. Rev. 108; ACES manual",
    default_system="SI",
)

_SOLVE = ("R_max (from A, B)", "A (from R_max, B)", "B (from R_max, A)")

INPUTS = (
    Field("solve_for", "Solve for", "choice", "", "", default="R_max (from A, B)",
          choices=_SOLVE, note="fix any two of {A, B, R_max}; the third is computed"),
    Field("A", "Scaling parameter A", "float", "m^B", "m^B", default=30000.0 ** 1.5,
          lo=1e-3, hi=1e12, note="length-scaling parameter; R_max = A^(1/B)"),
    Field("B", "Peakedness parameter B", "float", "", "", default=1.5, lo=0.5, hi=2.5,
          note="shape factor, typically ~1-2.5"),
    Field("R_max", "Radius of maximum wind", "float", "km", "nmi", default=30.0 * _KM,
          lo=1e3, hi=2e5, note="R_max = A^(1/B); used when solving for A or B"),
    Field("pc", "Central pressure", "float", "hPa", "hPa", default=940.0 * _HPA,
          lo=850.0 * _HPA, hi=1010.0 * _HPA, note="storm central pressure p_c"),
    Field("pn", "Peripheral pressure", "float", "hPa", "hPa", default=1013.0 * _HPA,
          lo=990.0 * _HPA, hi=1030.0 * _HPA, note="ambient pressure p_n; dP = p_n - p_c > 0"),
    Field("lat", "Latitude", "angle", "deg", "deg", default=20.0, lo=0.0, hi=80.0,
          note="for the Coriolis term in the gradient-wind balance"),
    Field("rho_air", "Air density", "float", "kg/m^3", "kg/m^3", default=RHO_AIR,
          lo=1.0, hi=1.3, note="ambient air density (standard value may be changed)"),
    Field("r_plot", "Maximum plot radius", "float", "km", "nmi", default=200.0 * _KM,
          lo=5e3, hi=1e6, note="radial extent of the output profiles"),
    Field("n_points", "Profile points", "int", "", "", default=200, lo=20, hi=2000,
          note="number of radial samples for the profiles"),
)

OUTPUTS = (
    Out("U_max", "Maximum wind speed (gradient)", "km/h", "kt", "scalar",
        note="peak gradient-balance wind speed of the storm, the maximum of V_gr(r) including the Coriolis term."),
    Out("r_at_max", "Radius of maximum wind", "km", "nmi", "scalar",
        note="radius at which the gradient wind peaks, shifted slightly inside R_max by the Coriolis term."),
    Out("U_max_cyclo", "Maximum wind (cyclostrophic)", "km/h", "kt", "scalar",
        note="closed-form cyclostrophic peak wind sqrt(B dP/(rho_a e)) at r = R_max, an upper bound on the gradient wind."),
    Out("R_max_out", "R_max (computed/echoed)", "km", "nmi", "scalar",
        note="radius of maximum wind R_max = A^(1/B), computed or echoed per the solve_for selector."),
    Out("A_out", "A (computed/echoed)", "m^B", "m^B", "scalar",
        note="Holland length-scaling parameter A (units m^B), computed or echoed per the solve_for selector."),
    Out("B_out", "B (computed/echoed)", "", "", "scalar",
        note="Holland peakedness/shape parameter B (~0.5-2.5) controlling how sharply the wind peaks at R_max."),
    Out("dP", "Pressure deficit", "hPa", "hPa", "scalar",
        note="central pressure deficit dP = p_n - p_c (positive) that drives the storm circulation."),
    Out("profile_r", "Profile: radial distance", "km", "nmi", "profile",
        note="radial distance r from the storm centre at which the profile quantities are sampled."),
    Out("profile_p", "Profile: pressure", "hPa", "hPa", "profile",
        note="Holland radial surface-pressure profile p(r), rising from p_c at the eye toward p_n far away."),
    Out("profile_Vgr", "Profile: gradient wind", "km/h", "kt", "profile",
        note="gradient-wind speed V_gr(r) including the Coriolis term, the balanced wind versus radius."),
    Out("profile_Vc", "Profile: cyclostrophic wind", "km/h", "kt", "profile",
        note="cyclostrophic wind speed V_c(r), the f -> 0 limit and an upper bound on the gradient wind."),
)


@dataclass
class Result:
    U_max: float
    r_at_max: float
    U_max_cyclo: float
    R_max_out: float
    A_out: float
    B_out: float
    dP: float
    profile_r: np.ndarray
    profile_p: np.ndarray
    profile_Vgr: np.ndarray
    profile_Vc: np.ndarray
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle") or f.key not in inp:
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo:g}, {f.hi:g}] ({f.note})")


def pressure(r, Rmax, pc, dP, B):
    """Holland radial pressure (Pa) at radius r (m)."""
    r = np.maximum(r, 1.0)
    return pc + dP * np.exp(-((Rmax / r) ** B))


def gradient_wind(r, Rmax, dP, B, f, rho_a):
    """Gradient-wind speed (m/s) at radius r (m); Holland (1980) eq with Coriolis."""
    r = np.maximum(r, 1.0)
    x = (Rmax / r) ** B
    inside = (B * dP / rho_a) * x * np.exp(-x) + (r * f / 2.0) ** 2
    return np.sqrt(np.maximum(inside, 0.0)) - r * f / 2.0


def cyclostrophic_wind(r, Rmax, dP, B, rho_a):
    """Cyclostrophic-wind speed (m/s) at radius r (m); the f -> 0 limit."""
    r = np.maximum(r, 1.0)
    x = (Rmax / r) ** B
    return np.sqrt((B * dP / rho_a) * x * np.exp(-x))


def _resolve_geometry(inp: dict) -> tuple:
    """Return (A, B, R_max) honoring the solve_for selector (R_max = A^(1/B))."""
    mode = str(inp.get("solve_for", _SOLVE[0]))
    A = float(inp["A"]); B = float(inp["B"]); Rmax = float(inp["R_max"])
    if mode.startswith("R_max"):
        Rmax = A ** (1.0 / B)
    elif mode.startswith("A"):
        A = Rmax ** B
    else:  # solve for B from R_max = A^(1/B)  ->  B = ln A / ln R_max
        if Rmax <= 1.0 or A <= 1.0:
            raise ValueError("solving for B needs A > 1 and R_max > 1 (in SI base units)")
        B = math.log(A) / math.log(Rmax)
        if not (0.3 <= B <= 3.0):
            raise ValueError(f"implied B = {B:.3f} outside the physical 0.5-2.5 range")
    return A, B, Rmax


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': "Reconstructs a tropical cyclone's radial structure from the Holland (1980) "
            'two-parameter model, returning the radial pressure profile, the gradient- and '
            'cyclostrophic-wind profiles, and the maximum wind speed.',
 'methods': [{'name': 'Holland (1980) two-parameter cyclone model',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'p(r) = p_c + \\Delta p \\, '
                                    '\\exp\\left[-\\left(\\frac{R_{max}}{r}\\right)^{B}\\right]',
                             'desc': 'Holland radial surface-pressure profile; p rises '
                                     'from p_c at the eye toward p_n far away.'},
                            {'tex': '\\Delta p = p_n - p_c',
                             'desc': 'Pressure deficit (must be positive) driving the '
                                     'storm.'},
                            {'tex': 'R_{max} = A^{1/B}',
                             'desc': 'Radius of maximum wind set by the scaling parameter '
                                     'A and peakedness B; fixing any two of {A, B, R_max} '
                                     'gives the third.'},
                            {'tex': 'V_{gr}(r) = \\sqrt{\\frac{B \\, \\Delta '
                                    'p}{\\rho_a}\\left(\\frac{R_{max}}{r}\\right)^{B} '
                                    '\\exp\\left[-\\left(\\frac{R_{max}}{r}\\right)^{B}\\right] '
                                    '+ \\left(\\frac{r f}{2}\\right)^{2}} - \\frac{r f}{2}',
                             'desc': 'Gradient-wind profile including the Coriolis term; f '
                                     '= 2 omega sin(lat).'},
                            {'tex': 'V_c(r) = \\sqrt{\\frac{B \\, \\Delta '
                                    'p}{\\rho_a}\\left(\\frac{R_{max}}{r}\\right)^{B} '
                                    '\\exp\\left[-\\left(\\frac{R_{max}}{r}\\right)^{B}\\right]}',
                             'desc': 'Cyclostrophic wind (the f -> 0 limit), an upper '
                                     'bound on the gradient wind.'},
                            {'tex': 'V_{max} = \\sqrt{\\frac{B \\, \\Delta p}{\\rho_a \\, '
                                    'e}}',
                             'desc': 'Closed-form cyclostrophic maximum wind, attained at '
                                     'r = R_max.'}]}],
 'symbols': [['p(r)', 'surface pressure at radius r'],
             ['p_c', 'storm central (eye) pressure'],
             ['p_n', 'peripheral / ambient pressure'],
             ['Delta p', 'pressure deficit, p_n minus p_c'],
             ['R_max', 'radius of maximum wind'],
             ['A', 'length-scaling parameter (units m^B)'],
             ['B', 'Holland peakedness/shape parameter (~0.5-2.5)'],
             ['V_gr', 'gradient wind speed at radius r'],
             ['V_c', 'cyclostrophic wind speed (f -> 0 limit)'],
             ['rho_a', 'ambient air density'],
             ['f', 'Coriolis parameter, 2 omega sin(latitude)']],
 'references': ['Holland (1980) Mon. Wea. Rev. 108, 1212-1218',
                'ACES help manual (Holland hurricane wind model)']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Holland (1980) wind/pressure profiles for SI inputs (pressures Pa, radii m)."""
    _validate(inp)
    A, B, Rmax = _resolve_geometry(inp)
    pc = float(inp["pc"]); pn = float(inp["pn"])
    dP = pn - pc
    if dP <= 0.0:
        raise ValueError("pressure deficit dP = p_n - p_c must be positive")
    lat = math.radians(float(inp["lat"]))
    f = 2.0 * OMEGA * math.sin(lat)
    rho_a = float(inp.get("rho_air", RHO_AIR))
    r_plot = float(inp["r_plot"]); n = int(inp["n_points"])

    r = np.linspace(max(r_plot / n, 100.0), r_plot, n)
    p = pressure(r, Rmax, pc, dP, B)
    Vgr = gradient_wind(r, Rmax, dP, B, f, rho_a)
    Vc = cyclostrophic_wind(r, Rmax, dP, B, rho_a)

    # gradient-wind maximum: cyclostrophic peak sits exactly at r = R_max; Coriolis
    # shifts it slightly inward, so locate it on a fine grid bracketing R_max.
    rr = np.linspace(0.3 * Rmax, 3.0 * Rmax, 4000)
    vgr_fine = gradient_wind(rr, Rmax, dP, B, f, rho_a)
    imax = int(np.argmax(vgr_fine))
    U_max = float(vgr_fine[imax]); r_at_max = float(rr[imax])
    U_max_cyclo = math.sqrt(B * dP / (rho_a * _E))  # analytic cyclostrophic peak

    notes = [f"Holland B={B:.3f}, R_max={Rmax/_KM:.1f} km, dP={dP/_HPA:.1f} hPa",
             f"V_max(gradient)={U_max*_MPS_TO_KT:.1f} kt, V_max(cyclostrophic)="
             f"{U_max_cyclo*_MPS_TO_KT:.1f} kt"]
    return Result(U_max=U_max, r_at_max=r_at_max, U_max_cyclo=U_max_cyclo,
                  R_max_out=Rmax, A_out=A, B_out=B, dP=dP,
                  profile_r=r, profile_p=p, profile_Vgr=Vgr, profile_Vc=Vc,
                  notes="; ".join(notes))


# --- self-tests (analytic: Holland limits + R_max identity) ---------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}
    r = compute(base)
    dP = 73.0 * _HPA

    # 1) cyclostrophic peak equals sqrt(B dP/(rho_a e)) exactly
    assert _approx(r.U_max_cyclo, math.sqrt(1.5 * dP / (RHO_AIR * _E)), 1e-9), r.U_max_cyclo
    # 2) the cyclostrophic profile peaks at r = R_max
    Vc = cyclostrophic_wind(r.profile_r, r.R_max_out, r.dP, r.B_out, RHO_AIR)
    r_peak = r.profile_r[int(np.argmax(Vc))]
    assert _approx(r_peak, r.R_max_out, 0.02), (r_peak, r.R_max_out)
    # 3) gradient wind is below cyclostrophic (Coriolis subtracts), and U_max < U_max_cyclo
    assert r.U_max < r.U_max_cyclo, (r.U_max, r.U_max_cyclo)
    assert np.all(r.profile_Vgr <= r.profile_Vc + 1e-9)
    # 4) pressure limits: p -> p_c at the eye, p -> p_n far away
    assert _approx(float(pressure(1.0, r.R_max_out, 940.0 * _HPA, r.dP, r.B_out)),
                   940.0 * _HPA, 1e-6)
    assert _approx(float(pressure(1e8, r.R_max_out, 940.0 * _HPA, r.dP, r.B_out)),
                   1013.0 * _HPA, 1e-6)
    # 5) R_max identity round-trips through all three solve modes
    rA = compute({**base, "solve_for": "A (from R_max, B)", "R_max": 30.0 * _KM, "B": 1.5})
    assert _approx(rA.A_out, (30.0 * _KM) ** 1.5, 1e-9), rA.A_out
    rB = compute({**base, "solve_for": "B (from R_max, A)", "R_max": 30.0 * _KM,
                  "A": (30.0 * _KM) ** 1.5})
    assert _approx(rB.B_out, 1.5, 1e-9), rB.B_out
    # consistency: all three modes give the same R_max for a consistent (A,B,R_max) set
    assert _approx(r.R_max_out, 30.0 * _KM, 1e-6), r.R_max_out
    # 6) deeper storm (larger dP) gives a stronger V_max (monotone)
    rdeep = compute({**base, "pc": 920.0 * _HPA})
    assert rdeep.U_max > r.U_max
    print("  self-tests: PASS (cyclostrophic V_max, pressure limits, R_max=A^(1/B) "
          "round-trip, gradient<cyclostrophic, dP monotonicity)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  storm: p_c=940 hPa, p_n=1013 hPa (dP={r.dP/_HPA:.0f} hPa), B={r.B_out:.2f}, "
          f"R_max={r.R_max_out/_KM:.0f} km, lat=20 deg")
    print("  OUTPUTS:")
    print(f"    Maximum wind (gradient)      U_max      = {r.U_max:6.2f} m/s "
          f"({r.U_max*_MPS_TO_KT:.1f} kt, {r.U_max*_MPS_TO_MPH:.0f} mph) at r={r.r_at_max/_KM:.1f} km")
    print(f"    Maximum wind (cyclostrophic) U_max_cyc  = {r.U_max_cyclo:6.2f} m/s "
          f"({r.U_max_cyclo*_MPS_TO_KT:.1f} kt)")
    print(f"    R_max = A^(1/B)              R_max      = {r.R_max_out/_KM:6.1f} km "
          f"(A={r.A_out:.3e}, B={r.B_out:.2f})")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
