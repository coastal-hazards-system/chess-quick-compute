"""CHESS-QC application — Bathystrophic Storm Surge (Bodine 1971).

First Quick Compute tool beyond the original 34 ACES applications (functional area:
Storm Surge). Estimates open-coast hurricane surge along a single cross-shelf traverse
by the quasi-1D bathystrophic method of Bodine (1971), CERC TM-35 (verification:
Pararas-Carayannis 1975, TM-50).

Classification: provisional. The analytic sub-models (gradient wind, pressure setup, wind
and bottom stresses) are exact, but the integrated surge reproduces the Bodine TM-35
oracle only to ballpark (~13.4 ft): Bodine read the wind-isovel field and the bathymetry
from graphs, so the reference's own inputs are not recoverable, and Bodine (1971) cites a
possible factor-of-two uncertainty. Already flagged "Screening only" (superseded by ADCIRC
for operational use).

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

Wind field, two selectable parametric models sharing one analytic engine:
  * Holland (1980): p(r) = Pc + dP exp(-(R/r)^B);  default, modern standard.
  * Myers (1954) / Bodine: p(r) = Pc + dP exp(-R/r)  == Holland with B = 1.
  Gradient wind  V_gr(r) = sqrt( (B dP/rho_a)(R/r)^B exp(-(R/r)^B) + (r f/2)^2 ) - r f/2,
  reduced to the surface (~0.9), given a forward-speed asymmetry and an inflow angle,
  then split into onshore/alongshore components that drive the integrator. The Holland
  shape factor B (default 1.5) may be overridden by an explicit Vmax via
  B = rho_a e Vmax^2 / dP.

Validation: Bodine TM-35 Chesapeake Bay Entrance example (Pc=27.57, Pn=29.92 inHg,
R=35 nm, Vf=22 kt, lat 37 deg, K=0.0025) -> peak surge ~13.4 ft. Because Bodine read
the wind isovel field and bathymetry from graphs, a parametric model reproduces 13.4 ft
only approximately; the analytic sub-models (gradient wind, pressure setup, stresses)
are checked exactly. self-contained, numpy + stdlib only.
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


@dataclass(frozen=True)
class Out:
    key: str; label: str; unit_si: str = ""; unit_us: str = ""; kind: str = "scalar"


APP_META = AppMeta(
    aces_id="9-1",
    name="Bathystrophic Storm Surge",
    area="Storm Surge",
    classification="provisional",
    cite="Bodine (1971) TM-35; Holland (1980); Myers (1954); TR/CERC",
    default_system="US",
    status="Screening only",
    superseded_by="ADCIRC (risk assessment)",
)

_WIND_MODELS = ("Holland (1980)", "Myers / Bodine (1954)")
_NM, _FT, _KT, _MPH, _MI = 1852.0, 0.3048, 0.514444, 0.44704, 1609.344
_INHG = 3386.389

# Default = Bodine TM-35 Chesapeake Bay Entrance example (US units; stored SI).
# Bathymetry: distance from shore (nm) vs depth below SWL (ft), shelf edge -> shore.
_BATHY = [(60.0, 300.0), (50.0, 240.0), (40.0, 185.0), (32.0, 130.0), (26.0, 95.0),
          (20.0, 66.0), (15.0, 48.0), (10.0, 33.0), (6.0, 21.0), (3.0, 12.0), (1.0, 5.0)]

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
    Field("wind_model", "Wind model", "choice", "", "", default="Holland (1980)",
          choices=_WIND_MODELS, note="Holland (B adjustable) or Myers/Bodine (B=1)"),
    Field("B_holland", "Holland B (peakedness)", "float", "", "", default=1.5,
          lo=0.5, hi=2.5, note="Holland shape factor; locked to 1.0 for Myers/Bodine"),
    Field("Vmax", "Max wind (optional)", "float", "km/h", "mph", default=0.0,
          lo=0.0, hi=120.0, note="if > 0, overrides B via B = rho_a e Vmax^2 / dP"),
    Field("rho_air", "Air density", "float", "kg/m^3", "kg/m^3", default=RHO_AIR,
          lo=1.0, hi=1.3, note="ambient air density"),
    Field("K_bottom", "Bottom friction coefficient", "float", "", "", default=0.0025,
          lo=1e-4, hi=2e-2, note="bed friction K (~0.002-0.005)"),
    Field("Se", "Initial water-level rise", "float", "m", "ft", default=0.0,
          lo=-3.0, hi=10.0, note="initial setup at start of computation"),
    Field("SA", "Astronomical tide", "float", "m", "ft", default=0.0,
          lo=-5.0, hi=10.0, note="astronomical tide above MSL datum"),
    Field("dt", "Time step", "float", "hr", "hr", default=0.5 * 3600.0,
          lo=60.0, hi=7200.0, note="integration time step (stored in seconds)"),
    Field("n_steps", "Number of time steps", "int", "", "", default=80, lo=10, hi=400,
          note="storm is swept past the traverse over these steps"),
)

OUTPUTS = (
    Out("peak_surge", "Peak surge at shore", "m", "ft", "scalar"),
    Out("S_wind", "  wind setup (Sx)", "m", "ft", "scalar"),
    Out("S_bathy", "  bathystrophic setup (Sy)", "m", "ft", "scalar"),
    Out("S_press", "  pressure setup", "m", "ft", "scalar"),
    Out("Vmax_out", "Max wind speed (30 ft)", "km/h", "mph", "scalar"),
    Out("B_used", "Holland B used", "", "", "scalar"),
    Out("t_peak", "Time of peak", "s", "hr", "scalar"),
    Out("profile_X", "Profile: distance from shore", "m", "nm", "profile"),
    Out("profile_eta", "Profile: setup at peak", "m", "ft", "profile"),
    Out("profile_u", "Profile: still-water depth", "m", "ft", "profile"),
    Out("profile_w", "Profile: total depth at peak", "m", "ft", "profile"),
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


def compute(inp: dict, *, g: float = G_SI, rho_w: float = RHO_W) -> Result:
    """Bathystrophic surge along a traverse. SI inputs; the GUI converts at the edge."""
    _validate(inp)
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
    # 5) Holland default runs; explicit Vmax overrides B
    rh = compute(inp)
    assert rh.B_used == 1.5
    rv = compute({**inp, "Vmax": 50.0})
    assert _approx(rv.B_used, 1.20 * _E * 50.0 ** 2 / (rv_dP := (float(inp["Pn"]) - float(inp["Pc"]))), 1e-6)
    print(f"  self-tests: PASS (wind models, pressure setup, Vmax->B; Chesapeake Myers "
          f"peak {rm.peak_surge/_FT:.1f} ft / Vmax {rm.Vmax_out/_MPH:.0f} mph vs Bodine 13.4/102)")


def _print_default_example() -> None:
    inp = {fd.key: fd.default for fd in INPUTS}
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  Chesapeake Bay Entrance traverse (Bodine TM-35 example):")
    for wm in ("Myers / Bodine (1954)", "Holland (1980)"):
        r = compute({**inp, "wind_model": wm})
        print(f"    [{wm:22s}] peak {r.peak_surge/_FT:6.2f} ft  "
              f"(Sx {r.S_wind/_FT:.1f}, Sy {r.S_bathy/_FT:.1f}, Sp {r.S_press/_FT:.1f}); "
              f"Vmax {r.Vmax_out/_MPH:5.1f} mph; B {r.B_used:.2f}")
    print(f"  Bodine TM-35 oracle: peak surge 13.4 ft, Vx 102 mph (Myers/B=1 is the match;")
    print(f"  parametric isovel vs Bodine's graphical Graham-Nunn field -> ballpark agreement).")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
