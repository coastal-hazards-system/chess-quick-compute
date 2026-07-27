"""CHESS-QC application 1-5 — Near-surface Wind Speeds.

Originating ACES application: 1-5 "Near-surface Wind Speeds" (functional area: Wave
Prediction; a later ACES addition). Given the geostrophic wind speed, the air-sea
temperature difference, latitude, and a height z, it returns the friction velocity,
the wind speed at height z, the drag coefficients (at z and at 10 m), the surface
roughness length, the Monin-Obukhov stability length, the stability function, and
the surface momentum flux.

Classification: standard. The full ACES planetary-boundary-layer (PBL) resistance law
(TR 1-1 eqs. 14--19) is solved simultaneously with the constant-stress equations. Silva
(2005, eqs. 7.16--7.19) records Sherlock's 1996 revision and constants: A0=0.8,
B0=B1=3.5, C=-7.0, and C2=0.0144/980.

Theory and references:
  * ACES full-PBL geostrophic drag law (TR 1-1 eqs. 14--18), including its
    stability-dependent A and B functions and cross-isobar angle.
  * Surface-layer Monin-Obukhov profile (ACES TR 1-1, eqs 5-13; shared with CHESS-QC
    1-1): U_z = (U_*/k) [ ln(z/z_0) - Psi(z/L') ], with the ACES sea-surface roughness
        z_0 = C_1/U_* + C_2 U_*^2 + C_3   (cgs: C_1=0.1525, C_2=0.019/980, C_3=-0.00371)
    the bulk Obukhov length L' = 1.79 (U_*^2/dT)[ln(z/z_0) - Psi] (TR eq 8), and the
        ACES KEYPS/Lumley-Panofsky stability function Psi.
  * Silva Casarin (2005), Analisis y descripcion estadistica del oleaje, eqs. 7.16--7.19
    (publishes the revised relation and values attributed to Sherlock, 1996).
  * ACES help manual (ACESManual.rtf, "Near-surface wind speeds").

Method. U_*, L', Psi, and the full-PBL stability functions are iterated together to the
ACES convergence criteria. The input convention is dT = T_air - T_sea, exactly as in
TR eq. 8, so L' and the PBL stability parameter mu=k U_*/(f L') retain their source signs.
There is no ACES numerical worked example for this late-added utility, so the
implementation is checked by its governing-equation residuals, neutral reduction, and
physical monotonicity.

Self-containment: zero sibling imports; embeds its own contract dataclasses. Runnable
standalone:
    python chessqc_1_5_near_surface_wind.py
which runs the analytic self-tests (neutral log-profile recovery, C_D band, U_* << |V_g|,
monotonicity) then prints the default example. stdlib only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

G_SI = 9.80665
RHO_AIR = 1.20            # kg/m^3
OMEGA = 7.2921159e-5      # rad/s (earth rotation)
K_VON_KARMAN = 0.40
# ACES full-PBL constants.  Silva (2005), eqs. 7.16--7.19, publishes the revised
# relation and values attributed to Sherlock (1996); ACES TR 1-1 supplies eqs. 5--19.
_A0, _B0, _B1 = 0.8, 3.5, 3.5
_PSI_STABLE_C = -7.0
# ACES sea-surface roughness constants (cgs: U_* in cm/s, z_0 in cm).  C2 is the
# full-PBL value specified by TR 1-1-6, rather than the constant-stress-only value.
_C1, _C2, _C3 = 0.1525, 0.0144 / 980.0, -0.00371
_KT = 0.514444
_MPS_TO_KT = 1.0 / _KT


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
    aces_id="1-5",
    name="Near-surface Wind Speeds",
    area="Wave Prediction",
    classification="standard",
    cite="ACES TR 1-1 eqs. 5-19; Silva (2005) eqs. 7.16-7.19; Lumley & Panofsky (1964)",
    default_system="SI",
)

INPUTS = (
    Field("Ug", "Geostrophic wind speed", "float", "m/s", "kt", default=30.0, lo=1.0, hi=120.0,
          note="free-atmosphere geostrophic wind |V_g| > 0"),
    Field("deltaT", "Air-sea temperature difference", "float", "deg C", "deg C", default=0.0,
          lo=-20.0, hi=20.0, note="dT = T_air - T_sea; <0 unstable (warm sea), >0 stable"),
    Field("lat", "Latitude", "angle", "deg", "deg", default=40.0, lo=1.0, hi=80.0,
          note="for the Coriolis parameter f = 2 Omega sin(lat) > 0"),
    Field("z", "Height above surface", "float", "m", "ft", default=10.0, lo=0.1, hi=300.0,
          note="elevation z at which the wind speed is reported"),
    Field("water", "Water type", "choice", "", "", default="Salt", choices=("Salt", "Fresh"),
          note="affects air density only weakly; retained for parity with ACES"),
    Field("rho_air", "Air density", "float", "kg/m^3", "kg/m^3", default=RHO_AIR,
          lo=1.0, hi=1.3, note="standard value may be changed"),
)

OUTPUTS = (
    Out("u_star", "Friction velocity U*", "m/s", "kt", "scalar",
        note="Friction velocity U_* = sqrt(tau/rho_a), the surface shear velocity scaling the turbulent momentum flux."),
    Out("Uz", "Wind speed at height z", "m/s", "kt", "scalar",
        note="Mean wind speed at the chosen height z from the stability-corrected logarithmic surface-layer profile."),
    Out("U10", "Wind speed at 10 m", "m/s", "kt", "scalar",
        note="Mean wind speed at the standard 10 m reference height from the same surface-layer profile."),
    Out("CDz", "Drag coefficient at z", "", "", "scalar",
        note="Drag coefficient referenced to height z, C_Dz = (U_*/U_z)^2."),
    Out("CD", "Drag coefficient at 10 m", "", "", "scalar",
        note="Drag coefficient referenced to 10 m, C_D = (U_*/U10)^2; typically in the 1.0-2.5e-3 band."),
    Out("z0", "Surface roughness length", "m", "ft", "scalar",
        note="Aerodynamic sea-surface roughness length z_0, the height where the extrapolated log wind profile reaches zero."),
    Out("L", "Monin-Obukhov length", "m", "ft", "scalar",
        note="Monin-Obukhov stability length L'; infinite for neutral, negative for unstable, positive for stable stratification."),
    Out("psi", "Stability function Psi(z/L)", "", "", "scalar",
        note="Businger-Dyer momentum stability correction Psi(z/L'); 0 neutral, >0 unstable, <0 stable."),
    Out("tau", "Surface momentum flux", "Pa", "Pa", "scalar",
        note="Surface momentum flux (wind stress) tau = rho_a U_*^2 transferred from air to the sea surface."),
    Out("alpha", "Cross-isobar angle", "deg", "deg", "scalar",
        note="Cross-isobar angle between the geostrophic wind aloft and the surface wind/stress, sin(alpha) = B U_*/(k |V_g|)."),
)


@dataclass
class Result:
    u_star: float
    Uz: float
    U10: float
    CDz: float
    CD: float
    z0: float
    L: float
    psi: float
    tau: float
    alpha: float
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle") or f.key not in inp:
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo:g}, {f.hi:g}] ({f.note})")


def _z0_cgs(u_star_cgs: float) -> float:
    """ACES sea-surface roughness z_0 (cm) from U_* (cm/s); TR 1-1 eq 6/7."""
    return max(_C1 / u_star_cgs + _C2 * u_star_cgs * u_star_cgs + _C3, 1e-4)


def _psi_m(zeta: float) -> float:
    """ACES KEYPS/Lumley-Panofsky momentum similarity function (TR 1-1 eqs. 9--11)."""
    if abs(zeta) < 1e-12:
        return 0.0
    if zeta > 0.0:
        return _PSI_STABLE_C * zeta
    rz = zeta
    for _ in range(100):
        new_rz = zeta * (1.0 - 18.0 * rz) ** 0.25
        if abs(new_rz - rz) < 1e-13:
            rz = new_rz
            break
        rz = new_rz
    phi = 1.0 / (1.0 - 18.0 * rz) ** 0.25
    return (1.0 - phi - 3.0 * math.log(phi) + 2.0 * math.log((1.0 + phi) / 2.0)
            + 2.0 * math.atan(phi) - math.pi / 2.0 + math.log((1.0 + phi * phi) / 2.0))


def _solve_u_star(G_cgs: float, f: float, A: float, B: float) -> float:
    """Solve ACES TR 1-1 eq. 14 for U* with prescribed PBL A and B."""
    def residual(us):
        z0 = _z0_cgs(us)
        q = K_VON_KARMAN * G_cgs / us
        return math.log(G_cgs / (f * z0)) - (A - math.log(us / G_cgs)
                                              + math.sqrt(q * q - B * B))
    lo = max(1e-4, K_VON_KARMAN * G_cgs / 1e5)
    hi = K_VON_KARMAN * G_cgs / B * (1.0 - 1e-12)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _surface_stability(u_star: float, dT: float, z: float, z0: float) -> tuple[float, float]:
    """Couple ACES eqs. 5, 8--11 for a specified friction velocity (cgs)."""
    if abs(dT) < 1e-12:
        return math.inf, 0.0
    psi = 0.0
    L = math.inf
    for _ in range(200):
        L = 1.79 * (u_star * u_star / dT) * (math.log(z / z0) - psi)
        new_psi = _psi_m(z / L)
        if abs(new_psi - psi) < 1e-11:
            return L, new_psi
        psi = 0.5 * (psi + new_psi)
    return L, psi


def _pbl_functions(u_star: float, f: float, L: float) -> tuple[float, float]:
    """Revised ACES PBL functions (Silva 2005 eqs. 7.16--7.18; Sherlock 1996)."""
    if not math.isfinite(L):
        return _A0, _B0
    mu = K_VON_KARMAN * u_star / (f * L)
    if mu <= 0.0:  # unstable/neutral: L' <= 0 for dT = T_air - T_sea <= 0
        a_fac = 1.0 - math.exp(0.015 * mu)
        b_fac = 1.0 - math.exp(0.03 * mu)
        return _A0 + (_PSI_STABLE_C - _A0) * a_fac, _B0 - _B1 * b_fac
    return _A0 - 0.96 * math.sqrt(mu) + math.log1p(mu), _B0 + 0.7 * math.sqrt(mu)


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Given the free-atmosphere geostrophic wind, air-sea temperature difference, '
            'latitude, and a height z, this app solves the planetary-boundary-layer '
            'resistance law for the friction velocity and then evaluates the '
            'stability-corrected surface-layer profile. It returns U*, the wind speed at z '
            'and at 10 m, the drag coefficients, the sea-surface roughness length, the '
            'Monin-Obukhov length and stability function, the surface momentum flux, and '
            'the cross-isobar angle.',
 'methods': [{'name': 'Geostrophic-drag resistance law with Monin-Obukhov surface-layer '
                      'profile',
              'when': None,
              'tag': '',
              'note': 'Friction velocity, the PBL resistance coefficients A and B, '
                      'the Obukhov length, and the ACES KEYPS stability function are '
                      'iterated together.  The revised neutral solution reduces to A=0.8 and '
                      'B=3.5.',
              'equations': [{'tex': '| V_g | = '
                                    '\\frac{U_*}{k}\\sqrt{\\left(\\ln\\frac{U_*}{f\\,z_0} '
                                    '- A\\right)^2 + B^2}',
                             'desc': 'ACES full-PBL geostrophic-drag (Rossby-number '
                                     'similarity) law.  It is solved with the '
                                     'stability-dependent PBL coefficients A and B; '
                                     'A=0.8 and B=3.5 for neutral conditions.'},
                            {'tex': 'z_0 = \\frac{C_1}{U_*} + C_2\\,U_*^2 + C_3',
                             'desc': 'ACES full-PBL sea-surface roughness length as a '
                                     'function of U* (cgs constants C_1=0.1525, C_2=0.0144/980, '
                                     'C_3=-0.00371).'},
                            {'tex': 'U_z = \\frac{U_*}{k}\\left[\\ln\\frac{z}{z_0} - '
                                    "\\Psi\\!\\left(\\frac{z}{L'}\\right)\\right]",
                             'desc': 'Stability-corrected logarithmic surface-layer wind '
                                     'profile; the same form gives U10 (z = 10 m).'},
                            {'tex': "L' = 1.79\\,\\frac{U_*^2}{\\Delta "
                                    'T}\\left[\\ln\\frac{z}{z_0} - \\Psi\\right]',
                             'desc': 'Bulk Monin-Obukhov stability length set by the '
                                     'air-sea temperature difference; solved '
                                     'self-consistently with Psi.'},
                            {'tex': '\\sin\\alpha = \\frac{B\\,U_*}{k\\,| V_g |}',
                             'desc': 'Cross-isobar angle between the geostrophic wind and '
                                     'the surface stress.'},
                            {'tex': '\\tau = \\rho_a\\,U_*^2',
                             'desc': 'Surface momentum flux; the drag coefficient at any '
                                     'level is C_D = (U_* / U)^2.'}]}],
 'symbols': [['V_g', 'Free-atmosphere geostrophic wind speed (magnitude)'],
             ['U_*', 'Friction velocity (surface shear velocity)'],
             ['U_z', 'Wind speed at height z (and U10 at 10 m)'],
             ['z_0', 'Sea-surface roughness length'],
             ['k', 'Von Karman constant (approx 0.40)'],
             ['f', 'Coriolis parameter, f = 2 Omega sin(lat)'],
             ["L'", 'Monin-Obukhov stability (Obukhov) length'],
             ['Psi', "ACES KEYPS/Lumley-Panofsky momentum stability function of z/L'"],
             ['Delta T', 'Air-sea temperature difference (T_air - T_sea)'],
             ['C_D',
              'Drag coefficient, (U_*/U)^2; tau is the momentum flux, rho_a air density, '
              'alpha cross-isobar angle']],
 'references': ['ACES Technical Reference 1-1 (eqs 5-19)',
                'Silva Casarin (2005), Analisis y descripcion estadistica del oleaje, eqs 7.16-7.19',
                'Lumley & Panofsky (1964)',
                'ACES help manual, Near-surface wind speeds']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Near-surface wind from the geostrophic wind, dT, latitude, and height (SI in/out)."""
    _validate(inp)
    G = float(inp["Ug"]) * 100.0                # cm/s (cgs internally, matching z_0 eq)
    dT = float(inp["deltaT"])
    f = 2.0 * OMEGA * math.sin(math.radians(float(inp["lat"])))
    z = float(inp["z"]) * 100.0                 # cm
    rho_a = float(inp.get("rho_air", RHO_AIR))

    # ACES TR 1-1 couples the surface layer to the full PBL resistance law.
    # Begin from its revised neutral reduction (A=A0, B=B0), then jointly iterate U*, L',
    # Psi, and the stability-dependent PBL coefficients (eqs. 14--19).
    u_star = _solve_u_star(G, f, _A0, _B0)      # cm/s
    L, psi_z, A, B = math.inf, 0.0, _A0, _B0
    for _ in range(200):
        z0 = _z0_cgs(u_star)
        L, psi_z = _surface_stability(u_star, dT, z, z0)
        A, B = _pbl_functions(u_star, f, L)
        solved = _solve_u_star(G, f, A, B)
        if abs(solved - u_star) < 1e-7 * max(1.0, u_star):
            u_star = solved
            z0 = _z0_cgs(u_star)
            L, psi_z = _surface_stability(u_star, dT, z, z0)
            A, B = _pbl_functions(u_star, f, L)
            break
        # Under-relaxation retains the source fixed-point method while making the
        # very stable/unstable limits robust for interactive input ranges.
        u_star = 0.5 * (u_star + solved)
    else:
        raise RuntimeError("ACES PBL iteration did not converge")

    psi_10 = _psi_m(1000.0 / L) if math.isfinite(L) else 0.0
    Uz = (u_star / K_VON_KARMAN) * (math.log(z / z0) - psi_z)              # cm/s
    U10 = (u_star / K_VON_KARMAN) * (math.log(1000.0 / z0) - psi_10)       # cm/s (10 m=1000 cm)
    CDz = (u_star / Uz) ** 2 if Uz > 0 else float("nan")
    CD = (u_star / U10) ** 2 if U10 > 0 else float("nan")
    sin_alpha = max(min(B * u_star / (K_VON_KARMAN * G), 1.0), -1.0)
    alpha = math.degrees(math.asin(sin_alpha))

    # convert to SI
    u_star_si = u_star / 100.0
    tau = rho_a * u_star_si ** 2                                          # Pa
    L_si = (L / 100.0) if math.isfinite(L) else float("inf")
    strat = "neutral" if abs(dT) < 1e-9 else ("unstable" if dT < 0 else "stable")
    notes = [f"{strat} (dT={dT:+.1f} C); U*={u_star_si:.3f} m/s, z0={z0/100.0*1000:.2f} mm",
             f"CD(10m)={CD*1e3:.2f}e-3, cross-isobar angle {abs(alpha):.1f} deg"]
    return Result(u_star=u_star_si, Uz=Uz / 100.0, U10=U10 / 100.0, CDz=CDz, CD=CD,
                  z0=z0 / 100.0, L=L_si, psi=psi_z, tau=tau, alpha=abs(alpha),
                  notes="; ".join(notes))


# --- self-tests (analytic: neutral log profile, C_D band, U* << G) --------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}
    r = compute(base)

    # 1) neutral: Psi = 0, the profile is the pure log law and L is infinite
    assert r.psi == 0.0 and not math.isfinite(r.L), (r.psi, r.L)
    z0 = r.z0
    Uz_log = (r.u_star / K_VON_KARMAN) * math.log(10.0 / z0)
    assert _approx(r.Uz, Uz_log, 1e-9), (r.Uz, Uz_log)
    assert _approx(r.U10, r.Uz, 1e-9)              # default z = 10 m

    # 2) drag coefficient lands in the observed 1.0-2.5e-3 band over the synoptic wind
    #    range (the ACES roughness relation is not intended to represent tropical-
    #    cyclone drag saturation at extreme wind speeds).
    for Ug in (10.0, 20.0, 30.0, 45.0):
        rr = compute({**base, "Ug": Ug})
        assert 1.0e-3 <= rr.CD <= 2.5e-3, (Ug, rr.CD)
        assert rr.u_star < 0.05 * Ug              # U* ~ 2.5-3% of geostrophic wind
        assert rr.u_star > 0.0

    # 3) momentum flux tau = rho_a U*^2
    assert _approx(r.tau, RHO_AIR * r.u_star ** 2, 1e-12)

    # 4) monotonicity: stronger geostrophic wind -> larger U*, z0, and 10-m wind
    r1 = compute({**base, "Ug": 20.0}); r2 = compute({**base, "Ug": 40.0})
    assert r2.u_star > r1.u_star and r2.z0 > r1.z0 and r2.U10 > r1.U10

    # 5) cross-isobar angle is positive and modest (a few to ~30 deg)
    assert 0.0 < r.alpha < 40.0, r.alpha

    # 6) stratification: a finite Obukhov length and a non-zero Psi appear for dT != 0,
    #    and the sign of Psi follows stability (unstable dT<0 -> Psi>0; stable dT>0 -> Psi<0)
    ru = compute({**base, "deltaT": -5.0}); rs = compute({**base, "deltaT": 5.0})
    assert math.isfinite(ru.L) and ru.psi > 0.0, (ru.L, ru.psi)
    assert math.isfinite(rs.L) and rs.psi < 0.0, (rs.L, rs.psi)

    # 7) exact revised PBL branches (Silva 2005 eqs. 7.16--7.18).  Use L values
    #    constructed to give prescribed mu and verify both unequal exponents and continuity.
    f = 1.0e-4; us = 50.0
    assert _pbl_functions(us, f, math.inf) == (_A0, _B0)
    mu_u = -10.0; Au, Bu = _pbl_functions(us, f, K_VON_KARMAN * us / (f * mu_u))
    assert _approx(Au, _A0 + (_PSI_STABLE_C - _A0) * (1.0 - math.exp(0.015 * mu_u)), 1e-12)
    assert _approx(Bu, _B0 - _B1 * (1.0 - math.exp(0.03 * mu_u)), 1e-12)
    mu_s = 10.0; As, Bs = _pbl_functions(us, f, K_VON_KARMAN * us / (f * mu_s))
    assert _approx(As, _A0 - 0.96 * math.sqrt(mu_s) + math.log1p(mu_s), 1e-12)
    assert _approx(Bs, _B0 + 0.7 * math.sqrt(mu_s), 1e-12)

    # 8) neutral geostrophic-drag equation (TR eq. 14) closes independently.
    G = float(base["Ug"]) * 100.0
    f = 2.0 * OMEGA * math.sin(math.radians(float(base["lat"])))
    us = r.u_star * 100.0; z0 = r.z0 * 100.0
    lhs = math.log(G / (f * z0))
    rhs = _A0 - math.log(us / G) + math.sqrt((K_VON_KARMAN * G / us) ** 2 - _B0 ** 2)
    assert _approx(lhs, rhs, 1e-10), (lhs, rhs)

    print("  self-tests: PASS (neutral log-profile recovery, C_D in 1.0-2.5e-3 band, "
          "U* ~ 3% of G, tau=rho U*^2, stability sign, exact PBL branches/residual)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  inputs: Ug=30 m/s, dT=0 C (neutral), lat=40 deg, z=10 m")
    print("  OUTPUTS:")
    print(f"    Friction velocity            U*    = {r.u_star:7.3f} m/s")
    print(f"    Wind at z                    Uz    = {r.Uz:7.2f} m/s ({r.Uz*_MPS_TO_KT:.1f} kt)")
    print(f"    Wind at 10 m                 U10   = {r.U10:7.2f} m/s")
    print(f"    Drag coefficient at z        CDz   = {r.CDz*1e3:7.3f} e-3")
    print(f"    Drag coefficient at 10 m     CD    = {r.CD*1e3:7.3f} e-3")
    print(f"    Roughness length             z0    = {r.z0*1000:7.3f} mm")
    print(f"    Monin-Obukhov length         L     = {'inf' if not math.isfinite(r.L) else f'{r.L:.1f} m'}")
    print(f"    Momentum flux                tau   = {r.tau:7.3f} Pa")
    print(f"    Cross-isobar angle           alpha = {r.alpha:7.1f} deg")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
