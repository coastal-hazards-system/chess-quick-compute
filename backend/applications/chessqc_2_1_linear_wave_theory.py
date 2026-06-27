"""CHESS-QC application 2-1 — Linear (Airy) Wave Theory.

Originating ACES application: 2-1 "Linear Wave Theory" (functional area: Wave Theory).
First-order (small-amplitude / sinusoidal / Airy) approximations of wave motion.

Classification: exact (full closed-form linear theory).
Theory and references: Airy (1845); explicit dispersion via Hunt (1979).
  Equations transcribed in docs/EQUATIONS.md, TR chapter 2-1 (eqs 1-24):
    dispersion (10), Hunt explicit c (11), L=cT (12), C_g (13), eta (14),
    E (15), P (16), pressure (17), displacements (18-19), velocities (20-21),
    accelerations (22-23), Ursell (24).  s = z + d ; theta = k*x - omega*t.

Self-containment: zero sibling imports; embeds its own Field/AppMeta/Result
dataclasses and dispersion solver. Runnable standalone:
    python chessqc_2_1_linear_wave_theory.py
which runs the analytic-limit self-tests, then prints an ACES-style tabulation of
the default example.  stdlib + numpy only.

I/O mirrors the ACES manual's *exact* Linear Wave Theory lists (the GUI must capture
all of them):
  Inputs : H (>0), T (>0), d (>0), z (no restriction), X/L (0-1)   [+ unit system]
  Outputs: eta, C, L, C_g, P, E, U_r, p, (xi,zeta), (u,w), (du/dt,dw/dt)
           + profile arrays (eta, u, w vs X over +/- one wavelength) for plots.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

# --- standard physical constants (overridable; SI internal) ---------------------
G_SI = 9.80665           # m/s^2
RHO_SALT = 1025.18       # kg/m^3 (sea water)
RHO_FRESH = 999.0        # kg/m^3 (fresh water)


# --- embedded contract dataclasses (self-contained; identical across all apps) --
@dataclass(frozen=True)
class AppMeta:
    aces_id: str
    name: str
    area: str
    classification: str          # "exact" | "standard" | "provisional"
    cite: str
    default_system: str = "SI"   # unit system the GUI opens in ("SI" | "US")


@dataclass(frozen=True)
class Field:
    """One GUI/contract input field."""
    key: str
    label: str
    kind: str = "float"            # float | int | choice | bool | angle | file
    unit_si: str = ""
    unit_us: str = ""
    default: object = 0.0
    lo: float = -math.inf
    hi: float = math.inf
    choices: tuple = ()
    note: str = ""


@dataclass(frozen=True)
class Out:
    """One output descriptor (metadata; values live in Result)."""
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"           # scalar | point | profile
    note: str = ""           # hover definition shown on the output label


# --- application metadata --------------------------------------------------------
APP_META = AppMeta(
    aces_id="2-1",
    name="Linear Wave Theory",
    area="Wave Theory",
    classification="exact",
    cite="Airy (1845); Hunt (1979); TR 2-1",
    default_system="US",     # opens on the User's Guide Example (US units)
)

# Complete input list (ACES manual "Input requirements" for Linear Wave Theory).
# Defaults are SI-internal; values shown are the ACES User's Guide Example 2-1
# (H=6.30 ft, T=8 s, d=20 ft, z=-12 ft, X/L=0.75 -> L=189.90 ft ...).
_FT = 0.3048
INPUTS = (
    Field("H",  "Wave height",        "float", "m", "ft",  default=6.30 * _FT,  lo=1e-6, hi=1e4,
          note="> 0"),
    Field("T",  "Wave period",        "float", "s", "s",   default=8.0,  lo=1e-3, hi=1e4,
          note="> 0"),
    Field("d",  "Water depth",        "float", "m", "ft",  default=20.0 * _FT, lo=1e-6, hi=1e5,
          note="> 0"),
    Field("z",  "Vertical coordinate","float", "m", "ft",  default=-12.0 * _FT,
          note="from SWL (z=0), +up; no restriction (clamped to [-d, eta])"),
    Field("xL", "Wavelength fraction (X/L)", "float", "", "", default=0.75, lo=0.0, hi=1.0,
          note="0.0 to 1.0 (phase position; 0 = crest)"),
    Field("water", "Water type", "choice", "", "", default="Salt", choices=("Salt", "Fresh"),
          note="sets water density used for energy/pressure"),
)

# Complete output list (ACES manual "Output" for Linear Wave Theory).
OUTPUTS = (
    Out("L",     "Wave length",                  "m",     "ft",    "scalar",
        note="Horizontal distance between successive crests, L = c*T from the linear dispersion relation."),
    Out("C",     "Celerity",                     "m/s",   "ft/s",  "scalar",
        note="Wave phase speed c at which the waveform propagates, from Hunt's explicit dispersion solution."),
    Out("Cg",    "Group velocity",               "m/s",   "ft/s",  "scalar",
        note="Speed at which wave energy travels, C_g = (c/2)(1 + 2kd/sinh(2kd))."),
    Out("E",     "Energy density",               "N/m",   "lb/ft", "scalar",
        note="Average total wave energy per unit surface area, E = rho*g*H^2/8."),
    Out("P",     "Energy flux (power)",          "N/s",   "lb/s",  "scalar",
        note="Average rate of wave-energy transport per unit crest width, P = E*C_g."),
    Out("Ur",    "Ursell parameter",             "",      "",      "scalar",
        note="Nonlinearity measure U_r = H*L^2/d^3; linear theory becomes questionable above about 26."),
    Out("kd",    "Relative depth k*d",           "",      "",      "scalar",
        note="Dimensionless relative depth k*d distinguishing shallow, intermediate, and deep water."),
    Out("eta",   "Surface elevation",            "m",     "ft",    "point",
        note="Water-surface displacement above SWL at the chosen phase, eta = (H/2)cos(theta), +up."),
    Out("p",     "Pressure",                     "Pa",    "psf",   "point",
        note="Total fluid pressure at depth z and phase: hydrostatic plus the dynamic wave contribution."),
    Out("xi",    "Horizontal displacement",      "m",     "ft",    "point",
        note="Horizontal excursion of a water particle from its mean position at depth z and phase."),
    Out("zeta",  "Vertical displacement",        "m",     "ft",    "point",
        note="Vertical excursion of a water particle from its mean position at depth z and phase, +up."),
    Out("u",     "Horizontal velocity",          "m/s",   "ft/s",  "point",
        note="Horizontal water-particle velocity at depth z and phase, positive in the direction of propagation."),
    Out("w",     "Vertical velocity",            "m/s",   "ft/s",  "point",
        note="Vertical water-particle velocity at depth z and phase, positive upward."),
    Out("dudt",  "Horizontal acceleration",      "m/s^2", "ft/s^2","point",
        note="Local horizontal acceleration of a water particle at depth z and phase."),
    Out("dwdt",  "Vertical acceleration",        "m/s^2", "ft/s^2","point",
        note="Local vertical acceleration of a water particle at depth z and phase, positive upward."),
    Out("profile_X",   "Profile: X (+/- one wavelength)", "m",   "ft",   "profile",
        note="Horizontal coordinate X spanning +/- one wavelength, the x-axis for the profile plots."),
    Out("profile_eta", "Profile: surface elevation",      "m",   "ft",   "profile",
        note="Surface elevation eta versus X over +/- one wavelength, (H/2)cos(kX)."),
    Out("profile_u",   "Profile: horizontal velocity",    "m/s", "ft/s", "profile",
        note="Horizontal particle velocity u versus X at the chosen depth over +/- one wavelength."),
    Out("profile_w",   "Profile: vertical velocity",      "m/s", "ft/s", "profile",
        note="Vertical particle velocity w versus X at the chosen depth over +/- one wavelength."),
)


@dataclass
class Result:
    # scalars
    L: float;  C: float;  Cg: float;  E: float;  P: float;  Ur: float;  kd: float
    # point kinematics at (z, X/L)
    eta: float; p: float; xi: float; zeta: float
    u: float;  w: float;  dudt: float; dwdt: float
    # profile arrays (for plots), over X in [-L, L]
    profile_X: np.ndarray
    profile_eta: np.ndarray
    profile_u: np.ndarray
    profile_w: np.ndarray
    notes: str = ""


# --- dispersion solver (Hunt 1979 explicit Pade; TR 2-1 eq 11) ------------------
_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_celerity(T: float, d: float, g: float = G_SI) -> float:
    """Explicit linear celerity c (m/s) via Hunt (1979), accuracy < 0.01%.
    c^2 = g*d * [ y + (1 + sum_{n=1}^9 d_n y^n)^-1 ]^-1 ,  y = omega^2 d / g."""
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    c2 = g * d / (y + 1.0 / denom)
    return math.sqrt(c2)


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue                                  # choice/bool fields carry no bounds
        v = inp[f.key]
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes first-order (Airy/linear) wave properties from wave height, period, '
            'and water depth: wavelength, celerity, group velocity, energy density and '
            'flux, the Ursell parameter, plus the pressure and particle kinematics '
            '(displacements, velocities, accelerations) at a chosen depth and phase.',
 'methods': [{'name': 'Linear (Airy) wave theory',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': '\\omega^{2} = g\\,k\\,\\tanh(kd)',
                             'desc': 'Linear dispersion relation linking radian frequency, '
                                     'wavenumber, and depth.'},
                            {'tex': 'c^{2} = g\\,d\\,[\\,y + (1 + \\sum_{n=1}^{9} '
                                    'd_{n}\\,y^{n})^{-1}\\,]^{-1}, \\quad y = \\omega^{2} '
                                    'd / g',
                             'desc': 'Hunt (1979) explicit Pade celerity solving the '
                                     'dispersion relation (accuracy < 0.01%).'},
                            {'tex': 'C_{g} = \\frac{c}{2}\\,(1 + \\frac{2kd}{\\sinh(2kd)})',
                             'desc': 'Group velocity (energy transport speed).'},
                            {'tex': 'E = \\frac{1}{8}\\,\\rho\\,g\\,H^{2}',
                             'desc': 'Average wave energy density per unit surface area.'},
                            {'tex': 'P = E\\,C_{g}',
                             'desc': 'Wave energy flux (power) per unit crest width.'},
                            {'tex': 'U_{r} = \\frac{H\\,L^{2}}{d^{3}}',
                             'desc': 'Ursell parameter; linear theory becomes questionable '
                                     'above about 26.'}]}],
 'symbols': [['H', 'Wave height (= 2a)'],
             ['T', 'Wave period'],
             ['d', 'Still-water depth'],
             ['L', 'Wavelength (L = cT)'],
             ['c', 'Wave celerity (phase speed)'],
             ['k', 'Wavenumber (2*pi/L)'],
             ['omega', 'Radian frequency (2*pi/T)'],
             ['C_g', 'Group velocity'],
             ['rho', 'Water density (salt or fresh)'],
             ['U_r', 'Ursell parameter']],
 'references': ['Airy (1845)',
                'Hunt (1979)',
                'SPM (1984)',
                'Dean & Dalrymple (1984)',
                'Stokes (1847); Ursell (1953)',
                'ACES TR 2-1']}


def compute(inp: dict, *, g: float = G_SI, rho: float | None = None,
            n_profile: int = 201) -> Result:
    """Linear-wave-theory results for SI inputs {H, T, d, z, xL, water}.

    `inp` values are in SI (the GUI converts US->SI at the edge). `z` is clamped to
    [-d, eta_crest] for kinematics. Water density follows the `water` field
    (Salt|Fresh) unless `rho` is given explicitly. Returns a fully-populated Result
    (all manual outputs + profile arrays)."""
    _validate(inp)
    if rho is None:
        rho = RHO_FRESH if str(inp.get("water", "Salt")) == "Fresh" else RHO_SALT
    H, T, d, z, xL = (float(inp[k]) for k in ("H", "T", "d", "z", "xL"))

    omega = 2.0 * math.pi / T
    C = wave_celerity(T, d, g)            # eq 11 (=> eq 10 dispersion)
    L = C * T                             # eq 12
    k = 2.0 * math.pi / L
    kd = k * d
    Cg = 0.5 * C * (1.0 + 2.0 * kd / math.sinh(2.0 * kd))   # eq 13

    E = rho * g * H * H / 8.0             # eq 15
    P = E * Cg                            # eq 16
    Ur = H * L * L / (d ** 3)             # eq 24

    theta = 2.0 * math.pi * xL            # phase angle (X/L fraction at a snapshot)

    # clamp z to the fluid column for kinematics (z from SWL, +up; bed at -d)
    z_use = min(max(z, -d), 0.5 * H)
    s = z_use + d                         # height above seabed
    sinh_kd = math.sinh(kd)
    cosh_kd = math.cosh(kd)
    cosh_ks = math.cosh(k * s)
    sinh_ks = math.sinh(k * s)
    ct, st = math.cos(theta), math.sin(theta)

    eta = 0.5 * H * ct                                            # eq 14
    p = -rho * g * z_use + 0.5 * rho * g * H * (cosh_ks / cosh_kd) * ct   # eq 17
    xi = -0.5 * H * (cosh_ks / sinh_kd) * st                     # eq 18
    zeta = 0.5 * H * (sinh_ks / sinh_kd) * ct                    # eq 19
    u = (math.pi * H / T) * (cosh_ks / sinh_kd) * ct             # eq 20
    w = (math.pi * H / T) * (sinh_ks / sinh_kd) * st             # eq 21
    dudt = (2.0 * math.pi ** 2 * H / T ** 2) * (cosh_ks / sinh_kd) * st   # eq 22
    dwdt = -(2.0 * math.pi ** 2 * H / T ** 2) * (sinh_ks / sinh_kd) * ct  # eq 23

    # profile arrays over +/- one wavelength (manual: plots vs X)
    X = np.linspace(-L, L, n_profile)
    th = k * X
    prof_eta = 0.5 * H * np.cos(th)
    prof_u = (math.pi * H / T) * (cosh_ks / sinh_kd) * np.cos(th)
    prof_w = (math.pi * H / T) * (sinh_ks / sinh_kd) * np.sin(th)

    notes = []
    rel = d / L
    if rel > 0.5:
        notes.append("deep water (d/L > 1/2)")
    elif rel < 0.05:
        notes.append("shallow water (d/L < 1/20)")
    else:
        notes.append("intermediate depth")
    if z + d < 0:
        notes.append("WARNING: point outside waveform (z below the bed, z + d < 0)")
    elif z > eta:
        notes.append("WARNING: point above the water surface (z > eta); kinematics at the surface")
    elif z != z_use:
        notes.append(f"z clamped to {z_use:.3f} m for kinematics")
    if Ur > 26:
        notes.append("Ursell > 26: linear theory questionable")

    return Result(L=L, C=C, Cg=Cg, E=E, P=P, Ur=Ur, kd=kd,
                  eta=eta, p=p, xi=xi, zeta=zeta, u=u, w=w, dudt=dudt, dwdt=dwdt,
                  profile_X=X, profile_eta=prof_eta, profile_u=prof_u, profile_w=prof_w,
                  notes="; ".join(notes))


# --- self-tests (analytic limits; "input data that produce correct answers") ----
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    g = G_SI
    # 1) Deep water: c -> gT/2pi, L -> gT^2/2pi, Cg -> c/2
    T, d, H = 8.0, 1000.0, 1.0
    r = compute({"H": H, "T": T, "d": d, "z": 0.0, "xL": 0.0}, g=g)
    c_deep = g * T / (2 * math.pi)
    assert _approx(r.C, c_deep, 2e-3), (r.C, c_deep)
    assert _approx(r.L, g * T * T / (2 * math.pi), 2e-3), r.L
    assert _approx(r.Cg, 0.5 * r.C, 2e-3), (r.Cg, r.C)

    # 2) Shallow water: c -> sqrt(g d), Cg -> c
    T, d = 60.0, 2.0
    r = compute({"H": 0.2, "T": T, "d": d, "z": 0.0, "xL": 0.0}, g=g)
    assert _approx(r.C, math.sqrt(g * d), 5e-3), (r.C, math.sqrt(g * d))
    assert _approx(r.Cg, r.C, 5e-3), (r.Cg, r.C)

    # 3) Dispersion consistency: omega^2 == g k tanh(kd)
    T, d = 10.0, 12.0
    r = compute({"H": 1.5, "T": T, "d": d, "z": -3.0, "xL": 0.0}, g=g)
    k = 2 * math.pi / r.L
    omega2 = (2 * math.pi / T) ** 2
    assert _approx(omega2, g * k * math.tanh(k * d), 1e-3), (omega2, g * k * math.tanh(k * d))

    # 4) Power identity P == E*Cg ; Ursell == H L^2/d^3
    assert _approx(r.P, r.E * r.Cg, 1e-9)
    assert _approx(r.Ur, 1.5 * r.L ** 2 / d ** 3, 1e-9)

    # 5) Crest kinematics: at X/L=0 (theta=0) horizontal velocity is max, w=0
    assert r.w == 0.0 or abs(r.w) < 1e-9
    assert r.u > 0.0

    print("  self-tests: PASS (deep/shallow limits, dispersion, P=E*Cg, Ursell, crest)")


def _print_default_example() -> None:
    """ACES-style tabulation of the default example (SI).
    Defaults are the ACES User's Guide Example 2-1 (US units; stored SI-internal)."""
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        vv = inp[f.key]
        sval = f"{vv:>10.4g}" if isinstance(vv, (int, float)) and f.kind != "choice" else f"{vv:>10}"
        print(f"    {f.label:28s} {f.key:4s} = {sval} {f.unit_si}")
    print("  OUTPUTS:")
    scal = {"L": r.L, "C": r.C, "Cg": r.Cg, "E": r.E, "P": r.P, "Ur": r.Ur, "kd": r.kd}
    pts = {"eta": r.eta, "p": r.p, "xi": r.xi, "zeta": r.zeta,
           "u": r.u, "w": r.w, "dudt": r.dudt, "dwdt": r.dwdt}
    by_key = {o.key: o for o in OUTPUTS}
    for kk, vv in {**scal, **pts}.items():
        o = by_key[kk]
        print(f"    {o.label:28s} {kk:5s} = {vv:>12.5g} {o.unit_si}")
    print(f"  notes: {r.notes}")
    print(f"  profile arrays: X/eta/u/w each length {len(r.profile_X)} (over +/- one L)")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
