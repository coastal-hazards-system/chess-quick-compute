"""CHESS-QC application 2-4 — Wave Parameters.

Originating ACES grouping: 2-4 "Wave Parameters" (functional area: Wave Theory).
A linear (Airy) wave-theory engine presented as the classic SPM dimensionless
"wave parameter" table, extended with two field-engineering conveniences:

  1. Period OR frequency input (a wave gauge often reports frequency).
  2. Pressure-transducer height inversion: given the dynamic pressure amplitude
     measured by a bottom- or mid-column gauge, recover the wave height H.

Classification: exact (closed-form linear theory; same core as 2-1).
Theory and references: Airy (1845); SPM (1984) Appendix C wave-parameter functions;
explicit dispersion via Hunt (1979). The kinematics equations are identical to
CHESS-QC 2-1 (docs/EQUATIONS.md, TR chapter 2-1); this application adds the
dimensionless shoaling/pressure-response/mass-transport quantities that the SPM
tabulates as functions of relative depth.

Self-containment: zero sibling imports; embeds its own Field/AppMeta/Out/Result
dataclasses and the Hunt dispersion solver. Runnable standalone:
    python chessqc_2_4_wave_parameters.py
which runs the analytic-limit self-tests, the pressure-inversion round-trip, and a
cross-check that the shared outputs reproduce the validated 2-1 example, then prints
an ACES-style tabulation. stdlib + numpy only.

Validation strategy (no dedicated User's Guide section for this app): the shared
linear-theory outputs (L, C, Cg, E, P, Ur, and the point kinematics) are required to
reproduce the User's Guide Example 2-1 numbers (so this app inherits the 2-1 oracle),
plus a pressure-inversion round-trip and the deep/shallow analytic limits.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

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
    status: str = "Current"      # operational currency: Current | Screening only | ...
    superseded_by: str = ""


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
    note: str = ""                 # hover definition shown on the output label


# --- application metadata --------------------------------------------------------
APP_META = AppMeta(
    aces_id="2-4",
    name="Wave Parameters",
    area="Wave Theory",
    classification="exact",
    cite="Airy (1845); SPM (1984) App. C; Hunt (1979)",
    default_system="US",     # opens on the 2-1 example inputs (US units)
)

_FT = 0.3048
# Default example reuses the validated 2-1 User's Guide case (H=6.30 ft, T=8 s,
# d=20 ft, z=-12 ft, phase 270 deg <=> X/L=0.75). The shared outputs must match 2-1.
INPUTS = (
    Field("mode", "Solve for", "choice", "", "", default="Forward (height given)",
          choices=("Forward (height given)", "Invert (pressure to height)"),
          note="Invert: recover H from a measured dynamic pressure amplitude at z"),
    Field("spec", "Wave specified by", "choice", "", "", default="Period (s)",
          choices=("Period (s)", "Frequency (1/s)"),
          note="interpret the value below as period T or frequency f = 1/T"),
    Field("Tf", "Period or frequency", "float", "s or 1/s", "s or 1/s", default=8.0,
          lo=1e-4, hi=1e4, note="> 0 (period in s, or frequency in 1/s per the choice)"),
    Field("H",  "Wave height", "float", "m", "ft", default=6.30 * _FT, lo=0.0, hi=1e4,
          note="forward mode only; > 0"),
    Field("p_gauge", "Measured dynamic pressure amplitude", "float", "Pa", "psf",
          default=0.0, lo=0.0, hi=1e9,
          note="invert mode only; amplitude (not incl. hydrostatic) of the wave pressure at z"),
    Field("d",  "Water depth", "float", "m", "ft", default=20.0 * _FT, lo=1e-6, hi=1e5,
          note="> 0"),
    Field("z",  "Vertical coordinate / gauge elevation", "float", "m", "ft",
          default=-12.0 * _FT,
          note="from SWL (z=0), +up; bed at -d. Also the pressure-gauge elevation in invert mode"),
    Field("theta_deg", "Phase angle", "angle", "deg", "deg", default=270.0, lo=0.0, hi=360.0,
          note="0 deg = crest; 270 deg <=> X/L = 0.75"),
    Field("water", "Water type", "choice", "", "", default="Salt", choices=("Salt", "Fresh"),
          note="sets water density used for energy / pressure"),
)

OUTPUTS = (
    # --- the SPM wave-parameter table (dimensionless + dimensional) ---
    Out("f",      "Frequency",                       "1/s",  "1/s",  "scalar",
        note="Wave frequency f = 1/T, the number of wave cycles per second."),
    Out("L",      "Wave length",                     "m",    "ft",   "scalar",
        note="Wavelength L = cT, the horizontal distance between successive wave crests at depth d."),
    Out("L0",     "Deepwater wave length",           "m",    "ft",   "scalar",
        note="Deepwater wavelength L0 = gT^2/(2 pi), the wavelength the wave would have in deep water."),
    Out("C",      "Celerity",                        "m/s",  "ft/s", "scalar",
        note="Phase celerity c = L/T, the speed at which an individual wave crest propagates at depth d."),
    Out("C0",     "Deepwater celerity",              "m/s",  "ft/s", "scalar",
        note="Deepwater phase celerity c0 = gT/(2 pi), the crest speed in deep water."),
    Out("Cg",     "Group velocity",                  "m/s",  "ft/s", "scalar",
        note="Group velocity Cg = n*c, the speed at which wave energy propagates."),
    Out("Cg0",    "Deepwater group velocity",        "m/s",  "ft/s", "scalar",
        note="Deepwater group velocity Cg0 = c0/2, the energy-transport speed in deep water."),
    Out("n",      "Group/phase ratio n = Cg/C",      "",     "",     "scalar",
        note="Ratio n = Cg/c of group to phase speed; ranges from 1/2 in deep water to 1 in shallow water."),
    Out("Cg_C0",  "Cg / C0",                         "",     "",     "scalar",
        note="Group velocity normalized by the deepwater celerity, Cg/c0."),
    Out("kd",     "Relative depth k*d",              "",     "",     "scalar",
        note="Relative depth kd = 2 pi d/L, the dimensionless depth governing the dispersion regime."),
    Out("d_L",    "d/L",                             "",     "",     "scalar",
        note="Depth-to-wavelength ratio d/L; >1/2 is deep water, <1/20 is shallow water."),
    Out("d_L0",   "d/L0",                            "",     "",     "scalar",
        note="Ratio of depth to deepwater wavelength, d/L0."),
    Out("tanh_kd","tanh(kd)",                        "",     "",     "scalar",
        note="Hyperbolic tangent of the relative depth kd, the dispersion factor (c^2 = (g/k) tanh(kd))."),
    Out("sinh_kd","sinh(kd)",                        "",     "",     "scalar",
        note="Hyperbolic sine of the relative depth kd, used in the kinematic depth functions."),
    Out("cosh_kd","cosh(kd)",                        "",     "",     "scalar",
        note="Hyperbolic cosine of the relative depth kd, used in the kinematic and pressure depth functions."),
    Out("csch_kd","1 / sinh(kd)",                    "",     "",     "scalar",
        note="Hyperbolic cosecant 1/sinh(kd) of the relative depth kd."),
    Out("sech_kd","1 / cosh(kd)",                    "",     "",     "scalar",
        note="Hyperbolic secant 1/cosh(kd) of the relative depth kd, equal to the bed pressure-response factor."),
    Out("Ks",     "Shoaling coefficient H/H0",       "",     "",     "scalar",
        note="Shoaling coefficient Ks = sqrt(Cg0/Cg) = H/H0 from conservation of energy flux."),
    Out("H0",     "Equivalent deepwater height H0",  "m",    "ft",   "scalar",
        note="Equivalent (unrefracted) deepwater wave height H0 = H/Ks."),
    Out("steep_ratio", "Steepness ratio (H/L)/(H0/L0)", "",  "",     "scalar",
        note="Ratio of local to deepwater wave steepness, (H/L)/(H0/L0) = Ks/tanh(kd)."),
    Out("Kp_z",   "Pressure response factor at z",   "",     "",     "scalar",
        note="Pressure-response (depth-attenuation) factor Kp(z) = cosh(k(z+d))/cosh(kd) at elevation z."),
    Out("Kp_bed", "Pressure response factor at bed", "",     "",     "scalar",
        note="Pressure-response factor at the bed, Kp(-d) = 1/cosh(kd)."),
    Out("E",      "Energy density",                  "N/m",  "lb/ft","scalar",
        note="Mean total wave energy density E = rho g H^2/8 per unit surface area."),
    Out("E0",     "Deepwater energy density",        "N/m",  "lb/ft","scalar",
        note="Deepwater mean energy density E0 = rho g H0^2/8 per unit surface area."),
    Out("Ek",     "Kinetic energy per crest length", "N-m/m","lb-ft/ft","scalar",
        note="Kinetic energy per unit crest length over one wavelength, equal to half the total (equipartition)."),
    Out("Etot",   "Total energy per crest length",   "N-m/m","lb-ft/ft","scalar",
        note="Total wave energy per unit crest length over one wavelength, Etot = E*L."),
    Out("P",      "Energy flux (wave power)",        "N/s",  "lb/s", "scalar",
        note="Energy flux (wave power) P = E*Cg per unit crest length."),
    Out("P0",     "Deepwater wave power",            "N/s",  "lb/s", "scalar",
        note="Deepwater energy flux P0 = E0*Cg0 per unit crest length."),
    Out("Ur",     "Ursell parameter",                "",     "",     "scalar",
        note="Ursell parameter Ur = H L^2/d^3; large values (> ~26) indicate linear theory is questionable."),
    # --- point quantities at (z, theta) ---
    Out("eta",    "Surface elevation",               "m",    "ft",   "point",
        note="Water-surface elevation eta = (H/2) cos(theta) above SWL at the given phase angle."),
    Out("p",      "Pressure (total)",                "Pa",   "psf",  "point",
        note="Total pressure at (z, theta): hydrostatic -rho g z plus the dynamic wave pressure."),
    Out("p_dyn",  "Dynamic pressure amplitude at z", "Pa",   "psf",  "point",
        note="Amplitude of the wave-induced (dynamic) pressure at z, p_d = rho g (H/2) Kp(z)."),
    Out("xi",     "Horizontal displacement",         "m",    "ft",   "point",
        note="Horizontal particle displacement from its mean position at (z, theta)."),
    Out("zeta",   "Vertical displacement",           "m",    "ft",   "point",
        note="Vertical particle displacement from its mean position at (z, theta)."),
    Out("u",      "Horizontal velocity",             "m/s",  "ft/s", "point",
        note="Horizontal water-particle velocity at (z, theta); positive in the wave-propagation direction."),
    Out("w",      "Vertical velocity",               "m/s",  "ft/s", "point",
        note="Vertical water-particle velocity at (z, theta); positive upward."),
    Out("dudt",   "Horizontal acceleration",         "m/s^2","ft/s^2","point",
        note="Local horizontal water-particle acceleration at (z, theta)."),
    Out("dwdt",   "Vertical acceleration",           "m/s^2","ft/s^2","point",
        note="Local vertical water-particle acceleration at (z, theta); positive upward."),
    Out("Us_z",   "Mass-transport velocity at z",    "m/s",  "ft/s", "point",
        note="Time-mean Stokes-drift (mass-transport) velocity at elevation z, directed downwave."),
    Out("Us_surf","Mass-transport velocity at SWL",  "m/s",  "ft/s", "point",
        note="Time-mean Stokes-drift velocity at the still-water level (z=0), the maximum in the column."),
    Out("Us_bed", "Mass-transport velocity at bed",  "m/s",  "ft/s", "point",
        note="Time-mean Stokes-drift velocity at the bed (z=-d), the minimum in the column."),
    Out("H_used", "Wave height used",                "m",    "ft",   "point",
        note="Wave height H used in the computation: the input value, or the value inverted from measured pressure."),
    # --- profile arrays (for plots), over X in [-L, L] ---
    Out("profile_X",   "Profile: X (+/- one wavelength)", "m",   "ft",   "profile",
        note="Horizontal coordinate X spanning plus/minus one wavelength, the abscissa of the profile plots."),
    Out("profile_eta", "Profile: surface elevation",      "m",   "ft",   "profile",
        note="Snapshot of surface elevation eta = (H/2) cos(kX) over one wavelength either side of the point."),
    Out("profile_u",   "Profile: horizontal velocity",    "m/s", "ft/s", "profile",
        note="Snapshot of horizontal particle velocity at the given z along the X profile."),
    Out("profile_w",   "Profile: vertical velocity",      "m/s", "ft/s", "profile",
        note="Snapshot of vertical particle velocity at the given z along the X profile."),
)


@dataclass
class Result:
    # SPM wave-parameter table
    f: float; L: float; L0: float; C: float; C0: float; Cg: float; Cg0: float
    n: float; Cg_C0: float
    kd: float; d_L: float; d_L0: float
    tanh_kd: float; sinh_kd: float; cosh_kd: float; csch_kd: float; sech_kd: float
    Ks: float; H0: float; steep_ratio: float; Kp_z: float; Kp_bed: float
    E: float; E0: float; Ek: float; Etot: float; P: float; P0: float; Ur: float
    # point quantities
    eta: float; p: float; p_dyn: float; xi: float; zeta: float
    u: float; w: float; dudt: float; dwdt: float
    Us_z: float; Us_surf: float; Us_bed: float; H_used: float
    # profile arrays
    profile_X: np.ndarray
    profile_eta: np.ndarray
    profile_u: np.ndarray
    profile_w: np.ndarray
    notes: str = ""


# --- dispersion solver (Hunt 1979 explicit Pade; TR 2-1 eq 11) ------------------
_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_celerity(T: float, d: float, g: float = G_SI) -> float:
    """Explicit linear celerity c (m/s) via Hunt (1979), accuracy < 0.01%."""
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    c2 = g * d / (y + 1.0 / denom)
    return math.sqrt(c2)


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue                                  # choice/bool fields carry no bounds
        v = float(inp.get(f.key, f.default))          # H/p_gauge are mode-specific (optional)
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes the full linear (Airy) wave-parameter table for a given '
            'period/frequency and depth: wavelength, celerity, group velocity, shoaling '
            'and pressure-response factors, energy and power, the Ursell number, and point '
            'kinematics. It can either take wave height as input or invert a measured '
            'dynamic-pressure amplitude to recover the height.',
 'methods': [{'name': 'Linear (Airy) wave parameters with shoaling and pressure response',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'c^{2} = \\frac{g}{k}\\tanh(kd)',
                             'desc': 'Linear dispersion relation; solved explicitly via '
                                     'the 9-term Hunt (1979) Pade form, then L = cT.'},
                            {'tex': 'C_{g} = \\frac{c}{2}\\left[1 + '
                                    '\\frac{2kd}{\\sinh(2kd)}\\right]',
                             'desc': 'Group velocity; the bracketed factor over 2 is the '
                                     'ratio n = C_g/c.'},
                            {'tex': 'K_{s} = \\sqrt{\\frac{C_{g0}}{C_{g}}}',
                             'desc': 'Shoaling coefficient H/H_0 from conservation of '
                                     'energy flux (deepwater C_{g0} = c_0/2).'},
                            {'tex': 'K_{p}(z) = \\frac{\\cosh(k(z+d))}{\\cosh(kd)}',
                             'desc': 'Pressure-response (depth-attenuation) factor at '
                                     'elevation z; at the bed it reduces to 1/cosh(kd).'},
                            {'tex': 'H = \\frac{2\\,p_{d}}{\\rho\\,g\\,K_{p}(z)}',
                             'desc': 'Invert mode: recover wave height from the measured '
                                     'dynamic-pressure amplitude p_d at the gauge '
                                     'elevation z.'}]}],
 'symbols': [['H', 'Wave height (= 2a)'],
             ['T', 'Wave period (or input frequency f = 1/T)'],
             ['d', 'Still-water depth'],
             ['z',
              'Vertical coordinate / gauge elevation from SWL, positive up (bed at -d)'],
             ['k', 'Wave number, k = 2 pi / L'],
             ['c', 'Phase celerity, c = L/T'],
             ['C_g', 'Group velocity; n = C_g/c'],
             ['K_s', 'Shoaling coefficient, H/H_0'],
             ['K_p', 'Pressure response (depth attenuation) factor'],
             ['p_d', 'Dynamic (wave-induced) pressure amplitude at z'],
             ['rho', 'Water density (salt or fresh)'],
             ['U_r', 'Ursell parameter, H L^2 / d^3']],
 'references': ['Airy (1845)',
                'SPM (1984) App. C',
                'Hunt (1979)',
                'Dean & Dalrymple (1984)',
                'Ursell (1953); Stokes (1847)']}


def compute(inp: dict, *, g: float = G_SI, rho: Optional[float] = None,
            n_profile: int = 201) -> Result:
    """Wave-parameter results for SI inputs.

    Required keys: mode, spec, Tf, d, z, theta_deg, water; plus H (forward) or
    p_gauge (invert). `inp` values are SI (the GUI converts US->SI at the edge)."""
    _validate(inp)
    if rho is None:
        rho = RHO_FRESH if str(inp.get("water", "Salt")) == "Fresh" else RHO_SALT

    spec = str(inp.get("spec", "Period (s)"))
    Tf = float(inp["Tf"])
    T = (1.0 / Tf) if spec.startswith("Frequency") else Tf
    if T <= 0.0:
        raise ValueError("period/frequency must give T > 0")
    d = float(inp["d"]); z = float(inp["z"])
    theta = math.radians(float(inp.get("theta_deg", 0.0)))

    omega = 2.0 * math.pi / T
    f = 1.0 / T
    C = wave_celerity(T, d, g)            # eq 11 (=> eq 10 dispersion)
    L = C * T                             # eq 12
    k = 2.0 * math.pi / L
    kd = k * d
    sinh_kd = math.sinh(kd); cosh_kd = math.cosh(kd); tanh_kd = math.tanh(kd)
    n = 0.5 * (1.0 + 2.0 * kd / math.sinh(2.0 * kd))
    Cg = n * C                            # eq 13

    # deepwater references and shoaling
    C0 = g * T / (2.0 * math.pi)
    L0 = g * T * T / (2.0 * math.pi)
    Cg0 = 0.5 * C0
    Cg_C0 = Cg / C0
    Ks = math.sqrt(Cg0 / Cg)              # shoaling coefficient H/H0
    d_L = d / L; d_L0 = d / L0
    csch_kd = 1.0 / sinh_kd
    sech_kd = 1.0 / cosh_kd

    # pressure response factor (depth attenuation) at z and at the bed
    z_use = min(max(z, -d), 0.0)          # gauge is in the water column
    Kp_z = math.cosh(k * (z_use + d)) / cosh_kd
    Kp_bed = 1.0 / cosh_kd

    # resolve the wave height: forward (given) or inverted from measured pressure
    mode = str(inp.get("mode", "Forward (height given)"))
    if mode.startswith("Invert"):
        p_gauge = float(inp.get("p_gauge", 0.0))
        # dynamic pressure amplitude p_d = rho g (H/2) Kp(z)  ->  H = 2 p_d /(rho g Kp)
        H = 2.0 * p_gauge / (rho * g * Kp_z)
    else:
        H = float(inp["H"])

    E = rho * g * H * H / 8.0             # eq 15: energy density (per unit surface area)
    P = E * Cg                            # eq 16: energy flux (wave power)
    Ur = H * L * L / (d ** 3)             # eq 24
    H0 = H / Ks                           # equivalent (unrefracted) deepwater height
    steep_ratio = Ks / tanh_kd           # (H/L)/(H0/L0)
    E0 = rho * g * H0 * H0 / 8.0          # deepwater energy density
    Etot = E * L                          # total energy per unit crest length over one L
    Ek = 0.5 * Etot                       # kinetic energy = half (equipartition)
    P0 = E0 * Cg0                         # deepwater wave power

    # point kinematics at (z, theta); clamp z to the column
    zk = min(max(z, -d), 0.5 * H)
    s = zk + d
    cosh_ks = math.cosh(k * s); sinh_ks = math.sinh(k * s)
    ct, st = math.cos(theta), math.sin(theta)

    eta = 0.5 * H * ct                                            # eq 14
    p_dyn = 0.5 * rho * g * H * (cosh_ks / cosh_kd)              # amplitude of wave pressure
    p = -rho * g * zk + p_dyn * ct                               # eq 17 (total)
    xi = -0.5 * H * (cosh_ks / sinh_kd) * st                     # eq 18
    zeta = 0.5 * H * (sinh_ks / sinh_kd) * ct                    # eq 19
    u = (math.pi * H / T) * (cosh_ks / sinh_kd) * ct             # eq 20
    w = (math.pi * H / T) * (sinh_ks / sinh_kd) * st             # eq 21
    dudt = (2.0 * math.pi ** 2 * H / T ** 2) * (cosh_ks / sinh_kd) * st   # eq 22
    dwdt = -(2.0 * math.pi ** 2 * H / T ** 2) * (sinh_ks / sinh_kd) * ct  # eq 23

    # mass-transport (Stokes drift) velocity, SPM (1984) eq 2-55:
    #   Us(z) = (pi H / L)^2 * C * cosh(2 k (z+d)) / (2 sinh^2(kd))
    def _drift(zz: float) -> float:
        ss = min(max(zz, -d), 0.0) + d
        return (math.pi * H / L) ** 2 * C * math.cosh(2.0 * k * ss) / (2.0 * sinh_kd ** 2)
    Us_z = _drift(zk); Us_surf = _drift(0.0); Us_bed = _drift(-d)

    # profile arrays over +/- one wavelength (for plots; snapshot at the given z)
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
    if mode.startswith("Invert"):
        notes.append(f"H recovered from measured pressure = {H / _FT:.3f} ft")
    if z + d < 0:
        notes.append("WARNING: point/gauge below the bed (z + d < 0)")
    elif z > 0.5 * H:
        notes.append("point above the surface; kinematics evaluated at the surface")
    if Ur > 26:
        notes.append("Ursell > 26: linear theory questionable (consider 2-2/2-3)")

    return Result(
        f=f, L=L, L0=L0, C=C, C0=C0, Cg=Cg, Cg0=Cg0, n=n, Cg_C0=Cg_C0,
        kd=kd, d_L=d_L, d_L0=d_L0,
        tanh_kd=tanh_kd, sinh_kd=sinh_kd, cosh_kd=cosh_kd, csch_kd=csch_kd, sech_kd=sech_kd,
        Ks=Ks, H0=H0, steep_ratio=steep_ratio, Kp_z=Kp_z, Kp_bed=Kp_bed,
        E=E, E0=E0, Ek=Ek, Etot=Etot, P=P, P0=P0, Ur=Ur,
        eta=eta, p=p, p_dyn=p_dyn, xi=xi, zeta=zeta, u=u, w=w, dudt=dudt, dwdt=dwdt,
        Us_z=Us_z, Us_surf=Us_surf, Us_bed=Us_bed, H_used=H,
        profile_X=X, profile_eta=prof_eta, profile_u=prof_u, profile_w=prof_w,
        notes="; ".join(notes))


# --- self-tests (analytic limits + round-trip + 2-1 cross-check) ----------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    g = G_SI
    base = dict(mode="Forward (height given)", spec="Period (s)", water="Salt")

    # 1) Deep water: C -> gT/2pi, L -> L0, Cg -> C/2, Ks -> 1
    r = compute({**base, "Tf": 8.0, "H": 1.0, "d": 1000.0, "z": 0.0, "theta_deg": 0.0}, g=g)
    assert _approx(r.C, g * 8.0 / (2 * math.pi), 2e-3), r.C
    assert _approx(r.L, r.L0, 2e-3), (r.L, r.L0)
    assert _approx(r.Cg, 0.5 * r.C, 2e-3), (r.Cg, r.C)
    assert _approx(r.Ks, 1.0, 5e-3), r.Ks
    assert _approx(r.n, 0.5, 5e-3), r.n

    # 2) Shallow water: C -> sqrt(g d), Cg -> C, n -> 1
    r = compute({**base, "Tf": 60.0, "H": 0.2, "d": 2.0, "z": 0.0, "theta_deg": 0.0}, g=g)
    assert _approx(r.C, math.sqrt(g * 2.0), 5e-3), r.C
    assert _approx(r.Cg, r.C, 5e-3), (r.Cg, r.C)
    assert _approx(r.n, 1.0, 1e-2), r.n

    # 3) Frequency input equals the equivalent period input
    rp = compute({**base, "spec": "Period (s)", "Tf": 10.0, "H": 1.5, "d": 12.0,
                  "z": -3.0, "theta_deg": 0.0}, g=g)
    rf = compute({**base, "spec": "Frequency (1/s)", "Tf": 0.1, "H": 1.5, "d": 12.0,
                  "z": -3.0, "theta_deg": 0.0}, g=g)
    assert _approx(rp.L, rf.L, 1e-9) and _approx(rp.f, rf.f, 1e-9), (rp.L, rf.L)

    # 4) Pressure-inversion round trip: forward p_dyn at z -> invert -> recover H
    fwd = compute({**base, "Tf": 8.0, "H": 6.30 * _FT, "d": 20.0 * _FT,
                   "z": -12.0 * _FT, "theta_deg": 0.0}, g=g)
    inv = compute({**base, "mode": "Invert (pressure to height)", "Tf": 8.0,
                   "p_gauge": fwd.p_dyn, "d": 20.0 * _FT, "z": -12.0 * _FT,
                   "theta_deg": 0.0}, g=g)
    assert _approx(inv.H_used, 6.30 * _FT, 1e-9), (inv.H_used, 6.30 * _FT)

    # 5) Power identity and Stokes-drift sign (drift is downwave-positive everywhere)
    assert _approx(fwd.P, fwd.E * fwd.Cg, 1e-9)
    assert fwd.Us_surf > fwd.Us_bed > 0.0

    # 6) extended SPM table identities: Etot = 2 Ek, sech(kd) = Kp(bed), Cg/C0 = n*(C/C0)
    assert _approx(fwd.Etot, 2.0 * fwd.Ek, 1e-12)
    assert _approx(fwd.sech_kd, fwd.Kp_bed, 1e-12)
    assert _approx(fwd.Cg_C0, fwd.n * fwd.C / fwd.C0, 1e-12)
    assert _approx(fwd.P0, fwd.E0 * fwd.Cg0, 1e-12)

    print("  self-tests: PASS (deep/shallow limits, f<->T, pressure round-trip, P=E*Cg, drift)")


def _cross_check_21() -> None:
    """Shared linear-theory outputs must reproduce the validated 2-1 example
    (User's Guide Example 2-1): phase 270 deg <=> X/L = 0.75."""
    r = compute(dict(mode="Forward (height given)", spec="Period (s)", Tf=8.0,
                     H=6.30 * _FT, p_gauge=0.0, d=20.0 * _FT, z=-12.0 * _FT,
                     theta_deg=270.0, water="Salt"))
    for attr, exp_us, fac, tol in (
            ("L", 189.90, _FT, 0.05), ("C", 23.74, _FT, 0.02), ("Cg", 20.87, _FT, 0.02),
            ("Ur", 28.40, 1.0, 0.05), ("w", -0.93, _FT, 0.01), ("xi", 4.59, _FT, 0.02),
            ("dudt", -2.83, _FT, 0.02)):
        got = getattr(r, attr) / fac
        assert abs(got - exp_us) <= tol, f"{attr}: got {got:.3f}, manual 2-1 {exp_us:.3f}"
    assert abs(r.p / 47.880259 - 767.83) <= 0.5, r.p / 47.880259
    print("  2-1 cross-check: PASS (inherits User's Guide Example 2-1 oracle)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        vv = inp[f.key]
        sval = f"{vv:>10.4g}" if isinstance(vv, (int, float)) else f"{vv:>10}"
        print(f"    {f.label:34s} {f.key:9s} = {sval} {f.unit_si}")
    print("  OUTPUTS:")
    by_key = {o.key: o for o in OUTPUTS}
    for kk in ("f", "L", "L0", "C", "C0", "Cg", "Cg0", "n", "Cg_C0", "kd", "d_L", "d_L0",
               "tanh_kd", "sinh_kd", "cosh_kd", "csch_kd", "sech_kd", "Ks", "H0", "steep_ratio",
               "Kp_z", "Kp_bed", "E", "E0", "Ek", "Etot", "P", "P0", "Ur",
               "eta", "p", "p_dyn", "u", "w", "Us_z", "Us_surf", "Us_bed"):
        o = by_key[kk]
        print(f"    {o.label:34s} {kk:11s} = {getattr(r, kk):>12.5g} {o.unit_si}")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _cross_check_21()
    _print_default_example()
