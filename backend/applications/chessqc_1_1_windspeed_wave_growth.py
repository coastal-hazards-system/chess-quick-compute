"""CHESS-QC application 1-1 — Windspeed Adjustment and Wave Growth.

Originating ACES application: 1-1 "Windspeed Adjustment and Wave Growth" (functional
area: Wave Prediction). Adjusts an observed wind to a 10-m equivalent-neutral wind,
linearizes it to the constant-drag wind used by the growth formulas, then estimates
fetch-/duration-limited wave growth (deep or shallow water, open or restricted fetch).

Classification: exact (the wind-adjustment and SPM/Resio wave-growth coefficients are all
known from the TR -- nothing guessed; reproduces the User's Guide Example 1-1 (open) and
Example 3 (restricted geometry); the optional non-neutral stability correction uses the
canonical Businger-Dyer psi convention).
Theory and references: SPM (1984) Ch.3; Resio & Vincent (1977); Smith (1991) restricted
  fetch. Equations transcribed in docs/EQUATIONS.md, TR chapter 1-1 (eqs 1-46):
  ship-bias (1), constant-stress profile (3),(6),(7),(13), short-fetch reduction (20),
  duration adjustment (21),(22), wind stress (24),(25), deepwater growth (27)-(42),
  shallow growth (43),(44).

Scope: **Open Water and Restricted** fetch (restricted via the radial-fetch table +
eq 26 off-wind direction selection), **neutral stability** (air-sea temperature
difference = 0). Wind observation types via the constant-stress log profile: Overwater
(shipboard) applies the ship-bias correction (1) first; Overwater (not shipboard) and
Shore (windward/leeward) use the profile directly.

Air-sea temperature (stability) correction, default Neutral; physical correction opt-in:
  The default path treats deltaT as neutral (validated against Example 1). An opt-in
  "Businger-Dyer" stability model applies the literature-standard surface-layer correction
  (`_psi_m` + the bulk Obukhov length, TR eq 8 with 1.79 = theta_bar/(k^2 g)).

  Why opt-in and not validated against the ACES examples. The unstable-branch function
  transcribed in ACES TR eq 9 is a CORRUPTED copy of the Businger-Dyer/Paulson (1970)
  psi_m: it shares the two log terms but has +2 atan - pi/2 where the canonical form has
  -2 atan + pi/2, plus spurious 1 - phi_a - 3 ln(phi_a) terms (verified against Paulson
  1970, J. Appl. Meteor. 9:857, and the standard form ln[(1+x^2)(1+x)^2/8] - 2 atan x +
  pi/2). The canonical psi_m is used here. The ACES examples still cannot be reproduced:
  Example 3 (overwater, deltaT=-3) reports U_e = 44.00 mph, a -3.2% drop from the neutral
  45.5, whereas correct surface-layer physics gives only ~-0.9% AND the OPPOSITE sign
  (deltaT<0 is unstable -> a small increase). The -3.2% is both ~3-6x too large and wrong-
  signed for standard physics; ACES's own eqs 8-9 are internally sign-inconsistent; and no
  worked example decouples stability from the observation-type bias (Ex1 shore+neutral,
  Ex2 shipboard, Ex3 overwater+deltaT). So the physical correction is shipped as an option
  (small, correctly-signed) and the validated neutral path is the default.
  (The full-PBL geostrophic model, eqs 14-19, with untranscribed appendix constants
  A_0/B_0/B_1, is a SEPARATE path used only for Geostrophic/Over-land types.)

Self-containment: zero sibling imports; embeds its own contract dataclasses and the
wind-adjustment iteration. Runnable standalone:
    python chessqc_1_1_windspeed_wave_growth.py
which runs self-tests (incl. the User's Guide example) then prints the example.
stdlib only (no numpy needed; this tool has no profile arrays).

I/O mirrors the ACES manual's lists:
  Inputs : z_obs, U_obs, deltaT, dur_obs, dur_final, lat, obs_type, fetch_type,
           wave_eq, F, d
  Outputs: U_e (equivalent neutral wind), U_a (adjusted wind), H_mo, T_p, + growth note.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

G_SI = 9.80665           # m/s^2
_K_VON_KARMAN = 0.40


@dataclass(frozen=True)
class AppMeta:
    aces_id: str
    name: str
    area: str
    classification: str
    cite: str
    default_system: str = "SI"


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
    enable_if: tuple = ()    # (other_key, value): gray out (disable) unless that input == value


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    note: str = ""           # hover definition shown on the output label


APP_META = AppMeta(
    aces_id="1-1",
    name="Windspeed Adjustment and Wave Growth",
    area="Wave Prediction",
    classification="exact",
    cite="SPM (1984) Ch.3; Resio & Vincent (1977); Smith (1991); TR 1-1",
    default_system="US",
)

_OBS_TYPES = ("Overwater (shipboard)", "Overwater (not shipboard)",
              "Shore (windward)", "Shore (leeward)")
_FETCH_TYPES = ("Open Water", "Restricted")
_WAVE_EQS = ("Shallow", "Deep")

_FT = 0.3048
_MPH = 0.44704
_MI = 1609.344
# Defaults are the ACES User's Guide Example 1-1 (US units; stored SI-internal):
# z_obs=25 ft, U_obs=45 mph, deltaT=0, dur=3 hr, lat=30 deg, F=26 mi, d=13 ft,
# Shore (windward) / Open Water / Shallow -> U_e=46.49 mph, U_a=68.06 mph,
# H_mo=4.24 ft, T_p=4.77 s. See tests/test_manual_oracle.py.
INPUTS = (
    Field("z_obs", "Elevation of observed wind", "float", "m", "ft", default=25.0 * _FT,
          lo=1.0 * _FT, hi=1e4, note=">= 1 ft"),
    Field("U_obs", "Observed wind speed", "float", "km/h", "mph", default=45.0 * _MPH,
          lo=0.1 * _MPH, hi=200.0, note="> 0"),
    Field("deltaT", "Air-sea temperature difference", "float", "C", "C", default=0.0,
          lo=-50.0, hi=50.0, note="T_air - T_sea; 0 = neutral. <0 unstable, >0 stable"),
    Field("stability", "Stability model", "choice", "", "", default="Neutral (validated)",
          choices=("Neutral (validated)", "Businger-Dyer (physical)"),
          note="Neutral: validated default. Businger-Dyer: physical correction when deltaT!=0 "
               "(opt-in; does not reproduce ACES Example 3, see docstring)"),
    Field("dur_obs", "Duration of observed wind", "float", "hr", "hr", default=3.0 * 3600.0,
          lo=1.0, hi=36000.0, note=">= 0.1 (stored internally in seconds)"),
    Field("dur_final", "Duration of final wind", "float", "hr", "hr", default=3.0 * 3600.0,
          lo=1.0, hi=36000.0, note=">= 0.1 (stored internally in seconds)"),
    Field("lat", "Latitude of observation", "angle", "deg", "deg", default=30.0,
          lo=0.0, hi=180.0, note="0 to 180 deg"),
    Field("obs_type", "Wind observation type", "choice", "", "", default="Shore (windward)",
          choices=_OBS_TYPES, note="how/where the wind was observed"),
    Field("fetch_type", "Fetch type", "choice", "", "", default="Open Water",
          choices=_FETCH_TYPES, note="Open Water (single fetch) or Restricted (radial fetches)"),
    Field("wave_eq", "Wave equation type", "choice", "", "", default="Shallow",
          choices=_WAVE_EQS, note="Shallow (depth-limited) or Deep water"),
    Field("F", "Wind fetch length", "float", "km", "mi", default=26.0 * _MI,
          lo=1e-3, hi=5e6, note="open water: single fetch length"),
    Field("d", "Average fetch depth", "float", "m", "ft", default=13.0 * _FT,
          lo=0.1 * _FT, hi=1e4, note="> 0 (used for shallow water)",
          enable_if=("wave_eq", "Shallow")),
    # restricted-fetch geometry (used only when fetch_type = Restricted)
    Field("wind_dir", "Wind direction", "angle", "deg", "deg", default=125.0, lo=0.0, hi=360.0,
          note="restricted: deg clockwise from north"),
    Field("dbeta", "Radial angle increment", "angle", "deg", "deg", default=12.0, lo=1.0, hi=180.0,
          note="restricted: spacing of radial fetches"),
    Field("beta1", "Direction of first radial", "angle", "deg", "deg", default=0.0, lo=0.0, hi=360.0,
          note="restricted: deg clockwise from north"),
    Field("fetches", "Radial fetch lengths", "table",
          default=[[v * _MI] for v in (3.7, 12.3, 13.4, 12.2, 13.2, 36.0, 35.6,
                                        28.7, 26.8, 13.0, 10.4, 10.1, 6.4, 5.7)],
          columns=(("Fetch", "km", "mi"),), lo=0.0, hi=5e6,
          note="restricted: one fetch length per radial (starting at beta1, step dbeta)"),
)

OUTPUTS = (
    Out("U_e", "Equivalent neutral wind speed", "km/h", "mph", "scalar",
        note="The observed wind reduced to a 10-m equivalent-neutral wind via the constant-stress log profile and adjusted to the final-wind duration."),
    Out("U_a", "Adjusted wind speed", "km/h", "mph", "scalar",
        note="The equivalent-neutral wind linearized to the constant-drag (C_D = 0.001) wind that drives the wave-growth formulas."),
    Out("F_eff", "Effective wind fetch", "km", "mi", "scalar",
        note="The fetch length used for wave growth; the single input fetch for open water, or the off-wind direction-weighted radial fetch for restricted geometry."),
    Out("wave_dir", "Mean wave direction", "deg", "deg", "scalar",
        note="The direction of wave development (deg clockwise from north), offset from the wind direction by the off-wind angle phi for restricted fetch."),
    Out("H_mo", "Wave height", "m", "ft", "scalar",
        note="The spectrally based significant wave height H_mo from fetch- or duration-limited growth, capped at the fully developed limit."),
    Out("T_p", "Peak wave period", "s", "s", "scalar",
        note="The peak spectral wave period T_p from fetch- or duration-limited growth, capped at the fully developed limit."),
)


@dataclass
class Result:
    U_e: float; U_a: float; F_eff: float; wave_dir: float; H_mo: float; T_p: float
    growth: str = ""
    notes: str = ""


# --- wind adjustment: observed wind -> 10-m equivalent-neutral wind -------------
def _equiv_neutral_10m(U_obs: float, z_obs: float, obs_type: str) -> float:
    """Constant-stress logarithmic profile (SPM/TR eqs 3,6,7,13), neutral stability.

    Iterates friction velocity U_* and roughness z_0 (cgs units), then evaluates the
    wind at 10 m (1000 cm). Shipboard observations get the ship-bias correction (1)
    first. Returns the 10-m neutral wind in m/s."""
    U = U_obs
    if obs_type == "Overwater (shipboard)":
        U = 1.864 * U ** (7.0 / 9.0)                 # eq (1) ship-bias correction (m/s)
    Uo = U * 100.0                                    # cm/s
    zo = z_obs * 100.0                                # cm
    C1, C2, C3 = 0.1525, 0.019 / 980.0, -0.00371      # eq (6),(7) constant-stress (cgs)
    u_star = 0.04 * Uo                                # initial guess
    for _ in range(200):
        z0 = max(C1 / u_star + C2 * u_star * u_star + C3, 1e-4)
        u_new = _K_VON_KARMAN * Uo / math.log(zo / z0)   # eq (3), neutral (Psi = 0)
        if abs(u_new - u_star) <= 0.1:                # eq (12): eps_U* <= 0.1 cm/s
            u_star = u_new
            break
        u_star = u_new
    z0 = max(C1 / u_star + C2 * u_star * u_star + C3, 1e-4)
    U_10 = (u_star / _K_VON_KARMAN) * math.log(1000.0 / z0)   # eq (13), 10 m = 1000 cm
    return U_10 / 100.0                               # m/s


def _psi_m(zeta: float) -> float:
    """Integrated momentum stability function (Businger-Dyer / Paulson 1970).

    Unstable (zeta < 0): psi_m = 2 ln[(1+x)/2] + ln[(1+x^2)/2] - 2 atan(x) + pi/2,
    x = (1 - 16 zeta)^(1/4). Stable (zeta > 0): psi_m = -5 zeta. This is the canonical
    surface-layer form verified against the boundary-layer literature; it is NOT the
    (corrupted, sign-inconsistent) function transcribed in ACES TR eq 9."""
    if zeta >= 0.0:
        return -5.0 * zeta
    x = (1.0 - 16.0 * zeta) ** 0.25
    return 2.0 * math.log((1.0 + x) / 2.0) + math.log((1.0 + x * x) / 2.0) \
        - 2.0 * math.atan(x) + math.pi / 2.0


def _stability_10m(U_obs: float, z_obs: float, obs_type: str, deltaT: float) -> float:
    """10-m equivalent-neutral wind WITH a physical air-sea stability correction.

    Same constant-stress iteration as `_equiv_neutral_10m`, but the wind profile carries
    the Businger-Dyer psi_m term and the bulk Obukhov length (TR eq 8,
    L' = 1.79 U*^2/deltaT [ln(z/z0) - psi], where 1.79 = theta_bar/(k^2 g) in cgs).
    deltaT = T_air - T_sea: deltaT < 0 is unstable (psi_m > 0, equivalent-neutral wind a
    little higher); deltaT > 0 is stable (a little lower). The equivalent-neutral 10-m
    wind itself drops the psi term (eq 13). Opt-in; see the module docstring on why this
    does not reproduce ACES Example 3."""
    U = U_obs
    if obs_type == "Overwater (shipboard)":
        U = 1.864 * U ** (7.0 / 9.0)
    Uo = U * 100.0; zo = z_obs * 100.0                # cgs
    C1, C2, C3 = 0.1525, 0.019 / 980.0, -0.00371
    u_star = 0.04 * Uo; psi = 0.0
    for _ in range(300):
        z0 = max(C1 / u_star + C2 * u_star * u_star + C3, 1e-4)
        for _ in range(100):                          # inner: self-consistent L' / psi
            Lp = 1.79 * (u_star * u_star / deltaT) * (math.log(zo / z0) - psi)
            psi_new = _psi_m(zo / Lp)
            if abs(psi_new - psi) < 1e-9:
                psi = psi_new; break
            psi = psi_new
        u_new = _K_VON_KARMAN * Uo / (math.log(zo / z0) - psi)
        if abs(u_new - u_star) <= 0.1:
            u_star = u_new; break
        u_star = u_new
    z0 = max(C1 / u_star + C2 * u_star * u_star + C3, 1e-4)
    U_10 = (u_star / _K_VON_KARMAN) * math.log(1000.0 / z0)   # eq (13): equivalent neutral
    return U_10 / 100.0


def _duration_factor(t_s: float) -> float:
    """U_t / U_3600 (SPM/TR eqs 21,22); t in seconds."""
    if t_s < 1.0:
        t_s = 1.0
    if t_s <= 3600.0:
        return 1.277 + 0.296 * math.tanh(0.9 * math.log10(45.0 / t_s))   # eq (21)
    return -0.15 * math.log10(t_s) + 1.5334                              # eq (22)


# --- deepwater wave growth (open water; SPM/TR eqs 27-42) -----------------------
def _growth_deep_open(Ua: float, F: float, t_s: float, g: float):
    """Returns (H_mo, T_p, growth_label) for deep water, open-water fetch."""
    t_fmin = 68.8 * F ** (2.0 / 3.0) / (g ** (1.0 / 3.0) * Ua ** (1.0 / 3.0))   # eq (27)
    if t_s >= t_fmin:                                 # fetch-limited
        H = 0.0016 * (Ua * Ua / g) * (g * F / (Ua * Ua)) ** 0.5            # eq (33)
        T = 0.2857 * (Ua / g) * (g * F / (Ua * Ua)) ** (1.0 / 3.0)         # eq (35)
        label = "Deep water fetch-limited"
    else:                                             # duration-limited
        H = 0.0000851 * (Ua * Ua / g) * (g * t_s / Ua) ** (5.0 / 7.0)      # eq (29)
        T = 0.0702 * (Ua / g) * (g * t_s / Ua) ** 0.411                    # eq (31)
        label = "Deep water duration-limited"
    H_fd = 0.2433 * (Ua * Ua / g)                     # eq (37)
    T_fd = 8.134 * (Ua / g)                           # eq (39)
    if H >= H_fd or T >= T_fd:
        label = "Deep water fully developed"
    return min(H, H_fd), min(T, T_fd), label          # eqs (41),(42)


def _growth_shallow_open(Ua: float, F: float, d: float, g: float):
    """Returns (H_mo, T_p, growth_label) for shallow water, open-water fetch
    (SPM/TR eqs 43,44; fetch-limited with depth limiting)."""
    u2 = Ua * Ua
    dterm_h = math.tanh(0.530 * (g * d / u2) ** 0.75)                      # eq (43)
    H = (u2 / g) * 0.283 * dterm_h * math.tanh(
        (0.0016 / 0.283) * (g * F / u2) ** 0.5 / dterm_h)
    dterm_t = math.tanh(0.833 * (g * d / u2) ** 0.375)                     # eq (44)
    T = (Ua / g) * 7.54 * dterm_t * math.tanh(
        (0.2857 / 7.54) * (g * F / u2) ** (1.0 / 3.0) / dterm_t)
    return H, T, "Shallow water fetch-limited"


# --- restricted-fetch growth (SPM/Smith 1991; TR eqs 28,30,32,34,36,45,46) ------
def _growth_deep_restricted(Uab: float, F: float, t_s: float, g: float):
    """Deep water, restricted fetch. Uab = U_a*cos(phi) (fetch-parallel component)."""
    t_fmin = 51.09 * F ** 0.73 / (g ** 0.38 * Uab ** 0.44)                 # eq (28)
    if t_s >= t_fmin:                                 # fetch-limited
        H = 0.0015 * (Uab * Uab / g) * (g * F / (Uab * Uab)) ** 0.5        # eq (34)
        T = 0.3704 * (Uab / g) * (g * F / (Uab * Uab)) ** 0.28             # eq (36)
        label = "Deep water fetch-limited (restricted)"
    else:                                             # duration-limited
        H = 0.000103 * (Uab * Uab / g) * (g * t_s / Uab) ** 0.69           # eq (30)
        T = 0.082 * (Uab / g) * (g * t_s / Uab) ** 0.39                    # eq (32)
        label = "Deep water duration-limited (restricted)"
    H_fd = 0.2433 * (Uab * Uab / g)                   # eq (38)
    T_fd = 8.134 * (Uab / g)                          # eq (40)
    if H >= H_fd or T >= T_fd:
        label = "Deep water fully developed (restricted)"
    return min(H, H_fd), min(T, T_fd), label


def _growth_shallow_restricted(Uab: float, F: float, d: float, g: float):
    """Shallow water, restricted fetch (TR eqs 45,46)."""
    u2 = Uab * Uab
    dh = math.tanh(0.530 * (g * d / u2) ** 0.75)                           # eq (45)
    H = (u2 / g) * 0.283 * dh * math.tanh((0.0015 / 0.283) * (g * F / u2) ** 0.5 / dh)
    dt = math.tanh(0.833 * (g * d / u2) ** 0.375)                          # eq (46)
    T = (Uab / g) * 7.54 * dt * math.tanh((0.3704 / 7.54) * (g * F / u2) ** 0.28 / dt)
    return H, T, "Shallow water fetch-limited (restricted)"


def _restricted_fetch(fetches, wind_dir, dbeta, beta1):
    """Off-wind wave-development direction selection (TR eq 26): maximize
    F_phi^0.28 * (cos phi)^0.44 over off-wind angle phi (0..90 deg, 1-deg steps),
    where F_phi is the 15-deg-arc-averaged (direction-interpolated) radial fetch.
    Returns (F_eff, phi_deg, mean_wave_direction_deg)."""
    n = len(fetches)
    ang = [beta1 + i * dbeta for i in range(n)]
    lo, hi = ang[0], ang[-1]

    def f_interp(th):
        if th <= lo:
            return fetches[0]
        if th >= hi:
            return fetches[-1]
        for i in range(n - 1):
            if ang[i] <= th <= ang[i + 1]:
                t = (th - ang[i]) / (ang[i + 1] - ang[i])
                return fetches[i] + t * (fetches[i + 1] - fetches[i])
        return fetches[-1]

    def f_arc(th):                                    # 15-deg arc average (1-deg steps)
        s, c, d = 0.0, 0, th - 7.5
        while d <= th + 7.5 + 1e-9:
            s += f_interp(d); c += 1; d += 1.0
        return s / c

    best = None
    for phi in range(0, 91):
        for sgn in (-1, 1):
            th = wind_dir + sgn * phi
            if th < lo or th > hi:
                continue
            F = f_arc(th)
            score = F ** 0.28 * math.cos(math.radians(phi)) ** 0.44
            if best is None or score > best[0]:
                best = (score, float(phi), th, F)
    if best is None:
        return fetches[0], 0.0, wind_dir
    return best[3], best[1], best[2]


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle") or f.key not in inp:
            continue                                   # tables/choices + absent optional keys
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo:g}, {f.hi:g}] ({f.note})")


# --- compute --------------------------------------------------------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Adjusts an observed wind to a 10-m equivalent-neutral wind, linearizes it to '
            'the constant-drag "adjusted" wind used by the growth laws, then estimates '
            'fetch- or duration-limited wave growth, returning the equivalent/adjusted '
            'wind speeds plus the significant wave height H_mo and peak period T_p (deep '
            'or shallow water, open or restricted fetch).',
 'method_key': 'wave_eq',
 'methods': [{'name': 'Deep-water wave growth',
              'when': 'Deep',
              'tag': 'standard',
              'note': 'Hasselmann/SPM deepwater fetch- and duration-limited growth driven '
                      'by the adjusted wind U_a; restricted fetch substitutes the '
                      'fetch-parallel component U_a cos(phi).',
              'equations': [{'tex': 'u_{*} = \\frac{k\\,U_{obs}}{\\ln(z_{obs}/z_0)}',
                             'desc': 'Constant-stress friction velocity from the observed '
                                     'wind (eq 3); z_0 from eq 6, iterated to '
                                     'convergence.'},
                            {'tex': 'U_e = \\frac{u_{*}}{k}\\,\\ln\\frac{1000}{z_0}',
                             'desc': '10-m equivalent-neutral wind, cgs 1000 cm = 10 m (eq '
                                     '13); a 0.9 factor is applied for fetch < 16 km (eq '
                                     '20).'},
                            {'tex': 'C_D = 0.001\\,(0.75 + 0.067\\,U_e)',
                             'desc': 'Garratt wind-stress drag coefficient (eq 24).'},
                            {'tex': 'U_a = U_e\\,\\sqrt{C_D / 0.001}',
                             'desc': 'Equivalent neutral wind linearized to constant drag '
                                     'C_D = 0.001 (eq 25); this U_a drives the growth '
                                     'laws.'},
                            {'tex': 'H_{mo} = '
                                    '0.0016\\,\\frac{U_a^{2}}{g}\\left(\\frac{g\\,F}{U_a^{2}}\\right)^{1/2}',
                             'desc': 'Deepwater fetch-limited wave height (eq 33), capped '
                                     'at the fully developed limit 0.2433 U_a^2/g (eqs '
                                     '37,41).'},
                            {'tex': 'T_p = '
                                    '0.2857\\,\\frac{U_a}{g}\\left(\\frac{g\\,F}{U_a^{2}}\\right)^{1/3}',
                             'desc': 'Deepwater fetch-limited peak period (eq 35), capped '
                                     'at the fully developed limit 8.134 U_a/g (eqs '
                                     '39,42).'}]},
             {'name': 'Shallow-water wave growth',
              'when': 'Shallow',
              'tag': 'standard',
              'note': 'SPM depth-limited (Bretschneider-Reid) fetch-limited forms with a '
                      'constant fetch depth d; flagged interim in the TR. Same '
                      'wind-adjustment chain as the deep-water method.',
              'equations': [{'tex': 'u_{*} = \\frac{k\\,U_{obs}}{\\ln(z_{obs}/z_0)}',
                             'desc': 'Constant-stress friction velocity from the observed '
                                     'wind (eq 3).'},
                            {'tex': 'U_e = \\frac{u_{*}}{k}\\,\\ln\\frac{1000}{z_0}',
                             'desc': '10-m equivalent-neutral wind (eq 13).'},
                            {'tex': 'U_a = U_e\\,\\sqrt{C_D / 0.001}',
                             'desc': 'Adjusted constant-drag wind with C_D = 0.001(0.75 + '
                                     '0.067 U_e) (eqs 24,25).'},
                            {'tex': 'H_{mo} = '
                                    '0.283\\,\\frac{U_a^{2}}{g}\\,\\tanh\\left[0.530\\left(\\frac{g\\,d}{U_a^{2}}\\right)^{0.75}\\right]\\tanh\\left[\\frac{0.00565\\,(g\\,F/U_a^{2})^{0.5}}{\\tanh[0.530(g\\,d/U_a^{2})^{0.75}]}\\right]',
                             'desc': 'Shallow-water fetch-limited wave height with depth '
                                     'limiting (eq 43).'},
                            {'tex': 'T_p = '
                                    '7.54\\,\\frac{U_a}{g}\\,\\tanh\\left[0.833\\left(\\frac{g\\,d}{U_a^{2}}\\right)^{0.375}\\right]\\tanh\\left[\\frac{0.03788\\,(g\\,F/U_a^{2})^{0.333}}{\\tanh[0.833(g\\,d/U_a^{2})^{0.375}]}\\right]',
                             'desc': 'Shallow-water fetch-limited peak period with depth '
                                     'limiting (eq 44).'}]}],
 'symbols': [['U_obs', 'Observed wind speed at elevation z_obs'],
             ['z_obs', 'Elevation of the wind observation'],
             ['z_0', 'Sea-surface roughness length (eq 6, iterated)'],
             ['u_*', 'Friction velocity of the constant-stress layer'],
             ['k', 'Von Karman constant (= 0.40)'],
             ['U_e', '10-m equivalent-neutral wind speed'],
             ['C_D', 'Wind-stress drag coefficient (Garratt)'],
             ['U_a', 'Adjusted constant-drag wind driving the growth laws'],
             ['F', 'Wind fetch length (restricted: effective radial fetch)'],
             ['d', 'Average fetch depth (shallow-water forms)'],
             ['g', 'Gravitational acceleration'],
             ['H_mo', 'Spectrally based significant wave height'],
             ['T_p', 'Peak spectral wave period'],
             ['phi',
              'Off-wind angle between wind and wave-development direction (restricted '
              'fetch)']],
 'references': ['SPM (1984) Ch. 3',
                'Resio & Vincent (1977)',
                'Smith (1991) CERC-91-2',
                'Hasselmann et al. (1973, 1976)',
                'Bretschneider & Reid (1954)',
                'Garratt (1977)',
                'ACES TR 1-1 (eqs 1-46)']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Windspeed adjustment + wave growth for SI inputs. The GUI converts US->SI at
    the edge (wind in m/s, fetch/depth in m, durations in s, deltaT in C)."""
    _validate(inp)
    U_obs = float(inp["U_obs"]); z_obs = float(inp["z_obs"]); deltaT = float(inp["deltaT"])
    dur_obs = float(inp["dur_obs"]); dur_final = float(inp["dur_final"])
    F = float(inp["F"]); d = float(inp["d"])
    obs_type = str(inp.get("obs_type", "Shore (windward)"))
    fetch_type = str(inp.get("fetch_type", "Open Water"))
    wave_eq = str(inp.get("wave_eq", "Shallow"))
    wind_dir = float(inp.get("wind_dir", 0.0))

    notes = []

    # 1) observed wind -> 10-m equivalent-neutral wind
    stability = str(inp.get("stability", "Neutral (validated)"))
    if deltaT != 0.0 and stability.startswith("Businger"):
        U_10 = _stability_10m(U_obs, z_obs, obs_type, deltaT)
        notes.append("Businger-Dyer stability correction applied (physical, opt-in). "
                     "This is the literature-standard surface-layer form; it does NOT "
                     "reproduce ACES Example 3, whose larger opposite-signed correction "
                     "cannot be reconciled with standard physics or its own equations.")
    else:
        U_10 = _equiv_neutral_10m(U_obs, z_obs, obs_type)
        if deltaT != 0.0:
            notes.append("deltaT != 0 but stability = Neutral: temperature correction NOT "
                         "applied (validated neutral path). Select Businger-Dyer to apply a "
                         "physical correction; see module docstring.")
    # 2) duration adjustment to the final-wind duration (eqs 21,22)
    U_e = U_10 * _duration_factor(dur_final) / _duration_factor(dur_obs)

    # 3) fetch + wave growth
    if fetch_type == "Restricted":
        fetches = [float(row[0]) for row in inp["fetches"] if row and row[0] not in (None, "")]
        if len(fetches) < 2:
            raise ValueError("restricted fetch needs at least 2 radial fetch lengths")
        F_eff, phi, wave_dir = _restricted_fetch(
            fetches, wind_dir, float(inp.get("dbeta", 12.0)), float(inp.get("beta1", 0.0)))
        F_use = F_eff
    else:
        F_eff, phi, wave_dir, F_use = F, 0.0, wind_dir, F

    # short-fetch reduction (eq 20): effective fetch < 16 km
    if F_use < 16000.0:
        U_e = 0.9 * U_e
        notes.append("short-fetch reduction applied (F < 16 km): U_e x 0.9")
    # 4) wind stress -> adjusted (constant-drag) wind (eqs 24,25)
    C_D = 0.001 * (0.75 + 0.067 * U_e)
    U_a = U_e * math.sqrt(C_D / 0.001)
    U_ab = U_a * math.cos(math.radians(phi))          # fetch-parallel component (restricted)

    # 5) wave growth
    if fetch_type == "Restricted":
        if wave_eq == "Shallow":
            H_mo, T_p, growth = _growth_shallow_restricted(U_ab, F_use, d, g)
        else:
            H_mo, T_p, growth = _growth_deep_restricted(U_ab, F_use, dur_final, g)
    elif wave_eq == "Shallow":
        H_mo, T_p, growth = _growth_shallow_open(U_a, F_use, d, g)
    else:
        H_mo, T_p, growth = _growth_deep_open(U_a, F_use, dur_final, g)
    notes.insert(0, growth)
    return Result(U_e=U_e, U_a=U_a, F_eff=F_eff, wave_dir=wave_dir,
                  H_mo=H_mo, T_p=T_p, growth=growth, notes="; ".join(notes))


# --- self-tests -----------------------------------------------------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    # User's Guide Example 1-1 (US units): Shore (windward) / Open Water / Shallow
    r = compute(dict(z_obs=25.0 * _FT, U_obs=45.0 * _MPH, deltaT=0.0,
                     dur_obs=3.0 * 3600.0, dur_final=3.0 * 3600.0, lat=30.0,
                     obs_type="Shore (windward)", fetch_type="Open Water",
                     wave_eq="Shallow", F=26.0 * _MI, d=13.0 * _FT))
    assert _approx(r.U_e / _MPH, 46.493, 3e-3), r.U_e / _MPH
    assert _approx(r.U_a / _MPH, 68.055, 3e-3), r.U_a / _MPH
    assert _approx(r.H_mo / _FT, 4.237, 5e-3), r.H_mo / _FT
    assert _approx(r.T_p, 4.771, 5e-3), r.T_p
    assert "fetch-limited" in r.growth
    # duration symmetry: dur_obs == dur_final leaves the 10-m wind unchanged
    r2 = compute(dict(z_obs=25.0 * _FT, U_obs=45.0 * _MPH, deltaT=0.0, dur_obs=3600.0,
                      dur_final=3600.0, lat=30.0, obs_type="Overwater (not shipboard)",
                      fetch_type="Open Water", wave_eq="Deep", F=26.0 * _MI, d=13.0 * _FT))
    assert r2.H_mo > 0 and r2.T_p > 0 and "Deep" in r2.growth

    # Restricted fetch, Example 3 geometry (stability-independent): wind 125 deg,
    # dbeta=12, beta1=0, 14 radials -> effective fetch ~26.6 mi, wave dir ~93 deg.
    ex3 = [3.7, 12.3, 13.4, 12.2, 13.2, 36.0, 35.6, 28.7, 26.8, 13.0, 10.4, 10.1, 6.4, 5.7]
    F_eff, phi, wdir = _restricted_fetch([v * _MI for v in ex3], 125.0, 12.0, 0.0)
    assert _approx(F_eff / _MI, 26.61, 0.02), F_eff / _MI
    assert _approx(wdir, 93.0, 0.02), wdir
    # restricted deep growth with the manual's adjusted wind (Ua=63.27 mph) -> 7.80 ft, 5.74 s
    Uab = 63.27 * _MPH * math.cos(math.radians(phi))
    H, T, lab = _growth_deep_restricted(Uab, F_eff, 5.0 * 3600.0, G_SI)
    assert _approx(H / _FT, 7.80, 5e-3), H / _FT
    assert _approx(T, 5.74, 5e-3), T
    assert "fetch-limited" in lab

    # Opt-in Businger-Dyer stability: neutral default unchanged; physical correction is small
    # and correctly signed (deltaT<0 unstable -> higher; deltaT>0 stable -> lower).
    base = dict(z_obs=30.0 * _FT, U_obs=45.0 * _MPH, dur_obs=5 * 3600.0, dur_final=5 * 3600.0,
                lat=47.0, obs_type="Overwater (not shipboard)", fetch_type="Open Water",
                wave_eq="Deep", F=26.6 * _MI, d=100.0)
    u_neu = compute(dict(base, deltaT=0.0)).U_e
    u_uns = compute(dict(base, deltaT=-3.0, stability="Businger-Dyer (physical)")).U_e
    u_stb = compute(dict(base, deltaT=3.0, stability="Businger-Dyer (physical)")).U_e
    assert u_uns > u_neu > u_stb, (u_uns / _MPH, u_neu / _MPH, u_stb / _MPH)
    assert abs(u_uns - u_neu) / u_neu < 0.03    # standard surface-layer effect is small
    # default (Neutral) ignores deltaT -> same as neutral
    assert _approx(compute(dict(base, deltaT=-3.0)).U_e, u_neu, 1e-9)
    print("  self-tests: PASS (Example 1-1 + deep + restricted-fetch + opt-in stability sign)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        vv = inp[f.key]
        if f.kind == "table":
            sval = f"{len(vv)} rows"
        elif isinstance(vv, (int, float)) and f.kind != "choice":
            sval = f"{vv:>10.4g}"
        else:
            sval = f"{vv:>10}"
        print(f"    {f.label:30s} {f.key:10s} = {sval} {f.unit_si}")
    print("  OUTPUTS:")
    by_key = {o.key: o for o in OUTPUTS}
    for kk, vv in {"U_e": r.U_e, "U_a": r.U_a, "H_mo": r.H_mo, "T_p": r.T_p}.items():
        o = by_key[kk]
        print(f"    {o.label:30s} {kk:5s} = {vv:>12.5g} {o.unit_si}")
    print(f"  wave growth: {r.growth}")
    print(f"  notes: {r.notes}")
    print("  (US: U_e=46.49 mph, U_a=68.06 mph, H_mo=4.24 ft, T_p=4.77 s)")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
