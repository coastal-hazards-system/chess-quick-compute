"""CHESS-QC application 2-2 — Cnoidal Wave Theory.

Originating ACES application: 2-2 "Cnoidal Wave Theory" (functional area: Wave Theory).
Finite-amplitude periodic long-wave theory in terms of Jacobian elliptic functions;
first-order (Isobe 1985) and second-order (Hardy & Kraus 1987) approximations.

Classification: exact (cnoidal theory of Isobe 1985 / Hardy & Kraus 1987 -- all series
coefficients known, nothing guessed; reproduces the User's Guide Example 2-2 (first order);
the order truncation is the documented method, not a gap).
Theory and references: Korteweg & de Vries (1895); Isobe (1985) 1st order; Hardy &
  Kraus (1987) 2nd order; elliptic identities Abramowitz & Stegun (1972).
  Equations transcribed in docs/EQUATIONS.md, TR chapter 2-2 (eqs 1-61):
    dispersion 1st (36) / 2nd (52), celerity (37/53), L=cT (38), eta (39/54),
    E (40/55), F (41/56), pressure (42/57), u (43/58), w (44/59), du/dt (45/60),
    dw/dt (46/61). epsilon=H/d ; theta = 2K[(x/L) - (t/T)] ; modulus kappa.

  NOTE on accelerations: the TR's du/dt (45,60) and dw/dt (46,61) are reproduced
  here by analytic theta-differentiation of the verified velocity fields (chain rule,
  d(theta)/dt = -2K/T). This regenerates eqs (45) and (60) *exactly*; for dw/dt it
  corrects a transcription artifact in (46)/(61) whose trig factor lost its squared
  products -- d(csd)/dtheta = dn^2(cn^2 - sn^2) - kappa^2 sn^2 cn^2 (csd = cn*sn*dn).

Self-containment: zero sibling imports; embeds its own Field/AppMeta/Result
dataclasses, the complete elliptic integrals K,E (AGM), and the Jacobian elliptic
functions sn,cn,dn (descending Landen / Numerical Recipes sncndn). Runnable
standalone:
    python chessqc_2_2_cnoidal_wave_theory.py
which runs analytic self-tests (elliptic-function identities, mean-zero surface),
then prints an ACES-style tabulation of the default example.  stdlib + numpy only.

I/O mirrors the ACES manual's *exact* Cnoidal Wave Theory lists:
  Inputs : H (>0), T (>0), d (>0), z (no restriction), X/L (0-1), Order (1|2)
  Outputs: eta, C, L, P, E, U_r, p, (u,w), (du/dt,dw/dt), modulus kappa
           + profile arrays (eta, u, w vs X over +/- one wavelength) for plots.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple

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
    aces_id="2-2",
    name="Cnoidal Wave Theory",
    area="Wave Theory",
    classification="exact",
    cite="Isobe (1985); Hardy & Kraus (1987); TR 2-2",
    default_system="US",     # opens on the User's Guide Example (US units)
)

# Complete input list (ACES manual "Input requirements" for Cnoidal Wave Theory).
# Defaults are SI-internal; values shown are the ACES User's Guide Example 2-2
# (H=10 ft, T=15 s, d=25 ft, z=-12.5 ft, X/L=0.5, Order=1 -> L=455.74 ft, C=30.38
# ft/s, eta=-2.14 ft ...). See tests/test_manual_oracle.py.
_FT = 0.3048
INPUTS = (
    Field("H",  "Wave height",        "float", "m", "ft",  default=10.0 * _FT, lo=1e-6, hi=1e4,
          note="> 0"),
    Field("T",  "Wave period",        "float", "s", "s",   default=15.0, lo=1e-3, hi=1e4,
          note="> 0"),
    Field("d",  "Water depth",        "float", "m", "ft",  default=25.0 * _FT, lo=1e-6, hi=1e5,
          note="> 0"),
    Field("z",  "Vertical coordinate","float", "m", "ft",  default=-12.5 * _FT,
          note="from SWL (z=0), +up; no restriction (clamped to [-d, eta])"),
    Field("xL", "Wavelength fraction (X/L)", "float", "", "", default=0.50, lo=0.0, hi=1.0,
          note="0.0 to 1.0 (phase position; 0 = crest)"),
    Field("order", "Order", "choice", "", "", default=1, choices=(1, 2),
          note="1 = Isobe (1985); 2 = Hardy & Kraus (1987)"),
    Field("water", "Water type", "choice", "", "", default="Salt", choices=("Salt", "Fresh"),
          note="sets water density used for energy/pressure"),
)

# Complete output list (ACES manual "Output" for Cnoidal Wave Theory).
OUTPUTS = (
    Out("L",     "Wave length",                  "m",     "ft",    "scalar",
        note="Cnoidal wavelength L = cT, the horizontal crest-to-crest distance of one wave."),
    Out("C",     "Celerity",                     "m/s",   "ft/s",  "scalar",
        note="Phase speed (celerity) c at which the wave form propagates."),
    Out("E",     "Energy density",               "N/m",   "lb/ft", "scalar",
        note="Mean total wave energy per unit surface area (sum of kinetic and potential energy)."),
    Out("P",     "Energy flux (power)",          "N/s",   "lb/s",  "scalar",
        note="Mean energy flux (wave power) transmitted per unit crest width."),
    Out("Ur",    "Ursell parameter (HL^2/d^3)",  "",      "",      "scalar",
        note="Ursell number U_r = HL^2/d^3 gauging nonlinearity; cnoidal validity is questionable below ~26."),
    Out("kappa", "Elliptic modulus",             "",      "",      "scalar",
        note="Jacobian elliptic modulus kappa in [0,1]; kappa -> 1 is the solitary-wave limit."),
    Out("eta",   "Surface elevation",            "m",     "ft",    "point",
        note="Water-surface elevation above still-water level (z=0) at the chosen phase X/L; + is above SWL."),
    Out("p",     "Pressure",                     "Pa",    "psf",   "point",
        note="Total fluid pressure at the point (z, X/L), including hydrostatic and dynamic contributions."),
    Out("u",     "Horizontal velocity",          "m/s",   "ft/s",  "point",
        note="Horizontal water-particle velocity at (z, X/L); + is in the direction of wave travel."),
    Out("w",     "Vertical velocity",            "m/s",   "ft/s",  "point",
        note="Vertical water-particle velocity at (z, X/L); + is upward."),
    Out("dudt",  "Horizontal acceleration",      "m/s^2", "ft/s^2","point",
        note="Local horizontal water-particle acceleration du/dt at (z, X/L)."),
    Out("dwdt",  "Vertical acceleration",        "m/s^2", "ft/s^2","point",
        note="Local vertical water-particle acceleration dw/dt at (z, X/L); + is upward."),
    Out("profile_X",   "Profile: X (+/- one wavelength)", "m",   "ft",   "profile",
        note="Horizontal coordinate X spanning +/- one wavelength, the abscissa for the profile plots."),
    Out("profile_eta", "Profile: surface elevation",      "m",   "ft",   "profile",
        note="Surface elevation above SWL versus X over +/- one wavelength (the cnoidal wave form)."),
    Out("profile_u",   "Profile: horizontal velocity",    "m/s", "ft/s", "profile",
        note="Horizontal water-particle velocity versus X at the chosen depth z, over +/- one wavelength."),
    Out("profile_w",   "Profile: vertical velocity",      "m/s", "ft/s", "profile",
        note="Vertical water-particle velocity versus X at the chosen depth z, over +/- one wavelength."),
)


@dataclass
class Result:
    # scalars
    L: float;  C: float;  E: float;  P: float;  Ur: float;  kappa: float
    # point kinematics at (z, X/L)
    eta: float; p: float; u: float;  w: float;  dudt: float; dwdt: float
    # profile arrays (for plots), over X in [-L, L]
    profile_X: np.ndarray
    profile_eta: np.ndarray
    profile_u: np.ndarray
    profile_w: np.ndarray
    notes: str = ""


# --- complete elliptic integrals K(m), E(m) via AGM (A&S 17.6; m = kappa^2) ------
def _cei(m: float) -> Tuple[float, float]:
    """Complete elliptic integrals (K, E), 1st and 2nd kind, parameter m = kappa^2.

    Arithmetic-geometric-mean iteration (quadratic convergence):
      a0=1, b0=sqrt(1-m), c0=sqrt(m); K = pi/(2 a_inf);
      E = K (1 - sum_{n=0}^inf 2^{n-1} c_n^2)   with c_{n+1}=(a_n-b_n)/2.
    """
    if m <= 0.0:
        return math.pi / 2.0, math.pi / 2.0
    if m >= 1.0:
        m = 1.0 - 1e-16                       # K diverges; clamp for a large finite value
    a, b = 1.0, math.sqrt(1.0 - m)
    s = 0.5 * m                               # n=0 term: 2^{-1} c0^2 = m/2
    p = 1.0                                    # 2^{n-1} for n=1 is 2^0 = 1
    for _ in range(60):
        c = 0.5 * (a - b)
        s += p * c * c
        p *= 2.0
        a, b = 0.5 * (a + b), math.sqrt(a * b)
        if c < 1e-16 * a:
            break
    K = math.pi / (2.0 * a)
    return K, K * (1.0 - s)


# --- Jacobian elliptic functions sn,cn,dn (descending Landen; NR "sncndn") -------
def _sncndn(uu: float, m: float) -> Tuple[float, float, float]:
    """Jacobian elliptic sn, cn, dn at argument uu, parameter m = kappa^2 in [0,1].

    Descending Landen transformation (Abramowitz & Stegun 16.4 / Numerical Recipes).
    For m=1 reduces to the solitary-wave limit (sn->tanh, cn=dn->sech)."""
    emc = 1.0 - m                              # complementary parameter kappa'^2
    u = uu
    if emc == 0.0:                             # m = 1: solitary-wave limit
        cn = 1.0 / math.cosh(u)
        return math.tanh(u), cn, cn
    em = [0.0] * 14
    en = [0.0] * 14
    a = 1.0
    dn = 1.0
    l = 0
    c = a
    for i in range(1, 14):
        l = i
        em[i] = a
        emc = math.sqrt(emc)
        en[i] = emc
        c = 0.5 * (a + emc)
        if abs(a - emc) <= 1.0e-12 * a:
            break
        emc = emc * a
        a = c
    u = c * u
    sn = math.sin(u)
    cn = math.cos(u)
    if sn != 0.0:
        a = cn / sn
        c = c * a
        for ii in range(l, 0, -1):
            b = em[ii]
            a = a * c
            c = c * dn
            dn = (en[ii] + a) / (b + a)
            a = c / b
        a = 1.0 / math.sqrt(c * c + 1.0)
        sn = a if sn >= 0.0 else -a
        cn = c * sn
    return sn, cn, dn


def _ellipj(u, m: float):
    """Vectorized sn,cn,dn over a scalar or array argument u (parameter m=kappa^2)."""
    arr = np.atleast_1d(np.asarray(u, dtype=float))
    sn = np.empty_like(arr); cn = np.empty_like(arr); dn = np.empty_like(arr)
    for i, ui in enumerate(arr.ravel()):
        s, c, d = _sncndn(float(ui), m)
        sn.ravel()[i], cn.ravel()[i], dn.ravel()[i] = s, c, d
    if np.isscalar(u) or np.ndim(u) == 0:
        return float(sn[0]), float(cn[0]), float(dn[0])
    return sn, cn, dn


# --- dispersion solve: find the elliptic modulus kappa --------------------------
def _solve_modulus(H: float, T: float, d: float, order: int, g: float) -> float:
    """Solve the cnoidal dispersion relation for the squared modulus m = kappa^2.

    1st order (eq 36):  16 kappa^2 K^2 / 3 = g H T^2 / d^2.
    2nd order (eq 52):  RHS scaled by [1 - eps (1 + 2 lambda)/4]  (lambda depends on
    kappa -> solved by fixed-point iteration seeded with the 1st-order root).
    m*K(m)^2 is strictly increasing on (0,1): 0 -> +inf, so the root is unique.
    """
    base = g * H * T * T / (d * d)             # = 16 m K^2 / 3   (1st order)
    eps = H / d

    def root_for(rhs: float) -> float:
        target = 3.0 * rhs / 16.0              # solve m * K(m)^2 = target
        lo, hi = 1e-14, 1.0 - 1e-15
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            K, _ = _cei(mid)
            if mid * K * K < target:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-15:
                break
        return 0.5 * (lo + hi)

    m = root_for(base)
    if order >= 2:                             # iterate the lambda-dependent RHS
        for _ in range(50):
            kp2 = 1.0 - m                       # kappa'^2
            lam = kp2 / m
            rhs = base * (1.0 - eps * (1.0 + 2.0 * lam) / 4.0)
            m_new = root_for(rhs)
            if abs(m_new - m) < 1e-14:
                m = m_new
                break
            m = m_new
    return m


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes finite-amplitude cnoidal (periodic long-wave) properties from wave '
            'height, period, and depth: wavelength, celerity, surface elevation, '
            'water-particle velocities and accelerations, pressure, energy density and '
            'flux, the Ursell number, and the elliptic modulus. The solution is expressed '
            'with Jacobian elliptic functions, using either first-order (Isobe 1985) or '
            'second-order (Hardy & Kraus 1987) cnoidal theory.',
 'method_key': 'order',
 'methods': [{'name': 'First-order cnoidal (Isobe 1985)',
              'when': '1',
              'tag': 'standard',
              'note': "Default; reproduces ACES User's Guide Example 2-2 and is adequate "
                      'for typical shallow-water design.',
              'equations': [{'tex': '\\frac{16\\,\\kappa^2 K^2}{3} = \\frac{g H T^2}{d^2}',
                             'desc': 'First-order dispersion relation; solved for the '
                                     'squared elliptic modulus m = kappa^2 (m K^2 is '
                                     'monotone increasing, so the root is unique).'},
                            {'tex': 'c = \\sqrt{g d}\\,\\left(1 + \\epsilon\\,\\frac{1 + '
                                    '2\\lambda - 3\\mu}{2}\\right)',
                             'desc': 'Celerity; wavelength follows from L = c T. Here '
                                     'epsilon = H/d.'},
                            {'tex': '\\eta = d\\,\\left(A_0 + '
                                    'A_1\\,\\mathrm{cn}^{2}\\theta\\right)',
                             'desc': 'Surface elevation, with A_0 = epsilon(lambda - mu) '
                                     'and A_1 = epsilon.'},
                            {'tex': '\\theta = 2K\\left(\\frac{x}{L} - '
                                    '\\frac{t}{T}\\right)',
                             'desc': 'Phase argument of the Jacobian elliptic functions '
                                     '(crest at theta = 0).'},
                            {'tex': 'E = \\rho g H^2 \\,\\frac{-\\lambda + 2\\mu + '
                                    '4\\lambda\\mu - \\lambda^2 - 3\\mu^2}{3}',
                             'desc': 'Mean energy density per unit surface area; the '
                                     'energy flux uses the same factor times sqrt(g d).'}]},
             {'name': 'Second-order cnoidal (Hardy & Kraus 1987)',
              'when': '2',
              'tag': 'preferred',
              'note': 'Higher-order truncation; adds the cn^4 surface term and '
                      'depth-dependent velocity corrections for steeper or more strongly '
                      'nonlinear waves.',
              'equations': [{'tex': '\\frac{16\\,\\kappa^2 K^2}{3} = \\frac{g H '
                                    'T^2}{d^2}\\left(1 - \\epsilon\\,\\frac{1 + '
                                    '2\\lambda}{4}\\right)',
                             'desc': 'Second-order dispersion relation; the modulus is '
                                     'found by fixed-point iteration seeded with the '
                                     'first-order root.'},
                            {'tex': 'c = \\sqrt{g d}\\,\\left(1 + \\epsilon C_1 + '
                                    '\\epsilon^2 C_2\\right)',
                             'desc': 'Celerity with C_1 = (1 + 2lambda - 3mu)/2 and the '
                                     'second-order coefficient C_2.'},
                            {'tex': '\\eta = d\\,\\left(A_0 + '
                                    'A_1\\,\\mathrm{cn}^{2}\\theta + '
                                    'A_2\\,\\mathrm{cn}^{4}\\theta\\right)',
                             'desc': 'Surface elevation, now including the second-order '
                                     'cn^4 term (A_2 = (3/4) epsilon^2).'},
                            {'tex': 'u = \\sqrt{g d}\\,\\left[\\left(B_{00} + '
                                    'B_{10}\\,\\mathrm{cn}^{2}\\theta + '
                                    'B_{20}\\,\\mathrm{cn}^{4}\\theta\\right) - '
                                    '\\frac{1}{2}\\left(\\frac{z+d}{d}\\right)^{2}\\left(B_{01} '
                                    '+ B_{11}\\,\\mathrm{cn}^{2}\\theta + '
                                    'B_{21}\\,\\mathrm{cn}^{4}\\theta\\right)\\right]',
                             'desc': 'Horizontal velocity with elevation-dependent (z+d)/d '
                                     'correction terms.'},
                            {'tex': 'L = c T',
                             'desc': 'Wavelength from celerity and period; the Ursell '
                                     'number U_r = H L^2 / d^3 gauges cnoidal validity '
                                     '(questionable below ~26).'}]}],
 'symbols': [['H', 'Wave height'],
             ['T', 'Wave period'],
             ['d', 'Still-water depth'],
             ['eta', 'Surface elevation above SWL (z = 0)'],
             ['epsilon', 'Perturbation parameter, epsilon = H/d'],
             ['kappa',
              'Elliptic modulus (auxiliary parameter; kappa -> 1 is the solitary limit)'],
             ['K', 'Complete elliptic integral of the first kind, K(kappa)'],
             ['lambda', "lambda = kappa'^2 / kappa^2 = (1 - kappa^2)/kappa^2"],
             ['mu',
              'mu = E(kappa) / (kappa^2 K), where E is the complete elliptic integral of '
              'the second kind'],
             ['theta', 'Phase argument of cn, sn, dn: theta = 2K(x/L - t/T)']],
 'references': ['Isobe (1985)',
                'Hardy & Kraus (1987)',
                'Korteweg & de Vries (1895)',
                'Abramowitz & Stegun (1972)',
                'ACES TR 2-2 (eqs 1-61)']}


def compute(inp: dict, *, g: float = G_SI, rho: float | None = None,
            n_profile: int = 201) -> Result:
    """Cnoidal-wave-theory results for SI inputs {H, T, d, z, xL, order, water}.

    `inp` values are in SI (the GUI converts US->SI at the edge). `order` is 1
    (Isobe 1985) or 2 (Hardy & Kraus 1987). Water density follows the `water` field
    (Salt|Fresh) unless `rho` is given. Returns a fully-populated Result (all manual
    outputs + profile arrays over +/- one wavelength)."""
    H, T, d, z, xL = (float(inp[k]) for k in ("H", "T", "d", "z", "xL"))
    order = int(inp.get("order", 1))
    if rho is None:
        rho = RHO_FRESH if str(inp.get("water", "Salt")) == "Fresh" else RHO_SALT
    if H <= 0 or T <= 0 or d <= 0:
        raise ValueError("H, T, d must all be > 0")
    if not (0.0 <= xL <= 1.0):
        raise ValueError(f"X/L = {xL} outside [0, 1]")
    if order not in (1, 2):
        raise ValueError(f"order = {order} must be 1 or 2")

    eps = H / d                                 # perturbation parameter (eq 22)
    sgd = math.sqrt(g * d)                       # sqrt(g d)

    # modulus + elliptic quantities
    m = _solve_modulus(H, T, d, order, g)        # m = kappa^2
    kappa = math.sqrt(m)
    K, E = _cei(m)
    kp2 = 1.0 - m                                # kappa'^2 (eq 47)
    lam = kp2 / m                                # lambda (eq 48)
    mu = E / (m * K)                             # mu     (eq 49)

    # series coefficients by order (eqs 37-46 first; 53-61 second) -------------
    e2 = eps * eps
    # surface elevation eta/d = A0 + A1 cn^2 + A2 cn^4   (eqs 39 / 54)
    A0 = eps * (lam - mu)
    A1 = eps
    A2 = 0.0
    # horizontal velocity u/sqrt(gd) = (B00 + B10 cn^2 + B20 cn^4)
    #                                  - 1/2 (s)^2 (B01 + B11 cn^2 + B21 cn^4)  (43/58)
    B00 = eps * (lam - mu)
    B10 = eps
    B20 = 0.0
    B01 = B11 = B21 = 0.0
    # celerity c = sqrt(gd)(C0 + eps C1 + eps^2 C2)   (eqs 37 / 53)
    C1 = (1.0 + 2.0 * lam - 3.0 * mu) / 2.0
    C2 = 0.0
    # Bernoulli/pressure constant p_b = rho g d (P0 + eps P1 + eps^2 P2)   (42/57)
    P0, P1, P2 = 1.5, C1, 0.0
    # energy density / flux factors (eqs 40-41 / 55-56)
    E0 = (-lam + 2.0 * mu + 4.0 * lam * mu - lam * lam - 3.0 * mu * mu) / 3.0
    E1 = 0.0
    F0 = E0
    F1 = 0.0

    if order >= 2:
        A0 += e2 * ((-2.0 * lam + mu - 2.0 * lam * lam + 2.0 * lam * mu) / 4.0)  # 54.1
        A1 += -0.75 * e2                                                          # 54.2
        A2 = 0.75 * e2                                                            # 54.3
        B00 += e2 * ((lam - mu - 2.0 * lam * lam + 2.0 * mu * mu) / 4.0)          # 58.1
        B10 += e2 * ((1.0 - 6.0 * lam + 2.0 * mu) / 4.0)                          # 58.2
        B20 = -e2                                                                  # 58.3
        B01 = (3.0 * lam / 2.0) * e2                                              # 58.4
        B11 = 3.0 * e2 * (1.0 - lam)                                             # 58.5
        B21 = -(9.0 / 2.0) * e2                                                  # 58.6
        C2 = (-6.0 - 16.0 * lam + 5.0 * mu - 16.0 * lam * lam
              + 10.0 * lam * mu + 15.0 * mu * mu) / 40.0                          # 53.3
        P2 = (-1.0 - 16.0 * lam + 15.0 * mu - 16.0 * lam * lam
              + 30.0 * lam * mu) / 40.0                                           # 57.4
        E1 = (lam - 2.0 * mu - 17.0 * lam * mu + 3.0 * lam * lam
              - 17.0 * lam * lam * mu + 2.0 * lam ** 3 + 15.0 * mu ** 3) / 30.0   # 55.2
        F1 = (-4.0 * lam + 8.0 * mu + 53.0 * lam * mu - 12.0 * lam * lam
              - 60.0 * mu * mu + 53.0 * lam * lam * mu - 120.0 * lam * mu * mu
              - 8.0 * lam ** 3 + 75.0 * mu ** 3) / 30.0                           # 56.2

    # celerity, wavelength, bulk integrals
    C = sgd * (1.0 + eps * C1 + e2 * C2)         # eq 37 / 53
    L = C * T                                    # eq 38
    Ur = H * L * L / (d ** 3)                     # Ursell number HL^2/d^3
    gamma = rho * g                               # specific weight
    Edens = gamma * H * H * (E0 + eps * E1)       # eq 40 / 55
    Flux = gamma * H * H * sgd * (F0 + eps * F1)  # eq 41 / 56
    p_b = gamma * d * (P0 + eps * P1 + e2 * P2)   # eq 42.1 / 57.1

    # ---- field evaluator: kinematics at phase theta and height s=(z+d)/d -------
    dtheta_dt = -2.0 * K / T                      # d(theta)/dt, theta = 2K(x/L - t/T)

    def fields(theta, s):
        sn, cn, dn = _ellipj(theta, m)
        cn2 = cn * cn
        cn4 = cn2 * cn2
        csd = cn * sn * dn                                  # eq 51
        dcsd = dn * dn * (cn2 - sn * sn) - m * sn * sn * cn2  # d(csd)/dtheta
        s2 = s * s
        eta = d * (A0 + A1 * cn2 + A2 * cn4)                # eq 39 / 54
        # horizontal velocity (eq 43 / 58)
        u = sgd * ((B00 + B10 * cn2 + B20 * cn4)
                   - 0.5 * s2 * (B01 + B11 * cn2 + B21 * cn4))
        # vertical velocity (eq 44 / 59): w = sqrt(gd)(4 K d / L) csd * W(s,cn^2)
        Wfac = s * (B10 + 2.0 * B20 * cn2) - (1.0 / 6.0) * s ** 3 * (B11 + 2.0 * B21 * cn2)
        w = sgd * (4.0 * K * d / L) * csd * Wfac
        # accelerations by chain rule d/dt = dtheta_dt * d/dtheta (regenerates 45/60,
        # and the squared-product corrected 46/61)
        dpoly_u = (B10 * (-2.0 * csd) + B20 * (-4.0 * cn2 * csd)
                   - 0.5 * s2 * (B11 * (-2.0 * csd) + B21 * (-4.0 * cn2 * csd)))
        dudt = sgd * dpoly_u * dtheta_dt
        dWfac = -4.0 * csd * (s * B20 - (1.0 / 6.0) * s ** 3 * B21)
        dwdt = sgd * (4.0 * K * d / L) * (dcsd * Wfac + csd * dWfac) * dtheta_dt
        return eta, u, w, dudt, dwdt

    # ---- point of interest (z, X/L) at t = 0 -----------------------------------
    eta_crest = d * (A0 + A1 + A2)                # cn^2 = 1
    z_use = min(max(z, -d), eta_crest)            # clamp into the fluid column
    s_pt = (z_use + d) / d
    theta_pt = 2.0 * K * xL
    eta_p, u_p, w_p, dudt_p, dwdt_p = fields(theta_pt, s_pt)
    # pressure (eq 42 / 57): p = p_b - (rho/2)[(u-c)^2 + w^2] - rho g (z+d)
    p_p = p_b - 0.5 * rho * ((u_p - C) ** 2 + w_p ** 2) - gamma * (z_use + d)

    # ---- profile arrays over +/- one wavelength --------------------------------
    X = np.linspace(-L, L, n_profile)
    theta = 2.0 * K * (X / L)
    prof_eta, prof_u, prof_w, _, _ = fields(theta, s_pt)

    notes = [f"order {order}: " + ("Isobe (1985)" if order == 1 else "Hardy & Kraus (1987)")]
    if Ur < 26.0:
        notes.append(f"Ursell {Ur:.1f} < 26: cnoidal theory questionable (use linear)")
    if kappa > 0.99999:
        notes.append("near-solitary limit (kappa -> 1)")
    if z + d < 0:
        notes.append("WARNING: point outside waveform (z below the bed, z + d < 0)")
    elif z != z_use:
        notes.append(f"z clamped to {z_use:.3f} m for kinematics")

    return Result(L=L, C=C, E=Edens, P=Flux, Ur=Ur, kappa=kappa,
                  eta=eta_p, p=p_p, u=u_p, w=w_p, dudt=dudt_p, dwdt=dwdt_p,
                  profile_X=X, profile_eta=prof_eta, profile_u=prof_u, profile_w=prof_w,
                  notes="; ".join(notes))


# --- self-tests (elliptic identities + internal consistency) --------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    # 1) complete elliptic integrals vs known values
    K0, E0 = _cei(0.0)
    assert _approx(K0, math.pi / 2, 1e-9) and _approx(E0, math.pi / 2, 1e-9)
    K5, E5 = _cei(0.5)
    assert _approx(K5, 1.8540746773, 1e-6), K5
    assert _approx(E5, 1.3506438810, 1e-6), E5

    # 2) Jacobian elliptic-function identities at several arguments / moduli
    for m in (0.1, 0.5, 0.9, 0.999):
        K, _ = _cei(m)
        for u in (0.0, 0.3 * K, K, 1.7 * K):
            sn, cn, dn = _sncndn(u, m)
            assert _approx(sn * sn + cn * cn, 1.0, 1e-7), (m, u)
            assert _approx(dn * dn + m * sn * sn, 1.0, 1e-7), (m, u)
        # boundary values: sn(K)=1, cn(K)=0, dn(K)=kappa'
        sn, cn, dn = _sncndn(K, m)
        assert _approx(sn, 1.0, 1e-6) and abs(cn) < 1e-6
        assert _approx(dn, math.sqrt(1.0 - m), 1e-6), (m, dn)

    # 3) example reproduces the manual (first order); mean surface ~ 0 over one L
    r = compute({"H": 10.0 * _FT, "T": 15.0, "d": 25.0 * _FT,
                 "z": -12.5 * _FT, "xL": 0.5, "order": 1})
    assert _approx(r.L / _FT, 455.74, 5e-3), r.L / _FT
    assert _approx(r.C / _FT, 30.38, 5e-3), r.C / _FT
    assert _approx(r.eta / _FT, -2.14, 5e-2), r.eta / _FT
    assert _approx(r.u / _FT, -2.43, 5e-2), r.u / _FT
    assert abs(r.w) < 1e-6, r.w                      # w = 0 at the trough
    # mean of the surface profile over +/- one wavelength ~ 0
    assert abs(np.mean(r.profile_eta)) < 0.02 * r.profile_eta.max(), np.mean(r.profile_eta)
    # crest horizontal velocity matches sqrt(gd)(B00+B10) (manual table 8.917 ft/s)
    assert _approx(r.profile_u.max() / _FT, 8.917, 1e-2), r.profile_u.max() / _FT

    print("  self-tests: PASS (elliptic K/E + sn/cn/dn identities, manual example, mean eta~0)")


def _print_default_example() -> None:
    """ACES-style tabulation of the default example (SI).
    Defaults are the ACES User's Guide Example 2-2 (US units; stored SI-internal)."""
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        v = inp[f.key]
        sval = f"{v:>10.4g}" if isinstance(v, (int, float)) and f.kind != "choice" else f"{v:>10}"
        print(f"    {f.label:28s} {f.key:5s} = {sval} {f.unit_si}")
    print("  OUTPUTS:")
    scal = {"L": r.L, "C": r.C, "E": r.E, "P": r.P, "Ur": r.Ur, "kappa": r.kappa}
    pts = {"eta": r.eta, "p": r.p, "u": r.u, "w": r.w, "dudt": r.dudt, "dwdt": r.dwdt}
    by_key = {o.key: o for o in OUTPUTS}
    for kk, vv in {**scal, **pts}.items():
        o = by_key[kk]
        print(f"    {o.label:30s} {kk:5s} = {vv:>12.5g} {o.unit_si}")
    print(f"  notes: {r.notes}")
    print(f"  profile arrays: X/eta/u/w each length {len(r.profile_X)} (over +/- one L)")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
