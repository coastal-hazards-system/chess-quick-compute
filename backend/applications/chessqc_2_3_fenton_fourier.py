"""CHESS-QC application 2-3 — Fourier Series Wave Theory (Fenton).

Originating ACES application: 2-3 "Fenton's Fourier Series Wave Theory" (functional
area: Wave Theory). Steady progressive wave of permanent form, solved by an N-term
stream-function Fourier series (Rienecker & Fenton 1981; Fenton 1988).

Classification: exact (numerical solution of the exact nonlinear steady-wave problem --
an N-term stream-function Fourier series Newton-solved to convergence; no empirical
coefficients; reproduces the User's Guide example to the digit).
Theory and references: Rienecker & Fenton (1981); Fenton (1988a, 1988b, 1990).
  Equations transcribed in docs/EQUATIONS.md, TR chapter 2-3 (eqs 1-27 + Tables):
    stream function (10), Newton system (Table 2-3-2: KFSBC f_{m+9}, DFSBC f_{N+10+m},
    relations f_1-f_8), Newton solve (11)-(14) with height ramping, kinematics
    (16)-(21), integral properties (22)-(27).

  Method here follows Fenton (1988) in a compact non-dimensional (depth d, sqrt(gd))
  form: 2N+5 unknowns [kd, ubar/sqrt(gd), B_1..B_N, eta_0..eta_N (/d), Q, r], solved
  by Newton with a finite-difference Jacobian (TR eq 13) and height ramping (TR
  eq 12) for robustness near breaking. Collocation at m=0..N over the half wave
  (phase j*m*pi/N). The two free-surface BCs + wave-height + mean-zero + the
  celerity/current condition close the square system.

  Datum note (resolved against the manual): the DFSBC Bernoulli constant is `r`
  (SWL datum); the reported "Bernoulli Constant R" = r + g*d (bed datum), and
  `r = R - g*d` (TR symbol list, eq 20 pressure uses r). TR eq (5) writes "= R"
  loosely; the physically-correct SWL-datum form `= r` is used here and matches
  the User's Guide example exactly.

Self-containment: zero sibling imports; embeds its own Field/AppMeta/Result
dataclasses and the full Newton solver (numpy.linalg.solve only). Runnable
standalone:
    python chessqc_2_3_fenton_fourier.py
which runs self-tests (linear-limit + manual example), then prints an ACES-style
tabulation of the default example.  stdlib + numpy only.

I/O mirrors the ACES manual's Fenton lists:
  Inputs : H (>0), T (>0), d (>0), celerity def (Euler|Stokes), mean velocity,
           N terms (1-25), height-ramp steps (1-10), z, X/L
  Outputs: C, L, ubar1, ubar2, q, R, Q, I, E_K, E_P, E density, Ub^2, S_xx, F,
           eta, (U,W), (a_x,a_z), p  + profile arrays (eta,U,W vs X) for plots.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# --- standard physical constants (overridable; SI internal) ---------------------
G_SI = 9.80665           # m/s^2
RHO_SALT = 1025.18       # kg/m^3 (sea water; gamma ~ 64 lb/ft^3 as in the manual)
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
    aces_id="2-3",
    name="Fourier Series Wave Theory",
    area="Wave Theory",
    classification="exact",
    cite="Rienecker & Fenton (1981); Fenton (1988); TR 2-3",
    default_system="US",     # opens on the User's Guide Example (US units)
)

_FT = 0.3048
# Complete input list (ACES manual "Input" for Fenton's Fourier Series).
# Defaults are the ACES User's Guide Example 2-3 (H=4.5 ft, T=9 s, d=22 ft, Euler,
# ubar=0, N=16, 5 ramp steps; point of interest x/L=0, Z=-5 ft -> C=25.62 ft/s,
# L=230.58 ft, eta_crest=2.803 ft ...). See tests/test_manual_oracle.py.
INPUTS = (
    Field("H",  "Wave height",        "float", "m", "ft",  default=4.50 * _FT, lo=1e-6, hi=200 * _FT,
          note="> 0"),
    Field("T",  "Wave period",        "float", "s", "s",   default=9.0, lo=1.0, hi=1000.0,
          note="> 0"),
    Field("d",  "Water depth",        "float", "m", "ft",  default=22.0 * _FT, lo=1.0 * _FT, hi=5000 * _FT,
          note="> 0"),
    Field("cdef", "Celerity definition", "choice", "", "", default="Euler",
          choices=("Euler", "Stokes"),
          note="Euler = mean Eulerian current; Stokes = mean mass-transport velocity"),
    Field("ubar", "Mean velocity", "float", "m/s", "ft/s", default=0.0, lo=-10.0, hi=10.0,
          note="current (Euler) or mass-transport (Stokes) velocity"),
    Field("N",  "Fourier terms",     "int", "", "", default=16, lo=1, hi=25,
          note="number of terms (1 to 25)"),
    Field("nramp", "Height-ramp steps", "int", "", "", default=5, lo=1, hi=10,
          note="wave-height ramping steps (1 to 10)"),
    Field("z",  "Vertical coordinate","float", "m", "ft",  default=-5.0 * _FT,
          note="from SWL (z=0), +up; clamped to [-d, eta]"),
    Field("xL", "Wavelength fraction (X/L)", "float", "", "", default=0.0, lo=0.0, hi=1.0,
          note="0.0 to 1.0 (phase position; 0 = crest)"),
    Field("water", "Water type", "choice", "", "", default="Salt", choices=("Salt", "Fresh"),
          note="sets water density used for energies/pressure"),
)

# Complete output list (ACES manual "Output" for Fenton's Fourier Series).
OUTPUTS = (
    Out("C",     "Celerity",                       "m/s",   "ft/s",   "scalar",
        note="Wave celerity (phase speed), c = L/T, the speed the wave form travels."),
    Out("L",     "Wave length",                    "m",     "ft",     "scalar",
        note="Wavelength L = 2*pi/k, the horizontal distance between successive crests."),
    Out("ubar1", "Mean Eulerian fluid velocity",   "m/s",   "ft/s",   "scalar",
        note="Mean Eulerian current u_1, the time-averaged horizontal velocity at a fixed point."),
    Out("ubar2", "Mean mass-transport velocity",   "m/s",   "ft/s",   "scalar",
        note="Depth-averaged mass-transport (Stokes drift) velocity u_2 = c - Q/d."),
    Out("q",     "Volume flux (wave)",             "m^2/s", "ft^2/s", "scalar",
        note="Wave-induced volume flux per unit width, q = u_bar*d - Q (Stokes transport)."),
    Out("Q",     "Volume flux (mean)",             "m^2/s", "ft^2/s", "scalar",
        note="Mean volume flow rate per unit width in the steady (wave-following) frame."),
    Out("R",     "Bernoulli constant",             "m^2/s^2", "ft^2/s^2", "scalar",
        note="Bernoulli constant (bed datum), R = r + g*d, the energy constant of the DFSBC."),
    Out("I",     "Impulse",                        "N*s/m^2", "lb*s/ft^2", "scalar",
        note="Wave impulse (momentum) per unit area, I = rho*(c*d - Q)."),
    Out("EK",    "Kinetic energy",                 "N/m",   "lb/ft",  "scalar",
        note="Mean kinetic energy of the wave per unit horizontal area."),
    Out("EP",    "Potential energy",               "N/m",   "lb/ft",  "scalar",
        note="Mean potential energy per unit horizontal area, EP = 0.5*gamma*mean(eta^2)."),
    Out("Edens", "Energy density",                 "N/m",   "lb/ft",  "scalar",
        note="Total mean wave energy density per unit horizontal area, EK + EP."),
    Out("Ub2",   "Mean square of bed velocity",    "m^2/s^2", "ft^2/s^2", "scalar",
        note="Time-mean of the squared horizontal fluid velocity at the bed."),
    Out("Sxx",   "Radiation stress",               "N/m",   "lb/ft",  "scalar",
        note="Radiation stress S_xx, the wave-induced excess flux of horizontal momentum."),
    Out("F",     "Wave power (energy flux)",       "N/s",   "lb/s",   "scalar",
        note="Mean wave energy flux (wave power) per unit crest width."),
    Out("eta",   "Surface elevation",              "m",     "ft",     "point",
        note="Free-surface elevation above SWL at the chosen phase X/L (+up)."),
    Out("U",     "Horizontal velocity",            "m/s",   "ft/s",   "point",
        note="Horizontal fluid velocity (fixed frame) at the point (z, X/L), +in wave direction."),
    Out("W",     "Vertical velocity",              "m/s",   "ft/s",   "point",
        note="Vertical fluid velocity at the point (z, X/L), +up."),
    Out("ax",    "Horizontal acceleration",        "m/s^2", "ft/s^2", "point",
        note="Horizontal fluid acceleration at the point (z, X/L)."),
    Out("az",    "Vertical acceleration",          "m/s^2", "ft/s^2", "point",
        note="Vertical fluid acceleration at the point (z, X/L), +up."),
    Out("p",     "Pressure",                       "Pa",    "psf",    "point",
        note="Total fluid pressure at the point (z, X/L), from the Bernoulli equation."),
    Out("profile_X",   "Profile: X (+/- one wavelength)", "m",   "ft",   "profile",
        note="Horizontal coordinate spanning +/- one wavelength for the profile plots."),
    Out("profile_eta", "Profile: surface elevation",      "m",   "ft",   "profile",
        note="Surface elevation eta vs X over the wave profile (+up from SWL)."),
    Out("profile_U",   "Profile: horizontal velocity",    "m/s", "ft/s", "profile",
        note="Horizontal velocity at depth z vs X along the wave profile."),
    Out("profile_W",   "Profile: vertical velocity",      "m/s", "ft/s", "profile",
        note="Vertical velocity at depth z vs X along the wave profile."),
)


@dataclass
class Result:
    # scalars (integral properties)
    C: float; L: float; ubar1: float; ubar2: float; q: float; Q: float; R: float
    I: float; EK: float; EP: float; Edens: float; Ub2: float; Sxx: float; F: float
    # point kinematics at (z, X/L)
    eta: float; U: float; W: float; ax: float; az: float; p: float
    # profile arrays (for plots), over X in [-L, L]
    profile_X: np.ndarray
    profile_eta: np.ndarray
    profile_U: np.ndarray
    profile_W: np.ndarray
    notes: str = ""


# --- linear dispersion (initial guess for the modulus kd) -----------------------
def _linear_kd(tau: float) -> float:
    """Solve linear dispersion kd*tanh(kd) = (2 pi / tau)^2 for kd
    (tau = T*sqrt(g/d) dimensionless period)."""
    rhs = (2.0 * math.pi / tau) ** 2
    lo, hi = 1e-6, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid * math.tanh(mid) < rhs:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --- Fenton stream-function Newton solver (non-dimensional: lengths/d, vel/sqrt(gd))
def _solve_fenton(Htil: float, tau: float, U1: float, U2: float,
                  N: int, cdef: str, nramp: int):
    """Return the converged non-dimensional unknown vector v (length 2N+5).

    v = [K=kd, U=ubar/sqrt(gd), b_1..b_N, e_0..e_N (=eta_m/d), Q, r].
    Htil = H/d (target).  tau = T sqrt(g/d).  U1,U2 = current/transport, /sqrt(gd).
    """
    M = 2 * N + 5
    jj = np.arange(1, N + 1, dtype=float)              # j = 1..N
    mm = np.arange(0, N + 1, dtype=float)              # m = 0..N
    phase = np.outer(mm, jj) * math.pi / N             # (N+1, N): j m pi/N
    cph, sph = np.cos(phase), np.sin(phase)

    def resid(v, Hcur):
        K = v[0]; U = v[1]
        b = v[2:2 + N]
        e = v[2 + N:2 * N + 3]                          # e_0..e_N (N+1)
        Q = v[2 * N + 3]; r = v[2 * N + 4]
        c = 2.0 * math.pi / (K * tau)                   # celerity / sqrt(gd)
        # stable cosh/sinh ratios at each collocation point (N+1, N)
        A = np.outer((1.0 + e), jj) * K                 # j*K*(1+e_m)
        B = jj * K                                      # j*K
        denom = 1.0 + np.exp(-2.0 * B)
        eApB = np.exp(A - B)
        emApB = np.exp(-A - B)
        Cj = (eApB + emApB) / denom                     # cosh(A)/cosh(B)
        Sj = (eApB - emApB) / denom                     # sinh(A)/cosh(B)
        sK = math.sqrt(K)
        psi = -U * (1.0 + e) + (1.0 / K ** 1.5) * np.sum(b * Sj * cph, axis=1)
        u = -U + (1.0 / sK) * np.sum(jj * b * Cj * cph, axis=1)
        w = (1.0 / sK) * np.sum(jj * b * Sj * sph, axis=1)
        F = np.empty(M)
        F[0:N + 1] = psi + Q                            # KFSBC: psi(x_m, eta_m) = -Q
        F[N + 1:2 * N + 2] = 0.5 * (u * u + w * w) + e - r   # DFSBC (SWL-datum r)
        F[2 * N + 2] = e[0] - e[N] - Hcur               # wave height
        F[2 * N + 3] = e[0] + e[N] + 2.0 * np.sum(e[1:N])    # mean elevation = 0
        if cdef == "Euler":
            F[2 * N + 4] = U - c + U1                   # c = ubar + ubar1
        else:
            F[2 * N + 4] = c - Q - U2                   # ubar2 = c - Q/d
        return F

    def jac(v, Hcur):
        F0 = resid(v, Hcur)
        J = np.empty((M, M))
        for j in range(M):
            dj = v[j] / 100.0 if abs(v[j]) > 1e-4 else 1e-3   # TR eq 13
            vp = v.copy(); vp[j] += dj
            J[:, j] = (resid(vp, Hcur) - F0) / dj
        return J

    # initial guess at the first ramp height (linear theory)
    K = _linear_kd(tau)
    c = 2.0 * math.pi / (K * tau)
    U = c - U1
    H1 = Htil / nramp
    v = np.zeros(M)
    v[0] = K
    v[1] = U
    v[2] = U * (H1 / 2.0) * K ** 1.5 / math.tanh(K)     # linear B_1
    v[2 + N:2 * N + 3] = (H1 / 2.0) * np.cos(mm * math.pi / N)
    v[2 * N + 3] = U                                    # Q ~ ubar
    v[2 * N + 4] = 0.5 * U * U                          # r ~ 1/2 ubar^2

    # Newton with height ramping (TR eq 12)
    for step in range(1, nramp + 1):
        Hcur = Htil * step / nramp
        for _ in range(30):
            F = resid(v, Hcur)
            J = jac(v, Hcur)
            dv = np.linalg.solve(J, -F)
            v = v + dv
            if np.sum(np.abs(dv)) < 1e-11:              # TR eq 14 convergence
                break
        else:
            raise RuntimeError(f"Fenton Newton did not converge at ramp step {step}")
    return v


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes a steady progressive wave of permanent form by Newton-solving an '
            "N-term stream-function Fourier series (Fenton's exact nonlinear theory), "
            'returning celerity, wavelength, the full velocity/acceleration/pressure field '
            'at a point and along the profile, plus integral properties (energies, '
            'impulse, radiation stress, energy flux).',
 'methods': [{'name': 'Fenton stream-function Fourier series',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': '\\psi(x,z) = -\\bar{u}\\,(d+z) + '
                                    '\\sqrt{\\frac{g}{k^3}}\\,\\sum_{j=1}^{N} '
                                    'B_{j}\\,\\frac{\\sinh[\\,j k (d+z)\\,]}{\\cosh(j k '
                                    'd)}\\,\\cos(j k x)',
                             'desc': 'N-term Fourier cosine series for the steady-frame '
                                     'stream function (eq 10); the B_j coefficients are '
                                     'the primary unknowns.'},
                            {'tex': '\\psi(x,\\eta) = -Q',
                             'desc': 'Kinematic free-surface boundary condition (KFSBC, eq '
                                     '4): the surface is a streamline carrying constant '
                                     'volume flux Q.'},
                            {'tex': '\\frac{1}{2}\\left[\\left(\\frac{\\partial '
                                    '\\psi}{\\partial x}\\right)^2 + '
                                    '\\left(\\frac{\\partial \\psi}{\\partial '
                                    'z}\\right)^2\\right] + g\\,\\eta(x) = R',
                             'desc': 'Dynamic free-surface boundary condition (DFSBC, eq '
                                     '5): Bernoulli at the surface, collocated at N+1 '
                                     'points to close the system.'},
                            {'tex': 'c = \\frac{L}{T} = \\frac{2\\pi}{k T} = \\bar{u} + '
                                    '\\bar{u}_{1}',
                             'desc': 'Celerity / dispersion closure (eq 9); Euler form '
                                     'uses mean current u_1 (shown), Stokes form replaces '
                                     'it with c = Q/d + u_2.'},
                            {'tex': 'u(x,z) = -\\bar{u} + '
                                    '\\sqrt{\\frac{g}{k}}\\,\\sum_{j=1}^{N} '
                                    'j\\,B_{j}\\,\\frac{\\cosh[\\,j k (d+z)\\,]}{\\cosh(j '
                                    'k d)}\\,\\cos(j k x)',
                             'desc': 'Horizontal velocity field from the converged '
                                     'coefficients (eq 16); w, accelerations and pressure '
                                     'follow by differentiation.'}]}],
 'symbols': [['psi', 'Stream function in the steady (wave-following) frame'],
             ['B_j', 'Dimensionless Fourier coefficients (the Newton unknowns), j = 1..N'],
             ['k', 'Wave number, k = 2*pi/L'],
             ['c', 'Wave celerity, c = L/T'],
             ['Q', 'Constant volume flow rate per unit width in the steady frame'],
             ['R', 'Bernoulli constant (bed datum); r = R - g*d is the SWL-datum form'],
             ['u_bar',
              'Mean fluid velocity; u_1 mean Eulerian current, u_2 depth-averaged '
              'mass-transport velocity'],
             ['eta', 'Free-surface elevation above SWL'],
             ['d', 'Still-water depth'],
             ['N', 'Number of Fourier terms (1 to 25)']],
 'references': ['Rienecker & Fenton (1981)',
                'Fenton (1988a, 1988b)',
                'Fenton (1990)',
                'ACES Technical Reference, Chapter 2-3 (Fourier Series Wave Theory)']}


def compute(inp: dict, *, g: float = G_SI, rho: float | None = None,
            n_profile: int = 201) -> Result:
    """Fenton Fourier-series results for SI inputs.

    `inp` keys: H, T, d (SI), cdef ("Euler"|"Stokes"), ubar (SI, m/s), N (int),
    nramp (int), z (SI), xL, water ("Salt"|"Fresh"). The GUI converts US->SI at the
    edge. Water density follows `water` unless `rho` is given."""
    if rho is None:
        rho = RHO_FRESH if str(inp.get("water", "Salt")) == "Fresh" else RHO_SALT
    H, T, d, z, xL = (float(inp[k]) for k in ("H", "T", "d", "z", "xL"))
    cdef = str(inp.get("cdef", "Euler"))
    ubar_in = float(inp.get("ubar", 0.0))
    N = int(inp.get("N", 16))
    nramp = int(inp.get("nramp", 5))
    if H <= 0 or T <= 0 or d <= 0:
        raise ValueError("H, T, d must all be > 0")
    if not (1 <= N <= 25):
        raise ValueError(f"N = {N} must be 1..25")
    if not (1 <= nramp <= 10):
        raise ValueError(f"ramp steps = {nramp} must be 1..10")
    if not (0.0 <= xL <= 1.0):
        raise ValueError(f"X/L = {xL} outside [0, 1]")

    sgd = math.sqrt(g * d)
    Htil = H / d
    tau = T * math.sqrt(g / d)
    U1 = ubar_in / sgd if cdef == "Euler" else 0.0       # Eulerian current /sqrt(gd)
    U2 = ubar_in / sgd if cdef == "Stokes" else 0.0      # mass transport /sqrt(gd)

    v = _solve_fenton(Htil, tau, U1, U2, N, cdef, nramp)

    # unpack + redimensionalize
    K = v[0]; Ud = v[1]
    b = v[2:2 + N]
    e = v[2 + N:2 * N + 3]
    Qd = v[2 * N + 3]; rd = v[2 * N + 4]
    k = K / d
    L = 2.0 * math.pi / k
    c = 2.0 * math.pi / (k * T)                          # = L/T
    ubar = Ud * sgd
    eta_m = e * d
    Q = Qd * d * sgd                                    # mean volume flux
    r = rd * g * d                                       # SWL-datum Bernoulli
    R = r + g * d                                       # bed-datum (reported)
    ubar1 = ubar_in if cdef == "Euler" else (c - ubar)
    ubar2 = c - Q / d
    q = ubar * d - Q                                    # wave-induced flux
    gamma = rho * g

    # integral properties (TR eqs 22-27) -----------------------------------------
    # mean of eta^2 over the wavelength (trapezoid over the half wave, symmetric)
    eta2_mean = (0.5 * eta_m[0] ** 2 + np.sum(eta_m[1:N] ** 2)
                 + 0.5 * eta_m[N] ** 2) / N
    EP = 0.5 * gamma * eta2_mean                         # eq 22
    Iimp = rho * (c * d - Q)                             # eq 23
    EK = 0.5 * (c * Iimp - rho * ubar1 * Q)              # eq 24
    Ub2 = 2.0 * r - c * c + 2.0 * ubar1 * c              # eq 25 (r = R - g d)
    F = ((3.0 * EK - 2.0 * EP - 2.0 * ubar1 * Iimp) * c
         + 0.5 * Ub2 * (Iimp + rho * c * d))             # eq 26
    Sxx = 4.0 * EK - 3.0 * EP + rho * d * Ub2 - 2.0 * ubar1 * Iimp   # eq 27
    Edens = EK + EP

    # ---- kinematics field at (x, z) (steady frame u,w; fixed frame U=u+c) -------
    jj = np.arange(1, N + 1, dtype=float)

    def kin(xphase, zz):
        """xphase = k*x (scalar or array); zz scalar depth. Returns U,W,ax,az,p,u,w."""
        xp = np.atleast_1d(np.asarray(xphase, dtype=float))
        A = jj * k * (d + zz)                            # j k (d+z)
        B = jj * k * d
        denom = 1.0 + np.exp(-2.0 * B)
        Cj = (np.exp(A - B) + np.exp(-A - B)) / denom     # cosh/cosh
        Sj = (np.exp(A - B) - np.exp(-A - B)) / denom     # sinh/cosh
        cjx = np.cos(np.outer(xp, jj))                    # (npts, N)
        sjx = np.sin(np.outer(xp, jj))
        coef = math.sqrt(g / k)
        u = -ubar + coef * (jj * b * Cj * cjx).sum(axis=1)            # eq 16
        w = coef * (jj * b * Sj * sjx).sum(axis=1)                     # eq 17
        # derivatives (eq 19): ux = -(gk)^.5 sum j^2 b Cj sin ; uz = (gk)^.5 sum j^2 b Sj cos
        cgk = math.sqrt(g * k)
        ux = -cgk * (jj ** 2 * b * Cj * sjx).sum(axis=1)
        uz = cgk * (jj ** 2 * b * Sj * cjx).sum(axis=1)
        wx = uz                                          # irrotational: w_x = u_z
        wz = -ux                                         # continuity: w_z = -u_x
        ax = u * ux + w * uz                             # eq 18
        az = u * wx + w * wz                             # eq 19
        p = rho * r - rho * g * zz - 0.5 * rho * (u * u + w * w)       # eq 20
        Ufix = u + c
        return Ufix, w, ax, az, p, u

    # ---- point of interest (z, X/L) --------------------------------------------
    eta_crest = eta_m[0]
    z_use = min(max(z, -d), eta_crest)
    xphase_pt = 2.0 * math.pi * xL                       # k*x = 2 pi (x/L)
    Up, Wp, axp, azp, pp, _ = kin(xphase_pt, z_use)
    Up, Wp, axp, azp, pp = (float(Up[0]), float(Wp[0]), float(axp[0]),
                            float(azp[0]), float(pp[0]))
    # surface elevation eta(x): Fourier coefficients f_j (TR eq 21), inverse
    # cosine transform of the nodal elevations
    mm = np.arange(0, N + 1)
    cos_jm = np.cos(np.outer(jj, mm) * math.pi / N)      # (N, N+1)
    weights = np.ones(N + 1); weights[0] = 0.5; weights[N] = 0.5
    fj = (2.0 / N) * (cos_jm * (eta_m * weights)).sum(axis=1)   # f_1..f_N
    eta_pt = float((fj * np.cos(jj * xphase_pt)).sum())

    # ---- profile arrays over +/- one wavelength --------------------------------
    X = np.linspace(-L, L, n_profile)
    xph = k * X
    prof_eta = (fj[None, :] * np.cos(np.outer(xph, jj))).sum(axis=1)
    prof_U, prof_W, _, _, _, _ = kin(xph, z_use)

    notes = [f"N={N} terms, {nramp} ramp steps, {cdef} celerity"]
    Hmax = d * ((0.141063 * (L / d) + 0.0095721 * (L / d) ** 2 + 0.0077829 * (L / d) ** 3)
                / (1 + 0.078834 * (L / d) + 0.0317567 * (L / d) ** 2 + 0.0093407 * (L / d) ** 3))
    if H > Hmax:
        notes.append(f"H exceeds Fenton H_max={Hmax:.3f} m (wave may be unstable)")
    if z + d < 0:
        notes.append("WARNING: point outside waveform (z below the bed, z + d < 0)")
    elif z != z_use:
        notes.append(f"z clamped to {z_use:.3f} m for kinematics")

    return Result(C=c, L=L, ubar1=ubar1, ubar2=ubar2, q=q, Q=Q, R=R, I=Iimp,
                  EK=EK, EP=EP, Edens=Edens, Ub2=Ub2, Sxx=Sxx, F=F,
                  eta=eta_pt, U=Up, W=Wp, ax=axp, az=azp, p=pp,
                  profile_X=X, profile_eta=prof_eta, profile_U=prof_U, profile_W=prof_W,
                  notes="; ".join(notes))


# --- self-tests -----------------------------------------------------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    # 1) linear limit: a very small wave (N=4) -> linear celerity c^2 = g/k tanh(kd)
    r = compute({"H": 0.01, "T": 8.0, "d": 10.0, "cdef": "Euler", "ubar": 0.0,
                 "N": 4, "nramp": 1, "z": 0.0, "xL": 0.0})
    k = 2 * math.pi / r.L
    c_lin = math.sqrt(G_SI / k * math.tanh(k * 10.0))
    assert _approx(r.C, c_lin, 1e-3), (r.C, c_lin)

    # 2) manual example (US units, first/second screens) reproduced
    r = compute({"H": 4.50 * _FT, "T": 9.0, "d": 22.0 * _FT, "cdef": "Euler",
                 "ubar": 0.0, "N": 16, "nramp": 5, "z": -5.0 * _FT, "xL": 0.0})
    assert _approx(r.C / _FT, 25.620, 5e-3), r.C / _FT
    assert _approx(r.L / _FT, 230.581, 5e-3), r.L / _FT
    assert _approx(r.eta / _FT, 2.803, 1e-2), r.eta / _FT
    assert _approx(r.U / _FT, 3.150, 2e-2), r.U / _FT
    assert abs(r.W) < 1e-6, r.W
    assert _approx(r.az / _FT, -1.302, 2e-2), r.az / _FT
    assert _approx(r.ubar2 / _FT, 0.140, 5e-2), r.ubar2 / _FT
    # bed velocity, energies (rough tol; manual rounds + uses gamma=64)
    assert _approx(r.Ub2 / (_FT ** 2), 2.6479, 2e-2), r.Ub2 / (_FT ** 2)

    print("  self-tests: PASS (linear limit + User's Guide Example 2-3)")


def _print_default_example() -> None:
    """ACES-style tabulation of the default example (SI)."""
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        vv = inp[f.key]
        sval = f"{vv:>10.4g}" if isinstance(vv, (int, float)) and f.kind not in ("choice",) else f"{vv:>10}"
        print(f"    {f.label:28s} {f.key:6s} = {sval} {f.unit_si}")
    print("  OUTPUTS:")
    scal = {"C": r.C, "L": r.L, "ubar1": r.ubar1, "ubar2": r.ubar2, "q": r.q,
            "Q": r.Q, "R": r.R, "I": r.I, "EK": r.EK, "EP": r.EP, "Edens": r.Edens,
            "Ub2": r.Ub2, "Sxx": r.Sxx, "F": r.F}
    pts = {"eta": r.eta, "U": r.U, "W": r.W, "ax": r.ax, "az": r.az, "p": r.p}
    by_key = {o.key: o for o in OUTPUTS}
    for kk, vv in {**scal, **pts}.items():
        o = by_key[kk]
        print(f"    {o.label:30s} {kk:6s} = {vv:>12.5g} {o.unit_si}")
    print(f"  notes: {r.notes}")
    print(f"  profile arrays: X/eta/U/W each length {len(r.profile_X)} (over +/- one L)")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
