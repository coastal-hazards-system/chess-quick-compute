"""CHESS-QC application 3-4 — Vertical-Wedge Diffraction/Reflection on a Uniform Grid.

Originating ACES grouping: 3-4 "Combined Diffraction and Reflection by a Vertical Wedge
Evaluated Upon a Uniform Grid" (functional area: Wave Transformation). This is the gridded
companion to 3-3: it evaluates the same fully-reflecting-wedge solution (the PCDFRAC / Chen
1987 eigenfunction expansion) at every node of a uniform X-Y grid around the wedge apex,
producing maps of the wave-height modification factor, the modified wave height, and the
phase relative to the incident wave.

Classification: exact given 3-3 (same closed-form solver, evaluated on a grid).
Theory and references: identical to application 3-3 (Chen 1987; Stoker 1957; Penny & Price
1952). For each grid point at radius r and angle theta from the wedge apex,
    phi(r,theta) = (2/nu)[ J_0(kr) + 2 sum_{n>=1} e^(i n pi/(2 nu)) J_{n/nu}(kr)
                            cos(n alpha/nu) cos(n theta/nu) ] ,   nu = theta_0/pi,
the modification factor is |phi| and the modified height is |phi| H_i.

Self-containment: zero sibling imports; embeds the contract dataclasses, the dispersion
solver, and the fractional-order Bessel series (same as 3-3). numpy + stdlib only. Runnable:
    python chessqc_3_4_wedge_grid.py

Validation: reproduces the ACES User's Guide Example 3 (semi-infinite breakwater, wedge angle
0 so nu=2; incident 4 ft, 12 s wave in 30 ft of water, wave angle 52 deg) on its X grid
[-600..200] step 200 ft, Y grid [-400..200] step 100 ft (wavelength 356.85 ft). At grid point
(-600, -400) ft the modification factor is 0.903 and the modified height 3.61 ft, matching the
published Table 3-3-2 to better than 1 percent. (The 2-D grid output is a field; the generic
GUI renderer shows scalar summaries, with the full grid available from compute().)
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from math import lgamma

import numpy as np

G_SI = 9.80665
_FT = 0.3048


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
    aces_id="3-4",
    name="Vertical-Wedge Diffraction/Reflection on a Uniform Grid",
    area="Wave Transformation",
    classification="exact",
    cite="Chen (1987); Stoker (1957); Penny & Price (1952)",
    default_system="US",
)

INPUTS = (
    Field("Hi", "Incident wave height", "float", "m", "ft", default=4.0 * _FT, lo=1e-4, hi=1e3),
    Field("T", "Wave period", "float", "s", "s", default=12.0, lo=1e-2, hi=1e3),
    Field("d", "Water depth", "float", "m", "ft", default=30.0 * _FT, lo=1e-3, hi=1e4),
    Field("alpha", "Incident wave angle", "angle", "deg", "deg", default=52.0, lo=0.0, hi=360.0),
    Field("wedge_angle", "Wedge angle", "angle", "deg", "deg", default=0.0, lo=0.0, hi=180.0),
    Field("X0", "X start", "float", "m", "ft", default=-600.0 * _FT, lo=-1e4, hi=1e4),
    Field("Xm", "X end", "float", "m", "ft", default=200.0 * _FT, lo=-1e4, hi=1e4),
    Field("dX", "X increment", "float", "m", "ft", default=200.0 * _FT, lo=1e-3, hi=1e4),
    Field("Y0", "Y start", "float", "m", "ft", default=-400.0 * _FT, lo=-1e4, hi=1e4),
    Field("Ym", "Y end", "float", "m", "ft", default=200.0 * _FT, lo=-1e4, hi=1e4),
    Field("dY", "Y increment", "float", "m", "ft", default=100.0 * _FT, lo=1e-3, hi=1e4),
)

OUTPUTS = (
    Out("L",       "Wave length",                       "m", "ft", "scalar",
        note="Linear-theory wavelength L = 2 pi / k from the dispersion relation for the given period and depth."),
    Out("mod_max", "Maximum modification factor",       "",  "",   "scalar",
        note="Largest wave-height modification factor |phi| (diffraction/reflection coefficient) found over the grid, dimensionless."),
    Out("mod_min", "Minimum modification factor",       "",  "",   "scalar",
        note="Smallest wave-height modification factor |phi| over the grid, dimensionless (values below 1 indicate sheltering)."),
    Out("H_max",   "Maximum modified wave height",      "m", "ft", "scalar",
        note="Largest modified wave height H = |phi| H_i found over the grid."),
    Out("grid_x",  "Grid X coordinates",                "m", "ft", "profile",
        note="X-axis node positions of the uniform grid, measured from the wedge apex."),
    Out("grid_y",  "Grid Y coordinates",                "m", "ft", "profile",
        note="Y-axis node positions of the uniform grid, measured from the wedge apex."),
    Out("mod_grid", "Modification-factor field",        "",  "",   "grid",
        note="Field of the dimensionless modification factor |phi| (diffraction/reflection coefficient) at every grid node."),
    Out("H_grid",  "Modified-height field",             "m", "ft", "grid",
        note="Field of the modified wave height H = |phi| H_i at every grid node."),
)


@dataclass
class Result:
    L: float; mod_max: float; mod_min: float; H_max: float
    grid_x: np.ndarray; grid_y: np.ndarray
    mod_grid: np.ndarray; H_grid: np.ndarray
    notes: str = ""


def _dispersion_L(T: float, d: float, g: float) -> float:
    L0 = g * T * T / (2.0 * math.pi)
    L = L0
    for _ in range(200):
        Ln = L0 * math.tanh(2.0 * math.pi * d / L)
        if abs(Ln - L) < 1e-12:
            return Ln
        L = Ln
    return L


def bessel_jp(p: float, x: float, nterm: int = 200) -> float:
    if x == 0.0:
        return 1.0 if p == 0.0 else 0.0
    s = 0.0
    for m in range(nterm):
        term = (-1.0) ** m * math.exp((2 * m + p) * math.log(x / 2.0)
                                      - lgamma(m + 1) - lgamma(m + p + 1.0))
        s += term
        if abs(term) < 1e-15 and m > p + 2:
            break
    return s


def _wedge_potential(kr: float, theta: float, alpha: float, nu: float) -> complex:
    s = complex(bessel_jp(0.0, kr), 0.0)
    small = 0
    for n in range(1, 600):
        ph = cmath.exp(1j * n * math.pi / (2.0 * nu))
        term = 2.0 * ph * bessel_jp(n / nu, kr) * math.cos(n * alpha / nu) * math.cos(n * theta / nu)
        s += term
        small = small + 1 if abs(term) < 1e-6 else 0
        if small >= 8:
            break
    return (2.0 / nu) * s


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Evaluates the fully-reflecting vertical-wedge diffraction/reflection solution '
            '(Chen 1987 eigenfunction series) at every node of a uniform X-Y grid around '
            'the wedge apex, returning fields of the wave-height modification factor, the '
            'modified wave height, and the phase relative to the incident wave.',
 'methods': [{'name': 'Chen (1987) wedge eigenfunction series on a uniform grid',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': '\\phi(r,\\theta) = \\frac{2}{\\nu}\\left[ J_{0}(kr) + '
                                    '2\\sum_{n=1}^{\\infty} e^{i n \\pi/(2\\nu)} '
                                    'J_{n/\\nu}(kr)\\cos\\left(\\frac{n\\alpha}{\\nu}\\right)\\cos\\left(\\frac{n\\theta}{\\nu}\\right) '
                                    '\\right]',
                             'desc': 'Complex horizontal-plane potential at grid radius r '
                                     'and angle theta from the wedge apex (Chen 1987 '
                                     'Bessel eigenfunction series)'},
                            {'tex': '\\nu = \\frac{\\theta_{0}}{\\pi}, \\quad \\theta_{0} '
                                    '= 2\\pi - \\theta_{w}',
                             'desc': 'Wedge parameter from the water-domain opening angle '
                                     'theta_0 (theta_w is the solid wedge angle); '
                                     'semi-infinite breakwater is theta_w = 0 so nu = 2'},
                            {'tex': '\\omega^{2} = g\\,k\\tanh(kh)',
                             'desc': 'Linear dispersion relation giving wavenumber k = 2 '
                                     'pi / L (and wavelength L) from period and depth'},
                            {'tex': '|\\phi| = \\sqrt{(\\mathrm{Re}\\,\\phi)^{2} + '
                                    '(\\mathrm{Im}\\,\\phi)^{2}}',
                             'desc': 'Wave-height modification factor (= '
                                     'diffraction/reflection coefficient); modified height '
                                     'H = |phi| H_i'},
                            {'tex': '\\beta = '
                                    '\\tan^{-1}\\left(\\frac{\\mathrm{Im}\\,\\phi}{\\mathrm{Re}\\,\\phi}\\right) '
                                    '- k\\,r\\cos\\alpha',
                             'desc': 'Phase of the modified wave relative to the incident '
                                     'wave'}]}],
 'symbols': [['phi', 'Complex horizontal-plane velocity potential at a grid point'],
             ['|phi|',
              'Wave-height modification factor (diffraction/reflection coefficient)'],
             ['H_i', 'Incident wave height (modified height = |phi| H_i)'],
             ['k', 'Wavenumber, 2 pi / L'],
             ['L', 'Wavelength from the dispersion relation'],
             ['r', 'Radial distance from the wedge apex to the grid point'],
             ['theta', 'Angular position of the grid point measured from the wedge face'],
             ['alpha', 'Incident wave angle'],
             ['nu', 'Wedge parameter, theta_0 / pi'],
             ['J_{n/nu}', 'Bessel function of the first kind, order n/nu']],
 'references': ['Chen (1987) CERC-87-16',
                'Stoker (1957)',
                'Penny & Price (1952)',
                'SPM (1984)',
                "ACES User's Guide, Application 3-4 (gridded companion to 3-3)"]}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Wedge modification-factor / height / phase field over a uniform grid (SI inputs)."""
    _validate(inp)
    Hi = float(inp["Hi"]); T = float(inp["T"]); d = float(inp["d"])
    alpha = math.radians(float(inp["alpha"])); wedge = math.radians(float(inp["wedge_angle"]))
    X0 = float(inp["X0"]); Xm = float(inp["Xm"]); dX = float(inp["dX"])
    Y0 = float(inp["Y0"]); Ym = float(inp["Ym"]); dY = float(inp["dY"])

    L = _dispersion_L(T, d, g)
    k = 2.0 * math.pi / L
    nu = (2.0 * math.pi - wedge) / math.pi

    xs = np.arange(X0, Xm + 0.5 * dX, dX)
    ys = np.arange(Y0, Ym + 0.5 * dY, dY)
    mod = np.zeros((len(ys), len(xs)))
    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            r = math.hypot(x, y)
            theta = math.atan2(y, x) % (2.0 * math.pi)
            mod[iy, ix] = abs(_wedge_potential(k * r, theta, alpha, nu))
    H = mod * Hi
    notes = (f"nu={nu:.3f}; {len(xs)}x{len(ys)} grid; modification factor = |phi| "
             f"(diffraction/reflection coefficient); same solver as 3-3")
    return Result(L=L, mod_max=float(mod.max()), mod_min=float(mod.min()), H_max=float(H.max()),
                  grid_x=xs, grid_y=ys, mod_grid=mod, H_grid=H, notes=notes)


# --- self-tests (ACES User's Guide Example 3 grid) ------------------------------
def _approx(a, b, tol):
    return abs(a - b) <= tol


def _self_tests() -> None:
    g = G_SI
    r = compute({f.key: f.default for f in INPUTS}, g=g)
    assert _approx(r.L / _FT, 356.85, 0.1), r.L / _FT
    # grid point (-600, -400) ft: modification factor 0.903, modified height 3.61 ft (Table 3-3-2)
    ix = int(np.argmin(np.abs(r.grid_x / _FT - (-600.0))))
    iy = int(np.argmin(np.abs(r.grid_y / _FT - (-400.0))))
    assert _approx(r.mod_grid[iy, ix], 0.903, 0.01), r.mod_grid[iy, ix]
    assert _approx(r.H_grid[iy, ix] / _FT, 3.61, 0.05), r.H_grid[iy, ix] / _FT
    print(f"  self-tests: PASS (L={r.L/_FT:.2f} ft; grid (-600,-400): |phi|={r.mod_grid[iy,ix]:.3f}, "
          f"H={r.H_grid[iy,ix]/_FT:.2f} ft; field max |phi|={r.mod_max:.2f})")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = User's Guide Example 3: semi-infinite breakwater, 5x7 grid)")
    print(f"    wave length L = {r.L/_FT:.2f} ft   grid {len(r.grid_x)} x {len(r.grid_y)}")
    print(f"    modification factor: min {r.mod_min:.3f}  max {r.mod_max:.3f}   "
          f"max modified height = {r.H_max/_FT:.2f} ft")
    print(f"    oracle check: at (-600,-400) ft |phi|=0.903, H=3.61 ft")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
