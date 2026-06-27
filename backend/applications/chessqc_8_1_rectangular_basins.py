"""CHESS-QC application 8-1 — Properties of Rectangular Basins.

Originating ACES grouping: 8-1 "Properties of Rectangular Basins" (functional area:
Harbor Design). Computes the natural (resonant) oscillation periods of a rectangular
harbor or basin, the seiche modes that a long wave, surge, or tsunami can excite, plus
the standing-wave water-particle kinematics at a node.

Classification: exact (closed-form shallow-water resonance theory).
Theory and references: Merian's formula for closed and open rectangular basins; the
two-dimensional extension; the Helmholtz (pumping) mode for a basin connected to the
sea by a channel. Standard references: SPM (1984); Wilson (1972); Sorensen (1993);
Ippen (1966). There is no ACES Technical Reference chapter and no User's Guide worked
example for this application (the Harbor Design tools were externally sourced), so the
validation oracle is the textbook closed forms themselves.

Resonant periods use the shallow-water celerity c = sqrt(g d), the appropriate limit
for seiche wavelengths that are long compared with the depth.

  Open basin (one end open, quarter-wave):  T_n = 4 lB / ((2n+1) sqrt(g d)),  n = 0,1,2,...
  Closed basin (both ends, half-wave):       T_n = 2 lB / (n sqrt(g d)),       n = 1,2,3,...
  Closed 2-D (Merian):  T_(n,m) = 2 / ( sqrt(g d) sqrt( (n/lB)^2 + (m/lC)^2 ) )
  Helmholtz (pumping):  T_H = 2 pi sqrt( (L_ch + L_corr) A_b / (g A_c) )

Self-containment: zero sibling imports; embeds its own contract dataclasses. Runnable
standalone:
    python chessqc_8_1_rectangular_basins.py
which runs the closed-form self-tests, then prints an ACES-style tabulation.
stdlib + numpy only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# --- standard physical constants (overridable; SI internal) ---------------------
G_SI = 9.80665           # m/s^2


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
    kind: str = "float"            # float | int | choice | bool | angle | file
    unit_si: str = ""
    unit_us: str = ""
    default: object = 0.0
    lo: float = -math.inf
    hi: float = math.inf
    choices: tuple = ()
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


# --- application metadata --------------------------------------------------------
APP_META = AppMeta(
    aces_id="8-1",
    name="Properties of Rectangular Basins",
    area="Harbor Design",
    classification="exact",
    cite="Merian's formula; Helmholtz mode; SPM (1984); Wilson (1972)",
    default_system="US",
)

_FT = 0.3048
BASIN_TYPES = (
    "Open (one end open)",
    "Closed (both ends)",
    "Closed 2-D (rectangular)",
    "Helmholtz (basin + channel)",
)
INPUTS = (
    Field("basin_type", "Basin type", "choice", "", "", default="Closed (both ends)",
          choices=BASIN_TYPES, note="open/closed 1-D, closed 2-D, or Helmholtz pumping"),
    Field("lB", "Basin length", "float", "m", "ft", default=3000.0 * _FT, lo=1e-3, hi=1e6,
          note="length along the resonant (longitudinal) axis"),
    Field("lC", "Basin width", "float", "m", "ft", default=2000.0 * _FT, lo=1e-3, hi=1e6,
          note="2-D only: transverse dimension",
          enable_if=("basin_type", "Closed 2-D (rectangular)")),
    Field("d", "Water depth", "float", "m", "ft", default=30.0 * _FT, lo=1e-3, hi=1e4,
          note="mean basin depth; resonance uses shallow-water c = sqrt(g d)"),
    Field("H", "Standing-wave height", "float", "m", "ft", default=1.0 * _FT, lo=0.0, hi=1e3,
          note="antinode crest-to-trough height; used for the node kinematics"),
    Field("n", "Longitudinal mode n", "int", "", "", default=1, lo=0, hi=50,
          note="closed/2-D: n>=1; open: n>=0 (n=0 is the fundamental)"),
    Field("m_mode", "Transverse mode m", "int", "", "", default=1, lo=0, hi=50,
          note="2-D only: m>=0 (m=0 reduces to the 1-D longitudinal mode)",
          enable_if=("basin_type", "Closed 2-D (rectangular)")),
    # Helmholtz-only geometry
    Field("Ab", "Basin surface area", "float", "m^2", "ft^2", default=3000.0 * 2000.0 * _FT * _FT,
          lo=1e-3, hi=1e12, note="Helmholtz only: plan area of the basin",
          enable_if=("basin_type", "Helmholtz (basin + channel)")),
    Field("Ac", "Channel cross-section area", "float", "m^2", "ft^2",
          default=200.0 * 30.0 * _FT * _FT, lo=1e-3, hi=1e10,
          note="Helmholtz only: flow area of the entrance channel",
          enable_if=("basin_type", "Helmholtz (basin + channel)")),
    Field("L_ch", "Channel length", "float", "m", "ft", default=500.0 * _FT, lo=0.0, hi=1e5,
          note="Helmholtz only: length of the entrance channel",
          enable_if=("basin_type", "Helmholtz (basin + channel)")),
    Field("L_corr", "Mouth length correction", "float", "m", "ft", default=100.0 * _FT,
          lo=0.0, hi=1e5, note="Helmholtz only: added-mass correction at the mouth",
          enable_if=("basin_type", "Helmholtz (basin + channel)")),
)

OUTPUTS = (
    Out("c_shallow", "Shallow-water celerity sqrt(g d)", "m/s", "ft/s", "scalar",
        note="shallow-water wave speed c = sqrt(g d) used to set every seiche period"),
    Out("T_mode", "Resonant period of selected mode", "s", "s", "scalar",
        note="natural oscillation period of the chosen seiche mode (n, or n,m for 2-D)"),
    Out("L_mode", "Modal wavelength (1-D)", "m", "ft", "scalar",
        note="wavelength of the selected 1-D longitudinal mode (zero for genuine 2-D modes)"),
    Out("f_mode", "Resonant frequency of selected mode", "1/s", "1/s", "scalar",
        note="resonant frequency 1/T_mode of the selected mode"),
    Out("T_helm", "Helmholtz (pumping) period", "s", "s", "scalar",
        note="Helmholtz pumping period for a basin connected to the sea by a channel; zero unless Helmholtz type"),
    Out("Vmax", "Max horizontal velocity at node (surface)", "m/s", "ft/s", "scalar",
        note="peak horizontal water-particle velocity at the standing-wave node, at the surface"),
    Out("excursion", "Horizontal particle semi-excursion at node", "m", "ft", "scalar",
        note="amplitude of horizontal particle motion at the node, Vmax/omega"),
    Out("Vbar", "Mean horizontal speed at node", "m/s", "ft/s", "scalar",
        note="time-mean horizontal speed over a cycle at the node, (2/pi) Vmax"),
    Out("mode_index", "Mode table: mode number", "", "", "profile",
        note="mode number for each row of the seiche spectrum table"),
    Out("mode_T", "Mode table: resonant period", "s", "s", "profile",
        note="resonant period of the first several modes (a small seiche spectrum)"),
)


@dataclass
class Result:
    c_shallow: float
    T_mode: float
    L_mode: float
    f_mode: float
    T_helm: float
    Vmax: float
    excursion: float
    Vbar: float
    mode_index: np.ndarray
    mode_T: np.ndarray
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp.get(f.key, f.default))         # Helmholtz fields are type-specific
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")


def _period_1d(basin_type: str, j: int, lB: float, c: float) -> float:
    """Resonant period of longitudinal mode j (open: quarter-wave; closed: half-wave)."""
    if basin_type.startswith("Open"):
        return 4.0 * lB / ((2 * j + 1) * c)          # j = 0,1,2,...
    return 2.0 * lB / (j * c)                          # closed, j = 1,2,3,...


# --- compute (the single entry point both front-ends call) ----------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes the natural (resonant) seiche periods of a rectangular harbor or '
            'basin and the standing-wave water-particle kinematics at a node. Supports '
            'open and closed 1-D basins, the 2-D (Merian) extension, and the Helmholtz '
            'pumping mode for a basin connected to the sea by an entrance channel.',
 'method_key': 'basin_type',
 'methods': [{'name': 'Open basin (one end open, quarter-wave)',
              'when': 'Open (one end open)',
              'tag': 'standard',
              'note': 'One open end; modes are odd quarter-wavelengths, n = 0,1,2,... with '
                      'n = 0 the fundamental.',
              'equations': [{'tex': 'c = \\sqrt{g\\,d}',
                             'desc': 'Shallow-water celerity used for all seiche periods '
                                     '(depth small vs wavelength).'},
                            {'tex': 'T_{n} = \\frac{4\\,l_{B}}{(2n+1)\\,\\sqrt{g\\,d}}',
                             'desc': 'Resonant period of longitudinal mode n for a basin '
                                     'open at one end.'},
                            {'tex': 'u_{max} = \\frac{\\pi H}{T}\\,\\frac{\\cosh(k '
                                    'd)}{\\sinh(k d)}',
                             'desc': 'Horizontal velocity amplitude at the node (surface) '
                                     'of the standing wave (clapotis).'}]},
             {'name': 'Closed basin (both ends, half-wave) — Merian',
              'when': 'Closed (both ends)',
              'tag': 'standard',
              'note': 'Both ends closed; modes are integer half-wavelengths, n = 1,2,3,...',
              'equations': [{'tex': 'c = \\sqrt{g\\,d}', 'desc': 'Shallow-water celerity.'},
                            {'tex': 'T_{n} = \\frac{2\\,l_{B}}{n\\,\\sqrt{g\\,d}}',
                             'desc': "Merian's formula: resonant period of longitudinal "
                                     'mode n for a closed basin.'},
                            {'tex': 'u_{max} = \\frac{\\pi H}{T}\\,\\frac{\\cosh(k '
                                    'd)}{\\sinh(k d)}',
                             'desc': 'Horizontal velocity amplitude at the node (surface) '
                                     'of the standing wave.'}]},
             {'name': 'Closed 2-D rectangular (Merian extension)',
              'when': 'Closed 2-D (rectangular)',
              'tag': 'standard',
              'note': 'Two-dimensional modes (n,m); setting one index to zero recovers the '
                      '1-D longitudinal mode.',
              'equations': [{'tex': 'T_{n,m} = '
                                    '\\frac{2}{\\sqrt{g\\,d}\\,\\sqrt{\\left(\\frac{n}{l_{B}}\\right)^{2} '
                                    '+ \\left(\\frac{m}{l_{C}}\\right)^{2}}}',
                             'desc': 'Resonant period of the (n,m) mode of a closed '
                                     'rectangular basin.'},
                            {'tex': 'c = \\sqrt{g\\,d}',
                             'desc': 'Shallow-water celerity.'}]},
             {'name': 'Helmholtz pumping mode (basin + channel)',
              'when': 'Helmholtz (basin + channel)',
              'tag': 'standard',
              'note': 'Whole basin rises and falls through the entrance channel; no '
                      'internal nodes. Length correction accounts for added mass at the '
                      'mouth.',
              'equations': [{'tex': 'T_{H} = 2\\pi\\,\\sqrt{\\frac{(L_{ch} + '
                                    'L_{corr})\\,A_{b}}{g\\,A_{c}}}',
                             'desc': 'Helmholtz (pumping) resonance period for a basin '
                                     'connected to the sea by a channel.'}]}],
 'symbols': [['T_n', 'Resonant period of longitudinal mode n (s)'],
             ['T_{n,m}', 'Resonant period of 2-D mode (n,m) (s)'],
             ['T_H', 'Helmholtz pumping-mode period (s)'],
             ['c', 'Shallow-water celerity, sqrt(g d) (m/s)'],
             ['l_B', 'Basin length along the resonant (longitudinal) axis (m)'],
             ['l_C', 'Basin width, transverse dimension, 2-D only (m)'],
             ['d', 'Mean basin water depth (m)'],
             ['n, m', 'Longitudinal and transverse mode numbers'],
             ['H', 'Standing-wave (antinode) crest-to-trough height (m)'],
             ['A_b, A_c',
              'Basin plan area and channel cross-section area; L_ch, L_corr are channel '
              'length and mouth length correction']],
 'references': ["Merian's formula",
                'Helmholtz mode',
                'SPM (1984)',
                'Wilson (1972)',
                'Sorensen (1993)',
                'Ippen (1966)']}


def compute(inp: dict, *, g: float = G_SI, n_modes: int = 6) -> Result:
    """Rectangular-basin resonance and node kinematics for SI inputs."""
    _validate(inp)
    basin_type = str(inp.get("basin_type", "Closed (both ends)"))
    lB = float(inp["lB"]); lC = float(inp["lC"]); d = float(inp["d"]); H = float(inp["H"])
    n = int(inp.get("n", 1)); m = int(inp.get("m_mode", 1))
    c = math.sqrt(g * d)

    T_helm = 0.0
    L_mode = 0.0
    notes = []

    if basin_type.startswith("Helmholtz"):
        Ab = float(inp["Ab"]); Ac = float(inp["Ac"])
        L_ch = float(inp["L_ch"]); L_corr = float(inp["L_corr"])
        T_helm = 2.0 * math.pi * math.sqrt((L_ch + L_corr) * Ab / (g * Ac))
        T_mode = T_helm
        Vmax = excursion = Vbar = 0.0
        mode_idx = np.array([0]); mode_T = np.array([T_helm])
        notes.append("Helmholtz pumping mode: whole basin rises and falls; no internal nodes")
        return Result(c_shallow=c, T_mode=T_mode, L_mode=0.0, f_mode=1.0 / T_mode,
                      T_helm=T_helm, Vmax=Vmax, excursion=excursion, Vbar=Vbar,
                      mode_index=mode_idx, mode_T=mode_T, notes="; ".join(notes))

    if basin_type.startswith("Closed 2-D"):
        if n == 0 and m == 0:
            raise ValueError("2-D mode (n, m) = (0, 0) is the still state (no oscillation)")
        kx = n / lB; ky = m / lC
        T_mode = 2.0 / (c * math.sqrt(kx * kx + ky * ky))
        # 1-D modal wavelength only meaningful when one index is zero
        if m == 0:
            L_mode = 2.0 * lB / n
        elif n == 0:
            L_mode = 2.0 * lC / m
        else:
            L_mode = 0.0
        if m == 0 or n == 0:
            notes.append("one transverse index is zero: reduces to a 1-D longitudinal mode")
    else:
        # 1-D open or closed
        is_open = basin_type.startswith("Open")
        if not is_open and n < 1:
            raise ValueError("closed basin requires longitudinal mode n >= 1 (n = 0 is static)")
        T_mode = _period_1d(basin_type, n, lB, c)
        L_mode = (4.0 * lB / (2 * n + 1)) if is_open else (2.0 * lB / n)

    # standing-wave (clapotis) node kinematics from H and the modal wavelength
    Vmax = excursion = Vbar = 0.0
    if L_mode > 0.0:
        k = 2.0 * math.pi / L_mode
        kd = k * d
        omega = 2.0 * math.pi / T_mode
        # horizontal velocity amplitude at the node, surface (clapotis, SPM 1984):
        #   u_max = (pi H / T) * cosh(kd)/sinh(kd)
        Vmax = (math.pi * H / T_mode) * (math.cosh(kd) / math.sinh(kd))
        excursion = Vmax / omega                       # semi-excursion (amplitude of motion)
        Vbar = (2.0 / math.pi) * Vmax                  # time-mean speed of a sinusoid

    # table of the first several resonant periods (a small seiche "spectrum")
    if basin_type.startswith("Closed 2-D"):
        idx = list(range(1, n_modes + 1))
        Ts = [2.0 / (c * math.sqrt((j / lB) ** 2 + (m / lC) ** 2)) if (j or m) else 0.0
              for j in idx]
    elif basin_type.startswith("Open"):
        idx = list(range(0, n_modes))
        Ts = [_period_1d(basin_type, j, lB, c) for j in idx]
    else:
        idx = list(range(1, n_modes + 1))
        Ts = [_period_1d(basin_type, j, lB, c) for j in idx]
    mode_idx = np.array(idx, dtype=float)
    mode_T = np.array(Ts, dtype=float)

    notes.append(f"fundamental period {mode_T[0]:.1f} s; shallow-water c = {c:.2f} m/s")
    if d / max(L_mode, 1e-9) > 0.05 and L_mode > 0:
        notes.append("note: depth is not small vs modal wavelength; shallow-water c is approximate")

    return Result(c_shallow=c, T_mode=T_mode, L_mode=L_mode, f_mode=1.0 / T_mode,
                  T_helm=0.0, Vmax=Vmax, excursion=excursion, Vbar=Vbar,
                  mode_index=mode_idx, mode_T=mode_T, notes="; ".join(notes))


# --- self-tests (closed forms + reductions) -------------------------------------
def _approx(a: float, b: float, tol: float = 1e-9) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    g = G_SI
    lB, lC, d = 914.4, 609.6, 9.144
    c = math.sqrt(g * d)

    # 1) Open-basin fundamental: T0 = 4 lB / sqrt(g d)
    r = compute({"basin_type": "Open (one end open)", "lB": lB, "lC": lC, "d": d,
                 "H": 0.3, "n": 0, "m_mode": 0}, g=g)
    assert _approx(r.T_mode, 4.0 * lB / c), (r.T_mode, 4.0 * lB / c)

    # 2) Closed-basin fundamental: T1 = 2 lB / sqrt(g d)
    r = compute({"basin_type": "Closed (both ends)", "lB": lB, "lC": lC, "d": d,
                 "H": 0.3, "n": 1, "m_mode": 0}, g=g)
    assert _approx(r.T_mode, 2.0 * lB / c), (r.T_mode, 2.0 * lB / c)

    # 3) 2-D reduces to 1-D closed when m = 0: T(n,0) == 2 lB / (n sqrt(g d))
    r2 = compute({"basin_type": "Closed 2-D (rectangular)", "lB": lB, "lC": lC, "d": d,
                  "H": 0.3, "n": 2, "m_mode": 0}, g=g)
    assert _approx(r2.T_mode, 2.0 * lB / (2 * c)), r2.T_mode

    # 4) 2-D symmetry: square basin -> T(n,m) == T(m,n)
    a = compute({"basin_type": "Closed 2-D (rectangular)", "lB": 1000.0, "lC": 1000.0,
                 "d": d, "H": 0.3, "n": 1, "m_mode": 2}, g=g)
    b = compute({"basin_type": "Closed 2-D (rectangular)", "lB": 1000.0, "lC": 1000.0,
                 "d": d, "H": 0.3, "n": 2, "m_mode": 1}, g=g)
    assert _approx(a.T_mode, b.T_mode), (a.T_mode, b.T_mode)

    # 5) Helmholtz period closed form
    Ab, Ac, Lch, Lcorr = 1.0e6, 6.0e3, 152.4, 30.48
    rh = compute({"basin_type": "Helmholtz (basin + channel)", "lB": lB, "lC": lC, "d": d,
                  "H": 0.3, "n": 1, "m_mode": 0, "Ab": Ab, "Ac": Ac,
                  "L_ch": Lch, "L_corr": Lcorr}, g=g)
    assert _approx(rh.T_helm, 2.0 * math.pi * math.sqrt((Lch + Lcorr) * Ab / (g * Ac))), rh.T_helm

    # 6) node kinematics: shallow-water Vmax -> (H/2) sqrt(g/d); excursion = Vmax/omega
    r = compute({"basin_type": "Closed (both ends)", "lB": lB, "lC": lC, "d": d,
                 "H": 0.3, "n": 1, "m_mode": 0}, g=g)
    omega = 2.0 * math.pi / r.T_mode
    assert _approx(r.excursion, r.Vmax / omega), (r.excursion, r.Vmax / omega)
    assert _approx(r.Vbar, (2.0 / math.pi) * r.Vmax)
    assert _approx(r.Vmax, 0.5 * 0.3 * math.sqrt(g / d), tol=2e-2), r.Vmax  # shallow limit

    print("  self-tests: PASS (open/closed/2-D/Helmholtz closed forms, reductions, kinematics)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  INPUTS (SI):")
    for f in INPUTS:
        vv = inp[f.key]
        sval = f"{vv:>12.5g}" if isinstance(vv, (int, float)) and f.kind != "choice" else f"{vv:>12}"
        print(f"    {f.label:28s} {f.key:11s} = {sval} {f.unit_si}")
    print("  OUTPUTS:")
    by_key = {o.key: o for o in OUTPUTS}
    for kk in ("c_shallow", "T_mode", "L_mode", "f_mode", "Vmax", "excursion", "Vbar"):
        o = by_key[kk]
        print(f"    {o.label:40s} {kk:11s} = {getattr(r, kk):>12.5g} {o.unit_si}")
    print(f"    first {len(r.mode_T)} resonant periods (s): "
          + ", ".join(f"{t:.1f}" for t in r.mode_T))
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
