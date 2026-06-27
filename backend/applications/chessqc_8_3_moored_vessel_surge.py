"""CHESS-QC application 8-3 — Surging of a Moored Vessel.

Originating ACES application: 8-3 "Surging of a Moored Vessel" (functional area: Harbor
Design; a later ACES addition). For a vessel held by a set of mooring lines, it resolves
the line stiffnesses onto the surge (longitudinal) axis, adds the hydrodynamic added
mass, and returns the natural surge period together with the per-line loading and the
forward / reverse / total surge spring constants.

Classification: exact (closed-form linear small-oscillation mooring mechanics -- line
spring rate, cos^2 surge projection, T_S = 2 pi sqrt(m_v/k); no empirical coefficient;
validated analytically to machine precision).
Theory and references: standard moored-ship dynamics (EM 1110-2-1100 Part II; PIANC
mooring guidelines; Bruun "Port Engineering"). Each line is an axial spring whose rate
follows its load-elongation property,
    k_axial = B / ( (e/100) * L )          [breaking strength B at elongation e% of L]
A line at horizontal angle alpha to the surge axis contributes k_axial * cos^2(alpha) to
the surge stiffness (small displacement projected onto the line). Pretensioned lines act
as parallel springs, so the total surge stiffness is the sum over all lines; the forward
and reverse spring constants partition that sum by whether the line's anchor lies forward
(cos alpha >= 0) or aft. The natural surge period follows from the virtual (ship + added)
mass:
    m_v = m (1 + C_a),     T_S = 2 pi sqrt( m_v / k_total ).
ACES help manual (ACESManual.rtf, "Surging of a Moored Vessel").

Self-containment: zero sibling imports; embeds its own contract dataclasses. Runnable
standalone:
    python chessqc_8_3_moored_vessel_surge.py
which runs the analytic self-tests (single-line spring period, symmetric config
fwd == rev, virtual mass, stiffness monotonicity, load and impact flags) then prints the
default example. stdlib only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

G_SI = 9.80665
_TONNE = 1000.0           # kg per metric tonne
_KN = 1000.0              # N per kN


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
    columns: tuple = ()
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
    aces_id="8-3",
    name="Surging of a Moored Vessel",
    area="Harbor Design",
    classification="exact",
    cite="EM 1110-2-1100 Part II; PIANC mooring guidelines; ACES manual",
    default_system="SI",
)

# default: a 5000-tonne vessel on a symmetric 4-line moor (two bow lines at +/-30 deg
# forward, two stern lines at +/-30 deg aft); L=50 m, pretension 100 kN, breaking 1000 kN,
# 3% elongation at break.
_DEF_LINES = [
    [30.0, 50.0, 100.0, 1000.0, 3.0],
    [-30.0, 50.0, 100.0, 1000.0, 3.0],
    [150.0, 50.0, 100.0, 1000.0, 3.0],
    [-150.0, 50.0, 100.0, 1000.0, 3.0],
]

INPUTS = (
    Field("m", "Vessel mass", "float", "tonne", "tonne", default=5000.0, lo=1.0, hi=5e6,
          note="displacement mass m > 0 (stored in kg internally is not used; tonnes here)"),
    Field("Ca", "Surge added-mass coefficient", "float", "", "", default=0.10, lo=0.0, hi=1.0,
          note="virtual mass m_v = m (1 + C_a); surge C_a ~ 0.05-0.25 for ships"),
    Field("lines", "Mooring lines", "table", default=_DEF_LINES,
          columns=(("Angle to surge axis", "deg", "deg"),
                   ("Line length", "m", "m"),
                   ("Pretension", "kN", "kN"),
                   ("Breaking strength", "kN", "kN"),
                   ("Elongation at break", "%", "%")),
          note="one row per line; angle from the +surge (forward) axis to the anchor"),
    Field("swl_fraction", "Safe working load fraction", "float", "", "", default=0.50,
          lo=0.1, hi=1.0, note="line-impact flag trips when load exceeds this x breaking"),
)

OUTPUTS = (
    Out("m_v", "Virtual (surge) mass", "kg", "kg", "scalar",
        note="Virtual surge mass m_v = m (1 + C_a): the vessel displacement mass plus its hydrodynamic added mass."),
    Out("k_fwd", "Forward surge spring constant", "N/m", "N/m", "scalar",
        note="Forward surge stiffness: sum of k_axial cos^2(alpha) over lines whose anchor points forward (cos alpha >= 0)."),
    Out("k_rev", "Reverse surge spring constant", "N/m", "N/m", "scalar",
        note="Reverse surge stiffness: sum of k_axial cos^2(alpha) over aft-pointing lines (cos alpha < 0)."),
    Out("k_total", "Total surge spring constant", "N/m", "N/m", "scalar",
        note="Total surge restoring stiffness k_total = k_fwd + k_rev, the parallel sum of all line stiffnesses projected onto the surge axis."),
    Out("T_S", "Natural surge period", "s", "s", "scalar",
        note="Natural surge oscillation period T_S = 2 pi sqrt(m_v / k_total) of the moored vessel."),
    Out("max_load", "Maximum line load (% breaking)", "%", "%", "scalar",
        note="Largest line pretension expressed as a percent of that line's breaking strength (load = T/B x 100)."),
    Out("n_overloaded", "Lines over safe working load", "", "", "scalar",
        note="Count of lines whose pretension exceeds the safe-working-load fraction times their breaking strength."),
    Out("impact", "Line-impact flag (1 = yes)", "", "", "scalar",
        note="Line-impact / failure-risk flag: 1 if any line exceeds its safe working load, 0 otherwise."),
)


@dataclass
class Result:
    m_v: float
    k_fwd: float
    k_rev: float
    k_total: float
    T_S: float
    max_load: float
    n_overloaded: int
    impact: int
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle") or f.key not in inp:
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo:g}, {f.hi:g}] ({f.note})")


def line_axial_stiffness(B_kN: float, e_pct: float, L: float) -> float:
    """Axial spring rate (N/m) of a line: breaking strength B at elongation e% of length L."""
    if e_pct <= 0.0 or L <= 0.0:
        raise ValueError("line elongation% and length must be positive")
    return (B_kN * _KN) / ((e_pct / 100.0) * L)


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Models a moored vessel as a mass on a set of axial-spring mooring lines: it '
            "resolves each line's stiffness onto the surge axis, adds hydrodynamic added "
            'mass, and returns the natural surge period along with the forward, reverse, '
            'and total surge spring constants and per-line loading.',
 'methods': [{'name': 'Linear moored-vessel surge mechanics',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'k_{axial} = \\frac{B}{(e/100)\\, L}',
                             'desc': 'Axial spring rate of a line with breaking strength B '
                                     'reached at elongation e (percent) of length L.'},
                            {'tex': 'k_{x} = k_{axial} \\cos^{2}\\alpha',
                             'desc': 'Surge-axis stiffness contributed by a line at angle '
                                     'alpha to the forward axis (small-displacement '
                                     'projection).'},
                            {'tex': 'k_{total} = k_{fwd} + k_{rev} = \\sum_{i} k_{axial,i} '
                                    '\\cos^{2}\\alpha_{i}',
                             'desc': 'Total surge stiffness summed over all pretensioned '
                                     'lines (parallel springs); partitioned into forward '
                                     'and reverse parts by the sign of cos alpha.'},
                            {'tex': 'm_{v} = m\\,(1 + C_{a})',
                             'desc': 'Virtual (surge) mass: ship displacement mass plus '
                                     'hydrodynamic added mass.'},
                            {'tex': 'T_{S} = 2\\pi \\sqrt{\\frac{m_{v}}{k_{total}}}',
                             'desc': 'Natural surge period of the moored vessel.'},
                            {'tex': '\\mathrm{load} = \\frac{T}{B} \\cdot 100',
                             'desc': 'Line load as a percent of breaking strength; the '
                                     'impact flag trips when T exceeds the '
                                     'safe-working-load fraction times B.'}]}],
 'symbols': [['B', 'Line breaking strength'],
             ['e', 'Line elongation at break, percent of length'],
             ['L', 'Line length'],
             ['k_{axial}', 'Axial spring rate of a single line'],
             ['alpha', 'Angle from the forward (surge) axis to the line anchor'],
             ['k_{fwd}',
              'Forward surge spring constant (sum over lines with cos alpha >= 0)'],
             ['k_{rev}', 'Reverse surge spring constant (sum over aft lines)'],
             ['k_{total}', 'Total surge spring constant'],
             ['m_{v}', 'Virtual surge mass (ship plus added mass)'],
             ['C_{a}', 'Surge added-mass coefficient (~0.05-0.25 for ships)'],
             ['T_{S}', 'Natural surge period'],
             ['T', 'Line pretension load']],
 'references': ['EM 1110-2-1100 (Coastal Engineering Manual) Part II',
                'PIANC mooring guidelines',
                'Bruun, Port Engineering',
                'ACES Help Manual, Surging of a Moored Vessel']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Moored-vessel surge stiffness, natural period, and line loading (SI inputs)."""
    _validate(inp)
    m = float(inp["m"]) * _TONNE                    # kg
    Ca = float(inp.get("Ca", 0.10))
    swl = float(inp.get("swl_fraction", 0.50))
    rows = [r for r in inp["lines"] if r and len(r) >= 5]
    if not rows:
        raise ValueError("at least one mooring line is required")

    m_v = m * (1.0 + Ca)
    k_fwd = 0.0
    k_rev = 0.0
    max_load = 0.0
    n_over = 0
    for ang, L, T, B, e in ((float(r[0]), float(r[1]), float(r[2]), float(r[3]), float(r[4]))
                            for r in rows):
        if B <= 0.0:
            raise ValueError("line breaking strength must be positive")
        ca = math.cos(math.radians(ang))
        kx = line_axial_stiffness(B, e, L) * ca * ca
        if ca >= 0.0:
            k_fwd += kx
        else:
            k_rev += kx
        load = T / B * 100.0
        max_load = max(max_load, load)
        if T >= swl * B:
            n_over += 1
    k_total = k_fwd + k_rev
    T_S = 2.0 * math.pi * math.sqrt(m_v / k_total) if k_total > 0.0 else float("inf")
    impact = 1 if n_over > 0 else 0

    notes = [f"{len(rows)} lines; m_v={m_v:.3e} kg (C_a={Ca:.2f}); "
             f"k_total={k_total:.3e} N/m -> T_S={T_S:.1f} s",
             f"max line load {max_load:.1f}% of breaking"]
    if impact:
        notes.append(f"WARNING: {n_over} line(s) over {swl*100:.0f}% safe working load "
                     "-> line-impact / failure risk")
    return Result(m_v=m_v, k_fwd=k_fwd, k_rev=k_rev, k_total=k_total, T_S=T_S,
                  max_load=max_load, n_overloaded=n_over, impact=impact,
                  notes="; ".join(notes))


# --- self-tests (analytic mooring mechanics) ------------------------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}
    r = compute(base)

    # 1) virtual mass m_v = m (1 + Ca)
    assert _approx(r.m_v, 5000.0 * _TONNE * 1.10, 1e-9), r.m_v

    # 2) single horizontal line: k_total = k_axial; T_S = 2 pi sqrt(m_v/k_axial)
    single = compute({**base, "lines": [[0.0, 50.0, 100.0, 1000.0, 3.0]]})
    ka = line_axial_stiffness(1000.0, 3.0, 50.0)
    assert _approx(single.k_total, ka, 1e-9), (single.k_total, ka)
    assert _approx(single.T_S, 2.0 * math.pi * math.sqrt(single.m_v / ka), 1e-9)
    assert single.k_rev == 0.0 and single.k_fwd > 0.0     # forward-pointing line only

    # 3) symmetric configuration: forward and reverse spring constants are equal
    assert _approx(r.k_fwd, r.k_rev, 1e-9), (r.k_fwd, r.k_rev)
    assert _approx(r.k_total, r.k_fwd + r.k_rev, 1e-12)

    # 4) a line at angle alpha contributes k_axial cos^2(alpha)
    ang = compute({**base, "lines": [[60.0, 50.0, 100.0, 1000.0, 3.0]]})
    assert _approx(ang.k_total, ka * math.cos(math.radians(60.0)) ** 2, 1e-9)

    # 5) stiffer lines (smaller elongation% or shorter length) shorten the period
    stiff = compute({**base, "lines": [[a, 25.0, T, B, 1.5] for a, L, T, B, e in _DEF_LINES]})
    assert stiff.T_S < r.T_S

    # 6) load (% breaking) and the impact flag
    assert _approx(r.max_load, 100.0 / 1000.0 * 100.0, 1e-9)   # 100 kN of 1000 kN = 10%
    hot = compute({**base, "lines": [[0.0, 50.0, 600.0, 1000.0, 3.0]]})  # 60% > 50% SWL
    assert hot.impact == 1 and hot.n_overloaded == 1 and _approx(hot.max_load, 60.0, 1e-9)

    # 7) heavier vessel -> longer surge period (T_S ~ sqrt(m_v))
    heavy = compute({**base, "m": 20000.0})
    assert _approx(heavy.T_S, r.T_S * math.sqrt(4.0), 1e-9)

    print("  self-tests: PASS (m_v=m(1+Ca), single-line spring period, symmetric fwd==rev, "
          "cos^2 projection, stiffness/mass monotonicity, load & impact flags)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  5000-tonne vessel, 4-line symmetric moor (+/-30 deg bow & stern), L=50 m, "
          f"B=1000 kN, e=3%, pretension 100 kN, C_a=0.10")
    print("  OUTPUTS:")
    print(f"    Virtual mass                 m_v   = {r.m_v:10.3e} kg")
    print(f"    Forward spring constant      k_fwd = {r.k_fwd:10.3e} N/m")
    print(f"    Reverse spring constant      k_rev = {r.k_rev:10.3e} N/m")
    print(f"    Total spring constant        k_tot = {r.k_total:10.3e} N/m")
    print(f"    Natural surge period         T_S   = {r.T_S:10.2f} s")
    print(f"    Maximum line load            load  = {r.max_load:10.1f} % of breaking")
    print(f"    Lines over safe working load       = {r.n_overloaded:10d}")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
