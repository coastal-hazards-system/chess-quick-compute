"""CHESS-QC application 8-2 — Vessel-Generated Waves.

Originating ACES application: 8-2 "Vessel-Generated Waves" (functional area: Harbor
Design; a later ACES addition). For a design vessel moving at a given speed in a
prismatic channel, it returns the depth Froude number, the celerity/period/propagation
direction of the generated wave system, and the Schijf (1949) one-dimensional
return-current and drawdown (water-level depression alongside the vessel).

Classification: exact (deterministic Schijf 1-D continuity+energy hydraulics + exact
Kelvin/Havelock ship-wave geometry; no empirical correlation; validated analytically to
machine precision).
Theory and references:
  * Return current & drawdown, Schijf (1949) one-dimensional canal theory (PIANC 1987;
    EM 1110-2-1100 Part II): continuity and energy between the undisturbed section and
    the section alongside the vessel,
        continuity:  V_s A_c = (V_s + V_r)(A_c - A_m - b dh)
        energy:      dh = [(V_s + V_r)^2 - V_s^2] / (2g)
    solved simultaneously for the return current V_r and the drawdown dh. A subcritical
    solution exists only below the Schijf limiting speed for the given blockage; above
    it the 1-D theory breaks down (flagged).
  * Wave kinematics, Kelvin / Havelock ship-wave theory: depth Froude number
    F = V_s / sqrt(g d); subcritical (F < 1) wakes carry a transverse wave that travels
    with the vessel (celerity C = V_s; period from the finite-depth dispersion) bounded
    by diverging waves whose crests make 35deg16' (= arctan(1/sqrt2)) with the sailing
    line in deep water; supercritical (F >= 1) wakes carry only diverging waves at the
    Mach crest angle arcsin(1/F), celerity sqrt(g d).
  ACES help manual (ACESManual.rtf, "Vessel-Generated Waves").

Self-containment: zero sibling imports; embeds its own contract dataclasses. Runnable
standalone:
    python chessqc_8_2_vessel_waves.py
which runs the analytic self-tests (S->0 -> no drawdown; Schijf continuity+energy
residual; deep-water subcritical crest angle 35.27deg; supercritical Mach angle) then
prints the default example. stdlib only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

G_SI = 9.80665
_KELVIN_CREST_DEG = math.degrees(math.atan(1.0 / math.sqrt(2.0)))  # 35.2644 deg = 35deg16'
_KT = 0.514444
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
    aces_id="8-2",
    name="Vessel-Generated Waves",
    area="Harbor Design",
    classification="exact",
    cite="Schijf (1949); PIANC (1987); EM 1110-2-1100; Kelvin/Havelock ship-wave theory",
    default_system="SI",
)

INPUTS = (
    Field("b", "Channel width", "float", "m", "ft", default=100.0, lo=1.0, hi=1e4,
          note="prismatic channel top width b > 0"),
    Field("d", "Channel depth", "float", "m", "ft", default=6.0, lo=0.5, hi=200.0,
          note="still-water channel depth d > 0"),
    Field("VS", "Vessel speed", "float", "km/h", "kt", default=3.0, lo=0.01, hi=30.0,
          note="vessel speed through the water V_s > 0"),
    Field("Am", "Wetted cross-sectional area", "float", "m^2", "ft^2", default=25.0,
          lo=0.0, hi=1e4, note="submerged midship section area A_m >= 0"),
)

OUTPUTS = (
    Out("F", "Depth Froude number", "", "", "scalar",
        note="Depth (shallow-water) Froude number F = V_s/sqrt(g d); F<1 is a subcritical wake, F>=1 supercritical."),
    Out("C", "Wave celerity", "m/s", "kt", "scalar",
        note="Phase speed of the generated wave system; equals V_s subcritically (transverse wave moves with the vessel) and sqrt(g d) supercritically."),
    Out("T", "Wave period", "s", "s", "scalar",
        note="Period of the generated wave, T = L/C, from the finite-depth dispersion relation."),
    Out("L", "Wavelength", "m", "ft", "scalar",
        note="Wavelength of the generated wave satisfying C^2 = (g/k) tanh(k d) with L = 2 pi/k."),
    Out("direction", "Wave crest angle to sailing line", "deg", "deg", "scalar",
        note="Angle of the diverging wave crests to the sailing line: deep-water Kelvin angle 35deg16' when F<1, Mach angle arcsin(1/F) when F>=1."),
    Out("drawdown", "Drawdown (water-level depression)", "m", "ft", "scalar",
        note="Schijf drawdown delta h, the depression of the water surface alongside the moving vessel."),
    Out("D", "Relative drawdown dh/d", "", "", "scalar",
        note="Drawdown normalized by channel depth, D = delta h / d (dimensionless)."),
    Out("S", "Blockage ratio A_m/(b d)", "", "", "scalar",
        note="Channel blockage ratio S = A_m/(b d), the fraction of the channel cross-section occupied by the vessel's submerged midship area."),
    Out("Vr", "Return current", "m/s", "kt", "scalar",
        note="Schijf return current V_r flowing backward past the vessel (positive opposes the sailing direction)."),
)


@dataclass
class Result:
    F: float
    C: float
    T: float
    L: float
    direction: float
    drawdown: float
    D: float
    S: float
    Vr: float
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle") or f.key not in inp:
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo:g}, {f.hi:g}] ({f.note})")


def schijf_return_flow(Vs: float, b: float, d: float, Am: float, g: float = G_SI):
    """Solve Schijf (1949) continuity+energy for (V_r, dh, u, over_limit).

    u = (V_s + V_r)/V_s. Returns the subcritical root (smallest drawdown). If the vessel
    exceeds the Schijf limiting speed for the blockage, returns the tangent (limiting)
    drawdown and over_limit=True."""
    Ac = b * d
    S = Am / Ac
    Fh2 = Vs * Vs / (g * d)
    if S < 1e-12 or Vs < 1e-9:
        return 0.0, 0.0, 1.0, False
    umin = 1.0 / (1.0 - S)

    def f(u):  # energy - continuity (root where both drawdowns agree)
        return Fh2 * (u * u - 1.0) / 2.0 - (1.0 - S - 1.0 / u)

    # scan upward from umin; f(umin) >= 0. The subcritical root is the first downward
    # zero-crossing; track the minimum to detect the limiting (no-crossing) case.
    N = 20000
    u_hi = umin + 2.0
    step = (u_hi - umin) / N
    prev = f(umin)
    u_at_min = umin
    f_min = prev
    root = None
    for i in range(1, N + 1):
        u = umin + i * step
        cur = f(u)
        if cur < f_min:
            f_min = cur
            u_at_min = u
        if prev >= 0.0 and cur < 0.0:
            a, c = u - step, u
            for _ in range(100):
                mu = 0.5 * (a + c)
                if f(mu) >= 0.0:
                    a = mu
                else:
                    c = mu
            root = 0.5 * (a + c)
            break
        prev = cur

    if root is not None:
        u = root
        over = False
    else:
        # no subcritical solution -> at/above the limiting speed; report the tangent point
        u = u_at_min
        over = True
    Vr = (u - 1.0) * Vs
    dh = max(1.0 - S - 1.0 / u, 0.0) * d
    return Vr, dh, u, over


def _dispersion_k(C: float, d: float, g: float) -> float:
    """Solve C^2 = (g/k) tanh(k d) for k (rad/m) by fixed-point iteration (C < sqrt(g d))."""
    k = g / (C * C)                       # deep-water guess
    for _ in range(300):
        k_new = g * math.tanh(k * d) / (C * C)
        if abs(k_new - k) <= 1e-12 * k_new:
            k = k_new
            break
        k = k_new
    return k


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'For a design vessel moving at a given speed through a prismatic channel, this '
            'app computes the depth Froude number, the generated wave system (celerity, '
            'period, wavelength and crest-angle to the sailing line) from Kelvin/Havelock '
            'ship-wave theory, and the Schijf (1949) one-dimensional return current and '
            'drawdown (water-level depression) alongside the vessel.',
 'methods': [{'name': 'Schijf 1-D return flow with Kelvin/Havelock wave kinematics',
              'when': None,
              'tag': '',
              'note': 'A subcritical Schijf solution exists only below the limiting speed '
                      'for the given blockage; above it the 1-D theory breaks down and the '
                      'tangent (limiting) drawdown is flagged.',
              'equations': [{'tex': 'F = \\frac{V_s}{\\sqrt{g\\,d}}',
                             'desc': 'Depth Froude number; selects subcritical (F<1) '
                                     'versus supercritical (F>=1) wake regime.'},
                            {'tex': 'V_s\\,A_c = (V_s + V_r)\\,(A_c - A_m - b\\,\\delta h)',
                             'desc': 'Schijf continuity between the undisturbed section '
                                     'and the section alongside the vessel (A_c = b d).'},
                            {'tex': '\\delta h = \\frac{(V_s + V_r)^2 - V_s^2}{2\\,g}',
                             'desc': 'Schijf energy (Bernoulli) relation giving the '
                                     'drawdown; solved simultaneously with continuity for '
                                     'V_r and delta h.'},
                            {'tex': 'C^2 = \\frac{g}{k}\\,\\tanh(k\\,d)',
                             'desc': 'Finite-depth dispersion; for F<1 the transverse wave '
                                     'travels with the vessel (C = V_s) setting period T = '
                                     'L/C.'},
                            {'tex': '\\theta = \\arcsin\\left(\\frac{1}{F}\\right)',
                             'desc': 'Supercritical (F>=1) Mach crest angle; for F<1 the '
                                     'diverging crests sit at the deep-water Kelvin angle '
                                     "arctan(1/sqrt2) = 35deg16'."}]}],
 'symbols': [['F', 'Depth Froude number V_s/sqrt(g d)'],
             ['V_s', 'Vessel speed through the water'],
             ['V_r', 'Return current alongside the vessel'],
             ['delta h', 'Drawdown (water-level depression alongside the vessel)'],
             ['A_c', 'Channel cross-sectional area, b d'],
             ['A_m', 'Submerged midship (wetted) cross-sectional area of the vessel'],
             ['b', 'Prismatic channel top width'],
             ['d', 'Still-water channel depth'],
             ['C', 'Wave celerity'],
             ['theta', 'Wave crest angle to the sailing line (Kelvin or Mach angle)']],
 'references': ['Schijf (1949) one-dimensional canal theory',
                'PIANC (1987)',
                'EM 1110-2-1100 Part II',
                'Kelvin/Havelock ship-wave theory',
                'ACES Manual, Vessel-Generated Waves']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Vessel-generated wave system and Schijf return flow for SI inputs."""
    _validate(inp)
    b = float(inp["b"]); d = float(inp["d"]); Vs = float(inp["VS"]); Am = float(inp["Am"])
    Ac = b * d
    if Am >= Ac:
        raise ValueError("wetted area A_m must be smaller than the channel area b*d")
    S = Am / Ac
    F = Vs / math.sqrt(g * d)

    Vr, dh, u, over = schijf_return_flow(Vs, b, d, Am, g)

    if F < 1.0:
        # subcritical: a transverse wave travels with the vessel (C = V_s); diverging
        # crests at the deep-water Kelvin angle 35deg16'.
        C = Vs
        k = _dispersion_k(C, d, g)
        L = 2.0 * math.pi / k
        T = L / C
        direction = _KELVIN_CREST_DEG
        regime = "subcritical"
    else:
        # supercritical: diverging waves only, crests at the Mach angle arcsin(1/F);
        # the fastest free wave (celerity sqrt(g d)) sets the celerity.
        C = math.sqrt(g * d)
        direction = math.degrees(math.asin(min(1.0 / F, 1.0)))
        L = 2.0 * math.pi * C * C / g          # equivalent free-wave length at this celerity
        T = L / C
        regime = "supercritical"

    notes = [f"{regime} (F={F:.3f}), blockage S={S:.4f}",
             f"V_r={Vr:.3f} m/s, drawdown={dh*100:.1f} cm (D={dh/d:.4f})"]
    if over:
        notes.append("AT/ABOVE Schijf limiting speed: 1-D drawdown is the limiting "
                     "(tangent) value; expect strong squat/wave-making")
    return Result(F=F, C=C, T=T, L=L, direction=direction, drawdown=dh, D=dh / d,
                  S=S, Vr=Vr, notes="; ".join(notes))


# --- self-tests (analytic: Schijf limits + Kelvin/Mach angles) ------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}
    r = compute(base)

    # 1) no blockage -> no return flow and no drawdown
    r0 = compute({**base, "Am": 0.0})
    assert r0.Vr == 0.0 and r0.drawdown == 0.0 and r0.D == 0.0

    # 2) vanishing speed -> drawdown -> 0 (it scales as V_s^2); the Schijf return-flow
    #    ratio V_r/V_s -> S/(1-S) (the blockage-set minimum), so V_r itself -> 0 with V_s
    rslow = compute({**base, "VS": 0.01})
    assert rslow.drawdown < 1e-3, rslow.drawdown
    assert _approx(rslow.Vr / 0.01, rslow.S / (1.0 - rslow.S), 1e-3), rslow.Vr / 0.01

    # 3) Schijf solution satisfies BOTH continuity and energy (residual ~ 0)
    g = G_SI
    b, d, Vs, Am = base["b"], base["d"], base["VS"], base["Am"]
    Vr, dh, u, over = schijf_return_flow(Vs, b, d, Am, g)
    assert not over
    cont = Vs * (b * d) - (Vs + Vr) * (b * d - Am - b * dh)        # = 0
    ener = dh - ((Vs + Vr) ** 2 - Vs ** 2) / (2.0 * g)            # = 0
    assert _approx(cont, 0.0, 1e-6) and _approx(ener, 0.0, 1e-9), (cont, ener)

    # 4) deep-water subcritical diverging crest angle = 35deg16' (arctan(1/sqrt2))
    rdeep = compute({**base, "d": 60.0, "VS": 5.0})
    assert rdeep.F < 1.0 and _approx(rdeep.direction, 35.2644, 1e-3), rdeep.direction
    # transverse wave travels with the vessel: C = V_s
    assert _approx(rdeep.C, 5.0, 1e-9)

    # 5) supercritical: Mach crest angle arcsin(1/F), celerity sqrt(g d)
    rsup = compute({**base, "d": 2.0, "VS": 8.0, "Am": 5.0})
    assert rsup.F > 1.0
    assert _approx(rsup.direction, math.degrees(math.asin(1.0 / rsup.F)), 1e-6)
    assert _approx(rsup.C, math.sqrt(G_SI * 2.0), 1e-9)

    # 6) drawdown increases with vessel speed and with blockage (monotone)
    r_fast = compute({**base, "VS": 4.0})
    r_block = compute({**base, "Am": 40.0})
    assert r_fast.drawdown > r.drawdown and r_block.drawdown > r.drawdown

    print("  self-tests: PASS (S->0 & V->0 no drawdown, Schijf continuity+energy residual ~0, "
          "deep crest 35.26deg, supercritical Mach angle, monotonicity)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  channel b=100 m, d=6 m; vessel V_s=3 m/s ({3.0/_KT:.1f} kt), A_m=25 m^2")
    print("  OUTPUTS:")
    print(f"    Depth Froude number          F     = {r.F:7.3f}")
    print(f"    Wave celerity                C     = {r.C:7.2f} m/s")
    print(f"    Wave period                  T     = {r.T:7.2f} s")
    print(f"    Wavelength                   L     = {r.L:7.2f} m")
    print(f"    Crest angle to sailing line  dir   = {r.direction:7.2f} deg")
    print(f"    Drawdown                     dh    = {r.drawdown*100:7.2f} cm")
    print(f"    Relative drawdown            D     = {r.D:7.4f}")
    print(f"    Blockage ratio               S     = {r.S:7.4f}")
    print(f"    Return current               Vr    = {r.Vr:7.3f} m/s")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
