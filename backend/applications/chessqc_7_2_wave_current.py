"""CHESS-QC application 7-2 — Wave-Current Interaction in Channels.

Originating ACES application: 7-2 "Wave-current Interaction in Channels" (functional
area: Inlet Processes; a later ACES addition). Computes how a wave train is modified
when it crosses a navigation channel carrying a steady current: the current Doppler-
shifts the dispersion relation (changing the wavelength) and changes the wave height
through conservation of wave action.

Classification: exact (deterministic linear wave-current theory).
Theory and references: linear wave-current interaction (Jonsson 1990; Jonsson, Skovgaard
& Wang 1970; Peregrine 1976). The absolute (ground-frame) angular frequency omega is
conserved across the steady current; the intrinsic frequency sigma satisfies the
still-water dispersion relation, and the two are linked by the Doppler shift along the
wave-propagation direction:
    omega = sigma + k U,   sigma = sqrt(g k tanh(k d)),   U = V cos(alpha)
with V the current speed, alpha the angle between the wave orthogonal and the current.
The wavelength with current follows from solving this for k. The wave height follows
from conservation of wave-action flux (E/sigma)(Cg + U):
    R_H = H/H0 = sqrt[ (sigma/omega) * Cg0 / (Cg + U) ]
where Cg is the intrinsic group velocity and the subscript 0 denotes the no-current
state. An opposing current blocks the waves when Cg + U -> 0 (effective Froude number
F = U/Cg -> -1); this is detected and flagged. ACES help manual (ACESManual.rtf,
"Wave-current Interaction in Channels").

Self-containment: zero sibling imports; embeds its own contract dataclasses. Runnable
standalone:
    python chessqc_7_2_wave_current.py
which runs the analytic self-tests (no current -> R_H = R_L = 1; following current
lengthens and lowers, opposing shortens and steepens; blocking flag) then prints the
default example. stdlib only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

G_SI = 9.80665


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
    aces_id="7-2",
    name="Wave-Current Interaction in Channels",
    area="Inlet Processes",
    classification="exact",
    cite="Jonsson (1990); Jonsson, Skovgaard & Wang (1970); Peregrine (1976)",
    default_system="SI",
)

INPUTS = (
    Field("T", "Wave period", "float", "s", "s", default=8.0, lo=0.5, hi=30.0,
          note="absolute wave period (ground frame) T > 0"),
    Field("alpha", "Angle of wave orthogonal to current", "angle", "deg", "deg",
          default=0.0, lo=0.0, hi=180.0,
          note="0 = following current, 180 = directly opposing"),
    Field("V", "Channel current velocity", "float", "m/s", "kt", default=1.5, lo=-10.0, hi=10.0,
          note="current speed; combined with the angle (positive magnitude)"),
    Field("dT", "Channel depth", "float", "m", "ft", default=10.0, lo=0.3, hi=300.0,
          note="still-water channel depth d > 0"),
)

OUTPUTS = (
    Out("omega_star", "Dimensionless frequency (omega sqrt(d/g))", "", "", "scalar",
        note="Nondimensional absolute wave frequency, omega sqrt(d/g), with omega = 2 pi / T; characterizes the relative depth."),
    Out("F", "Effective Froude number U/Cg", "", "", "scalar",
        note="Effective Froude number U/Cg comparing the along-wave current to the intrinsic group velocity; approaches -1 as an opposing current blocks the waves."),
    Out("R_H", "Wave height factor H/H0", "", "", "scalar",
        note="Wave-height modification factor H/H0 from action conservation; < 1 for a following current (lowering), > 1 for an opposing current (steepening)."),
    Out("R_L", "Wavelength factor L/L0", "", "", "scalar",
        note="Wavelength modification factor L/L0; > 1 for a following current (lengthening), < 1 for an opposing current (shortening)."),
    Out("L0", "Wavelength (no current)", "m", "ft", "scalar",
        note="Wavelength in still water (no current) from the intrinsic dispersion relation at depth d."),
    Out("L", "Wavelength (with current)", "m", "ft", "scalar",
        note="Wavelength in the presence of the current, from the Doppler-shifted dispersion relation omega = sigma + k U."),
)


@dataclass
class Result:
    omega_star: float
    F: float
    R_H: float
    R_L: float
    L0: float
    L: float
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle") or f.key not in inp:
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo:g}, {f.hi:g}] ({f.note})")


def _sigma(k: float, d: float, g: float) -> float:
    """Intrinsic angular frequency sqrt(g k tanh(k d))."""
    return math.sqrt(g * k * math.tanh(k * d))


def _cg_intrinsic(k: float, d: float, g: float) -> float:
    """Intrinsic group velocity dsigma/dk = 0.5 (sigma/k)(1 + 2kd/sinh(2kd))."""
    s = _sigma(k, d, g)
    return 0.5 * (s / k) * (1.0 + 2.0 * k * d / math.sinh(2.0 * k * d))


def _k_no_current(omega: float, d: float, g: float) -> float:
    """Solve sigma(k) = omega for k (no current) by fixed-point iteration."""
    k = omega * omega / g
    for _ in range(500):
        k_new = omega * omega / (g * math.tanh(k * d))
        if abs(k_new - k) <= 1e-14 * k_new:
            k = k_new
            break
        k = k_new
    return k


def _k_with_current(omega: float, d: float, U: float, g: float):
    """Solve omega = sigma(k) + k U for k. Returns (k, blocked).

    The physical root is the one continuously connected to the no-current wavenumber.
    For an opposing current (U < 0) the wave is blocked when no such root exists (the
    maximum of sigma(k)+kU falls below omega, i.e. Cg + U -> 0)."""
    def f(k):
        return _sigma(k, d, g) + k * U - omega

    if U >= 0.0:
        lo, hi = 1e-7, 1.0
        while f(hi) < 0.0 and hi < 1e6:
            hi *= 2.0
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            if f(mid) < 0.0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi), False

    # opposing current: locate the maximum of f (where Cg_intrinsic + U = 0)
    lo, hi = 1e-7, 1.0
    while _cg_intrinsic(hi, d, g) + U > 0.0 and hi < 1e6:
        hi *= 2.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if _cg_intrinsic(mid, d, g) + U > 0.0:
            lo = mid
        else:
            hi = mid
    k_star = 0.5 * (lo + hi)
    if f(k_star) < 0.0:
        return None, True                       # blocking: no real wave
    a, b = 1e-7, k_star                          # smaller root (connected to k0)
    for _ in range(300):
        mid = 0.5 * (a + b)
        if f(mid) < 0.0:
            a = mid
        else:
            b = mid
    return 0.5 * (a + b), False


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Computes how a wave train is modified when it crosses a navigation channel '
            'carrying a steady current. Returns the wavelength and wave-height '
            'modification factors (relative to the no-current state), the effective Froude '
            'number, and flags opposing-current wave blocking.',
 'methods': [{'name': 'Linear wave-current interaction (action conservation)',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': '\\omega = \\sigma + k U',
                             'desc': 'Absolute (ground-frame) frequency is conserved; '
                                     'intrinsic frequency sigma is Doppler-shifted by the '
                                     'current component k U.'},
                            {'tex': '\\sigma = \\sqrt{g k \\tanh(k d)}',
                             'desc': 'Intrinsic (still-water) linear dispersion relation, '
                                     'solved for the wavenumber k with and without '
                                     'current.'},
                            {'tex': 'U = V \\cos(\\alpha)',
                             'desc': 'Current component along the wave-propagation '
                                     'direction (alpha measured from the wave orthogonal '
                                     'to the current).'},
                            {'tex': 'C_{g} = \\frac{\\sigma}{2 k}\\left(1 + \\frac{2 k '
                                    'd}{\\sinh(2 k d)}\\right)',
                             'desc': 'Intrinsic group velocity d sigma / dk used in the '
                                     'wave-action flux.'},
                            {'tex': 'R_{H} = \\frac{H}{H_{0}} = '
                                    '\\sqrt{\\frac{\\sigma}{\\omega}\\,\\frac{C_{g0}}{C_{g} '
                                    '+ U}}',
                             'desc': 'Wave-height factor from conservation of wave-action '
                                     'flux (E/sigma)(Cg + U); subscript 0 is the '
                                     'no-current state.'},
                            {'tex': 'F = \\frac{U}{C_{g}}',
                             'desc': 'Effective Froude number; an opposing current blocks '
                                     'the waves as C_g + U -> 0 (F -> -1).'}]}],
 'symbols': [['omega',
              'Absolute (ground-frame) angular frequency, 2 pi / T, conserved across the '
              'current'],
             ['sigma',
              'Intrinsic (relative-to-current) angular frequency satisfying still-water '
              'dispersion'],
             ['k', 'Wavenumber, 2 pi / L (solved with and without current)'],
             ['U', 'Current component along wave propagation, V cos(alpha)'],
             ['V', 'Channel current speed (magnitude)'],
             ['alpha',
              'Angle between the wave orthogonal and the current (0 following, 180 '
              'opposing)'],
             ['d', 'Still-water channel depth'],
             ['C_g', 'Intrinsic group velocity; C_{g0} is the no-current value'],
             ['R_H', 'Wave-height modification factor H/H0'],
             ['F', 'Effective Froude number U/Cg (-> -1 at blocking)']],
 'references': ['Jonsson (1990)',
                'Jonsson, Skovgaard & Wang (1970)',
                'Peregrine (1976)',
                "ACES Technical Reference / User's Guide, App. 7-2 Wave-Current "
                'Interaction in Channels']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Wave-current interaction across a channel for SI inputs."""
    _validate(inp)
    T = float(inp["T"]); alpha = math.radians(float(inp["alpha"]))
    V = float(inp["V"]); d = float(inp["dT"])
    omega = 2.0 * math.pi / T
    U = V * math.cos(alpha)                       # current component along wave propagation

    k0 = _k_no_current(omega, d, g)
    L0 = 2.0 * math.pi / k0
    Cg0 = _cg_intrinsic(k0, d, g)
    omega_star = omega * math.sqrt(d / g)

    k, blocked = _k_with_current(omega, d, U, g)
    if blocked:
        notes = (f"OPPOSING-CURRENT WAVE BLOCKING (U={U:.2f} m/s; the current speed reaches "
                 f"the group velocity, Cg0={Cg0:.2f} m/s): no steady wave crosses the channel")
        return Result(omega_star=omega_star, F=-1.0, R_H=float("inf"), R_L=float("inf"),
                      L0=L0, L=float("inf"), notes=notes)

    L = 2.0 * math.pi / k
    sigma = _sigma(k, d, g)
    Cg = _cg_intrinsic(k, d, g)
    R_L = L / L0
    R_H = math.sqrt((sigma / omega) * Cg0 / (Cg + U))
    F = U / Cg                                    # effective Froude number; -> -1 at blocking

    sense = "no current" if abs(U) < 1e-9 else ("following" if U > 0 else "opposing")
    notes = [f"{sense} current U={U:.2f} m/s; L {L0:.1f}->{L:.1f} m (R_L={R_L:.3f}), "
             f"R_H={R_H:.3f}, F={F:.3f}"]
    if F < -0.85:
        notes.append("near opposing-current blocking (F -> -1)")
    return Result(omega_star=omega_star, F=F, R_H=R_H, R_L=R_L, L0=L0, L=L,
                  notes="; ".join(notes))


# --- self-tests (analytic wave-current limits) ----------------------------------
def _approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def _self_tests() -> None:
    base = {f.key: f.default for f in INPUTS}

    # 1) no current -> R_H = R_L = 1, L = L0, F = 0
    r0 = compute({**base, "V": 0.0})
    assert _approx(r0.R_H, 1.0, 1e-9) and _approx(r0.R_L, 1.0, 1e-9), (r0.R_H, r0.R_L)
    assert _approx(r0.L, r0.L0, 1e-9) and _approx(r0.F, 0.0, 1e-12)

    # 2) following current (alpha=0, V>0): waves lengthen and lower (R_L>1, R_H<1)
    rf = compute({**base, "alpha": 0.0, "V": 1.5})
    assert rf.R_L > 1.0 and rf.R_H < 1.0, (rf.R_L, rf.R_H)
    assert rf.F > 0.0

    # 3) opposing current (alpha=180): waves shorten and steepen (R_L<1, R_H>1)
    ro = compute({**base, "alpha": 180.0, "V": 1.5})
    assert ro.R_L < 1.0 and ro.R_H > 1.0, (ro.R_L, ro.R_H)
    assert ro.F < 0.0

    # 4) the with-current dispersion is satisfied: omega = sigma(k) + k U
    g = G_SI
    omega = 2.0 * math.pi / base["T"]
    U = 1.5 * math.cos(0.0)
    k, blk = _k_with_current(omega, base["dT"], U, g)
    assert not blk and _approx(_sigma(k, base["dT"], g) + k * U, omega, 1e-9)

    # 5) action conservation identity: R_H^2 (Cg+U)/sigma == Cg0/omega
    r = compute({**base, "alpha": 0.0, "V": 1.5})
    k0 = _k_no_current(omega, base["dT"], g); Cg0 = _cg_intrinsic(k0, base["dT"], g)
    k, _ = _k_with_current(omega, base["dT"], U, g)
    Cg = _cg_intrinsic(k, base["dT"], g); sig = _sigma(k, base["dT"], g)
    assert _approx(r.R_H ** 2 * (Cg + U) / sig, Cg0 / omega, 1e-9)

    # 6) 90-degree crossing (alpha=90): U=0, no modification
    r90 = compute({**base, "alpha": 90.0, "V": 5.0})
    assert _approx(r90.R_H, 1.0, 1e-9) and _approx(r90.R_L, 1.0, 1e-9)

    # 7) strong opposing current -> blocking flagged (R_H -> inf, F = -1)
    rb = compute({**base, "alpha": 180.0, "V": 6.0})
    assert not math.isfinite(rb.R_H) and rb.F == -1.0, (rb.R_H, rb.F)

    print("  self-tests: PASS (no current R_H=R_L=1, following lengthens/lowers, opposing "
          "shortens/steepens, dispersion & action identities, 90-deg no-op, blocking flag)")


def _print_default_example() -> None:
    inp = {f.key: f.default for f in INPUTS}
    r = compute(inp)
    print(f"\nCHESS-QC {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print(f"  T=8 s, alpha=0 deg (following), V=1.5 m/s, d=10 m")
    print("  OUTPUTS:")
    print(f"    Dimensionless frequency      w*    = {r.omega_star:7.4f}")
    print(f"    Effective Froude number      F     = {r.F:7.4f}")
    print(f"    Wave height factor           R_H   = {r.R_H:7.4f}")
    print(f"    Wavelength factor            R_L   = {r.R_L:7.4f}")
    print(f"    Wavelength (no current)      L0    = {r.L0:7.2f} m")
    print(f"    Wavelength (with current)    L     = {r.L:7.2f} m")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
