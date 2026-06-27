"""CHESS-QC application 4-2 — Toe Protection Design.

Designs the toe-apron width and toe-stone weight for a vertical wall / bulkhead /
revetment. Apron width is the larger of a geotechnical (Rankine passive) width and
two hydraulic minima; toe-stone weight uses the Tanimoto, Yagyu & Goda (1982)
stability number N_s in a Hudson-form sizing equation.

Classification: exact (Tanimoto, Yagyu & Goda 1982 toe-stability number + EM 1110-2-1614
apron-width rules; all empirical coefficients known from the sources, nothing guessed;
reproduces the User's Guide Example 4-2 #1).

Self-contained (zero sibling imports): embeds the AppMeta/Field/Out/Result
dataclasses and the Hunt (1979) dispersion solver. No wave profile -> value rows only.

Theory: TR 4-2 (eqs 1-6); EM 1110-2-1614; Tanimoto, Yagyu & Goda (1982). Validated
against the ACES User's Guide Example 4-2 #1 (see tests/test_manual_oracle.py).

Run:
    python chessqc_4_2_toe_protection.py    # self-tests + tabulate the manual example
"""
from __future__ import annotations

import math
from dataclasses import dataclass

# --- constants ------------------------------------------------------------------
G_SI = 9.80665
_FT = 0.3048               # ft -> m
_LBF_PER_FT3 = 157.08746   # lb/ft^3 -> N/m^3
_LBF = 4.4482216           # lbf -> N


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
    aces_id="4-2",
    name="Toe Protection Design",
    area="Structural Design",
    classification="exact",
    cite="EM 1110-2-1614; Tanimoto, Yagyu & Goda (1982); Hunt (1979); TR 4-2",
    default_system="US",
)

# Complete input list (ACES User's Guide 4-2). Defaults = User's Guide Example 4-2 #1.
INPUTS = (
    Field("H_i", "Incident wave height", "float", "m", "ft", default=5.00 * _FT, lo=1e-6, hi=1e4,
          note="> 0 (incident/design wave height at the structure)"),
    Field("T", "Wave period", "float", "s", "s", default=12.0, lo=1e-3, hi=1e4, note="> 0"),
    Field("d_s", "Water depth at structure", "float", "m", "ft", default=20.00 * _FT, lo=1e-6, hi=1e5,
          note="> 0"),
    Field("cot_phi", "Cotangent of nearshore slope", "float", "", "", default=100.0, lo=1e-6, hi=1e6,
          note="cot(beach slope); collected for context (not used in toe sizing)"),
    Field("K_p", "Passive earth-pressure coefficient", "float", "", "", default=1.50, lo=0.0, hi=100.0,
          note="Rankine passive coefficient K_p (0 if no geotechnical width)"),
    Field("d_e", "Sheet-pile penetration depth", "float", "m", "ft", default=10.00 * _FT, lo=0.0, hi=1e4,
          note=">= 0 (0 if no sheet pile)"),
    Field("h_b", "Height of toe layer above mudline", "float", "m", "ft", default=4.50 * _FT, lo=0.0, hi=1e4,
          note="0 <= h_b < d_s; d1 = d_s - h_b"),
    Field("w_r", "Unit weight of rock", "float", "kN/m^3", "lb/ft^3", default=165.0 * _LBF_PER_FT3,
          lo=1.0, hi=1e6, note="unit weight of toe-stone; must exceed water unit weight"),
    Field("w_w", "Water unit weight", "float", "kN/m^3", "lb/ft^3", default=64.0 * _LBF_PER_FT3,
          lo=1.0, hi=1e6, note="64 lb/ft^3 seawater, 62.4 lb/ft^3 fresh"),
)

# Complete output list (ACES User's Guide 4-2) + N_s (key design parameter).
OUTPUTS = (
    Out("B",   "Width of toe protection apron",          "m", "ft", "scalar",
        note="Design width of the toe-protection apron, the larger of the Rankine geotechnical width K_p*d_e and the two hydraulic minima 2*H_i and 0.4*d_s."),
    Out("W",   "Weight of individual armor unit",        "N", "lb", "scalar",
        note="Weight of an individual toe-stone armor unit from the Hudson-form sizing equation W = w_r*H_i^3 / (N_s^3*(S_r-1)^3)."),
    Out("d_1", "Water depth at top of toe layer",        "m", "ft", "scalar",
        note="Water depth at the top of the toe layer, d_1 = d_s - h_b (still-water depth above the toe stone)."),
    Out("N_s", "Stability number (Tanimoto-Yagyu-Goda)", "",  "",   "scalar",
        note="Tanimoto-Yagyu-Goda (1982) toe-stability number, floored at 1.8; larger N_s yields lighter required toe stone."),
)


@dataclass
class Result:
    B: float; W: float; d_1: float; N_s: float
    notes: str = ""


# --- Hunt (1979) dispersion (identical kernel to 2-1/3-1) -----------------------
_HUNT_D = (0.66667, 0.35550, 0.16084, 0.06320, 0.02174,
           0.00654, 0.00171, 0.00039, 0.00011)


def wave_celerity(T: float, d: float, g: float = G_SI) -> float:
    omega = 2.0 * math.pi / T
    y = omega * omega * d / g
    denom = 1.0 + sum(dn * y ** (n + 1) for n, dn in enumerate(_HUNT_D))
    return math.sqrt(g * d / (y + 1.0 / denom))


def _validate(inp: dict) -> None:
    for f in INPUTS:
        v = inp[f.key]
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")
    if inp["h_b"] >= inp["d_s"]:
        raise ValueError("Toe-layer height h_b must be less than water depth d_s")
    if inp["w_r"] <= inp["w_w"]:
        raise ValueError("Rock unit weight must exceed water unit weight (S_r > 1)")


# --- compute --------------------------------------------------------------------
# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {'summary': 'Sizes the toe-protection apron for a vertical wall, bulkhead, or revetment: '
            'the apron width is the larger of a geotechnical (Rankine passive) width and '
            'two hydraulic minima, and the toe-stone weight follows a Hudson-form sizing '
            'equation using the Tanimoto-Yagyu-Goda (1982) stability number.',
 'methods': [{'name': 'Toe apron width and toe-stone weight (Tanimoto-Yagyu-Goda)',
              'when': None,
              'tag': '',
              'note': None,
              'equations': [{'tex': 'B = \\max(K_p\\,d_e,\\, 2\\,H_i,\\, 0.4\\,d_s)',
                             'desc': 'Design apron width: largest of the Rankine '
                                     'geotechnical width and the two EM 1110-2-1614 '
                                     'hydraulic minima (eqs 1-4).'},
                            {'tex': 'W = \\frac{w_r\\,H_i^{3}}{N_s^{3}\\,(S_r - 1)^{3}}',
                             'desc': 'Toe-stone weight, Hudson/SPM (1984) form, with S_r = '
                                     'w_r/w_w (eq 5).'},
                            {'tex': 'N_s = \\max\\left[ '
                                    '1.3\\,\\frac{1-K}{K^{1/3}}\\,\\frac{d_1}{H_i} + '
                                    '1.8\\,\\exp\\left(-1.5\\,\\frac{(1-K)^{2}}{K^{1/3}}\\,\\frac{d_1}{H_i}\\right),\\, '
                                    '1.8 \\right]',
                             'desc': 'Tanimoto-Yagyu-Goda (1982) stability number, floored '
                                     'at 1.8 (eq 6).'},
                            {'tex': 'K = '
                                    '\\frac{4\\pi\\,d_1/L}{\\sinh(4\\pi\\,d_1/L)}\\,\\sin^{2}\\!\\left(\\frac{2\\pi\\,B}{L}\\right)',
                             'desc': 'Velocity parameter at the apron edge from '
                                     'standing-wave linear theory; L is the wavelength at '
                                     'depth d_1 (eq 6).'}]}],
 'symbols': [['B', 'Design width of the toe-protection apron'],
             ['W', 'Weight of an individual toe-stone armor unit'],
             ['N_s', 'Tanimoto-Yagyu-Goda stability number (floored at 1.8)'],
             ['K', 'Velocity parameter at the apron edge from standing-wave theory'],
             ['K_p', 'Rankine passive earth-pressure coefficient'],
             ['d_e', 'Sheet-pile penetration depth (0 if no pile)'],
             ['H_i', 'Incident (design) wave height at the structure'],
             ['d_s', 'Water depth at the structure'],
             ['d_1', 'Water depth at top of the toe layer (d_s - h_b)'],
             ['S_r', 'Specific gravity of armor stone, w_r/w_w']],
 'references': ['EM 1110-2-1614',
                'Tanimoto, Yagyu & Goda (1982)',
                'Hudson (1959)',
                'SPM (1984)',
                'Eckert (1983)',
                'Eckert & Callender (1987)',
                'Hunt (1979)',
                'ACES TR 4-2']}


def compute(inp: dict, *, g: float = G_SI) -> Result:
    """Toe-protection design for SI inputs
    {H_i, T, d_s, cot_phi, K_p, d_e, h_b, w_r[N/m^3], w_w[N/m^3]}."""
    _validate(inp)
    H_i = float(inp["H_i"]); T = float(inp["T"]); d_s = float(inp["d_s"])
    K_p = float(inp["K_p"]); d_e = float(inp["d_e"]); h_b = float(inp["h_b"])
    w_r = float(inp["w_r"]); w_w = float(inp["w_w"])

    d_1 = d_s - h_b                                   # depth at top of toe layer
    # (1)-(4) apron width = max(geotechnical, two hydraulic minima)
    B1 = K_p * d_e
    B2 = 2.0 * H_i
    B3 = 0.4 * d_s
    B = max(B1, B2, B3)

    # wavelength at depth d_1 (linear theory)
    C = wave_celerity(T, d_1, g)
    L = C * T
    # (6) Tanimoto-Yagyu-Goda stability number
    x = 4.0 * math.pi * d_1 / L
    K = (x / math.sinh(x)) * math.sin(2.0 * math.pi * B / L) ** 2
    ratio = d_1 / H_i
    k13 = K ** (1.0 / 3.0)
    N_s = 1.3 * ((1.0 - K) / k13) * ratio + 1.8 * math.exp(-1.5 * ((1.0 - K) ** 2 / k13) * ratio)
    N_s = max(N_s, 1.8)

    # (5) toe-stone weight (W in N: w_r in N/m^3, H_i in m)
    S_r = w_r / w_w
    W = w_r * H_i ** 3 / (N_s ** 3 * (S_r - 1.0) ** 3)

    which = "geotechnical K_p*d_e" if B == B1 else ("hydraulic 2*H_i" if B == B2 else "hydraulic 0.4*d_s")
    notes = f"B governed by {which}; S_r={S_r:.3f}; N_s={N_s:.2f}"
    return Result(B=B, W=W, d_1=d_1, N_s=N_s, notes=notes)


# --- self-tests + manual-example tabulation -------------------------------------
def _self_tests() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    assert abs(r.B / _FT - 15.00) < 0.05, r.B / _FT
    assert abs(r.W / _LBF - 12.99) < 0.10, r.W / _LBF
    assert abs(r.d_1 / _FT - 15.50) < 0.05, r.d_1 / _FT
    # sanity: bigger waves -> heavier toe stone
    r2 = compute({**{f.key: f.default for f in INPUTS}, "H_i": 8.0 * _FT})
    assert r2.W > r.W
    print("  self-tests: PASS (matches User's Guide Example 4-2 #1)")


def _tab() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (values in US units; matches User's Guide Example 4-2 #1)")
    print(f"  Apron width B          = {r.B / _FT:8.2f} ft   [manual 15.00]")
    print(f"  Toe-stone weight W     = {r.W / _LBF:8.2f} lb   [manual 12.99]")
    print(f"  Depth at toe top d_1   = {r.d_1 / _FT:8.2f} ft   [manual 15.50]")
    print(f"  Stability number N_s   = {r.N_s:8.2f}")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _tab()
