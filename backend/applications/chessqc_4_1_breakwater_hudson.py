"""CHESS-QC application 4-1 — Breakwater Design Using Hudson and Related Equations.

Sizes the primary armor units of a rubble-mound breakwater/revetment from the Hudson
(1953-61) stability equation, and reports crest width, cover-layer thickness, and
armor-unit placement density (SPM 1984 Ch. 7; EM 1110-2-2904).

Classification: exact (closed-form Hudson stability equation and the related SPM crest-
width / thickness / placement-density formulas; reproduces the User's Guide Example 4-1 to
the digit. K_D, k_delta and porosity are user-supplied table coefficients, as in ACES.)
A selectable Van der Meer (1988) rock-armor stability method is also provided (the modern
CIRIA Rock Manual standard, accounting for wave period, storm duration N, notional
permeability P and damage level S). Hudson stays the default so the Example 4-1
reproduction is preserved; the Van der Meer path is validated for internal consistency and
physical limits (plunging/surging branch, stability-number range), as ACES has no oracle
for it.

Self-contained (zero sibling imports): embeds the AppMeta/Field/Out/Result dataclasses.
No wave kinematics -> value-rows-only (no profile plot).

Theory: TR 4-1 (eqs 1-4). Validated against the ACES User's Guide Example 4-1
(see tests/test_manual_oracle.py).

Run:
    python chessqc_4_1_breakwater_hudson.py    # self-tests + tabulate the manual example
"""
from __future__ import annotations

import math
from dataclasses import dataclass

# --- constants ------------------------------------------------------------------
_FT = 0.3048               # ft -> m
_LBF_PER_FT3 = 157.08746   # lb/ft^3 -> N/m^3
_LBF = 4.4482216           # lbf -> N
_TON = 8896.4432           # US short ton-force (2000 lbf) -> N
G = 9.80665                # m/s^2 (Van der Meer wave-steepness)


# --- contract dataclasses -------------------------------------------------------
@dataclass(frozen=True)
class AppMeta:
    aces_id: str
    name: str
    area: str
    classification: str
    cite: str
    default_system: str = "SI"
    status: str = "Current"          # Current | Screening only | Superseded
    superseded_by: str = ""          # newer method, if any (surfaced in the docs)


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
    enable_if: tuple = ()    # (other_key, value): gray out (disable) unless that input == value


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    note: str = ""           # hover definition shown on the output label


# armor-unit types (informational; help the user choose K_D, k_delta, P from SPM tables)
_ARMOR_TYPES = (
    "Quarrystone (smooth, rounded)", "Quarrystone (rough, angular)", "Graded riprap",
    "Tribar (trunk, nonbreaking)", "Tribar (trunk, breaking)", "Tetrapod", "Quadripod",
    "Dolos", "Modified cube", "Hexapod", "Toskane", "Other",
)

APP_META = AppMeta(
    aces_id="4-1",
    name="Breakwater Design (Hudson)",
    area="Structural Design",
    classification="exact",
    cite="Hudson (1953-61); SPM (1984) Ch.7; EM 1110-2-2904; TR 4-1",
    default_system="US",     # opens on the User's Guide Example (US units)
    superseded_by="Van der Meer (1988) stability formulae (preferred for many cases)",
)

# Complete input list (ACES User's Guide 4-1). Defaults = User's Guide Example 4-1.
INPUTS = (
    Field("armor_type", "Type of armor unit", "choice", "", "",
          default="Tribar (trunk, nonbreaking)", choices=_ARMOR_TYPES,
          note="optional/informational; pick K_D, k_delta, P from SPM tables accordingly"),
    Field("method", "Sizing method", "choice", "", "", default="Hudson",
          choices=("Hudson", "Van der Meer"),
          note="Hudson (ACES; reproduces Example 4-1) or Van der Meer (1988) rock-armor stability"),
    Field("w_r", "Armor unit weight", "float", "kN/m^3", "lb/ft^3", default=165.0 * _LBF_PER_FT3,
          lo=1.0, hi=1e6, note="unit weight of armor material; must exceed water unit weight"),
    Field("H", "Wave height", "float", "m", "ft", default=11.50 * _FT, lo=1e-6, hi=1e4,
          note="design wave height (H or H_i)"),
    Field("w_w", "Water unit weight", "float", "kN/m^3", "lb/ft^3", default=64.0 * _LBF_PER_FT3,
          lo=1.0, hi=1e6, note="64 lb/ft^3 seawater, 62.4 lb/ft^3 fresh"),
    Field("K_D", "Stability coefficient", "float", "", "", default=10.0, lo=1e-3, hi=1e4,
          note="K_D from SPM Table 7-8 (depends on armor type / slope / wave condition)",
          enable_if=("method", "Hudson")),
    Field("k_delta", "Layer coefficient", "float", "", "", default=1.02, lo=1e-3, hi=10.0,
          note="layer coefficient k_delta (SPM Table 7-13)"),
    Field("P", "Average porosity of armor layer", "float", "%", "%", default=54.0, lo=0.0, hi=99.0,
          note="cover-layer porosity, percent (SPM Table 7-13)"),
    Field("cot_theta", "Cotangent of structure slope", "float", "", "", default=2.00, lo=1e-3, hi=1e3,
          note="cot(theta); theta = seaward slope angle"),
    Field("n", "Number of armor units (layer thickness)", "int", "", "", default=2, lo=1, hi=10,
          note="number of armor-unit layers (>= 2 typical)"),
    # --- Van der Meer (1988) parameters (used only when method = Van der Meer) ---
    Field("Tm", "Mean wave period (Van der Meer)", "float", "s", "s", default=8.0, lo=1e-2, hi=1e3,
          note="Van der Meer only: mean period T_m for the surf-similarity parameter",
          enable_if=("method", "Van der Meer")),
    Field("N_waves", "Number of waves (Van der Meer)", "int", "", "", default=7500, lo=1, hi=1000000,
          note="Van der Meer only: storm duration in waves (typ. <= 7500)",
          enable_if=("method", "Van der Meer")),
    Field("perm", "Notional permeability P (Van der Meer)", "float", "", "", default=0.4, lo=0.1, hi=0.6,
          note="Van der Meer only: 0.1 impermeable core ... 0.5-0.6 homogeneous mound",
          enable_if=("method", "Van der Meer")),
    Field("S_damage", "Damage level S (Van der Meer)", "float", "", "", default=2.0, lo=1.0, hi=30.0,
          note="Van der Meer only: 2 = start of damage; higher = more allowed damage",
          enable_if=("method", "Van der Meer")),
)

# Complete output list (ACES User's Guide 4-1).
OUTPUTS = (
    Out("W",   "Weight of individual armor unit", "kN", "tons", "scalar",
        note="Required weight of a single primary armor unit (Hudson, or via D_n50 for Van der Meer)."),
    Out("B",   "Crest width of breakwater",       "m",  "ft",   "scalar",
        note="Minimum crest width, using a 3-unit-wide cap (B = 3 k_delta (W/gamma_r)^(1/3))."),
    Out("r",   "Average cover layer thickness",   "m",  "ft",   "scalar",
        note="Average two-layer cover thickness over the armor (r = n k_delta (W/gamma_r)^(1/3))."),
    Out("N_r", "Armor units per 1000 ft^2",       "",   "",     "scalar",
        note="Placement density: number of armor units per 1000 ft^2 of structure surface."),
    Out("Ns",  "Stability number Hs/(d_n50)",     "",   "",     "scalar",
        note="Dimensionless stability number H_s/(Delta d_n50); compare against the method's allowable."),
)

# --- "Method / Equations" panel content (rendered on the tool page in both front-ends).
# `method_key` names the choice input that selects the formulation; each method's `when`
# matches that input's value so the active method is highlighted. Equation `tex` strings
# stay in the KaTeX-intersect-matplotlib-mathtext subset so both renderers agree.
ABOUT = {
    "summary": (
        "Sizes the primary armor units of a rubble-mound breakwater or revetment and "
        "reports crest width, cover-layer thickness, and armor-unit placement density. "
        "Choose the classic Hudson stability equation or the modern Van der Meer (1988) "
        "rock-armor formulae."
    ),
    "method_key": "method",
    "methods": [
        {
            "name": "Hudson (1953-61)",
            "when": "Hudson",
            "tag": "legacy",
            "note": "Superseded for new USACE design, but retained to back-check legacy "
                    "designs and prior analyses. Reproduces ACES User's Guide Example 4-1.",
            "equations": [
                {"tex": r"W = \frac{\gamma_r\,H^3}{K_D\,(S_r-1)^3\,\cot\theta}",
                 "desc": "armor-unit weight"},
                {"tex": r"B = 3\,k_\Delta\,(W/\gamma_r)^{1/3}",
                 "desc": "crest width (minimum 3 units)"},
                {"tex": r"r = n\,k_\Delta\,(W/\gamma_r)^{1/3}",
                 "desc": "average cover-layer thickness (n layers)"},
                {"tex": r"N_r = 1000\,n\,k_\Delta\,(1-P/100)\,(\gamma_r/W)^{2/3}",
                 "desc": "armor units per 1000 ft^2"},
            ],
        },
        {
            "name": "Van der Meer (1988)",
            "when": "Van der Meer",
            "tag": "preferred",
            "note": "Modern CIRIA Rock Manual standard; accounts for wave period, storm "
                    "duration N, notional permeability P, and allowed damage level S.",
            "equations": [
                {"tex": r"\xi_m = \frac{\tan\alpha}{\sqrt{2\pi H / (g\,T_m^2)}}",
                 "desc": "surf-similarity parameter (mean period)"},
                {"tex": r"\frac{H}{\Delta D_{n50}} = 6.2\,P^{0.18}\,(S/\sqrt{N})^{0.2}\,\xi_m^{-0.5}",
                 "desc": "plunging waves"},
                {"tex": r"\frac{H}{\Delta D_{n50}} = P^{-0.13}\,(S/\sqrt{N})^{0.2}\,\sqrt{\cot\theta}\,\xi_m^{\,P}",
                 "desc": "surging waves"},
                {"tex": r"W = \gamma_r\,D_{n50}^{\,3}",
                 "desc": "armor-unit weight from nominal diameter"},
            ],
        },
    ],
    "symbols": [
        ["W", "armor-unit weight"],
        ["H", "design wave height"],
        ["K_D", "stability coefficient (SPM Table 7-8)"],
        ["S_r", "specific gravity of armor = gamma_r / gamma_w"],
        ["Delta", "relative buoyant density, S_r - 1"],
        ["k_delta", "layer coefficient (SPM Table 7-13)"],
        ["theta", "seaward structure slope angle (cot theta input)"],
        ["D_n50", "nominal median armor diameter, (W/gamma_r)^(1/3)"],
    ],
    "references": [
        "Hudson (1953-61)", "SPM (1984) Ch. 7", "EM 1110-2-2904",
        "Van der Meer (1988)", "CIRIA/CUR Rock Manual (2007)",
    ],
}

_N_CREST = 3      # ACES/SPM: crest width uses a minimum of 3 armor units


@dataclass
class Result:
    W: float; B: float; r: float; N_r: float; Ns: float
    notes: str = ""


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind in ("float", "int", "angle"):
            v = inp.get(f.key, f.default)     # optional Van der Meer inputs fall back to defaults
            if not (f.lo <= v <= f.hi):
                raise ValueError(f"{f.label} ({f.key}) = {v} outside [{f.lo}, {f.hi}] ({f.note})")
    if inp["w_r"] <= inp["w_w"]:
        raise ValueError("Armor unit weight must exceed water unit weight (S_r > 1)")


# --- compute (single entry point both front-ends call) --------------------------
def compute(inp: dict) -> Result:
    """Hudson breakwater design for SI inputs
    {armor_type, w_r[N/m^3], H[m], w_w[N/m^3], K_D, k_delta, P[%], cot_theta, n}."""
    _validate(inp)
    w_r = float(inp["w_r"]); H = float(inp["H"]); w_w = float(inp["w_w"])
    K_D = float(inp["K_D"]); k_delta = float(inp["k_delta"]); P = float(inp["P"])
    cot_theta = float(inp["cot_theta"]); n = int(inp["n"])
    method = str(inp.get("method", "Hudson"))

    S_r = w_r / w_w                                         # specific gravity of armor
    Delta = S_r - 1.0                                       # relative buoyant density
    if method == "Van der Meer":
        # Van der Meer (1988) deep-water rock-armor stability (W in N, via D_n50 [m]).
        Tm = float(inp.get("Tm", 8.0)); N_w = float(inp.get("N_waves", 7500))
        perm = float(inp.get("perm", 0.4)); S = float(inp.get("S_damage", 2.0))
        tan_a = 1.0 / cot_theta
        s_om = 2.0 * math.pi * H / (G * Tm * Tm)            # mean wave steepness Hs/L_om
        xi_m = tan_a / math.sqrt(s_om)                      # surf-similarity parameter
        xi_mc = (6.2 * perm ** 0.31 * math.sqrt(tan_a)) ** (1.0 / (perm + 0.5))
        SN = (S / math.sqrt(N_w)) ** 0.2
        if xi_m < xi_mc:                                    # plunging waves
            Ns = 6.2 * perm ** 0.18 * SN * xi_m ** -0.5
            branch = "plunging"
        else:                                              # surging waves
            Ns = perm ** -0.13 * SN * math.sqrt(cot_theta) * xi_m ** perm
            branch = "surging"
        Dn50 = H / (Delta * Ns)                            # nominal median diameter [m]
        W = w_r * Dn50 ** 3                                # weight (N)
        method_note = (f"Van der Meer 1988 ({branch}); xi_m={xi_m:.2f} (xi_mc={xi_mc:.2f}), "
                       f"N_s={Ns:.2f}, P={perm}, S={S:.0f}, N={int(N_w)}")
    else:
        # (1) Hudson armor-unit weight (W in N, since w_r in N/m^3 and H in m)
        W = w_r * H ** 3 / (K_D * (S_r - 1.0) ** 3 * cot_theta)
        method_note = f"Hudson; K_D={K_D}"

    cube = (W / w_r) ** (1.0 / 3.0)                         # (W/w_r)^(1/3) = D_n50  [m]
    # (2) crest width (ACES uses n = 3 armor units for the crest)
    B = _N_CREST * k_delta * cube
    # (3) average cover-layer thickness (uses the input number of layers n)
    r = n * k_delta * cube
    # (4) placement density: number of units per 1000 ft^2 (US-convention count;
    #     (w_r/W)^(2/3) is dimensional, so evaluate on US-unit volume to match ACES)
    w_r_us = w_r / _LBF_PER_FT3
    W_us = W / _LBF
    N_r = 1000.0 * n * k_delta * (1.0 - P / 100.0) * (w_r_us / W_us) ** (2.0 / 3.0)
    Ns_out = H / (Delta * cube)                            # stability number Hs/(Delta*Dn50)

    notes = f"{inp.get('armor_type', 'armor')}; {method_note}; S_r={S_r:.3f}; n={n} (crest uses {_N_CREST})"
    return Result(W=W, B=B, r=r, N_r=N_r, Ns=Ns_out, notes=notes)


# --- self-tests + manual-example tabulation -------------------------------------
def _self_tests() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    assert abs(r.W / _TON - 1.59) < 0.02, r.W / _TON       # tons
    assert abs(r.B / _FT - 8.21) < 0.05, r.B / _FT         # ft
    assert abs(r.r / _FT - 5.47) < 0.05, r.r / _FT         # ft
    assert abs(r.N_r - 130.30) < 0.5, r.N_r
    # sanity: heavier armor (larger K_D) -> lighter unit
    r2 = compute({**{f.key: f.default for f in INPUTS}, "K_D": 20.0})
    assert r2.W < r.W
    # Van der Meer (1988): selectable; validate internal consistency + physical ranges
    base = {f.key: f.default for f in INPUTS}
    vd = compute({**base, "method": "Van der Meer"})
    assert 1.0 < vd.Ns < 4.0, vd.Ns                        # rock stability number, sane range
    assert vd.W > 0 and vd.B > 0 and vd.r > 0
    vd_S = compute({**base, "method": "Van der Meer", "S_damage": 8.0})
    assert vd_S.W < vd.W                                   # more allowed damage -> smaller unit
    vd_H = compute({**base, "method": "Van der Meer", "H": 16.0 * _FT})
    assert vd_H.W > vd.W                                   # larger wave -> larger unit
    print("  self-tests: PASS (Hudson matches Example 4-1; Van der Meer consistent)")


def _tab() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (values in US units; matches User's Guide Example 4-1)")
    print(f"  Armor unit weight W   = {r.W / _TON:8.2f} tons   [manual 1.59]")
    print(f"  Crest width B         = {r.B / _FT:8.2f} ft     [manual 8.21]")
    print(f"  Cover layer thickness = {r.r / _FT:8.2f} ft     [manual 5.47]")
    print(f"  Units per 1000 ft^2   = {r.N_r:8.2f}        [manual 130.30]")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _tab()
