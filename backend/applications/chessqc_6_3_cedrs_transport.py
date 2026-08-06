"""CHESS-QC application 6-3 — Longshore Transport from a CEDRS Wave Climate.

Originating ACES application: 6-1 "Longshore Sediment Transport", its third capability -
"Transport using CEDRS statistical data: Percent Occurrence of wave height and period by
direction" (`Lstran.for`, subroutines `lscband`, `lscomp` and `LSCTRN`). CHESS-QC carries
it as its own application because its inputs and outputs are nothing like the
single-condition case, which is CHESS-QC 6-1.

Classification: exact. Reproduces the ACES source, recompiled and run on the CEDRS file
ACES ships (`g1033.810`), to better than a tenth of a percent.

What this computes. A CEDRS (Coastal Engineering Data Retrieval System) station record
gives, for each of sixteen 22.5-degree direction bands, the percent occurrence of every
combination of wave height class and peak-period class over the record. For a shoreline
of a given azimuth, ACES selects the bands whose waves can reach it, refracts and shoals
each height/period combination from the gauge depth in to breaking, applies the CERC
energy-flux transport formula at the breaker line, and weights the result by that
combination's occurrence. The net is the signed sum over all of it.

    for each direction band from the shore-left to the shore-right limit:
      for each height class h and period class t with non-zero occurrence:
        solve for the breaking depth d_b and wavelength L_b            [LSCTRN]
        H_b = 0.78 d_b ,  alpha_b = asin(L_b sin(alpha_0) / L_0)
        P_ls = c_b rho g^1.5 H_b^2.5 sin(2 alpha_b)
        Q    = K P_ls / [(rho_s - rho) g a'] * (occurrence / 100)
      weight the band by its contributing fraction and accumulate

The breaking condition comes from a two-variable Newton solve on depth and wavelength
that conserves the shoaled, refracted energy flux, floored at 3 ft wavelength and 0.5 ft
depth and capped at the gauge values, exactly as the source writes it.

This differs from what CHESS-QC previously did here, which applied the *deepwater* CERC
formula directly to each band's height with no period, depth or refraction at all. That
is a different calculation and it needed a fitted sediment density to match the published
example; see the note below and FINDINGS.md E12.

The two coefficient vintages. ACES ships its longshore source twice, and the 2/98
revision replaced solitary-wave energy-flux coefficients with linear-wave ones exactly
1.25 times smaller (see the 6-1 docstring, which sets this out in full). The published
User's Guide Example 3 net of -854,849 yd^3/yr is pre-1998 output: running the 3/92
source on `g1033.810` at a 40-degree shore azimuth gives -854,883, while the current 2/98
source gives -683,902. Both are offered here, defaulting to the linear pair because that
is what the current source computes.

An earlier reading recorded that the example needed an effective sediment density near
2319 kg/m^3, "below quartz and not stated in the Technical Reference". That was wrong on
both counts: the density stood in for the coefficient vintage (transport goes as
1/(rho_s - rho), and 2325 buys exactly the 1.25), and ACES's own sediment density is the
5.14 slug/ft^3 at `Lstran.for:1152`, which is 2649 kg/m^3 - quartz. No fitted quantity
remains.

Units: the application contract is SI - inputs arrive in SI and results are returned in
SI. LSCTRN is stated in US customary units throughout and carries hard limits in feet, so
compute() converts on entry and back on exit, using ACES's own g = 32.17 ft/s^2.

Self-containment: zero sibling imports; embeds the contract dataclasses, the eight
regional height and period tables, the band geometry, the breaking solver and the shipped
G1033 station record. stdlib only. Runnable:
    python chessqc_6_3_cedrs_transport.py

Validation, against the ACES source recompiled and run (tests/aces_oracle/fortran,
`sh build.sh lst`, which builds both vintages) on `g1033.810` at 40 degrees, K = 0.39:

    solitary (3/92)   -854,883 yd^3/yr      published Example 3: -854,849
    linear   (2/98)   -683,902 yd^3/yr

reproduced here to better than 0.1 percent on the net and on each of the nine
contributing bands.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

# ACES's own constants (Lstran.for:1770-1779), English units throughout.
G_US = 32.17              # ft/s^2
_CONV = 1168800.0         # ft^3/s -> yd^3/yr; equals 31,557,600 s/yr over 27 ft^3/yd^3
_SLUG_FT3 = 515.378818    # kg/m^3 per slug/ft^3
_FT = 0.3048
_M3_TO_YD3 = 1.0 / 0.764554858

# Energy-flux coefficient for the breaking form, both vintages (see the 6-1 module).
_LINEAR = "Linear (2/98 revision, SMS/Genesis consistent)"
_SOLITARY = "Solitary (original ACES, reproduces the User's Guide example)"
_THEORIES = (_LINEAR, _SOLITARY)
_COEF_BREAK = {_LINEAR: 0.07071,        # Lstran.for:1143
               _SOLITARY: 0.088388}     # lstran.forg:1127


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
    enable_if: tuple = ()


@dataclass(frozen=True)
class Out:
    key: str
    label: str
    unit_si: str = ""
    unit_us: str = ""
    kind: str = "scalar"
    note: str = ""


APP_META = AppMeta(
    aces_id="6-3",
    name="Longshore Transport from a CEDRS Wave Climate",
    area="Littoral Processes",
    classification="exact",
    cite="SPM (1984) Ch.4; Galvin (1979, 1980); Gravens (1988); "
         "Wang, Kraus & Davis (1998); ACES Technical Reference Ch. 6-1",
    default_system="SI",
)

# --- regional class tables (subroutine lscmb) -----------------------------------
# Mid-range wave heights in FEET and mid-band wave periods in seconds, by region.
# Atlantic, Pacific and the Great Lakes tables are carried so any CEDRS station can
# be used, not only the shipped Gulf one.
_REGIONS = ("Atlantic", "Gulf", "Pacific", "Lake Erie", "Lake Huron",
            "Lake Michigan", "Lake Ontario", "Lake Superior")
_RHEIGHT = {
    "Atlantic":      (1.64, 4.92, 8.20, 11.48, 14.76, 18.04, 21.33, 24.61, 27.89, 31.17, 34.45),
    "Gulf":          (0.82, 2.46, 4.10, 5.74, 7.38, 9.02, 10.66, 12.30, 13.94, 15.58, 17.22),
    "Pacific":       (1.64, 4.92, 8.20, 11.48, 14.76, 18.04, 21.33, 24.61, 27.89, 31.17, 34.45),
    "Lake Erie":     (0.39, 1.21, 2.03, 2.85, 3.67, 4.49, 5.31, 6.14, 6.96, 7.78, 8.60, 9.42, 10.24, 11.06, 11.88),
    "Lake Huron":    (0.39, 1.21, 2.03, 2.85, 3.67, 4.49, 5.31, 6.14, 6.96, 7.78, 8.60, 9.42, 10.24, 11.06, 11.88),
    "Lake Michigan": (0.39, 1.21, 2.03, 2.85, 3.67, 4.49, 5.31, 6.14, 6.96, 7.78, 8.60, 9.42, 10.24, 11.06, 11.88),
    "Lake Ontario":  (0.39, 1.21, 2.03, 2.85, 3.67, 4.49, 5.31, 6.14, 6.96, 7.78, 8.60, 9.42, 10.24, 11.06, 11.88),
    "Lake Superior": (0.39, 1.21, 2.03, 2.85, 3.67, 4.49, 5.31, 6.14, 6.96, 7.78, 8.60, 9.42, 10.24, 11.06, 11.88),
}
_RPERIOD = {
    "Atlantic":      (3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0, 11.0, 12.0),
    "Gulf":          (4.00, 4.50, 5.88, 6.90, 8.00, 9.10, 10.0, 11.1, 12.5, 14.3),
    "Pacific":       (5.13, 6.90, 8.70, 10.0, 11.1, 12.5, 14.3, 16.7, 20.0, 25.0),
    "Lake Erie":     (2.50, 3.23, 4.30, 5.41, 6.45, 7.41, 8.30, 9.10, 10.0, 11.77),
    "Lake Huron":    (2.00, 3.00, 4.35, 5.41, 6.45, 7.41, 8.30, 9.10, 10.0, 12.5),
    "Lake Michigan": (2.00, 3.00, 4.30, 5.41, 6.45, 7.41, 8.30, 9.10, 10.0, 12.5),
    "Lake Ontario":  (2.50, 3.23, 4.39, 5.41, 6.45, 7.41, 8.30, 9.10, 10.0, 11.77),
    "Lake Superior": (2.00, 3.00, 4.35, 5.41, 6.45, 7.41, 8.30, 9.10, 10.0, 13.33),
}

# Direction-band geometry (subroutine lscband). Sixteen 22.5-degree bands, extended to
# twenty-five so a shoreline whose window wraps past due north can be handled by index.
_CADB = tuple(22.5 * i for i in range(25))
_LLIMIT = tuple(-11.25 + 22.5 * i for i in range(25))
_RLIMIT = tuple(11.25 + 22.5 * i for i in range(25))
_DIRBND = tuple(list(range(1, 17)) + list(range(1, 10)))


# ACES CEDRS station G1033 (Gulf of Mexico, 29.0 N 85.5 W, 68.0 m), the file
# ACES ships as g1033.810 and names as its own default. Percent occurrence x1000
# of wave height by peak period, for each of the 16 direction bands: 16 blocks of
# 11 height classes, each row holding the 10 period columns. Band-major order.
_G1033 = (
    # band 1: central azimuth 0.0 deg
    (34, 30, 5, 0, 0, 0, 0, 0, 0, 0),
    (378, 532, 30, 0, 0, 0, 0, 0, 0, 0),
    (0, 780, 99, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 345, 3, 0, 0, 0, 0, 0, 0),
    (0, 0, 42, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 2: central azimuth 22.5 deg
    (66, 41, 15, 0, 0, 0, 0, 0, 0, 0),
    (321, 602, 56, 0, 0, 0, 0, 0, 0, 0),
    (1, 888, 104, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 359, 1, 0, 0, 0, 0, 0, 0),
    (0, 0, 58, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 3: central azimuth 45.0 deg
    (75, 18, 23, 0, 0, 0, 0, 0, 0, 0),
    (629, 915, 39, 0, 0, 0, 0, 0, 0, 0),
    (1, 1302, 135, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 592, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 126, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 4: central azimuth 67.5 deg
    (213, 53, 135, 0, 0, 0, 0, 0, 0, 0),
    (675, 2250, 181, 0, 0, 0, 0, 0, 0, 0),
    (0, 3151, 532, 3, 0, 0, 0, 0, 0, 0),
    (0, 3, 1052, 6, 0, 0, 0, 0, 0, 0),
    (0, 0, 147, 3, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 5: central azimuth 90.0 deg
    (472, 296, 75, 0, 0, 0, 0, 0, 0, 0),
    (980, 4433, 374, 0, 0, 0, 0, 0, 0, 0),
    (0, 2320, 4712, 159, 0, 0, 0, 0, 0, 0),
    (0, 18, 874, 831, 44, 0, 0, 0, 0, 0),
    (0, 0, 71, 39, 100, 0, 0, 0, 0, 0),
    (0, 0, 5, 0, 5, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 6: central azimuth 112.5 deg
    (412, 345, 112, 0, 0, 0, 0, 0, 0, 0),
    (925, 4404, 545, 0, 0, 0, 0, 0, 0, 0),
    (0, 716, 4609, 547, 0, 0, 0, 0, 0, 0),
    (0, 0, 246, 828, 58, 0, 0, 0, 0, 0),
    (0, 0, 6, 23, 102, 1, 0, 0, 0, 0),
    (0, 0, 0, 3, 27, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 1, 5, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 7: central azimuth 135.0 deg
    (496, 381, 39, 0, 0, 0, 0, 0, 0, 0),
    (853, 3417, 662, 1, 0, 0, 0, 0, 0, 0),
    (1, 638, 3307, 513, 3, 0, 0, 0, 0, 0),
    (0, 0, 119, 961, 136, 0, 0, 0, 0, 0),
    (0, 0, 0, 23, 160, 3, 0, 0, 0, 0),
    (0, 0, 0, 0, 15, 15, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 11, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 8: central azimuth 157.5 deg
    (133, 128, 77, 0, 0, 0, 0, 0, 0, 0),
    (627, 1875, 477, 6, 0, 0, 0, 0, 0, 0),
    (1, 487, 2345, 319, 1, 0, 0, 0, 0, 0),
    (0, 0, 119, 766, 203, 0, 0, 0, 0, 0),
    (0, 0, 0, 25, 249, 6, 0, 0, 0, 0),
    (0, 0, 0, 0, 39, 5, 0, 0, 0, 0),
    (0, 0, 0, 0, 3, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 9: central azimuth 180.0 deg
    (131, 71, 90, 0, 0, 0, 0, 0, 0, 0),
    (453, 1483, 277, 0, 0, 0, 0, 0, 0, 0),
    (0, 369, 1185, 191, 8, 0, 0, 0, 0, 0),
    (0, 0, 100, 503, 97, 0, 0, 0, 0, 0),
    (0, 0, 1, 34, 159, 1, 0, 0, 0, 0),
    (0, 0, 0, 0, 6, 5, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 10: central azimuth 202.5 deg
    (75, 66, 66, 0, 0, 0, 0, 0, 0, 0),
    (458, 1401, 241, 11, 0, 0, 0, 0, 0, 0),
    (0, 241, 922, 160, 17, 0, 0, 0, 0, 0),
    (0, 0, 87, 302, 53, 6, 1, 0, 0, 0),
    (0, 0, 1, 15, 71, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 15, 3, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 3, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 11: central azimuth 225.0 deg
    (78, 210, 124, 0, 0, 0, 0, 0, 0, 0),
    (557, 1810, 272, 25, 1, 0, 0, 0, 0, 0),
    (0, 361, 1630, 188, 22, 6, 0, 0, 0, 0),
    (0, 0, 83, 311, 59, 5, 0, 0, 0, 0),
    (0, 0, 0, 20, 85, 3, 0, 0, 0, 0),
    (0, 0, 0, 1, 35, 6, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 13, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 12: central azimuth 247.5 deg
    (131, 191, 73, 0, 0, 0, 0, 0, 0, 0),
    (557, 1813, 260, 6, 0, 0, 0, 0, 0, 0),
    (1, 278, 968, 181, 6, 0, 0, 0, 0, 0),
    (0, 0, 80, 244, 37, 0, 0, 0, 0, 0),
    (0, 0, 0, 23, 63, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 6, 0, 1, 0, 0, 0),
    (0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 13: central azimuth 270.0 deg
    (121, 97, 46, 0, 0, 0, 0, 0, 0, 0),
    (391, 1031, 361, 10, 0, 0, 0, 0, 0, 0),
    (1, 140, 740, 189, 10, 0, 0, 0, 0, 0),
    (0, 0, 106, 290, 68, 0, 0, 0, 0, 0),
    (0, 0, 8, 46, 90, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 25, 1, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 14: central azimuth 292.5 deg
    (44, 56, 58, 0, 0, 0, 0, 0, 0, 0),
    (301, 939, 330, 1, 0, 0, 0, 0, 0, 0),
    (0, 210, 869, 116, 1, 0, 0, 0, 0, 0),
    (0, 1, 179, 453, 6, 0, 0, 0, 0, 0),
    (0, 0, 8, 119, 42, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 3, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 15: central azimuth 315.0 deg
    (49, 71, 35, 0, 0, 0, 0, 0, 0, 0),
    (256, 872, 106, 0, 0, 0, 0, 0, 0, 0),
    (0, 669, 588, 27, 0, 0, 0, 0, 0, 0),
    (0, 1, 285, 200, 1, 0, 0, 0, 0, 0),
    (0, 0, 10, 37, 1, 0, 0, 0, 0, 0),
    (0, 0, 1, 5, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    # band 16: central azimuth 337.5 deg
    (63, 41, 10, 0, 0, 0, 0, 0, 0, 0),
    (369, 619, 100, 0, 0, 0, 0, 0, 0, 0),
    (1, 795, 123, 0, 0, 0, 0, 0, 0, 0),
    (0, 5, 306, 3, 0, 0, 0, 0, 0, 0),
    (0, 0, 58, 3, 0, 0, 0, 0, 0, 0),
    (0, 0, 3, 1, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
)


INPUTS = (
    Field("shore_azimuth", "Shore-normal azimuth", "angle", "deg", "deg", default=40.0,
          lo=0.0, hi=360.0,
          note="compass bearing of the outward shore normal; sets which direction "
               "bands reach the beach and at what angle"),
    Field("K", "Empirical coefficient K", "float", "", "", default=0.39, lo=0.0, hi=2.0,
          note="CERC coefficient; 0.39 for field data with significant wave height"),
    Field("gauge_depth", "Water depth at the gauge", "float", "m", "ft",
          default=68.0, lo=1.0, hi=5000.0,
          note="depth the CEDRS wave record applies at; the station header carries it "
               "(68.0 m for the shipped G1033 record)"),
    Field("region", "Coastline region", "choice", "", "", default="Gulf",
          choices=_REGIONS,
          note="selects the height and period class tables, which differ by region"),
    Field("wave_theory", "Energy-flux wave theory", "choice", "", "", default=_LINEAR,
          choices=_THEORIES,
          note="ACES shipped both. The solitary coefficient is 1.25x the linear one and "
               "is what the published example was computed with"),
    Field("rho_water", "Water density", "float", "kg/m^3", "kg/m^3", default=1025.18,
          lo=900.0, hi=1100.0, note="seawater ~1025, fresh ~1000"),
    Field("rho_sand", "Sediment density", "float", "kg/m^3", "kg/m^3", default=2650.0,
          lo=1500.0, hi=3500.0,
          note="quartz sand ~2650, which is what ACES uses: its 5.14 slug/ft^3 is "
               "2649 kg/m^3"),
    Field("porosity", "Sediment porosity", "float", "", "", default=0.40, lo=0.0, hi=0.7,
          note="pore fraction; solids fraction a' = 1 - porosity (ACES uses a' = 0.6)"),
    Field("occ", "CEDRS percent occurrence (x1000)", "matrix", "", "", default=_G1033,
          note="one row per (direction band, wave height class), band-major: the first "
               "11 rows are band 1's height classes, the next 11 band 2's, and so on for "
               "16 bands. Each row holds the 10 peak-period columns. This is the CEDRS "
               "station file laid out flat"),
)

OUTPUTS = (
    Out("Q_net", "Net longshore transport", "m^3/yr", "yd^3/yr", "scalar",
        note="Signed sum of the transport from every contributing direction band; negative is transport toward shore-left."),
    Out("Q_gross", "Gross longshore transport", "m^3/yr", "yd^3/yr", "scalar",
        note="Sum of the magnitudes of the per-band rates, the total volume moved in both directions. A CHESS-QC addition; ACES reports only the net and the per-band values."),
    Out("Q_right", "Transport to shore-right", "m^3/yr", "yd^3/yr", "scalar",
        note="Sum of the positive per-band rates."),
    Out("Q_left", "Transport to shore-left", "m^3/yr", "yd^3/yr", "scalar",
        note="Sum of the negative per-band rates."),
    Out("profile_X", "Wave approach angle", "deg", "deg", "profile",
        note="Angle of each contributing band's waves to the shore normal, positive toward shore-right."),
    Out("band_pct", "Contributing fraction", "%", "%", "profile",
        note="Fraction of the band that lies within the shoreline's 180-degree window. The two end bands are only partly exposed; the rest contribute in full."),
    Out("band_Q", "Transport by band", "m^3/yr", "yd^3/yr", "profile",
        note="Transport contributed by each direction band, after its contributing fraction is applied."),
    Out("n_bands", "Contributing bands", "", "", "scalar",
        note="Number of direction bands that reach this shoreline; nine for any azimuth, being a 180-degree window across 22.5-degree bands."),
)


@dataclass
class Result:
    Q_net: float; Q_gross: float; Q_right: float; Q_left: float
    profile_X: list; band_pct: list; band_Q: list
    n_bands: int
    notes: str = ""


def band_geometry(snaz: float):
    """Which direction bands reach a shoreline, and at what angle (subroutine lscband).

    Returns (lband, rband, perleft, perrght, zzin) with 1-based band indices into the
    25-element extended arrays."""
    ssnaz = snaz + 360.0 if snaz < 78.75 else snaz
    left = ssnaz - 90.0
    right = ssnaz + 90.0
    lband = rband = 0
    for i in range(25):
        if _LLIMIT[i] <= left < _RLIMIT[i]:
            lband = i + 1
            break
    for i in range(25):
        if _LLIMIT[i] < right <= _RLIMIT[i]:
            rband = i + 1
            break
    if lband == 0 or rband == 0:
        raise ValueError(f"shore azimuth {snaz} does not fall in the band table")
    perleft = (_CADB[lband - 1] - left + 11.25) / 22.5
    perrght = (right - _CADB[rband - 1] + 11.25) / 22.5
    zzin = tuple(ssnaz - _CADB[i] for i in range(25))
    return lband, rband, perleft, perrght, zzin


def breaking_transport(hin: float, tper: float, zin: float, din: float,
                       rper: float, k: float, coef: float,
                       rho: float, rho_s: float, a_solid: float,
                       g: float = G_US) -> float:
    """Transport from one height/period/direction cell (subroutine LSCTRN).

    All arguments in US customary units: height and depth in feet, period in seconds,
    angle in degrees, densities in slug/ft^3. `rper` is the percent occurrence. Returns
    yd^3/yr.

    Shoals and refracts the wave from the gauge depth in to breaking by a two-variable
    Newton solve on breaking depth and wavelength, then applies the CERC formula at the
    breaker line."""
    zin = zin * math.pi / 180.0

    # wavelength at the gauge, by Newton on the dispersion relation
    c1 = (g * tper ** 2.0) / (math.pi * 2.0)
    l = c1 * 2.0
    l2 = l / 2.0
    c2 = 2.0 * math.pi * din
    c3 = (g * tper ** 2.0) * din
    for _ in range(200):
        if abs(l2 - l) / l <= 0.001:
            break
        l = l2
        f1 = l - c1 * math.tanh(c2 / l)
        df1 = 1.0 + c3 / l ** 2.0 * (1.0 / math.cosh(c2 / l)) ** 2.0
        l2 = l - (f1 / df1)

    l1 = l2
    l2 = l2 / 2.5
    d2 = hin
    c2a1 = l1 * math.sinh(4.0 * math.pi * din / l1)
    c2 = hin ** 2.0 * l1 * math.cos(zin) / 0.6048 * (1.0 + 4.0 * math.pi * din / c2a1)
    c3 = math.sin(zin) / l1
    c4 = c1 * 2.0 * math.pi
    c5 = math.sin(zin) * math.sin(zin)
    c6 = l1 ** 2.0

    d2n = d2
    l2n = l2
    it = 0
    # the source loops on GOTO with no iteration cap; bounded here
    for _ in range(500):
        pid4 = 4.0 * math.pi * d2
        pidl = pid4 / l2
        pidl2 = pidl / 2.0
        root = 1.0 - abs(l2 * c3) ** 2.0
        root = math.sqrt(root) if root > 0.0 else 0.0
        f1dl = c2 - (d2 ** 2.0 * l2 * (1.0 + pid4 / (l2 * math.sinh(pidl)) * root))
        f2dl = l2 - c1 * math.tanh(pidl2)
        df1d = -d2 * root * (2.0 * l2 + (2.0 * pid4 / math.sinh(pidl))
                             * (1.5 - pid4 / (2.0 * l2 * math.tanh(pidl))))
        df1l1 = root * d2 ** 2.0 * l2
        df1l2 = 1.0 / l2 + l2 * c5 / (c6 - l2 ** 2.0 * c5)
        df1l3 = pid4 * c5 / (math.sinh(pidl) * c6 - l2 ** 2.0 * c5)
        df1l4 = pid4 ** 2.0 / (l2 ** 3.0 * math.sinh(pidl) * math.tanh(pidl))
        df1l = df1l1 * (df1l2 + df1l3 + df1l4)
        df2d = -c4 / l2 * (1.0 / math.cosh(pidl2)) ** 2.0
        df2l = 1.0 + c4 * d2 / l2 ** 2.0 * (1.0 / math.cosh(pidl2)) ** 2.0
        det = df1d * df2l - (df1l * df2d)
        if det == 0.0:
            break
        d = 1.0 / det
        d2n = d2 - (d * df2l * f1dl + d * (-df1l) * f2dl)
        l2n = l2 - (d * (-df2d) * f1dl + d * df1d * f2dl)
        it += 1
        if l2n < 3.0:
            l2n = 3.0            # breaking wavelength floor, 3 ft
        if d2n < 0.5:
            d2n = 0.5            # breaking depth floor, 0.5 ft
        if l2n > l1:
            l2n = l1
        if d2n > din:
            d2n = din
        if abs(d2n - d2) / d2 <= 0.001:
            break
        if it > 1:
            d2 = (d2 + d2n) / 2.0
            l2 = (l2 + l2n) / 2.0
            if l2 >= l1:
                l2 = l1 / 2.5
        else:
            d2 = d2n
            l2 = l2n

    hbr = 0.78 * d2n                                   # wave height at breaking
    arg = c3 * l2n
    arg = max(-1.0, min(1.0, arg))                     # the source does not guard asin
    zbr = math.asin(arg)                               # wave angle at breaking

    pls = coef * rho * (g ** 1.5) * (hbr ** 2.5) * math.sin(2.0 * zbr)
    q = pls * (k * _CONV) / ((rho_s - rho) * g * a_solid)
    return q * rper / 100.0


def _validate(inp: dict) -> None:
    for f in INPUTS:
        if f.kind not in ("float", "int", "angle"):
            continue
        v = float(inp[f.key])
        if not (f.lo <= v <= f.hi):
            raise ValueError(f"{f.label} out of range: {v} not in [{f.lo}, {f.hi}]")
    rows = inp.get("occ") or ()
    if len(rows) % 16 != 0:
        raise ValueError("the occurrence table must hold 16 equal direction bands, "
                         f"got {len(rows)} rows")


def compute(inp: dict, *, g: float = G_US) -> Result:
    """Net and per-band longshore transport from a CEDRS climate. SI in, SI out."""
    _validate(inp)
    snaz = float(inp["shore_azimuth"])
    k = float(inp["K"])
    din = float(inp["gauge_depth"]) / _FT                 # gauge depth, feet
    region = str(inp.get("region", "Gulf"))
    theory = str(inp.get("wave_theory", _LINEAR))
    if theory not in _COEF_BREAK:
        raise ValueError(f"unknown wave theory: {theory!r}")
    if region not in _RHEIGHT:
        raise ValueError(f"unknown region: {region!r}")
    coef = _COEF_BREAK[theory]
    rho = float(inp["rho_water"]) / _SLUG_FT3              # slug/ft^3
    rho_s = float(inp["rho_sand"]) / _SLUG_FT3
    if rho_s <= rho:
        raise ValueError("sediment density must exceed water density")
    a_solid = 1.0 - float(inp["porosity"])

    hin = _RHEIGHT[region]
    tper = _RPERIOD[region]
    rows = [list(r) for r in inp["occ"]]
    per_band = len(rows) // 16                             # height classes per band
    nh = min(per_band, len(hin))

    lband, rband, perleft, perrght, zzin = band_geometry(snaz)

    angles = []; pcts = []; qs = []
    for bandid in range(lband, rband + 1):
        band = _DIRBND[bandid - 1]
        if band == _DIRBND[lband - 1]:
            frac = perleft
        elif band == _DIRBND[rband - 1]:
            frac = perrght
        else:
            frac = 1.0
        z = zzin[bandid - 1]
        block = rows[(band - 1) * per_band:band * per_band]
        qband = 0.0
        for i1 in range(nh):
            row = block[i1]
            for i2 in range(min(10, len(row))):
                raw = float(row[i2])
                if raw > 0.0:
                    qband += breaking_transport(hin[i1], tper[i2], z, din,
                                                raw / 1000.0, k, coef,
                                                rho, rho_s, a_solid, g) * frac
        angles.append(z)
        pcts.append(frac * 100.0)
        qs.append(qband)

    # emit in ascending angle so the plotted series is single-valued
    order = sorted(range(len(angles)), key=lambda i: angles[i])
    angles = [angles[i] for i in order]
    pcts = [pcts[i] for i in order]
    qs = [qs[i] for i in order]

    yd3_to_m3 = 0.764554858
    q_net = sum(qs)
    q_right = sum(q for q in qs if q > 0.0)
    q_left = sum(q for q in qs if q < 0.0)
    direction = "shore-right" if q_net > 0 else ("shore-left" if q_net < 0 else "balanced")
    notes = (f"{region} class tables; {len(qs)} bands from {lband} to {rband} "
             f"(end fractions {perleft:.3f} / {perrght:.3f}); "
             f"{'linear' if theory == _LINEAR else 'solitary'} coefficient; "
             f"net {q_net:,.0f} yd^3/yr toward {direction}")
    return Result(
        Q_net=q_net * yd3_to_m3, Q_gross=(q_right - q_left) * yd3_to_m3,
        Q_right=q_right * yd3_to_m3, Q_left=q_left * yd3_to_m3,
        profile_X=angles, band_pct=pcts,
        band_Q=[q * yd3_to_m3 for q in qs],
        n_bands=len(qs), notes=notes)


# --- 'Method & equations' panel content (see chessqc_4_1 for the schema). ---
ABOUT = {
    'summary':
        'Estimates net and gross potential longshore sediment transport at a site from a '
        'CEDRS directional wave climate: the percent occurrence of every wave height and '
        'peak-period combination in each of sixteen 22.5-degree direction bands. For a '
        'shoreline of a given azimuth the bands that can reach it are selected, each '
        'height/period combination is refracted and shoaled from the gauge depth in to '
        'breaking, and the CERC energy-flux formula is applied at the breaker line and '
        'weighted by that combination\'s occurrence. Which of the two coefficient '
        'vintages ACES shipped is selectable: the linear-wave values its 2/98 revision '
        'adopted, or the original solitary-wave values 1.25 times larger that the '
        'published example was computed with.',
    'methods': [{
        'name': 'CEDRS climate, refracted and shoaled to breaking',
        'when': None,
        'tag': '',
        'note': None,
        'equations': [
            {'tex': r'\alpha_{i} = \theta_{n} - \phi_{i}, \quad '
                    r'\phi_{i} = 22.5^{\circ}(i-1)',
             'desc': 'Approach angle of band i to the shore normal, from the shore-normal '
                     'azimuth and the band central azimuth. Bands whose angle falls '
                     'outside +/-90 degrees do not reach the beach.'},
            {'tex': r'f_{L} = \frac{\phi_{L} - (\theta_{n}-90^{\circ}) + 11.25^{\circ}}'
                    r'{22.5^{\circ}}',
             'desc': 'Contributing fraction of the shore-left end band, which the '
                     '180-degree window only partly covers; the shore-right end band is '
                     'the mirror of this and the bands between contribute in full.'},
            {'tex': r'E_{f} = \frac{H^{2} L_{0} \cos\alpha_{0}}{0.6048}'
                    r'\left(1 + \frac{4\pi d_{0}}{L_{0}\sinh(4\pi d_{0}/L_{0})}\right)',
             'desc': 'The shoaled and refracted energy flux carried in from the gauge, '
                     'held constant to the breaker line.'},
            {'tex': r'd_{b}^{2} L_{b}\left(1 + \frac{4\pi d_{b}}'
                    r'{L_{b}\sinh(4\pi d_{b}/L_{b})}\sqrt{1-(L_{b}\sin\alpha_{0}/L_{0})^{2}}'
                    r'\right) = E_{f}, \quad L_{b} = L_{0}\tanh\frac{2\pi d_{b}}{L_{b}}',
             'desc': 'The coupled pair solved by Newton iteration for the breaking depth '
                     'and wavelength, floored at 0.5 ft depth and 3 ft wavelength and '
                     'capped at the gauge values.'},
            {'tex': r'H_{b} = 0.78\,d_{b}, \quad '
                    r'\alpha_{b} = \arcsin\!\left(L_{b}\sin\alpha_{0}/L_{0}\right)',
             'desc': 'Depth-limited breaker height and the refracted breaker angle from '
                     'Snell\'s law.'},
            {'tex': r'P_{ls} = c_{b}\,\rho\,g^{1.5}\,H_{b}^{2.5}\sin(2\alpha_{b})',
             'desc': 'Longshore energy-flux factor at the breaker line. c_b = 0.07071 '
                     'for linear wave theory, 0.088388 for the original solitary '
                     'approximation.'},
            {'tex': r'Q = \frac{K\,P_{ls}}{(\rho_{s}-\rho)\,g\,a}\cdot\frac{p}{100}',
             'desc': 'CERC transport rate for the cell, weighted by its percent '
                     'occurrence p. a = 1 - porosity. Summing over every cell of every '
                     'contributing band, each scaled by its band fraction, gives the net.'},
        ]}],
    'symbols': [
        ['Q', 'Potential volumetric longshore transport rate'],
        ['\\theta_n', 'Shore-normal azimuth, the outward compass bearing of the beach'],
        ['\\phi_i', 'Central azimuth of direction band i'],
        ['\\alpha_0', 'Wave approach angle at the gauge, relative to the shore normal'],
        ['\\alpha_b', 'Wave crest angle to the shoreline at breaking'],
        ['H_b', 'Breaker height, depth-limited at 0.78 of the breaking depth'],
        ['d_b, L_b', 'Water depth and wavelength at breaking'],
        ['d_0, L_0', 'Water depth and wavelength at the gauge'],
        ['E_f', 'Shoaled and refracted energy flux carried in from the gauge'],
        ['P_{ls}', 'Longshore energy-flux factor'],
        ['p', 'Percent occurrence of this height/period/direction combination'],
        ['f_L', 'Contributing fraction of an end band'],
        ['K', 'Empirical CERC coefficient (0.39 for field data with significant height)'],
        ['c_b', 'Breaking energy-flux coefficient; see the note on vintages'],
        ['a', 'Solids fraction of the bed, 1 - porosity'],
    ],
    'references': ['SPM (1984) Ch. 4', 'Galvin (1979, 1980)', 'Gravens (1988)',
                   'Wang, Kraus & Davis (1998)',
                   'ACES Technical Reference Ch. 6-1 (CEDRS capability)'],
}


# --- self-tests, against the recompiled ACES source -----------------------------
def _approx(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


# tests/aces_oracle/fortran, `sh build.sh lst`, driven over g1033.810 at 40 degrees
# with K = 0.39 and ACES's own constants. Per contributing band, in ACES's band order
# (shore-left first), yd^3/yr.
_SRC_BANDS_LINEAR = (7858.2, 54339.5, 71378.8, 48554.5, -24439.3,
                     -223345.9, -445583.6, -172637.6, -27.2)
_SRC_NET_LINEAR = -683902.5
_SRC_NET_SOLITARY = -854883.0
_PUBLISHED_EXAMPLE = -854849.0


def _aces_inputs(**over):
    """The shipped case with ACES's own constants, so the comparison is of the method."""
    inp = {f.key: f.default for f in INPUTS}
    inp.update(rho_water=1.989 * _SLUG_FT3, rho_sand=5.14 * _SLUG_FT3)
    inp.update(over)
    return inp


def _self_tests() -> None:
    # 1) band geometry for the example shoreline (subroutine lscband)
    lband, rband, perleft, perrght, zzin = band_geometry(40.0)
    assert (lband, rband) == (15, 23), (lband, rband)
    assert _approx(perleft, 0.72222, 1e-4), perleft
    assert _approx(perrght, 0.27778, 1e-4), perrght
    # zzin = azimuth - band centre, so the shore-LEFT band (smaller centre) is +85
    assert _approx(zzin[lband - 1], 85.0, 1e-9), zzin[lband - 1]
    assert _approx(zzin[rband - 1], -95.0, 1e-9), zzin[rband - 1]

    # 2) every azimuth sees exactly nine bands: a 180-degree window across 22.5-degree
    #    bands always spans nine of them, and the two end fractions sum to one
    for az in (0.0, 40.0, 78.6, 90.0, 180.0, 270.0, 359.0):
        lb, rb, pl, pr, _z = band_geometry(az)
        assert rb - lb == 8, (az, lb, rb)
        assert _approx(pl + pr, 1.0, 1e-9), (az, pl, pr)

    # 3) the march, against the source, band by band
    r = compute(_aces_inputs(wave_theory=_LINEAR))
    assert r.n_bands == 9, r.n_bands
    got = list(r.band_Q)[::-1]          # emitted by ascending angle; source is by band
    # combined absolute and relative tolerance: the source is single precision, and the
    # smallest contributing band here is -27 yd^3/yr, five orders below the net, where a
    # purely relative bound would be measuring rounding noise
    worst = worst_rel = 0.0
    for a, b in zip(_SRC_BANDS_LINEAR, got):
        d = abs(b * _M3_TO_YD3 - a)
        assert d <= 0.05 + 2e-5 * abs(a), (a, b * _M3_TO_YD3, d)
        worst = max(worst, d)
        worst_rel = max(worst_rel, d / abs(a))
    net = r.Q_net * _M3_TO_YD3
    assert _approx(net, _SRC_NET_LINEAR, 1.0), net

    # 4) the solitary vintage is exactly 1.25 times the linear one, and reproduces the
    #    published example
    rs = compute(_aces_inputs(wave_theory=_SOLITARY))
    nets = rs.Q_net * _M3_TO_YD3
    # the shipped coefficients are 1.25 to five figures (0.088388/0.07071 = 1.2500071)
    assert _approx(nets / net, 1.25, 1e-4), nets / net
    assert _approx(nets, _SRC_NET_SOLITARY, 5.0), nets
    assert abs(nets - _PUBLISHED_EXAMPLE) / abs(_PUBLISHED_EXAMPLE) < 1e-3, nets

    # 5) net is the signed sum of the bands; gross the sum of magnitudes
    assert _approx(r.Q_net, sum(r.band_Q), 1e-6)
    assert _approx(r.Q_gross, sum(abs(q) for q in r.band_Q), 1e-6)
    assert r.Q_gross >= abs(r.Q_net)
    assert _approx(r.Q_right + r.Q_left, r.Q_net, 1e-6)

    # 6) transport is linear in K, and the bands come out in ascending angle so the
    #    plotted series is single-valued
    rk = compute(_aces_inputs(wave_theory=_LINEAR, K=0.78))
    assert _approx(rk.Q_net, 2.0 * r.Q_net, abs(r.Q_net) * 1e-9), (rk.Q_net, r.Q_net)
    assert r.profile_X == sorted(r.profile_X), r.profile_X

    # 7) quartz and standard seawater, the shipped defaults, sit within a tenth of a
    #    percent of the source: the only difference is ACES's slightly different densities
    rq = compute({f.key: f.default for f in INPUTS})
    dq = abs(rq.Q_net * _M3_TO_YD3 - _SRC_NET_LINEAR) / abs(_SRC_NET_LINEAR)
    assert dq < 1e-3, dq

    print(f"  self-tests: PASS (bands 15-23 with fractions {perleft:.4f}/{perrght:.4f}; "
          f"9 bands match the source within {worst:.3f} yd^3/yr; net {net:,.0f} vs source "
          f"{_SRC_NET_LINEAR:,.0f}; solitary {nets:,.0f} vs published "
          f"{_PUBLISHED_EXAMPLE:,.0f})")


def _print_default_example() -> None:
    r = compute({f.key: f.default for f in INPUTS})
    print(f"\nACES application {APP_META.aces_id} - {APP_META.name}  [{APP_META.classification}]")
    print(f"  cite: {APP_META.cite}")
    print("  (default = the CEDRS record ACES ships, station G1033 in the Gulf of Mexico")
    print("   at 68 m, on a shoreline of azimuth 40 deg, K = 0.39, quartz sand)")
    print(f"    {'angle':>9} {'fraction':>9} {'transport':>16}")
    print(f"    {'(deg)':>9} {'(%)':>9} {'(yd^3/yr)':>16}")
    for a, p, q in zip(r.profile_X, r.band_pct, r.band_Q):
        print(f"    {a:9.2f} {p:9.2f} {q * _M3_TO_YD3:16,.0f}")
    print(f"    net   = {r.Q_net * _M3_TO_YD3:>14,.0f} yd^3/yr "
          f"({r.Q_net:,.0f} m^3/yr)")
    print(f"    gross = {r.Q_gross * _M3_TO_YD3:>14,.0f} yd^3/yr")
    print("  note: the User's Guide prints -854,849 yd^3/yr for this case. That is the")
    print("        pre-1998 solitary coefficient, 1.25x the linear one used here;")
    print("        select the solitary theory to reproduce it.")
    print(f"  notes: {r.notes}")


if __name__ == "__main__":
    print(f"CHESS-QC {APP_META.aces_id} {APP_META.name} - running self-tests...")
    _self_tests()
    _print_default_example()
