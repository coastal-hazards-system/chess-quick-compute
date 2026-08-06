# CHESS-QC — application status

What each application is, what source material exists for it, and what it has
been checked against. Generated from the repository by `common/gen_app_status.py`;
re-run it after any change so this cannot drift from the code.


40 applications. 25 of them implement the 23 applications of the original
ACES suite: ACES 3-3 and 6-1 each serve two here, because CHESS-QC splits the
wedge kernel into point and grid forms, and the longshore method into
single-condition and wave-climate forms. The remaining 15 are Quick Compute
additions in the same style, with no ACES counterpart.

**Classification** is the module's own: *exact* means the method is reproduced
with nothing fitted or guessed; *standard* means a documented approximation or an
unresolved residual remains. Both are stated per application below.

## Summary

| ID | Application | Class | ACES | FORTRAN | MATLAB | Python | DOS cases | Source run |
|---|---|---|---|---|---|---|---|---|
| 1-1 | Windspeed Adjustment and Wave Growth | exact | 1-1 | yes | yes | yes | 1017 | — |
| 1-2 | Beta-Rayleigh Distribution | exact | 1-2 | yes | yes | yes | 315 | — |
| 1-3 | Extremal Significant Wave Height Analysis | exact | 1-3 | yes | yes | yes | — | — |
| 1-4 | Constituent Tide Record Generation | exact | 1-4 | yes | yes | yes | — | — |
| 1-5 | Near-surface Wind Speeds | standard | - | — | — | — | — | — |
| 1-6 | Holland Hurricane Wind Model | exact | - | — | — | — | — | — |
| 2-1 | Linear Wave Theory | exact | 2-1 | yes | yes | yes | 516 | — |
| 2-2 | Cnoidal Wave Theory | exact | 2-2 | yes | yes | yes | 236 | — |
| 2-3 | Fourier Series Wave Theory | exact | 2-3 | yes | yes | — | — | fwt |
| 2-4 | Wave Parameters | exact | - | — | — | — | — | — |
| 2-5 | Solitary Wave Theory | exact | - | — | — | — | — | — |
| 3-1 | Linear Wave Theory with Snell's Law | exact | 3-1 | yes | yes | yes | 4579 | — |
| 3-2 | Irregular Wave Transformation (Goda's Method) | standard | 3-2 | yes | yes | yes | 794 | goda |
| 3-3 | Combined Diffraction and Reflection by a Vertical Wedge | exact | 3-3 | yes | yes | yes | 214 | — |
| 3-4 | Vertical-Wedge Diffraction/Reflection on a Uniform Grid | exact | 3-3 | yes | — | — | — | — |
| 4-1 | Breakwater Design (Hudson) | exact | 4-1 | yes | yes | yes | 1530 | — |
| 4-2 | Toe Protection Design | exact | 4-2 | yes | yes | yes | 676 | — |
| 4-3 | Nonbreaking Wave Forces at Vertical Walls | exact | 4-3 | yes | yes | yes | 2512 | — |
| 4-4 | Rubble-Mound Revetment Design | exact | 4-4 | yes | yes | yes | 312 | — |
| 5-1 | Irregular Wave Runup on Beaches | exact | 5-1 | yes | yes | yes | 2146 | — |
| 5-2 | Wave Runup and Overtopping on Impermeable Structures | exact | 5-2 | yes | yes | yes | 3200 | — |
| 5-3 | Wave Transmission on Impermeable Structures | exact | 5-3 | yes | yes | yes | 8982 | — |
| 5-4 | Wave Transmission through Permeable Structures | exact | 5-4 | yes | yes | yes | 510 | — |
| 5-5 | Wave Setup | exact | - | — | — | — | — | — |
| 6-1 | Longshore Sediment Transport | exact | 6-1 | yes | yes | — | — | lst / lst92 |
| 6-2 | Time-Dependent Beach and Dune Erosion | exact | 6-2 | yes | — | — | — | kd |
| 6-3 | Longshore Transport from a CEDRS Wave Climate | exact | 6-1 | yes | — | — | — | lst / lst92 |
| 6-4 | Beach Nourishment Overfill Ratio and Volume | exact | 6-4 | yes | yes | yes | 13 | — |
| 6-5 | Composite Grain-Size Distribution | exact | 6-3 | yes | yes | — | — | cgs |
| 7-1 | Spatially Integrated Numerical Model for Inlet Hydraulics | exact | 7-1 | yes | — | — | — | inlet |
| 7-2 | Wave-Current Interaction in Channels | exact | - | — | — | — | — | — |
| 8-1 | Properties of Rectangular Basins | exact | - | — | — | — | — | — |
| 8-2 | Vessel-Generated Waves | exact | - | — | — | — | — | — |
| 8-3 | Surging of a Moored Vessel | exact | - | — | — | — | — | — |
| 9-1 | Bathystrophic Storm Surge | standard | - | — | — | — | — | — |
| 10-1 | Water Level Detrending | standard | - | — | — | — | — | — |
| 10-2 | Non-Tidal Residual | standard | - | — | — | — | — | — |
| 10-3 | Peaks Over Threshold | standard | - | — | — | — | — | — |
| 10-4 | Probabilistic Simulation Technique | standard | - | — | — | — | — | — |
| M-1 | Miscellaneous Breaker and Steepness Routines | standard | - | yes | — | — | — | — |

*ACES* is the application's number in the original suite, which differs from
CHESS-QC's in two places: CHESS-QC 6-3 is a capability of ACES 6-1, and CHESS-QC
6-5 is ACES 6-3. *DOS cases* counts recorded cases from the Hawaii comparison
workbooks. *Source run* names the driver that compiles and runs the original
FORTRAN (`tests/aces_oracle/fortran`).

## Wave Prediction

### 1-1 — Windspeed Adjustment and Wave Growth

*exact*. Adjusts an observed wind to a 10 m, 1-hour, overwater equivalent, then grows a wave field from it over open water or a restricted fetch.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 1-1 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `WAWG.FOR`, driver `WAWG`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/wind_adj.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/wind_adj.py`
- Recorded DOS output: 1017 cases over 2 sheet(s), 2 of 4 fields comparing clean
- Literature: SPM (1984) Ch.3; Resio & Vincent (1977); Smith (1991); TR 1-1

**Verification and validation**

- 1,017 recorded DOS cases; wind adjustment residual 0.79% tracks the 0.60% in the wind itself (FINDINGS D5). Growth outputs open at 1.4x tolerance
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 1-2 — Beta-Rayleigh Distribution

*exact*. Gives the full shallow-water wave-height distribution from a significant height, period and depth, using the beta-Rayleigh form rather than assuming Rayleigh.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 1-2 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `BETAR.FOR`, driver `BETRDRV`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/beta_rayleigh.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/beta_rayleigh.py`
- Recorded DOS output: 315 cases over 2 sheet(s), 10 of 10 fields comparing clean
- Literature: Hughes & Borgman (1987); Thompson & Vincent (1985); TR 1-2

**Verification and validation**

- 315 DOS cases, all 10 outputs clean. The H1/10 6.55-vs-6.30 dispute was settled from BETAR.FOR: ACES's own 100-bin quadrature (A1)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 1-3 — Extremal Significant Wave Height Analysis

*exact*. Fits extremal distributions to a series of annual maxima and returns design wave heights with confidence intervals for chosen return periods.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 1-3 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `EXTGODA.FOR`, driver `EXTWDRV`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/ext_Hs_analysis.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/ext_Hs_analysis.py`
- Literature: Goda (1988); Gringorten (1963); EM 1110-2-1414; TR 1-3

**Verification and validation**

- two recorded runs of the original program, all fitted distributions
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 1-4 — Constituent Tide Record Generation

*exact*. Predicts a water-level time series from tidal harmonic constituents, with the Schureman node factors and equilibrium arguments for the epoch.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 1-4 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `TIDES.FOR, CNS2LIB.FOR`, driver `TIDRV`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/tide_generation.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/tide_generation.py`
- Literature: Schureman (1971); Table A-5; TR 1-4

**Verification and validation**

- 481 points of a full printed record; rms 0.078 ft, max 0.160 ft. The astronomy is ruled out - node factors agree to 2.4e-6 (E6). Open (B9)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 1-5 — Near-surface Wind Speeds

*standard*. Converts a geostrophic wind to the near-surface wind over water through the Ekman-layer relations, accounting for air-sea temperature difference.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: ACES TR 1-1 eqs. 5-19; Silva (2005) eqs. 7.16-7.19; Lumley & Panofsky (1964)

**Verification and validation**

- analytic: the Ekman-layer relations and their limits
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 1-6 — Holland Hurricane Wind Model

*exact*. Builds the radial wind and pressure profile of a tropical cyclone from Holland's two-parameter model.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Holland (1980) Mon. Wea. Rev. 108; ACES manual

**Verification and validation**

- analytic: Holland's profile, its gradient-wind balance and known limits
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

## Wave Theory

### 2-1 — Linear Wave Theory

*exact*. First-order (Airy) wave theory: celerity, length, kinematics, pressure and energy for a wave of given height, period and depth.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 2-1 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `LWT.FOR`, driver `LWT`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/linear_wave_theory.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/linear_wave_theory.py`
- Recorded DOS output: 516 cases over 4 sheet(s), 48 of 56 fields comparing clean
- Literature: Airy (1845); Hunt (1979); TR 2-1

**Verification and validation**

- 516 DOS cases over 4 sheets; 48 of 56 fields clean
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 2-2 — Cnoidal Wave Theory

*exact*. Finite-amplitude periodic long-wave theory in terms of the cnoidal functions, for waves too long and steep for linear theory.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 2-2 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `CWT.FOR`, driver `CWT`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/cnoidal_wave_theory.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/cnoidal_wave_theory.py`
- Recorded DOS output: 236 cases over 2 sheet(s), 21 of 22 fields comparing clean
- Literature: Isobe (1985); Hardy & Kraus (1987); TR 2-2

**Verification and validation**

- DOS cases at both first and second order (the workbook stacks them, C5)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 2-3 — Fourier Series Wave Theory

*exact*. Steady progressive wave of permanent form solved as an N-term stream-function Fourier series, valid to near breaking where the analytic theories fail.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 2-3 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `FWT.FOR`, driver `FWT`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/fourier.m`
- **The original FORTRAN, compiled and runnable**: `sh tests/aces_oracle/fortran/build.sh fwt`
- Literature: Rienecker & Fenton (1981); Fenton (1988); TR 2-3

**Verification and validation**

- the ACES FORTRAN recompiled and run: 8 cases spanning deep and shallow water, near-breaking steepness, both celerity definitions and both current directions, to 0.008% on wavelength and celerity (E13)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 2-4 — Wave Parameters

*exact*. The classical wave-parameter presentation of linear theory: the dimensionless ratios and the deep, transitional and shallow-water limits.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Airy (1845); SPM (1984) App. C; Hunt (1979)

**Verification and validation**

- analytic identities against 2-1, of which it is a presentation
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 2-5 — Solitary Wave Theory

*exact*. Solitary wave theory: a single wave of translation lying entirely above the still-water level, with its celerity, profile and kinematics.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: McCowan (1891); Munk (1949); SPM (1984); CEM (EM 1110-2-1100); Zaroodny (1972)

**Verification and validation**

- analytic: the solitary-wave relations and their conservation properties
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

## Wave Transformation

### 3-1 — Linear Wave Theory with Snell's Law

*exact*. Transforms a wave from one depth to another, and to deep water, by linear theory with Snell's law refraction and energy-flux shoaling.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 3-1 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `LWTS.FOR`, driver `LWTS`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/snells_law.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/snells_law.py`
- Recorded DOS output: 4579 cases over 4 sheet(s), 92 of 92 fields comparing clean
- Literature: O'Brien (1942); Hunt (1979); Weggel (1972); TR 3-1

**Verification and validation**

- 4,579 DOS cases. Three workbook columns were found to be mislabelled and were re-identified on 100% of rows (C1)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 3-2 — Irregular Wave Transformation (Goda's Method)

*standard*. Transforms an irregular wave train shoreward by Goda's method, marching the spectrum through shoaling, refraction, breaking and wave setup.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 3-2 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `WSU.FOR`, driver `WSU`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/irr_wave_trans.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/irr_wave_trans.py`
- Recorded DOS output: 794 cases over 1 sheet(s), 4 of 11 fields comparing clean
- **The original FORTRAN, compiled and runnable**: `sh tests/aces_oracle/fortran/build.sh goda`
- Literature: Goda (1975, 1985); Mitsuyasu (1975); Shuto (1974)

**Verification and validation**

- the FORTRAN recompiled: all 11 outputs of a 133-case sweep to print rounding. Two Hawaii mistranscriptions corrected from it (E1, E2)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 3-3 — Combined Diffraction and Reflection by a Vertical Wedge

*exact*. Wave height, phase and pressure in the lee of a semi-infinite vertical wedge, combining diffraction and reflection at a single point.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 3-3 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `DFRAC.FOR`, driver `DFRAC`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/refdiff_vert_wedge.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/refdiff_vert_wedge.py`
- Recorded DOS output: 214 cases over 2 sheet(s), 6 of 8 fields comparing clean
- Literature: Chen (1987); Stoker (1957); Penny & Price (1952)

**Verification and validation**

- 149 DOS cases; the far-field metric rows sit at ~2% where ACES sums up to 200 alternating single-precision terms (B1, B2 fixed)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 3-4 — Vertical-Wedge Diffraction/Reflection on a Uniform Grid

*exact*. The same wedge kernel evaluated over a uniform grid, giving a field rather than a single point.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 3-3 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `DFRAC.FOR (grid evaluation of the same kernel)`, driver `DFRAC`
- Literature: Chen (1987); Stoker (1957); Penny & Price (1952)

**Verification and validation**

- the same kernel as 3-3 on a grid; checked against it pointwise
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

## Structural Design

### 4-1 — Breakwater Design (Hudson)

*exact*. Sizes primary armour units for a rubble-mound breakwater or revetment from the Hudson stability equation, with crest width, layer thickness and placement density.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 4-1 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `HUDSON.FOR`, driver `HUDSON`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/breakwater_Hudson.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/breakwater_Hudson.py`
- Recorded DOS output: 1530 cases over 2 sheet(s), 8 of 8 fields comparing clean
- Literature: Hudson (1953-61); SPM (1984) Ch.7; EM 1110-2-2904; TR 4-1

**Verification and validation**

- DOS cases across both unit systems. The armour-weight column was found to switch units mid-column and is handled (C2)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 4-2 — Toe Protection Design

*exact*. Designs the toe apron width and toe-stone weight for a vertical wall, bulkhead or revetment.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 4-2 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `TOEPRO.FOR`, driver `TOEPRO`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/toe_design.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/toe_design.py`
- Recorded DOS output: 676 cases over 3 sheet(s), 7 of 7 fields comparing clean
- Literature: EM 1110-2-1614; Tanimoto, Yagyu & Goda (1982); Hunt (1979); TR 4-2

**Verification and validation**

- DOS cases; the wide-apron stability number is chaotic by construction and is documented rather than forced (B7)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 4-3 — Nonbreaking Wave Forces at Vertical Walls

*exact*. Standing-wave (clapotis) forces and overturning moments on a vertical wall, by both the Sainflou and Miche-Rundgren methods.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 4-3 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `WFVW.FOR`, driver `WFVW`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/wave_forces.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/wave_forces.py`
- Recorded DOS output: 2512 cases over 2 sheet(s), 24 of 24 fields comparing clean
- Literature: Sainflou (1928); Miche (1944); Rundgren (1958); Hunt (1979); TR 4-3

**Verification and validation**

- DOS cases, Sainflou and Miche-Rundgren, crest and trough
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 4-4 — Rubble-Mound Revetment Design

*exact*. Sizes the armour and filter layers of a rubble-mound revetment, with runup and the resulting crest elevation.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 4-4 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `RUBBLE.FOR`, driver `RUBBLE`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/rubble_mound.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/rubble_mound.py`
- Recorded DOS output: 312 cases over 2 sheet(s), 20 of 20 fields comparing clean
- Literature: Ahrens (1981); van der Meer (1988); Hudson (1958); Ahrens & Heimbaugh (1988)

**Verification and validation**

- DOS cases; runup exact
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

## Wave Runup, Transmission, and Overtopping

### 5-1 — Irregular Wave Runup on Beaches

*exact*. Estimates irregular-wave runup on a beach or a rough impermeable slope.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 5-1 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `RUBCH.FOR`, driver `RUBCH`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/irregular_runup.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/irregular_runup.py`
- Recorded DOS output: 2146 cases over 2 sheet(s), 10 of 10 fields comparing clean
- Literature: Mase (1989); Hunt (1959); Walton & Ahrens (1989)

**Verification and validation**

- DOS cases against User's Guide Example 5-1-4
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 5-2 — Wave Runup and Overtopping on Impermeable Structures

*exact*. Runup and overtopping rate on an impermeable structure, including the effect of an onshore wind.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 5-2 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `RUNOVT.FOR`, driver `RUNOVT`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/runup_overtopping.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/runup_overtopping.py`
- Recorded DOS output: 3200 cases over 4 sheet(s), 12 of 16 fields comparing clean
- Literature: Ahrens & McCartney (1975); Ahrens & Titus (1985); Weggel (1976); SPM (1984)

**Verification and validation**

- DOS cases, User's Guide Examples 1-7
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 5-3 — Wave Transmission on Impermeable Structures

*exact*. Wave transmission over an impermeable structure by overtopping, and the transmitted height behind it.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 5-3 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `RUNTRN.FOR`, driver `RUNTRN`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/wavetrans_imperm.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/wavetrans_imperm.py`
- Recorded DOS output: 8982 cases over 4 sheet(s), 3 of 3 fields comparing clean
- Literature: Seelig (1980); Seelig (1976); Ahrens & McCartney (1975); Ahrens & Titus (1985)

**Verification and validation**

- DOS cases, User's Guide Examples 1-4
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 5-4 — Wave Transmission through Permeable Structures

*exact*. Wave transmission through and over a permeable rubble-mound structure, by Madsen's method with the internal flow resistance.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 5-4 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `MADTRANS.FOR`, driver `MADTRANS`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/wavetrans_perm.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/wavetrans_perm.py`
- Recorded DOS output: 510 cases over 8 sheet(s), 28 of 36 fields comparing clean
- Literature: Madsen & White (1976); Seelig (1980); Ahrens & McCartney (1975)

**Verification and validation**

- DOS cases over three breakwater geometries; five source corrections made from MADTRANS.FOR. Three isolated rows remain open (B8)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 5-5 — Wave Setup

*exact*. Wave-induced change in mean water level across the surf zone: the set-down seaward of breaking and the set-up shoreward of it.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Longuet-Higgins & Stewart (1963); Weggel (1972); Singamsetti & Wind (1980)

**Verification and validation**

- closed-form Longuet-Higgins/Weggel. NOT an ACES application: the setup ACES reports comes from WSU, which serves 3-2. Earlier mapping was wrong
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

## Littoral Processes

### 6-1 — Longshore Sediment Transport

*exact*. Potential longshore sand transport rate from a single wave condition, by the CERC energy-flux method, from either deepwater or breaking inputs.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 6-1 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `Lstran.for (2/98) and lstran.forg (3/92)`, driver `LSTDRV`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/longshore_trans.m`
- **The original FORTRAN, compiled and runnable**: `sh tests/aces_oracle/fortran/build.sh lst`
- Literature: SPM (1984) Ch.4; Galvin (1979); TR 6-1

**Verification and validation**

- the FORTRAN recompiled, both shipped vintages: the published examples to 6 figures on the 3/92 source, and 5e-5 on both with quartz (E8)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 6-2 — Time-Dependent Beach and Dune Erosion

*exact*. Time-marches beach and dune profile change through a storm by the XSHORE explicit finite-difference scheme, reporting eroded volume and contour recession.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 6-2 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `KDMAIN.FOR, KD2-KD18`, driver `KDDRV1`
- **The original FORTRAN, compiled and runnable**: `sh tests/aces_oracle/fortran/build.sh kd`
- Literature: Kriebel (1984b, EBEACH); Kriebel & Dean (1985); Moore (1982); Dean (1977); ACES Technical Reference Ch. 6-2

**Verification and validation**

- the FORTRAN recompiled: 4 of the 5 shipped XSHORE decks, volume and all four contour changes to better than half the source's printed digit (E11)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 6-3 — Longshore Transport from a CEDRS Wave Climate

*exact*. Net and gross longshore transport from a full CEDRS directional wave climate, refracting and shoaling each height/period/direction combination in to breaking.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 6-1 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `Lstran.for, CEDRS capability (lscband, lscomp, LSCTRN)`, driver `LSTDRV`
- **The original FORTRAN, compiled and runnable**: `sh tests/aces_oracle/fortran/build.sh lst`
- Literature: SPM (1984) Ch.4; Galvin (1979, 1980); Gravens (1988); Wang, Kraus & Davis (1998); ACES Technical Reference Ch. 6-1

**Verification and validation**

- the FORTRAN recompiled over the shipped CEDRS station file: all 9 contributing bands and the net to within 0.6 yd^3/yr (E12)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 6-4 — Beach Nourishment Overfill Ratio and Volume

*exact*. Tells a nourishment designer how much borrow material a given beach fill needs, from the grain-size distributions of the native and borrow sand.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 6-4 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `BEACH.FOR`, driver `BFDRV`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/beach_nourishment.m`
- Hawaii 2017 Python conversion: `ACES_PYTHON/drivers/beach_nourishment.py`
- Recorded DOS output: 13 cases over 2 sheet(s), 6 of 6 fields comparing clean
- Literature: James (1975) TM-60; SPM (1984); TR 6-4

**Verification and validation**

- 13 usable DOS rows - all the workbook holds; ACES has no multi-case mode here. ACES's normal CDF is invalid for negative arguments (B6)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`
- In the recorded-output regression gate (`tests/test_aces_dos_oracle.py`)

### 6-5 — Composite Grain-Size Distribution

*exact*. Grain-size statistics of a composited sediment sample: the method of moments and the Folk graphic measures.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 6-3 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `CGS.FOR, CGSANL.FOR`, driver `CGS`
- Hawaii 2017 MATLAB conversion: `ACES_MATLAB/drivers/compositeGrain.m`
- **The original FORTRAN, compiled and runnable**: `sh tests/aces_oracle/fortran/build.sh cgs`
- Literature: Folk (1974); Krumbein (1934, 1938); SPM (1984); TR 6-3

**Verification and validation**

- the FORTRAN recompiled: all 9 statistics. It found a half-class bias in CHESS-QC's moment mean that the old self-referential test could not (E5)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

## Inlet Processes

### 7-1 — Spatially Integrated Numerical Model for Inlet Hydraulics

*exact*. Time-marches inlet discharge and bay water level through a tidal cycle for a multi-section inlet, resolving the channel flow net across each cross-section.

**Sources available**

- ACES Technical Reference and User's Guide, chapter 7-1 (`private/ACESTR.PDF`, `private/hawaii/.../Aces_UsersGuide`)
- ACES FORTRAN source: `INLMAIN.FOR, IH1-IH14`, driver `IHDRV`
- **The original FORTRAN, compiled and runnable**: `sh tests/aces_oracle/fortran/build.sh inlet`
- Literature: Seelig (1977); Seelig, Harris & Herchenroder (1977); Harris & Bodine (1977); Keulegan (1967); Schureman (1971)

**Verification and validation**

- the FORTRAN recompiled on the shipped deck: the full equal-discharge flow net, peak discharges to 0.01% and exchange volumes to 0.04% (E6)
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 7-2 — Wave-Current Interaction in Channels

*exact*. Wave-current interaction in a channel: how a following or opposing current changes wave height and length, and when the wave is blocked.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Jonsson (1990); Jonsson, Skovgaard & Wang (1970); Peregrine (1976)

**Verification and validation**

- analytic: dispersion and action identities, blocking, no-current unity
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

## Harbor Design

### 8-1 — Properties of Rectangular Basins

*exact*. Natural resonant oscillation periods of a rectangular harbour basin, open or closed, including the Helmholtz mode.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Merian's formula; Helmholtz mode; SPM (1984); Wilson (1972)

**Verification and validation**

- analytic: Merian, 2-D and Helmholtz closed forms and their reductions
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 8-2 — Vessel-Generated Waves

*exact*. Waves generated by a vessel moving in a channel, and the resulting drawdown and return current.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Schijf (1949); PIANC (1987); EM 1110-2-1100; Kelvin/Havelock ship-wave theory

**Verification and validation**

- analytic: Schijf continuity and energy, the 35.26 deg deep-water crest
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

### 8-3 — Surging of a Moored Vessel

*exact*. Surge motion of a vessel held at a berth by mooring lines, and the resulting line loads.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: EM 1110-2-1100 Part II; PIANC mooring guidelines; ACES manual

**Verification and validation**

- analytic: single-line period, symmetry, projection and scaling laws
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

## Storm Surge

### 9-1 — Bathystrophic Storm Surge

*standard*. Open-coast hurricane surge along a single cross-shelf traverse by the quasi-1D bathystrophic method.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Bodine (1971) TM-35; Holland (1980); Myers (1954); TR/CERC

**Verification and validation**

- analytic against Bodine (1971) TM-35; beyond the ACES set
- Self-tests in the module; covered by the shared suites

## Coastal Hazards

### 10-1 — Water Level Detrending

*standard*. Removes the long-term linear sea-level trend from a water-level record by least-squares regression.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Zervas (2009) NOAA CO-OPS 053; NTDE datum convention

**Verification and validation**

- analytic: least-squares trend removal and its invariants
- Self-tests in the module; covered by the shared suites

### 10-2 — Non-Tidal Residual

*standard*. Separates the non-tidal residual from a water-level record by removing the predicted astronomical tide.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: NOAA CO-OPS tide/residual practice; PyStorm NTR

**Verification and validation**

- analytic; chained from 10-1
- Self-tests in the module; covered by the shared suites

### 10-3 — Peaks Over Threshold

*standard*. Extracts independent storm peaks from a residual series by a peaks-over-threshold declustering.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Coles (2001); USACE coastal-hazards practice; PyStorm POT

**Verification and validation**

- analytic; chained from 10-2
- Self-tests in the module; covered by the shared suites

### 10-4 — Probabilistic Simulation Technique

*standard*. Builds a synthetic storm population from a fitted distribution and propagates it to a hazard curve.

**Sources available**

- No numbered ACES application; built from the cited literature below
- Literature: Coles (2001); Nadal-Caraballo et al. PST; PyStorm PST

**Verification and validation**

- analytic
- Self-tests in the module; covered by the shared suites

## Miscellaneous Routines

### M-1 — Miscellaneous Breaker and Steepness Routines

*standard*. The shared breaker and steepness relations the other applications call: the Miche steepness limit, the McCowan and Weggel breaking criteria, and the breaker height at a structure.

**Sources available**

- No numbered ACES application; built from the cited literature below
- ACES FORTRAN source: `HLIMIT.FOR (shared breaker utilities, not a menu entry)`, driver `HLIMIT`
- Literature: Weggel (1972); Miche (1944); McCowan (1894); Singamsetti & Wind (1980)

**Verification and validation**

- the breaker and steepness limits against their published forms, and for consistency with 2-1 and 3-1 which use them
- Self-tests in the module; a dedicated case in `tests/test_manual_oracle.py`

---

## How to read the verification column

Three kinds of evidence appear, and they are not equally strong:

1. **The original FORTRAN, recompiled and run.** The strongest available: the
   ACES computational routines separate cleanly from their screen layer and build
   with the compiler shipped alongside them, so results can be compared at full
   precision on any input. Eight applications now have this.
2. **Recorded DOS output.** 28,742 cases from the Hawaii workbooks, compared with
   absolute tolerances derived from the printed field widths rather than
   percentages. Sixteen applications have this.
3. **Analytic checks.** For applications with no ACES counterpart: closed-form
   identities, limiting cases and conservation properties. Weaker than either of
   the above, and labelled as such.

Where a published example and the source disagree, the source wins, and the
disagreement is recorded. Where the source and a recorded DOS listing disagree,
the listing wins. Open items and the full evidence are in
`tests/aces_oracle/FINDINGS.md`.
