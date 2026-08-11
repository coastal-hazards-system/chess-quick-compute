# CHESS-QC — User Manual
*Coastal Hazards, Engineering, and Structures System (CHESS) — Quick Compute (QC)*

Each tool computes in SI internally and displays in the selected unit system; ranges/defaults below are shown in display units. *Auto-generated from the application contracts by `common/gen_docs.py`; do not edit by hand.*

Each app header carries a fidelity class:

- **(I) exact** — every coefficient and variable-relationship is known from the source (nothing guessed) and the results are validated.
- **(II) standard** — a named method that involves a self-made convention/inference or only partial validation.
- **(III) provisional** — a needed coefficient or relationship had to be guessed or is not recoverable from the sources (or there is no numeric oracle).

## Contents
- **Wave Prediction**
  - 1-1 Windspeed Adjustment and Wave Growth
  - 1-2 Beta-Rayleigh Distribution
  - 1-3 Extremal Significant Wave Height Analysis
  - 1-4 Constituent Tide Record Generation
  - 1-5 Near-surface Wind Speeds
  - 1-6 Holland Hurricane Wind Model
- **Wave Theory**
  - 2-1 Linear Wave Theory
  - 2-2 Cnoidal Wave Theory
  - 2-3 Fourier Series Wave Theory
  - 2-4 Wave Parameters
  - 2-5 Solitary Wave Theory
- **Wave Transformation**
  - 3-1 Linear Wave Theory with Snell's Law
  - 3-2 Irregular Wave Transformation (Goda's Method)
  - 3-3 Combined Diffraction and Reflection by a Vertical Wedge
  - 3-4 Vertical-Wedge Diffraction/Reflection on a Uniform Grid
- **Structural Design**
  - 4-1 Breakwater Design (Hudson)
  - 4-2 Toe Protection Design
  - 4-3 Nonbreaking Wave Forces at Vertical Walls
  - 4-4 Rubble-Mound Revetment Design
- **Wave Runup, Transmission, and Overtopping**
  - 5-1 Irregular Wave Runup on Beaches
  - 5-2 Wave Runup and Overtopping on Impermeable Structures
  - 5-3 Wave Transmission on Impermeable Structures
  - 5-4 Wave Transmission through Permeable Structures
  - 5-5 Wave Setup
- **Littoral Processes**
  - 6-1 Longshore Sediment Transport
  - 6-2 Time-Dependent Beach and Dune Erosion
  - 6-3 Longshore Transport from a CEDRS Wave Climate
  - 6-4 Beach Nourishment Overfill Ratio and Volume
  - 6-5 Composite Grain-Size Distribution
- **Inlet Processes**
  - 7-1 Spatially Integrated Numerical Model for Inlet Hydraulics
  - 7-2 Wave-Current Interaction in Channels
- **Harbor Design**
  - 8-1 Properties of Rectangular Basins
  - 8-2 Vessel-Generated Waves
  - 8-3 Surging of a Moored Vessel
- **Storm Surge**
  - 9-1 Bathystrophic Storm Surge
- **Coastal Hazards**
  - 10-1 Water Level Detrending
  - 10-2 Non-Tidal Residual
  - 10-3 Peaks Over Threshold
  - 10-4 Probabilistic Simulation Technique
- **Miscellaneous Routines**
  - M-1 Miscellaneous Breaker and Steepness Routines


## Wave Prediction

### 1-1 — Windspeed Adjustment and Wave Growth  `[I]`

**Status:** Current.

Originating ACES application: 1-1 "Windspeed Adjustment and Wave Growth" (functional area: Wave Prediction). Adjusts an observed wind to a 10-m equivalent-neutral wind, linearizes it to the constant-drag wind used by the growth formulas, then estimates fetch-/duration-limited wave growth (deep or shallow water, open or restricted fetch).

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Elevation of observed wind | z_obs | ft / m | 1 to 32808.4 | 25 | >= 1 ft |
| Observed wind speed | U_obs | mph / km/h | 0.1 to 447.387 | 45 | > 0 |
| Air-sea temperature difference | deltaT | C | -50 to 50 | 0 | T_air - T_sea; 0 = neutral. <0 unstable, >0 stable |
| Stability model | stability | (none) | choices: Neutral, Businger-Dyer (ACES) | Neutral | Neutral reproduces the User's Guide worked examples. Businger-Dyer is what the ACES executable's own output prefers by a wide margin over the 1,017-case DOS sweep (median error 1.9% vs 10.2% in U_e) -- the two oracles disagree; see the docstring. |
| Duration of observed wind | dur_obs | hr | 0.000277778 to 10 | 3 | >= 0.1 (stored internally in seconds) |
| Duration of final wind | dur_final | hr | 0.000277778 to 10 | 3 | >= 0.1 (stored internally in seconds) |
| Latitude of observation | lat | deg | 0 to 180 | 30 | 0 to 180 deg |
| Wind observation type | obs_type | (none) | choices: Overwater (shipboard), Overwater (not shipboard), Shore (windward), Shore (leeward) | Shore (windward) | how/where the wind was observed |
| Fetch type | fetch_type | (none) | choices: Open Water, Restricted | Open Water | Open Water (single fetch) or Restricted (radial fetches) |
| Wave equation type | wave_eq | (none) | choices: Shallow, Deep | Shallow | Shallow (depth-limited) or Deep water |
| Wind fetch length | F | mi / km | 6.21371e-07 to 3106.86 | 26 | open water: single fetch length |
| Average fetch depth | d | ft / m | 0.1 to 32808.4 | 13 | > 0 (used for shallow water) |
| Wind direction | wind_dir | deg | 0 to 360 | 125 | restricted: deg clockwise from north |
| Radial angle increment | dbeta | deg | 1 to 180 | 12 | restricted: spacing of radial fetches |
| Direction of first radial | beta1 | deg | 0 to 360 | 0 | restricted: deg clockwise from north |
| Radial fetch lengths | fetches | (none) | table: Fetch | 14 default rows | restricted: one fetch length per radial (starting at beta1, step dbeta) |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Equivalent neutral wind speed | U_e | mph / km/h | scalar |
| Adjusted wind speed | U_a | mph / km/h | scalar |
| Effective wind fetch | F_eff | mi / km | scalar |
| Mean wave direction | wave_dir | deg | scalar |
| Wave height | H_mo | ft / m | scalar |
| Peak wave period | T_p | s | scalar |

*Reference:* SPM (1984) Ch.3; Resio & Vincent (1977); Smith (1991); TR 1-1

*Module:* `backend/applications/chessqc_1_1_windspeed_wave_growth.py`

### 1-2 — Beta-Rayleigh Distribution  `[I]`

**Status:** Current.

Originating ACES application: 1-2 "Beta-Rayleigh Distribution" (functional area: Wave Prediction). Given an energy-based significant wave height, a peak period, and a depth, it returns the characteristic individual wave heights of the sea state (root-mean-square, median, and the means of the highest third, tenth, and hundredth) together with the probability-density curve.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Energy-based wave height (Hmo) | Hmo | ft / m | 0.000328084 to 3280.84 | 5 | > 0 (zero-moment / significant height of the sea state) |
| Peak spectral period (Tp) | Tp | s | 0.01 to 1000 | 6.3 | > 0 |
| Water depth | d | ft / m | 0.000328084 to 32808.4 | 10.2 | > 0; the distribution reverts to Rayleigh when d/(g Tp^2) >= 0.01 |
| Breaking-height coefficient Hb/d | Hb_coef | (none) | choices: 1.0 (ACES), 0.9, 0.78 (SPM) | 1.0 (ACES) | upper bound that truncates the distribution, as a fraction of depth; ACES BETAR uses the depth itself |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Root-mean-square height | Hrms | ft / m | scalar |
| Median height | Hmed | ft / m | scalar |
| Mean of highest 1/3 (H1/3) | H13 | ft / m | scalar |
| Mean of highest 1/10 (H1/10) | H110 | ft / m | scalar |
| Mean of highest 1/100 (H1/100) | H1100 | ft / m | scalar |
| Breaking (upper-bound) height | Hb | ft / m | scalar |
| Root-mean-quad height (length^2) | Hrmq | ft^2 / m^2 | scalar |
| Beta-Rayleigh shape alpha | alpha | (none) | scalar |
| Beta-Rayleigh shape beta | beta | (none) | scalar |
| Relative depth d/(g Tp^2) | rel_depth | (none) | scalar |
| Distribution used | regime | (none) | scalar |
| Profile: wave height | profile_H | ft / m | profile |
| Profile: probability density | profile_pdf | 1/ft / 1/m | profile |

*Reference:* Hughes & Borgman (1987); Thompson & Vincent (1985); TR 1-2

*Module:* `backend/applications/chessqc_1_2_beta_rayleigh.py`

### 1-3 — Extremal Significant Wave Height Analysis  `[I]`

**Status:** Current.

Originating ACES application: 1-3 "Extremal Significant Wave Height Analysis" (functional area: Wave Prediction). Fits extremal probability distributions to a sample of storm significant wave heights and estimates design wave heights at specified return periods, with confidence intervals.

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Significant wave heights | heights | (none) | table: H_s | 15 default rows | one significant height per storm/event |
| Total number of events | N_T | (none) | 1 to 100000 | 20 | events during the record (>= number of heights) |
| Record length | K | yr | 0.001 to 10000 | 20 | years of record |
| Water depth | depth | ft / m | 0.001 to 100000 | 500 | for depth-limited breaking cap (inactive in deep water) |
| Confidence level | conf | (none) | choices: 80, 85, 90, 95, 99 | 90 | confidence interval, % |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Best-fit correlation | corr | (none) | scalar |
| H_s, return period 2 yr | Hs2 | ft / m | scalar |
| H_s, return period 5 yr | Hs5 | ft / m | scalar |
| H_s, return period 10 yr | Hs10 | ft / m | scalar |
| H_s, return period 25 yr | Hs25 | ft / m | scalar |
| H_s, return period 50 yr | Hs50 | ft / m | scalar |
| H_s, return period 100 yr | Hs100 | ft / m | scalar |
| Design (100 yr) lower bound | Hs100_lo | ft / m | scalar |
| Design (100 yr) upper bound | Hs100_hi | ft / m | scalar |

*Reference:* Goda (1988); Gringorten (1963); EM 1110-2-1414; TR 1-3

*Module:* `backend/applications/chessqc_1_3_extremal_hs.py`

### 1-4 — Constituent Tide Record Generation  `[I]`

**Status:** Current.

Originating ACES application: 1-4 "Constituent Tide Record Generation" (functional area: Wave Prediction). Predicts a water-level time series at a gage from the harmonic tidal constituents (amplitude and epoch per constituent) by the classical harmonic method.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Start year | year | (none) | 1800 to 2200 | 1989 |  |
| Start month | month | (none) | 1 to 12 | 1 |  |
| Start day | day | (none) | 1 to 31 | 10 |  |
| Start hour (local) | hour | hr | 0 to 24 | 10 | hour of day at the gage (0 to 24) |
| Record length | length_hr | hr | 0.1 to 100000 | 120 |  |
| Output interval | interval_min | min | 0.1 to 1440 | 15 |  |
| Mean water level above datum | H0 | ft / m | -3280.84 to 3280.84 | 1.79 | datum offset (e.g. MLLW) |
| Gage longitude | gage_lon | deg | -180 to 360 | 70.62 | degrees West (positive) |
| Constituents | constituents | (none) | table: Table A-5 #, Amplitude, Epoch | 25 default rows | one row per constituent: Table A-5 index (1=M2,2=S2,3=N2,4=K1,5=M4,6=O1,...), amplitude, epoch in degrees. See CANON / module docstring for the full index list. |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Maximum elevation | h_max | ft / m | scalar |
| Minimum elevation | h_min | ft / m | scalar |
| Tidal range (max - min) | range | ft / m | scalar |
| Elevation at start | h_start | ft / m | point |
| Number of active constituents | n_constituents | (none) | scalar |
| Profile: time | profile_t | hr | profile |
| Profile: tide elevation | profile_h | ft / m | profile |

*Reference:* Schureman (1971); Table A-5; TR 1-4

*Module:* `backend/applications/chessqc_1_4_tide_record.py`

### 1-5 — Near-surface Wind Speeds  `[II]`

**Status:** Current.

Originating ACES application: 1-5 "Near-surface Wind Speeds" (functional area: Wave Prediction; a later ACES addition). Given the geostrophic wind speed, the air-sea temperature difference, latitude, and a height z, it returns the friction velocity, the wind speed at height z, the drag coefficients (at z and at 10 m), the surface roughness length, the Monin-Obukhov stability length, the stability function, and the surface momentum flux.

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Geostrophic wind speed | Ug | kt / m/s | 1 to 120 | 30 | free-atmosphere geostrophic wind \|V_g\| > 0 |
| Air-sea temperature difference | deltaT | deg C | -20 to 20 | 0 | dT = T_air - T_sea; <0 unstable (warm sea), >0 stable |
| Latitude | lat | deg | 1 to 80 | 40 | for the Coriolis parameter f = 2 Omega sin(lat) > 0 |
| Height above surface | z | ft / m | 0.1 to 300 | 10 | elevation z at which the wind speed is reported |
| Water type | water | (none) | choices: Salt, Fresh | Salt | affects air density only weakly; retained for parity with ACES |
| Air density | rho_air | kg/m^3 | 1 to 1.3 | 1.2 | standard value may be changed |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Friction velocity U* | u_star | kt / m/s | scalar |
| Wind speed at height z | Uz | kt / m/s | scalar |
| Wind speed at 10 m | U10 | kt / m/s | scalar |
| Drag coefficient at z | CDz | (none) | scalar |
| Drag coefficient at 10 m | CD | (none) | scalar |
| Surface roughness length | z0 | ft / m | scalar |
| Monin-Obukhov length | L | ft / m | scalar |
| Stability function Psi(z/L) | psi | (none) | scalar |
| Surface momentum flux | tau | Pa | scalar |
| Cross-isobar angle | alpha | deg | scalar |

*Reference:* ACES TR 1-1 eqs. 5-19; Silva (2005) eqs. 7.16-7.19; Lumley & Panofsky (1964)

*Module:* `backend/applications/chessqc_1_5_near_surface_wind.py`

### 1-6 — Holland Hurricane Wind Model  `[I]`

**Status:** Current.

Originating ACES application: 1-6 "Holland Hurricane Wind Model" (functional area: Wave Prediction; a later ACES addition). Given the two Holland (1980) profile parameters and the storm's central / peripheral pressures, it reconstructs the radial pressure profile and the gradient- and cyclostrophic-wind profiles of a tropical cyclone, and reports the maximum wind speed.

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Solve for | solve_for | (none) | choices: R_max (from A, B), A (from R_max, B), B (from R_max, A) | R_max (from A, B) | fix any two of {A, B, R_max}; the third is computed |
| Scaling parameter A | A | m^B | 0.001 to 1e+12 | 5.19615e+06 | length-scaling parameter; R_max = A^(1/B) |
| Peakedness parameter B | B | (none) | 0.5 to 2.5 | 1.5 | shape factor, typically ~1-2.5 |
| Radius of maximum wind | R_max | nmi / km | 1 to 200 | 30 | R_max = A^(1/B); used when solving for A or B |
| Central pressure | pc | hPa | 850 to 1010 | 940 | storm central pressure p_c |
| Peripheral pressure | pn | hPa | 990 to 1030 | 1013 | ambient pressure p_n; dP = p_n - p_c > 0 |
| Latitude | lat | deg | 0 to 80 | 20 | for the Coriolis term in the gradient-wind balance |
| Air density | rho_air | kg/m^3 | 1 to 1.3 | 1.15 | ambient air density (standard value may be changed) |
| Maximum plot radius | r_plot | nmi / km | 5 to 1000 | 200 | radial extent of the output profiles |
| Profile points | n_points | (none) | 20 to 2000 | 200 | number of radial samples for the profiles |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Maximum wind speed (gradient) | U_max | kt / km/h | scalar |
| Radius of maximum wind | r_at_max | nmi / km | scalar |
| Maximum wind (cyclostrophic) | U_max_cyclo | kt / km/h | scalar |
| R_max (computed/echoed) | R_max_out | nmi / km | scalar |
| A (computed/echoed) | A_out | m^B | scalar |
| B (computed/echoed) | B_out | (none) | scalar |
| Pressure deficit | dP | hPa | scalar |
| Profile: radial distance | profile_r | nmi / km | profile |
| Profile: pressure | profile_p | hPa | profile |
| Profile: gradient wind | profile_Vgr | kt / km/h | profile |
| Profile: cyclostrophic wind | profile_Vc | kt / km/h | profile |

*Reference:* Holland (1980) Mon. Wea. Rev. 108; ACES manual

*Module:* `backend/applications/chessqc_1_6_holland_hurricane.py`


## Wave Theory

### 2-1 — Linear Wave Theory  `[I]`

**Status:** Current.

Originating ACES application: 2-1 "Linear Wave Theory" (functional area: Wave Theory). First-order (small-amplitude / sinusoidal / Airy) approximations of wave motion.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Wave height | H | ft / m | 3.28084e-06 to 32808.4 | 6.3 | > 0 |
| Wave period | T | s | 0.001 to 10000 | 8 | > 0 |
| Water depth | d | ft / m | 3.28084e-06 to 328084 | 20 | > 0 |
| Vertical coordinate | z | ft / m | any | -12 | from SWL (z=0), +up; no restriction (clamped to [-d, eta]) |
| Wavelength fraction (X/L) | xL | (none) | 0 to 1 | 0.75 | 0.0 to 1.0 (phase position; 0 = crest) |
| Water type | water | (none) | choices: Salt, Fresh | Salt | sets water density used for energy/pressure |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Wave length | L | ft / m | scalar |
| Celerity | C | ft/s / m/s | scalar |
| Group velocity | Cg | ft/s / m/s | scalar |
| Energy density | E | lb/ft / N/m | scalar |
| Energy flux (power) | P | lb/s / N/s | scalar |
| Ursell parameter | Ur | (none) | scalar |
| Relative depth k*d | kd | (none) | scalar |
| Surface elevation | eta | ft / m | point |
| Pressure | p | psf / Pa | point |
| Horizontal displacement | xi | ft / m | point |
| Vertical displacement | zeta | ft / m | point |
| Horizontal velocity | u | ft/s / m/s | point |
| Vertical velocity | w | ft/s / m/s | point |
| Horizontal acceleration | dudt | ft/s^2 / m/s^2 | point |
| Vertical acceleration | dwdt | ft/s^2 / m/s^2 | point |
| Profile: X (+/- one wavelength) | profile_X | ft / m | profile |
| Profile: surface elevation | profile_eta | ft / m | profile |
| Profile: horizontal velocity | profile_u | ft/s / m/s | profile |
| Profile: vertical velocity | profile_w | ft/s / m/s | profile |

*Reference:* Airy (1845); Hunt (1979); TR 2-1

*Module:* `backend/applications/chessqc_2_1_linear_wave_theory.py`

### 2-2 — Cnoidal Wave Theory  `[I]`

**Status:** Current.

Originating ACES application: 2-2 "Cnoidal Wave Theory" (functional area: Wave Theory). Finite-amplitude periodic long-wave theory in terms of Jacobian elliptic functions; first-order (Isobe 1985) and second-order (Hardy & Kraus 1987) approximations.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Wave height | H | ft / m | 3.28084e-06 to 32808.4 | 10 | > 0 |
| Wave period | T | s | 0.001 to 10000 | 15 | > 0 |
| Water depth | d | ft / m | 3.28084e-06 to 328084 | 25 | > 0 |
| Vertical coordinate | z | ft / m | any | -12.5 | from SWL (z=0), +up; no restriction (clamped to [-d, eta]) |
| Wavelength fraction (X/L) | xL | (none) | 0 to 1 | 0.5 | 0.0 to 1.0 (phase position; 0 = crest) |
| Order | order | (none) | choices: 1, 2 | 1 | 1 = Isobe (1985); 2 = Hardy & Kraus (1987) |
| Water type | water | (none) | choices: Salt, Fresh | Salt | sets water density used for energy/pressure |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Wave length | L | ft / m | scalar |
| Celerity | C | ft/s / m/s | scalar |
| Energy density | E | lb/ft / N/m | scalar |
| Energy flux (power) | P | lb/s / N/s | scalar |
| Ursell parameter (HL^2/d^3) | Ur | (none) | scalar |
| Elliptic modulus | kappa | (none) | scalar |
| Surface elevation | eta | ft / m | point |
| Pressure | p | psf / Pa | point |
| Horizontal velocity | u | ft/s / m/s | point |
| Vertical velocity | w | ft/s / m/s | point |
| Horizontal acceleration | dudt | ft/s^2 / m/s^2 | point |
| Vertical acceleration | dwdt | ft/s^2 / m/s^2 | point |
| Profile: X (+/- one wavelength) | profile_X | ft / m | profile |
| Profile: surface elevation | profile_eta | ft / m | profile |
| Profile: horizontal velocity | profile_u | ft/s / m/s | profile |
| Profile: vertical velocity | profile_w | ft/s / m/s | profile |

*Reference:* Isobe (1985); Hardy & Kraus (1987); TR 2-2

*Module:* `backend/applications/chessqc_2_2_cnoidal_wave_theory.py`

### 2-3 — Fourier Series Wave Theory  `[I]`

**Status:** Current.

Originating ACES application: 2-3 "Fenton's Fourier Series Wave Theory" (functional area: Wave Theory). Steady progressive wave of permanent form, solved by an N-term stream-function Fourier series (Rienecker & Fenton 1981; Fenton 1988).

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Wave height | H | ft / m | 3.28084e-06 to 200 | 4.5 | > 0 |
| Wave period | T | s | 1 to 1000 | 9 | > 0 |
| Water depth | d | ft / m | 1 to 5000 | 22 | > 0 |
| Celerity definition | cdef | (none) | choices: Euler, Stokes | Euler | Euler = mean Eulerian current; Stokes = mean mass-transport velocity |
| Mean velocity | ubar | ft/s / m/s | -32.8084 to 32.8084 | 0 | current (Euler) or mass-transport (Stokes) velocity |
| Fourier terms | N | (none) | 1 to 25 | 16 | number of terms (1 to 25) |
| Height-ramp steps | nramp | (none) | 1 to 10 | 5 | wave-height ramping steps (1 to 10) |
| Vertical coordinate | z | ft / m | any | -5 | from SWL (z=0), +up; clamped to [-d, eta] |
| Wavelength fraction (X/L) | xL | (none) | 0 to 1 | 0 | 0.0 to 1.0 (phase position; 0 = crest) |
| Water type | water | (none) | choices: Salt, Fresh | Salt | sets water density used for energies/pressure |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Celerity | C | ft/s / m/s | scalar |
| Wave length | L | ft / m | scalar |
| Mean Eulerian fluid velocity | ubar1 | ft/s / m/s | scalar |
| Mean mass-transport velocity | ubar2 | ft/s / m/s | scalar |
| Volume flux (wave) | q | ft^2/s / m^2/s | scalar |
| Volume flux (mean) | Q | ft^2/s / m^2/s | scalar |
| Bernoulli constant | R | ft^2/s^2 / m^2/s^2 | scalar |
| Impulse | I | lb*s/ft^2 / N*s/m^2 | scalar |
| Kinetic energy | EK | lb/ft / N/m | scalar |
| Potential energy | EP | lb/ft / N/m | scalar |
| Energy density | Edens | lb/ft / N/m | scalar |
| Mean square of bed velocity | Ub2 | ft^2/s^2 / m^2/s^2 | scalar |
| Radiation stress | Sxx | lb/ft / N/m | scalar |
| Wave power (energy flux) | F | lb/s / N/s | scalar |
| Surface elevation | eta | ft / m | point |
| Horizontal velocity | U | ft/s / m/s | point |
| Vertical velocity | W | ft/s / m/s | point |
| Horizontal acceleration | ax | ft/s^2 / m/s^2 | point |
| Vertical acceleration | az | ft/s^2 / m/s^2 | point |
| Pressure | p | psf / Pa | point |
| Profile: X (+/- one wavelength) | profile_X | ft / m | profile |
| Profile: surface elevation | profile_eta | ft / m | profile |
| Profile: horizontal velocity | profile_U | ft/s / m/s | profile |
| Profile: vertical velocity | profile_W | ft/s / m/s | profile |

*Reference:* Rienecker & Fenton (1981); Fenton (1988); TR 2-3

*Module:* `backend/applications/chessqc_2_3_fenton_fourier.py`

### 2-4 — Wave Parameters  `[I]`

**Status:** Current.

Originating ACES grouping: 2-4 "Wave Parameters" (functional area: Wave Theory). A linear (Airy) wave-theory engine presented as the classic SPM dimensionless "wave parameter" table, extended with two field-engineering conveniences:

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Solve for | mode | (none) | choices: Forward (height given), Invert (pressure to height) | Forward (height given) | Invert: recover H from a measured dynamic pressure amplitude at z |
| Wave specified by | spec | (none) | choices: Period (s), Frequency (1/s) | Period (s) | interpret the value below as period T or frequency f = 1/T |
| Period or frequency | Tf | s or 1/s | 0.0001 to 10000 | 8 | > 0 (period in s, or frequency in 1/s per the choice) |
| Wave height | H | ft / m | 0 to 32808.4 | 6.3 | forward mode only; > 0 |
| Measured dynamic pressure amplitude | p_gauge | psf / Pa | 0 to 2.08854e+07 | 0 | invert mode only; amplitude (not incl. hydrostatic) of the wave pressure at z |
| Water depth | d | ft / m | 3.28084e-06 to 328084 | 20 | > 0 |
| Vertical coordinate / gauge elevation | z | ft / m | any | -12 | from SWL (z=0), +up; bed at -d. Also the pressure-gauge elevation in invert mode |
| Phase angle | theta_deg | deg | 0 to 360 | 270 | 0 deg = crest; 270 deg <=> X/L = 0.75 |
| Water type | water | (none) | choices: Salt, Fresh | Salt | sets water density used for energy / pressure |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Frequency | f | 1/s | scalar |
| Wave length | L | ft / m | scalar |
| Deepwater wave length | L0 | ft / m | scalar |
| Celerity | C | ft/s / m/s | scalar |
| Deepwater celerity | C0 | ft/s / m/s | scalar |
| Group velocity | Cg | ft/s / m/s | scalar |
| Deepwater group velocity | Cg0 | ft/s / m/s | scalar |
| Group/phase ratio n = Cg/C | n | (none) | scalar |
| Cg / C0 | Cg_C0 | (none) | scalar |
| Relative depth k*d | kd | (none) | scalar |
| d/L | d_L | (none) | scalar |
| d/L0 | d_L0 | (none) | scalar |
| tanh(kd) | tanh_kd | (none) | scalar |
| sinh(kd) | sinh_kd | (none) | scalar |
| cosh(kd) | cosh_kd | (none) | scalar |
| 1 / sinh(kd) | csch_kd | (none) | scalar |
| 1 / cosh(kd) | sech_kd | (none) | scalar |
| Shoaling coefficient H/H0 | Ks | (none) | scalar |
| Equivalent deepwater height H0 | H0 | ft / m | scalar |
| Steepness ratio (H/L)/(H0/L0) | steep_ratio | (none) | scalar |
| Pressure response factor at z | Kp_z | (none) | scalar |
| Pressure response factor at bed | Kp_bed | (none) | scalar |
| Energy density | E | lb/ft / N/m | scalar |
| Deepwater energy density | E0 | lb/ft / N/m | scalar |
| Kinetic energy per crest length | Ek | lb-ft/ft / N-m/m | scalar |
| Total energy per crest length | Etot | lb-ft/ft / N-m/m | scalar |
| Energy flux (wave power) | P | lb/s / N/s | scalar |
| Deepwater wave power | P0 | lb/s / N/s | scalar |
| Ursell parameter | Ur | (none) | scalar |
| Surface elevation | eta | ft / m | point |
| Pressure (total) | p | psf / Pa | point |
| Dynamic pressure amplitude at z | p_dyn | psf / Pa | point |
| Horizontal displacement | xi | ft / m | point |
| Vertical displacement | zeta | ft / m | point |
| Horizontal velocity | u | ft/s / m/s | point |
| Vertical velocity | w | ft/s / m/s | point |
| Horizontal acceleration | dudt | ft/s^2 / m/s^2 | point |
| Vertical acceleration | dwdt | ft/s^2 / m/s^2 | point |
| Mass-transport velocity at z | Us_z | ft/s / m/s | point |
| Mass-transport velocity at SWL | Us_surf | ft/s / m/s | point |
| Mass-transport velocity at bed | Us_bed | ft/s / m/s | point |
| Wave height used | H_used | ft / m | point |
| Profile: X (+/- one wavelength) | profile_X | ft / m | profile |
| Profile: surface elevation | profile_eta | ft / m | profile |
| Profile: horizontal velocity | profile_u | ft/s / m/s | profile |
| Profile: vertical velocity | profile_w | ft/s / m/s | profile |

*Reference:* Airy (1845); SPM (1984) App. C; Hunt (1979)

*Module:* `backend/applications/chessqc_2_4_wave_parameters.py`

### 2-5 — Solitary Wave Theory  `[I]`

**Status:** Current.

Originating ACES grouping: 2-5 "Solitary Wave Theory" (functional area: Wave Theory). A solitary wave is a single wave of translation lying entirely above the still-water level, with no trough; long waves such as tsunamis and surge-driven bores approximate it. The app returns the wave kinematics and integral properties for a wave of height H in depth d.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Wave height | H | ft / m | 0.000328084 to 3280.84 | 3 |  |
| Water depth | d | ft / m | 0.00328084 to 32808.4 | 10 |  |
| Vertical coordinate (from bottom) | z | ft / m | 0 to 32808.4 | 10 | 0 at the bed; up to the surface |
| Horizontal distance from crest | x | ft / m | -32808.4 to 32808.4 | 0 |  |
| Beach slope (tan theta) | m | (none) | 0 to 0.2 | 0.02 | for the empirical breaking criterion |
| Water | water | (none) | choices: Salt, Fresh | Salt |  |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Wave celerity | C | ft/s / m/s | scalar |
| Surface elevation at x (above SWL) | eta | ft / m | scalar |
| Horizontal particle velocity at (x,z) | u | ft/s / m/s | scalar |
| Vertical particle velocity at (x,z) | w | ft/s / m/s | scalar |
| Dynamic pressure at bed under crest | dp_crest | psf / Pa | scalar |
| Total energy per unit crest width | E | lb/ft / N | scalar |
| McCowan-Munk coefficient M | M | (none) | scalar |
| McCowan-Munk coefficient N | N | (none) | scalar |
| Breaking height (flat bed, McCowan) | Hb_flat | ft / m | scalar |
| Breaking height (sloping bed) | Hb_slope | ft / m | scalar |
| Relative height H/d | relative_height | (none) | scalar |

*Reference:* McCowan (1891); Munk (1949); SPM (1984); CEM (EM 1110-2-1100); Zaroodny (1972)

*Module:* `backend/applications/chessqc_2_5_solitary_wave.py`


## Wave Transformation

### 3-1 — Linear Wave Theory with Snell's Law  `[I]`

**Status:** Current.

Transforms a wave of known height/period/direction at one depth to deep water and to a "subject" depth, using linear wave theory (Hunt 1979 dispersion), Snell's law (O'Brien 1942) for refraction, and energy-flux conservation for shoaling. Reports wave height / crest angle / length / celerity / group velocity / energy density / energy flux at three locations (Known, Deep water, Subject), Ursell numbers, deepwater steepness, and the Weggel (1972) breaker height/depth.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Wave height (known) | H1 | ft / m | 3.28084e-06 to 32808.4 | 10 | > 0 (height at the known depth) |
| Wave period | T | s | 0.001 to 10000 | 7.5 | > 0 |
| Water depth (known) | d1 | ft / m | 3.28084e-06 to 328084 | 25 | > 0 |
| Wave crest angle (known) | alpha1 | deg | 0 to 90 | 10 | angle between wave crest and depth contour, 0-90 deg |
| Cotan of nearshore slope | cot_phi | (none) | 1e-06 to 1e+06 | 100 | cot(beach slope); e.g. 100 = 1:100 slope |
| Water depth (subject) | d2 | ft / m | 3.28084e-06 to 328084 | 20 | > 0 (depth to transform the wave to) |
| Water type | water | (none) | choices: Salt, Fresh | Salt | sets water density used for energy density/flux |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Wave height (known) | H_k | ft / m | scalar |
| Wave crest angle (known) | a_k | deg | scalar |
| Wavelength (known) | L_k | ft / m | scalar |
| Wave celerity (known) | C_k | ft/s / m/s | scalar |
| Group velocity (known) | Cg_k | ft/s / m/s | scalar |
| Energy density (known) | E_k | lb/ft / N/m | scalar |
| Energy flux (known) | P_k | lb/s / N/s | scalar |
| Ursell number (known) | Ur_k | (none) | scalar |
| Wave height (deep) | H_0 | ft / m | scalar |
| Wave crest angle (deep) | a_0 | deg | scalar |
| Wavelength (deep) | L_0 | ft / m | scalar |
| Wave celerity (deep) | C_0 | ft/s / m/s | scalar |
| Group velocity (deep) | Cg_0 | ft/s / m/s | scalar |
| Energy density (deep) | E_0 | lb/ft / N/m | scalar |
| Energy flux (deep) | P_0 | lb/s / N/s | scalar |
| Deepwater wave steepness | steep_0 | (none) | scalar |
| Wave height (subject) | H_2 | ft / m | scalar |
| Wave crest angle (subject) | a_2 | deg | scalar |
| Wavelength (subject) | L_2 | ft / m | scalar |
| Wave celerity (subject) | C_2 | ft/s / m/s | scalar |
| Group velocity (subject) | Cg_2 | ft/s / m/s | scalar |
| Energy density (subject) | E_2 | lb/ft / N/m | scalar |
| Energy flux (subject) | P_2 | lb/s / N/s | scalar |
| Ursell number (subject) | Ur_2 | (none) | scalar |
| Breaker height | Hb | ft / m | scalar |
| Breaker depth | db | ft / m | scalar |

*Reference:* O'Brien (1942); Hunt (1979); Weggel (1972); TR 3-1

*Module:* `backend/applications/chessqc_3_1_snell.py`

### 3-2 — Irregular Wave Transformation (Goda's Method)  `[II]`

**Status:** Current.

Originating ACES grouping: 3-2 "Irregular Wave Transformation (Goda's method)" (functional area: Wave Transformation). Transforms an irregular (spectral) deepwater sea state to a nearshore depth over straight, parallel bottom contours, accounting for refraction, shoaling, and depth-limited breaking, and reports the transformed wave-height statistics plus shoaling and effective-refraction coefficients, surf beat, and wave setup.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Significant deepwater wave height | H0 | ft / m | 0.000328084 to 3280.84 | 20 |  |
| Water depth | d | ft / m | 9.97375 to 32808.4 | 50 | min ~10 ft / 3.04 m |
| Significant wave period | Ts | s | 1 to 16 | 8 | <= 16 s |
| Cotangent of nearshore slope | cot_phi | (none) | 1 to 10000 | 100 |  |
| Principal incident direction | theta | deg | -75 to 75 | 10 | from shore normal; \|theta\| <= 75 deg |
| Directional spreading | s_max | (none) | choices: Auto (ACES), 10 (wind waves), 25 (steep swell), 75 (flat swell) | Auto (ACES) | ACES selects the band from period and steepness (T_s <= 10 s -> 10; longer with H/L > 0.02 -> 25; otherwise 75). Override to pin it. |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Significant wave height (at depth) | Hs | ft / m | scalar |
| Mean wave height | Hmean | ft / m | scalar |
| Root-mean-square wave height | Hrms | ft / m | scalar |
| Average of highest 1/10 | H10 | ft / m | scalar |
| Average of highest 2% | H2 | ft / m | scalar |
| Maximum wave height | Hmax | ft / m | scalar |
| Shoaling coefficient | Ks | (none) | scalar |
| Effective refraction coefficient | Kr | (none) | scalar |
| RMS surf beat | surf_beat | ft / m | scalar |
| Wave setup at depth | setup | ft / m | scalar |
| Deepwater steepness H0/L0 | steepness | (none) | scalar |

*Reference:* Goda (1975, 1985); Mitsuyasu (1975); Shuto (1974)

*Module:* `backend/applications/chessqc_3_2_goda_transformation.py`

### 3-3 — Combined Diffraction and Reflection by a Vertical Wedge  `[I]`

**Status:** Current.

Originating ACES grouping: 3-3 "Combined Diffraction and Reflection by a Vertical Wedge" (functional area: Wave Transformation). Computes the wave-height modification factor and phase at a point near a fully-reflecting vertical wedge (a breakwater tip, a corner, or a semi-infinite breakwater), where an incident monochromatic wave is simultaneously diffracted and reflected. This is the PCDFRAC solver (Chen 1987).

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Incident wave height | Hi | ft / m | 0.000328084 to 3280.84 | 2 |  |
| Wave period | T | s | 0.01 to 1000 | 6 |  |
| Water depth | d | ft / m | 0.00328084 to 32808.4 | 12 |  |
| Incident wave angle | alpha | deg | 0 to 360 | 133 |  |
| Wedge angle | wedge_angle | deg | 0 to 180 | 0 | 0 = semi-infinite breakwater; solid wedge angle |
| X coordinate of point | X | ft / m | -32808.4 to 32808.4 | 33 |  |
| Y coordinate of point | Y | ft / m | -32808.4 to 32808.4 | -17 |  |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Wave length | L | ft / m | scalar |
| Modification factor \|phi\| | mod_factor | (none) | scalar |
| Wave phase | phase | rad | scalar |
| Modified wave height | H | ft / m | scalar |

*Reference:* Chen (1987); Stoker (1957); Penny & Price (1952)

*Module:* `backend/applications/chessqc_3_3_wedge_diffraction.py`

### 3-4 — Vertical-Wedge Diffraction/Reflection on a Uniform Grid  `[I]`

**Status:** Current.

Originating ACES grouping: 3-4 "Combined Diffraction and Reflection by a Vertical Wedge Evaluated Upon a Uniform Grid" (functional area: Wave Transformation). This is the gridded companion to 3-3: it evaluates the same fully-reflecting-wedge solution (the PCDFRAC / Chen 1987 eigenfunction expansion) at every node of a uniform X-Y grid around the wedge apex, producing maps of the wave-height modification factor, the modified wave height, and the phase relative to the incident wave.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Incident wave height | Hi | ft / m | 0.000328084 to 3280.84 | 4 |  |
| Wave period | T | s | 0.01 to 1000 | 12 |  |
| Water depth | d | ft / m | 0.00328084 to 32808.4 | 30 |  |
| Incident wave angle | alpha | deg | 0 to 360 | 52 |  |
| Wedge angle | wedge_angle | deg | 0 to 180 | 0 |  |
| X start | X0 | ft / m | -32808.4 to 32808.4 | -600 |  |
| X end | Xm | ft / m | -32808.4 to 32808.4 | 200 |  |
| X increment | dX | ft / m | 0.00328084 to 32808.4 | 200 |  |
| Y start | Y0 | ft / m | -32808.4 to 32808.4 | -400 |  |
| Y end | Ym | ft / m | -32808.4 to 32808.4 | 200 |  |
| Y increment | dY | ft / m | 0.00328084 to 32808.4 | 100 |  |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Wave length | L | ft / m | scalar |
| Maximum modification factor | mod_max | (none) | scalar |
| Minimum modification factor | mod_min | (none) | scalar |
| Maximum modified wave height | H_max | ft / m | scalar |
| Grid X coordinates | grid_x | ft / m | profile |
| Grid Y coordinates | grid_y | ft / m | profile |
| Modification-factor field | mod_grid | (none) | grid |
| Modified-height field | H_grid | ft / m | grid |

*Reference:* Chen (1987); Stoker (1957); Penny & Price (1952)

*Module:* `backend/applications/chessqc_3_4_wedge_grid.py`


## Structural Design

### 4-1 — Breakwater Design (Hudson)  `[I]`

**Status:** Current. Newer method: Van der Meer (1988) stability formulae (preferred for many cases).

Sizes the primary armor units of a rubble-mound breakwater/revetment from the Hudson (1953-61) stability equation, and reports crest width, cover-layer thickness, and armor-unit placement density (SPM 1984 Ch. 7; EM 1110-2-2904).

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Type of armor unit | armor_type | (none) | choices: Quarrystone (smooth, rounded), Quarrystone (rough, angular), Graded riprap, Tribar (trunk, nonbreaking), Tribar (trunk, breaking), Tetrapod, Quadripod, Dolos, Modified cube, Hexapod, Toskane, Other | Tribar (trunk, nonbreaking) | optional/informational; pick K_D, k_delta, P from SPM tables accordingly |
| Sizing method | method | (none) | choices: Hudson, Van der Meer | Hudson | Hudson (ACES; reproduces Example 4-1) or Van der Meer (1988) rock-armor stability |
| Armor unit weight | w_r | lb/ft^3 / kN/m^3 | 0.00636588 to 6365.88 | 165 | unit weight of armor material; must exceed water unit weight |
| Wave height | H | ft / m | 3.28084e-06 to 32808.4 | 11.5 | design wave height (H or H_i) |
| Water unit weight | w_w | lb/ft^3 / kN/m^3 | 0.00636588 to 6365.88 | 64 | 64 lb/ft^3 seawater, 62.4 lb/ft^3 fresh |
| Stability coefficient | K_D | (none) | 0.001 to 10000 | 10 | K_D from SPM Table 7-8 (depends on armor type / slope / wave condition) |
| Layer coefficient | k_delta | (none) | 0.001 to 10 | 1.02 | layer coefficient k_delta (SPM Table 7-13) |
| Average porosity of armor layer | P | % | 0 to 99 | 54 | cover-layer porosity, percent (SPM Table 7-13) |
| Cotangent of structure slope | cot_theta | (none) | 0.001 to 1000 | 2 | cot(theta); theta = seaward slope angle |
| Number of armor units (layer thickness) | n | (none) | 1 to 10 | 2 | number of armor-unit layers (>= 2 typical) |
| Mean wave period (Van der Meer) | Tm | s | 0.01 to 1000 | 8 | Van der Meer only: mean period T_m for the surf-similarity parameter |
| Number of waves (Van der Meer) | N_waves | (none) | 1 to 1e+06 | 7500 | Van der Meer only: storm duration in waves (typ. <= 7500) |
| Notional permeability P (Van der Meer) | perm | (none) | 0.1 to 0.6 | 0.4 | Van der Meer only: 0.1 impermeable core ... 0.5-0.6 homogeneous mound |
| Damage level S (Van der Meer) | S_damage | (none) | 1 to 30 | 2 | Van der Meer only: 2 = start of damage; higher = more allowed damage |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Weight of individual armor unit | W | tons / kN | scalar |
| Crest width of breakwater | B | ft / m | scalar |
| Average cover layer thickness | r | ft / m | scalar |
| Armor units per 1000 ft^2 | N_r | (none) | scalar |
| Stability number Hs/(d_n50) | Ns | (none) | scalar |

*Reference:* Hudson (1953-61); SPM (1984) Ch.7; EM 1110-2-2904; TR 4-1

*Module:* `backend/applications/chessqc_4_1_breakwater_hudson.py`

### 4-2 — Toe Protection Design  `[I]`

**Status:** Current.

Designs the toe-apron width and toe-stone weight for a vertical wall / bulkhead / revetment. Apron width is the larger of a geotechnical (Rankine passive) width and two hydraulic minima; toe-stone weight uses the Tanimoto, Yagyu & Goda (1982) stability number N_s in a Hudson-form sizing equation.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Incident wave height | H_i | ft / m | 3.28084e-06 to 32808.4 | 5 | > 0 (incident/design wave height at the structure) |
| Wave period | T | s | 0.001 to 10000 | 12 | > 0 |
| Water depth at structure | d_s | ft / m | 3.28084e-06 to 328084 | 20 | > 0 |
| Cotangent of nearshore slope | cot_phi | (none) | 1e-06 to 1e+06 | 100 | cot(beach slope); collected for context (not used in toe sizing) |
| Passive earth-pressure coefficient | K_p | (none) | 0 to 100 | 1.5 | Rankine passive coefficient K_p (0 if no geotechnical width) |
| Sheet-pile penetration depth | d_e | ft / m | 0 to 32808.4 | 10 | >= 0 (0 if no sheet pile) |
| Height of toe layer above mudline | h_b | ft / m | 0 to 32808.4 | 4.5 | 0 <= h_b < d_s; d1 = d_s - h_b |
| Unit weight of rock | w_r | lb/ft^3 / kN/m^3 | 0.00636588 to 6365.88 | 165 | unit weight of toe-stone; must exceed water unit weight |
| Water unit weight | w_w | lb/ft^3 / kN/m^3 | 0.00636588 to 6365.88 | 64 | 64 lb/ft^3 seawater, 62.4 lb/ft^3 fresh |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Width of toe protection apron | B | ft / m | scalar |
| Weight of individual armor unit | W | lb / N | scalar |
| Water depth at top of toe layer | d_1 | ft / m | scalar |
| Stability number (Tanimoto-Yagyu-Goda) | N_s | (none) | scalar |

*Reference:* EM 1110-2-1614; Tanimoto, Yagyu & Goda (1982); Hunt (1979); TR 4-2

*Module:* `backend/applications/chessqc_4_2_toe_protection.py`

### 4-3 — Nonbreaking Wave Forces at Vertical Walls  `[I]`

**Status:** Current.

Standing-wave (clapotis) forces and overturning moments on a vertical wall, by two methods: Sainflou (1928) and Miche-Rundgren [Miche (1944), Rundgren (1958)], for the crest and the trough at the wall. The wall surface elevation (height above bottom) is the Miche-Rundgren second-order result (used for both methods per the TR); force and moment come from numerically integrating the Lagrangian pressure over 90 increments, weighting each at-rest increment by its stretched (elevated) thickness.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Depth at SWL | d | ft / m | 3.28084e-06 to 328084 | 15 | still-water depth at the wall |
| Incident wave height | H_i | ft / m | 3.28084e-06 to 32808.4 | 8 | incident wave height (the clapotis reaches ~(1+chi)*H_i at the wall) |
| Wave period | T | s | 0.001 to 10000 | 10 | > 0 |
| Wave reflection coefficient | chi | (none) | 0 to 1 | 1 | 1.0 = full reflection (smooth wall); do not use < 0.9 for design |
| Cotangent of nearshore slope | cot_phi | (none) | 1e-06 to 1e+06 | 100 | collected for context (not used in the wall-force calculation) |
| Water unit weight | gamma | lb/ft^3 / kN/m^3 | 0.00636588 to 6365.88 | 64 | 64 lb/ft^3 seawater, 62.4 lb/ft^3 fresh |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Crest height above bottom | hgt_crest | ft / m | scalar |
| Trough height above bottom | hgt_trough | ft / m | scalar |
| Miche-Rundgren crest force | mr_F_crest | lb/ft / N/m | scalar |
| Miche-Rundgren crest moment | mr_M_crest | lb-ft/ft / N-m/m | scalar |
| Miche-Rundgren trough force | mr_F_trough | lb/ft / N/m | scalar |
| Miche-Rundgren trough moment | mr_M_trough | lb-ft/ft / N-m/m | scalar |
| Sainflou crest force | sf_F_crest | lb/ft / N/m | scalar |
| Sainflou crest moment | sf_M_crest | lb-ft/ft / N-m/m | scalar |
| Sainflou trough force | sf_F_trough | lb/ft / N/m | scalar |
| Sainflou trough moment | sf_M_trough | lb-ft/ft / N-m/m | scalar |

*Reference:* Sainflou (1928); Miche (1944); Rundgren (1958); Hunt (1979); TR 4-3

*Module:* `backend/applications/chessqc_4_3_vertical_wall_forces.py`

### 4-4 — Rubble-Mound Revetment Design  `[I]`

**Status:** Current.

Originating ACES grouping: 4-4 "Rubble-Mound Revetment Design" (functional area: Structural Design). Sizes the armor and filter (bedding) stone for a riprap revetment under irregular waves, reporting median stone weight and the full weight/size gradation of each layer, the layer thicknesses, and the expected and conservative wave runup.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Significant wave height | Hs | ft / m | 0.000328084 to 3280.84 | 5 |  |
| Significant wave period | Ts | s | 0.01 to 1000 | 10 |  |
| Water depth at toe of revetment | ds | ft / m | 0.00328084 to 32808.4 | 9 |  |
| Cotangent of structure slope | cot_theta | (none) | 1 to 10 | 2 |  |
| Unit weight of rock | wr | lb/ft^3 / N/m^3 | 6.36588 to 318.294 | 165 |  |
| Water unit weight | w_w | lb/ft^3 / N/m^3 | 6.36588 to 318.294 | 64 | 64 lb/ft^3 seawater, 62.4 lb/ft^3 fresh (ACES RUBBLE.FOR:949) |
| Permeability coefficient | P | (none) | 0.1 to 0.6 | 0.1 | 0.1 impermeable core, 0.4-0.5 permeable, 0.6 homogeneous (Fig 4-4-2) |
| Damage level | S | (none) | 1 to 20 | 2 | van der Meer damage S (Table 4-4-1) |
| Runup method | runup_method | (none) | choices: Ahrens-Heimbaugh, EurOtop | Ahrens-Heimbaugh | Ahrens & Heimbaugh 1988 (ACES) or EurOtop 2018 Ru2% (modern standard) |
| Roughness factor (EurOtop) | gamma_f | (none) | 0.3 to 1 | 0.55 | EurOtop only: ~0.40 permeable rock, 0.55 impermeable rock, 1.0 smooth |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Governing stability number | N_s | (none) | scalar |
| Median armor weight | W50 | lb / N | scalar |
| Median armor dimension | D50 | ft / m | scalar |
| Armor-layer thickness | r_armor | ft / m | scalar |
| Filter-layer thickness | r_filter | ft / m | scalar |
| Armor W15 | W15 | lb / N | scalar |
| Armor W85 | W85 | lb / N | scalar |
| Armor Wmax | Wmax | lb / N | scalar |
| Armor Wmin | Wmin | lb / N | scalar |
| Expected maximum runup | R_expected | ft / m | scalar |
| Conservative runup | R_conservative | ft / m | scalar |

*Reference:* Ahrens (1981); van der Meer (1988); Hudson (1958); Ahrens & Heimbaugh (1988)

*Module:* `backend/applications/chessqc_4_4_revetment_design.py`


## Wave Runup, Transmission, and Overtopping

### 5-1 — Irregular Wave Runup on Beaches  `[I]`

**Status:** Current.

Originating ACES grouping: 5-1 "Irregular Wave Runup on Beaches" (functional area: Wave Runup, Transmission, and Overtopping). Estimates how high irregular (real, many-period) waves run up a smooth, uniform, impermeable beach slope, reporting five runup statistics: the maximum, the level exceeded by 2 percent of runups, and the averages of the highest tenth, highest third, and all runups.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Deepwater significant wave height | Hs0 | ft / m | 0.000328084 to 3280.84 | 4.6 | > 0 (H_s0) |
| Peak energy wave period | Tp | s | 0.01 to 1000 | 9.5 | > 0 |
| Cotangent of foreshore slope | cot_theta | (none) | 0.001 to 10000 | 13 | cot(theta) > 0; tan(theta) = 1/cot(theta) |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Deepwater wave length | L0 | ft / m | scalar |
| Iribarren (surf-similarity) number | xi | (none) | scalar |
| Maximum runup | R_max | ft / m | scalar |
| Runup exceeded by 2% of runups | R_2 | ft / m | scalar |
| Average of highest 1/10 of runups | R_1_10 | ft / m | scalar |
| Average of highest 1/3 of runups | R_1_3 | ft / m | scalar |
| Average runup | R_mean | ft / m | scalar |

*Reference:* Mase (1989); Hunt (1959); Walton & Ahrens (1989)

*Module:* `backend/applications/chessqc_5_1_irregular_runup_beaches.py`

### 5-2 — Wave Runup and Overtopping on Impermeable Structures  `[I]`

**Status:** Current.

Originating ACES grouping: 5-2 "Wave Runup and Overtopping on Impermeable Structures" (functional area: Wave Runup, Transmission, and Overtopping). Computes wave runup on a smooth or rough impermeable slope (a seawall or revetment face) and, when the runup exceeds the crest, the resulting overtopping rate, for both monochromatic and irregular waves, with an optional onshore-wind correction.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Incident wave height at toe | Hi | ft / m | 0.000328084 to 3280.84 | 7.5 | significant height for irregular waves |
| Wave period | T | s | 0.01 to 1000 | 10 |  |
| Water depth at structure toe | ds | ft / m | 0.00328084 to 32808.4 | 12.5 |  |
| Cotangent of structure slope | cot_theta | (none) | 0.001 to 1000 | 3 |  |
| Structure height above toe | hs | ft / m | 0.00328084 to 32808.4 | 20 |  |
| Wave type | wave_type | (none) | choices: Monochromatic, Irregular | Monochromatic |  |
| Slope type | slope_type | (none) | choices: Rough (riprap), Smooth | Rough (riprap) |  |
| Compute overtopping | want_overtopping | (none) | yes / no | yes |  |
| Known runup (0 = compute) | R_known | ft / m | 0 to 32808.4 | 0 | if > 0, used directly instead of the runup formula |
| Rough-slope coefficient a | a | (none) | 0 to 10 | 0.956 | Ahrens & McCartney; per armor type (Appendix A) |
| Rough-slope coefficient b | b | (none) | 0 to 10 | 0.398 |  |
| Overtopping coefficient alpha | alpha | (none) | 0.0001 to 1 | 0.076463 | SPM figures; or set alpha_from_slope |
| Overtopping coefficient Q*0 | Qstar0 | (none) | 0 to 10 | 0.025 |  |
| Use alpha = 0.06 - 0.01431 sin(theta) | alpha_from_slope | (none) | yes / no | no |  |
| Onshore wind velocity | U | kt / km/h | 0 to 388.769 | 35 | 0 = no wind correction |
| Refraction coefficient | KR | (none) | 0 to 1 | 1 | H'0 = KR * H0 |
| Overtopping method | overtopping_method | (none) | choices: Weggel, EurOtop | Weggel | Weggel 1976 (ACES) or EurOtop 2018 mean discharge (modern standard) |
| Roughness factor (EurOtop) | gamma_f | (none) | 0.3 to 1 | 0.55 | EurOtop only: ~0.55 rough rock, 1.0 smooth |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Deepwater wave height | H0 | ft / m | scalar |
| Relative height ds/H0 | ds_H0 | (none) | scalar |
| Wave steepness H0/gT^2 | steepness | (none) | scalar |
| Deepwater wave length | L0 | ft / m | scalar |
| Surf-similarity (Iribarren) number | xi | (none) | scalar |
| Wave runup | R | ft / m | scalar |
| Crest freeboard (hs - ds) | F | ft / m | scalar |
| Wind correction factor | Cw | (none) | scalar |
| Overtopping rate per unit length | Q | ft^3/s/ft / m^3/s/m | scalar |

*Reference:* Ahrens & McCartney (1975); Ahrens & Titus (1985); Weggel (1976); SPM (1984)

*Module:* `backend/applications/chessqc_5_2_runup_overtopping_impermeable.py`

### 5-3 — Wave Transmission on Impermeable Structures  `[I]`

**Status:** Current.

Originating ACES grouping: 5-3 "Wave Transmission on Impermeable Structures" (functional area: Wave Runup, Transmission, and Overtopping). Estimates the height of the wave transmitted past an impermeable structure, either a sloped structure overtopped by waves or a vertical/composite (caisson-on-berm) structure, as a transmission coefficient K_T = H_T / H_i.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Incident wave height | Hi | ft / m | 0.000328084 to 3280.84 | 7.5 |  |
| Wave period | T | s | 0.01 to 1000 | 10 |  |
| Water depth at structure toe | ds | ft / m | 0.00328084 to 32808.4 | 10 |  |
| Structure height above toe | hs | ft / m | 0.00328084 to 32808.4 | 15 |  |
| Structure crest width | B | ft / m | 0 to 32808.4 | 7.5 |  |
| Structure type | structure_type | (none) | choices: Sloped, Vertical or composite | Sloped |  |
| Transmission method (sloped) | transmission_method | (none) | choices: Seelig, d'Angremond | Seelig | Seelig 1980 (ACES) or d'Angremond 1996 (CEM standard); sloped structures only |
| Cotangent of structure slope | cot_theta | (none) | 0.001 to 1000 | 3 | sloped structures |
| Slope type | slope_type | (none) | choices: Rough (riprap), Smooth | Rough (riprap) | sloped structures |
| Rough-slope coefficient a | a | (none) | 0 to 10 | 0.956 |  |
| Rough-slope coefficient b | b | (none) | 0 to 10 | 0.398 |  |
| Known runup (0 = compute) | R_known | ft / m | 0 to 32808.4 | 0 |  |
| Berm height above toe | berm_height | ft / m | 0 to 32808.4 | 0 | vertical/composite; 0 = no berm (pure vertical wall) |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Wave runup (sloped structures) | R | ft / m | scalar |
| Crest freeboard (hs - ds) | F | ft / m | scalar |
| Transmission coefficient | K_TO | (none) | scalar |
| Transmitted wave height | H_T | ft / m | scalar |

*Reference:* Seelig (1980); Seelig (1976); Ahrens & McCartney (1975); Ahrens & Titus (1985)

*Module:* `backend/applications/chessqc_5_3_transmission_impermeable.py`

### 5-4 — Wave Transmission through Permeable Structures  `[I]`

**Status:** Current.

Originating ACES grouping: 5-4 "Wave Transmission through Permeable Structures" (functional area: Wave Runup, Transmission, and Overtopping). Estimates the height of the wave transmitted past a permeable, multilayered, trapezoidal rubble-mound breakwater, combining transmission by overtopping with transmission through the porous structure.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Incident wave height | Hi | ft / m | 0.000328084 to 3280.84 | 6.56 |  |
| Wave period | T | s | 0.01 to 1000 | 20 |  |
| Water depth at structure | ds | ft / m | 0.00328084 to 32808.4 | 15.75 |  |
| Structure height above toe | hs | ft / m | 0.00328084 to 32808.4 | 19.69 |  |
| Cotangent of structure slope | cot_theta | (none) | 0.001 to 1000 | 1.5 |  |
| Structure crest width | B | ft / m | 0 to 32808.4 | 8.27 |  |
| Method | method | (none) | choices: Madsen-White, d'Angremond + Zanuttigh | Madsen-White | Madsen-White (ACES layered model) or d'Angremond 1996 transmission + Zanuttigh-van der Meer 2008 reflection (empirical; ignores the layer geometry) |
| Material median diameters (list) | d50 | ft / m | any | None | armor, underlayer, core, ... (one per material) |
| Material porosities (list) | porosity | (none) | any | None |  |
| Layer thicknesses (list, bottom->top) | TH | ft / m | any | None |  |
| Material length per layer LL[material][layer] | LL | ft / m | any | None |  |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Equivalent rectangle width | le | ft / m | scalar |
| Seaward-slope reflection | R_si | (none) | scalar |
| Internal reflection | R_ti | (none) | scalar |
| Internal transmission | T_ti | (none) | scalar |
| Through-transmission coefficient | K_Tt | (none) | scalar |
| Overtopping transmission coeff. | K_To | (none) | scalar |
| Total transmission coefficient | K_T | (none) | scalar |
| Reflection coefficient | K_R | (none) | scalar |
| Transmitted wave height | H_T | ft / m | scalar |

*Reference:* Madsen & White (1976); Seelig (1980); Ahrens & McCartney (1975)

*Module:* `backend/applications/chessqc_5_4_transmission_permeable.py`

### 5-5 — Wave Setup  `[I]`

**Status:** Current.

Originating ACES grouping: 5-5 "Wave Setup" (functional area: Wave Runup, Transmission, and Overtopping). Predicts the wave-induced change in mean water level across the surf zone: the small set-down just seaward of breaking and the larger set-up that raises the mean water line toward and onto the beach.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Deepwater wave height | H0 | ft / m | 0.000328084 to 3280.84 | 6 | > 0 (significant or monochromatic deepwater height) |
| Wave period | T | s | 0.01 to 1000 | 10 | > 0 |
| Beach slope (tan beta) | m | (none) | 0.0001 to 1 | 0.02 | bottom slope; Singamsetti & Wind valid roughly 0.02 to 0.2 |
| Refraction coefficient | KR | (none) | 0 to 1 | 1 | K_R = sqrt(b0/b); 1.0 = no refraction. H'0 = K_R * H0 |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Unrefracted deepwater height H'0 | H0p | ft / m | scalar |
| Deepwater wave length | L0 | ft / m | scalar |
| Weggel breaker coefficient a(m) | a_index | (none) | scalar |
| Weggel breaker coefficient b(m) | b_index | (none) | scalar |
| Breaker height | Hb | ft / m | scalar |
| Breaker depth | db | ft / m | scalar |
| Breaker index Hb/db | gamma_b | (none) | scalar |
| Set-down at breaking | setdown_b | ft / m | scalar |
| Set-up gradient d(eta)/dx | setup_slope | (none) | scalar |
| Surf-zone width (breaker to SWL shoreline) | surf_width | ft / m | scalar |
| Set-up at still-water shoreline | setup_swl | ft / m | scalar |
| Maximum set-up (at waterline) | setup_max | ft / m | scalar |
| Shoreline displacement | dx_shore | ft / m | scalar |
| Profile: distance shoreward of breaker | profile_x | ft / m | profile |
| Profile: bed elevation | profile_bed | ft / m | profile |
| Profile: mean water level (set-up) | profile_mwl | ft / m | profile |

*Reference:* Longuet-Higgins & Stewart (1963); Weggel (1972); Singamsetti & Wind (1980)

*Module:* `backend/applications/chessqc_5_5_wave_setup.py`


## Littoral Processes

### 6-1 — Longshore Sediment Transport  `[I]`

**Status:** Current.

Originating ACES application: 6-1 "Longshore Sediment Transport" (functional area: Littoral Processes). Estimates the potential volumetric longshore sand transport rate from wave conditions, using the CERC energy-flux method, from breaking-wave or deepwater inputs.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Wave data type | method | (none) | choices: Deepwater wave conditions, Breaking wave conditions | Deepwater wave conditions | deepwater (H_s0, alpha_0) or breaking (H_b, alpha_b) |
| Wave height | H | ft / m | 0.000328084 to 3280.84 | 1.75 | deepwater significant height H_s0, or breaker height H_b |
| Wave angle to shoreline | angle | deg | 0 to 90 | 15 | deepwater crest angle alpha_0, or breaker angle alpha_b |
| Empirical coefficient K | K | (none) | 0 to 2 | 0.39 | CERC coefficient; 0.39 for field data with significant wave height |
| Water density | rho_water | kg/m^3 | 900 to 1100 | 1025.18 | seawater ~1025, fresh ~1000 |
| Sediment density | rho_sand | kg/m^3 | 1500 to 3500 | 2650 | quartz sand ~2650, which is what ACES uses: its 5.14 slug/ft^3 is 2649 kg/m^3 |
| Energy-flux wave theory | wave_theory | (none) | choices: Linear (2/98 revision, SMS/Genesis consistent), Solitary (original ACES, reproduces the User's Guide examples) | Linear (2/98 revision, SMS/Genesis consistent) | ACES shipped both. The solitary coefficients are 1.25x the linear ones and are what the published User's Guide examples were computed with; the 2/98 revision replaced them with the linear ones |
| Sediment porosity | porosity | (none) | 0 to 0.7 | 0.4 | pore fraction; solids fraction a' = 1 - porosity (TR uses a' = 0.6) |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Longshore transport rate | Q | yd^3/yr / m^3/yr | scalar |
| Transport rate (volume/sec) | Q_m3s | m^3/s | scalar |
| Longshore energy flux factor | P_ls | lb/s / N/s | scalar |
| Q / P_ls | cerc_factor | (none) | scalar |

*Reference:* SPM (1984) Ch.4; Galvin (1979); TR 6-1

*Module:* `backend/applications/chessqc_6_1_longshore_transport.py`

### 6-2 — Time-Dependent Beach and Dune Erosion  `[I]`

**Status:** Current.

Originating ACES application: 6-2 "Numerical Simulation of Time-Dependent Beach and Dune Erosion" (functional area: Littoral Processes). A transcription of the XSHORE explicit finite-difference scheme ACES runs (`KDMAIN.FOR` subroutine `KD1`, with `KD2`-`KD5B`), after Kriebel (1984b, EBEACH) and Kriebel & Dean (1985).

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Profile definition | profile_type | (none) | choices: Generic (dune, berm and slopes), Surveyed profile | Generic (dune, berm and slopes) | build a generic profile from dune, berm and slope descriptors, or supply a surveyed one |
| Dune crest elevation | dune_elev | ft / m | 0 to 200 | 20 | elevation of the dune top above the datum |
| Dune base elevation | dune_base_elev | ft / m | 0 to 200 | 6 | elevation at the foot of the dune face; must be at least 5 ft below the crest |
| Dune crest width | dune_width | ft / m | 0 to 1000 | 50 |  |
| Berm elevation | berm_elev | ft / m | 0 to 200 | 6 |  |
| Berm width | berm_width | ft / m | 0 to 2000 | 100 |  |
| Dune face slope (run per rise) | slope_dune | (none) | 0.1 to 1000 | 2 | 1 on this many; 2 means 1:2 |
| Beach face slope (run per rise) | slope_beach | (none) | 0.1 to 1000 | 10 |  |
| Nearshore slope (run per rise) | slope_nearshore | (none) | 0.1 to 1000 | 20 |  |
| Surveyed profile | survey | m | any | ((0.0, 4.29768), (1.2192, 4.08432), (3.41376, 3.99288), (7.65048, 3.23088), (13.716, 4.572), (16.58112, 4.29768), (23.07336, 3.81), (32.06496, 3.93192), (42.55008, 4.1148), (49.95672, 3.81), (57.72912, 3.2004), (62.6364, 2.62128), (73.914, 1.31064), (85.83168, 0.70104), (97.71888, 0.33528), (114.20856, 0.12192), (119.99976, 0.06096), (128.41224, -0.1524), (138.25728, -0.94488), (151.57704, -2.10312), (164.34816, -2.1336), (175.99152, -2.01168), (191.01816, -2.31648), (194.67576, -2.65176), (205.00848, -2.98704), (220.03512, -2.95656), (224.24136, -2.68224), (235.57992, -2.7432), (250.48464, -2.31648), (269.41272, -2.286), (291.6936, -3.048), (297.36288, -3.44424), (304.1904, -3.6576), (313.3344, -4.08432), (327.9648, -4.90728), (341.376, -5.51688), (351.4344, -5.85216), (362.712, -6.2484), (373.6848, -6.5532), (391.668, -6.97992), (401.1168, -7.0104), (418.1856, -7.43712), (433.1208, -7.7724), (452.628, -8.0772), (466.9536, -8.29056), (483.108, -8.62584), (495.3, -8.90016), (512.6736, -9.17448), (525.1704, -9.26592), (541.6296, -9.47928), (555.0408, -9.6012), (569.976, -9.81456), (583.9968, -9.87552)) | rows of (distance seaward from the baseline, bed elevation above datum); resampled onto the 4 ft grid by natural cubic spline |
| Mean grain size | grain_size | mm | 0.01 to 10 | 0.22 | sets the equilibrium profile shape factor A (Moore 1982) |
| Depth of the wave records | gauge_depth | ft / m | 5 to 9999 | 60 | still-water depth at which the wave conditions below were measured |
| Wave records | waves | m | any | ((2.4384, 8.0, 10.0),) | rows of (height, period s, crest angle to the shoreline deg), one per wave interval; heights are SI on the contract |
| Wave record interval | wave_dt_hr | h | 1 to 1000 | 20 |  |
| Length of simulation | length_hr | h | 1 to 1000 | 20 |  |
| Tabular output interval | out_interval_hr | h | 1 to 1000 | 2 | ignored when explicit output times are given below |
| Explicit output times | out_times | h | any | () | report at these hours instead of at a fixed interval; empty to use the interval. ACES takes one or the other, never both |
| Still-water level series | water_level | ft / m | any | () | surge or tide elevation above the datum, at the interval below; empty for none. ACES applies it as a uniform shift of the whole profile |
| Water-level interval | wl_dt_hr | h | 1 to 1000 | 1 |  |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Distance seaward from baseline | profile_X | ft / m | profile |
| Initial bed elevation | z_initial | ft / m | profile |
| Final bed elevation | z_final | ft / m | profile |
| Change in subaerial volume | vol_total | yd^3/ft / m^3/m | scalar |
| Shoreline recession | recession | ft / m | scalar |
| Contour change, +5 ft | contour_5_final | ft / m | scalar |
| Contour change, +10 ft | contour_10_final | ft / m | scalar |
| Contour change, +15 ft | contour_15_final | ft / m | scalar |
| Equilibrium profile shape factor | A_shape | ft^(1/3) / m^(1/3) | scalar |
| Equilibrium energy dissipation | D_eq | (none) | scalar |

*Reference:* Kriebel (1984b, EBEACH); Kriebel & Dean (1985); Moore (1982); Dean (1977); ACES Technical Reference Ch. 6-2

*Module:* `backend/applications/chessqc_6_2_dune_erosion.py`

### 6-3 — Longshore Transport from a CEDRS Wave Climate  `[I]`

**Status:** Current.

Originating ACES application: 6-1 "Longshore Sediment Transport", its third capability - "Transport using CEDRS statistical data: Percent Occurrence of wave height and period by direction" (`Lstran.for`, subroutines `lscband`, `lscomp` and `LSCTRN`). CHESS-QC carries it as its own application because its inputs and outputs are nothing like the single-condition case, which is CHESS-QC 6-1.

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Shore-normal azimuth | shore_azimuth | deg | 0 to 360 | 40 | compass bearing of the outward shore normal; sets which direction bands reach the beach and at what angle |
| Empirical coefficient K | K | (none) | 0 to 2 | 0.39 | CERC coefficient; 0.39 for field data with significant wave height |
| Water depth at the gauge | gauge_depth | ft / m | 1 to 5000 | 68 | depth the CEDRS wave record applies at; the station header carries it (68.0 m for the shipped G1033 record) |
| Coastline region | region | (none) | choices: Atlantic, Gulf, Pacific, Lake Erie, Lake Huron, Lake Michigan, Lake Ontario, Lake Superior | Gulf | selects the height and period class tables, which differ by region |
| Energy-flux wave theory | wave_theory | (none) | choices: Linear (2/98 revision, SMS/Genesis consistent), Solitary (original ACES, reproduces the User's Guide example) | Linear (2/98 revision, SMS/Genesis consistent) | ACES shipped both. The solitary coefficient is 1.25x the linear one and is what the published example was computed with |
| Water density | rho_water | kg/m^3 | 900 to 1100 | 1025.18 | seawater ~1025, fresh ~1000 |
| Sediment density | rho_sand | kg/m^3 | 1500 to 3500 | 2650 | quartz sand ~2650, which is what ACES uses: its 5.14 slug/ft^3 is 2649 kg/m^3 |
| Sediment porosity | porosity | (none) | 0 to 0.7 | 0.4 | pore fraction; solids fraction a' = 1 - porosity (ACES uses a' = 0.6) |
| CEDRS percent occurrence (x1000) | occ | (none) | any | ((34, 30, 5, 0, 0, 0, 0, 0, 0, 0), (378, 532, 30, 0, 0, 0, 0, 0, 0, 0), (0, 780, 99, 0, 0, 0, 0, 0, 0, 0), (0, 0, 345, 3, 0, 0, 0, 0, 0, 0), (0, 0, 42, 0, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (66, 41, 15, 0, 0, 0, 0, 0, 0, 0), (321, 602, 56, 0, 0, 0, 0, 0, 0, 0), (1, 888, 104, 0, 0, 0, 0, 0, 0, 0), (0, 0, 359, 1, 0, 0, 0, 0, 0, 0), (0, 0, 58, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (75, 18, 23, 0, 0, 0, 0, 0, 0, 0), (629, 915, 39, 0, 0, 0, 0, 0, 0, 0), (1, 1302, 135, 0, 0, 0, 0, 0, 0, 0), (0, 0, 592, 0, 0, 0, 0, 0, 0, 0), (0, 0, 126, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (213, 53, 135, 0, 0, 0, 0, 0, 0, 0), (675, 2250, 181, 0, 0, 0, 0, 0, 0, 0), (0, 3151, 532, 3, 0, 0, 0, 0, 0, 0), (0, 3, 1052, 6, 0, 0, 0, 0, 0, 0), (0, 0, 147, 3, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (472, 296, 75, 0, 0, 0, 0, 0, 0, 0), (980, 4433, 374, 0, 0, 0, 0, 0, 0, 0), (0, 2320, 4712, 159, 0, 0, 0, 0, 0, 0), (0, 18, 874, 831, 44, 0, 0, 0, 0, 0), (0, 0, 71, 39, 100, 0, 0, 0, 0, 0), (0, 0, 5, 0, 5, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (412, 345, 112, 0, 0, 0, 0, 0, 0, 0), (925, 4404, 545, 0, 0, 0, 0, 0, 0, 0), (0, 716, 4609, 547, 0, 0, 0, 0, 0, 0), (0, 0, 246, 828, 58, 0, 0, 0, 0, 0), (0, 0, 6, 23, 102, 1, 0, 0, 0, 0), (0, 0, 0, 3, 27, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 5, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (496, 381, 39, 0, 0, 0, 0, 0, 0, 0), (853, 3417, 662, 1, 0, 0, 0, 0, 0, 0), (1, 638, 3307, 513, 3, 0, 0, 0, 0, 0), (0, 0, 119, 961, 136, 0, 0, 0, 0, 0), (0, 0, 0, 23, 160, 3, 0, 0, 0, 0), (0, 0, 0, 0, 15, 15, 0, 0, 0, 0), (0, 0, 0, 0, 0, 11, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (133, 128, 77, 0, 0, 0, 0, 0, 0, 0), (627, 1875, 477, 6, 0, 0, 0, 0, 0, 0), (1, 487, 2345, 319, 1, 0, 0, 0, 0, 0), (0, 0, 119, 766, 203, 0, 0, 0, 0, 0), (0, 0, 0, 25, 249, 6, 0, 0, 0, 0), (0, 0, 0, 0, 39, 5, 0, 0, 0, 0), (0, 0, 0, 0, 3, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (131, 71, 90, 0, 0, 0, 0, 0, 0, 0), (453, 1483, 277, 0, 0, 0, 0, 0, 0, 0), (0, 369, 1185, 191, 8, 0, 0, 0, 0, 0), (0, 0, 100, 503, 97, 0, 0, 0, 0, 0), (0, 0, 1, 34, 159, 1, 0, 0, 0, 0), (0, 0, 0, 0, 6, 5, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (75, 66, 66, 0, 0, 0, 0, 0, 0, 0), (458, 1401, 241, 11, 0, 0, 0, 0, 0, 0), (0, 241, 922, 160, 17, 0, 0, 0, 0, 0), (0, 0, 87, 302, 53, 6, 1, 0, 0, 0), (0, 0, 1, 15, 71, 0, 0, 0, 0, 0), (0, 0, 0, 0, 15, 3, 0, 0, 0, 0), (0, 0, 0, 0, 0, 3, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (78, 210, 124, 0, 0, 0, 0, 0, 0, 0), (557, 1810, 272, 25, 1, 0, 0, 0, 0, 0), (0, 361, 1630, 188, 22, 6, 0, 0, 0, 0), (0, 0, 83, 311, 59, 5, 0, 0, 0, 0), (0, 0, 0, 20, 85, 3, 0, 0, 0, 0), (0, 0, 0, 1, 35, 6, 0, 0, 0, 0), (0, 0, 0, 0, 0, 13, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (131, 191, 73, 0, 0, 0, 0, 0, 0, 0), (557, 1813, 260, 6, 0, 0, 0, 0, 0, 0), (1, 278, 968, 181, 6, 0, 0, 0, 0, 0), (0, 0, 80, 244, 37, 0, 0, 0, 0, 0), (0, 0, 0, 23, 63, 0, 0, 0, 0, 0), (0, 0, 0, 0, 6, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 1, 3, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (121, 97, 46, 0, 0, 0, 0, 0, 0, 0), (391, 1031, 361, 10, 0, 0, 0, 0, 0, 0), (1, 140, 740, 189, 10, 0, 0, 0, 0, 0), (0, 0, 106, 290, 68, 0, 0, 0, 0, 0), (0, 0, 8, 46, 90, 0, 0, 0, 0, 0), (0, 0, 0, 0, 25, 1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (44, 56, 58, 0, 0, 0, 0, 0, 0, 0), (301, 939, 330, 1, 0, 0, 0, 0, 0, 0), (0, 210, 869, 116, 1, 0, 0, 0, 0, 0), (0, 1, 179, 453, 6, 0, 0, 0, 0, 0), (0, 0, 8, 119, 42, 0, 0, 0, 0, 0), (0, 0, 0, 0, 3, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (49, 71, 35, 0, 0, 0, 0, 0, 0, 0), (256, 872, 106, 0, 0, 0, 0, 0, 0, 0), (0, 669, 588, 27, 0, 0, 0, 0, 0, 0), (0, 1, 285, 200, 1, 0, 0, 0, 0, 0), (0, 0, 10, 37, 1, 0, 0, 0, 0, 0), (0, 0, 1, 5, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (63, 41, 10, 0, 0, 0, 0, 0, 0, 0), (369, 619, 100, 0, 0, 0, 0, 0, 0, 0), (1, 795, 123, 0, 0, 0, 0, 0, 0, 0), (0, 5, 306, 3, 0, 0, 0, 0, 0, 0), (0, 0, 58, 3, 0, 0, 0, 0, 0, 0), (0, 0, 3, 1, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) | one row per (direction band, wave height class), band-major: the first 11 rows are band 1's height classes, the next 11 band 2's, and so on for 16 bands. Each row holds the 10 peak-period columns. This is the CEDRS station file laid out flat |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Net longshore transport | Q_net | yd^3/yr / m^3/yr | scalar |
| Gross longshore transport | Q_gross | yd^3/yr / m^3/yr | scalar |
| Transport to shore-right | Q_right | yd^3/yr / m^3/yr | scalar |
| Transport to shore-left | Q_left | yd^3/yr / m^3/yr | scalar |
| Wave approach angle | profile_X | deg | profile |
| Contributing fraction | band_pct | % | profile |
| Transport by band | band_Q | yd^3/yr / m^3/yr | profile |
| Contributing bands | n_bands | (none) | scalar |

*Reference:* SPM (1984) Ch.4; Galvin (1979, 1980); Gravens (1988); Wang, Kraus & Davis (1998); ACES Technical Reference Ch. 6-1

*Module:* `backend/applications/chessqc_6_3_cedrs_transport.py`

### 6-4 — Beach Nourishment Overfill Ratio and Volume  `[I]`

**Status:** Current.

Originating ACES application: 6-4 "Beach Nourishment Overfill Ratio and Volume" (functional area: Littoral Processes). Tells a nourishment designer how much borrow sand to place to obtain a given volume of usable beach (the overfill ratio), and how much faster or slower the borrow sand erodes compared with the native sand (the renourishment factor).

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Initial (usable) volume | VOL_I | yd^3 / m^3 | 0 to 1.30795e+12 | 800000 | target volume of usable beach fill |
| Native mean | M_n | phi | -5 to 15 | 1.8 | phi mean grain size of the native beach sand |
| Native standard deviation | sigma_n | phi | 0.0001 to 10 | 0.45 | phi sorting (standard deviation) of the native sand |
| Borrow mean | M_b | phi | -5 to 15 | 2.25 | phi mean grain size of the borrow sand |
| Borrow standard deviation | sigma_b | phi | 0.0001 to 10 | 0.76 | phi sorting (standard deviation) of the borrow sand |
| Winnowing function | W | (none) | 0 to 5 | 1 | James (1975) renourishment winnowing parameter (recommended 1.0) |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Overfill ratio | R_A | (none) | scalar |
| Renourishment factor | R_J | (none) | scalar |
| Design (borrow) volume | VOL_D | yd^3 / m^3 | scalar |
| Phi-mean difference (scaled) | delta | (none) | scalar |
| Sorting ratio | sigma_ratio | (none) | scalar |
| James category | category | (none) | scalar |

*Reference:* James (1975) TM-60; SPM (1984); TR 6-4

*Module:* `backend/applications/chessqc_6_4_beach_nourishment.py`

### 6-5 — Composite Grain-Size Distribution  `[I]`

**Status:** Current.

Originating ACES application: 6-5 "Composite Grain Size" (functional area: Littoral Processes; TR chapter 6-3). Computes grain-size statistics for a sediment sample (or a composite of samples) from sieve data, by both the Folk graphic method and the method of moments.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Sieve data | sieve | (none) | table: Grain size, Weight | 14 default rows | one row per sieve: phi size and weight retained (any consistent weight unit). For a composite, enter the combined (averaged) distribution. |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Mean (moments) | mom_mean | phi | scalar |
| Sorting (moments) | mom_sigma | phi | scalar |
| Skewness (moments) | mom_skew | (none) | scalar |
| Kurtosis (moments) | mom_kurt | (none) | scalar |
| Median (Folk, phi50) | folk_median | phi | scalar |
| Mean (Folk graphic) | folk_mean | phi | scalar |
| Sorting (Folk graphic) | folk_sigma | phi | scalar |
| Skewness (Folk graphic) | folk_skew | (none) | scalar |
| Kurtosis (Folk graphic) | folk_kurt | (none) | scalar |
| Median diameter | d50_mm | mm | scalar |
| Profile: phi size | profile_phi | phi | profile |
| Profile: cumulative percent | profile_cum | % | profile |

*Reference:* Folk (1974); Krumbein (1934, 1938); SPM (1984); TR 6-3

*Module:* `backend/applications/chessqc_6_5_composite_grain_size.py`


## Inlet Processes

### 7-1 — Spatially Integrated Numerical Model for Inlet Hydraulics  `[I]`

**Status:** Current.

Originating ACES grouping: 7-1 "A Spatially Integrated Numerical Model for Inlet Hydraulics" (functional area: Inlet Processes). This is the single largest, most input-intensive ACES application: a 1-D continuity + momentum model (Seelig 1977; Seelig, Harris & Herchenroder 1977; momentum after Harris & Bodine 1977) that time-marches the coupled inlet discharge Q(t) and bay water level h_b(t) under a constituent sea tide, through a multi-cross-section inlet, by 4th-order Runge-Kutta. It is NOT the Keulegan lumped-parameter model.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| M2 tide amplitude | tide_amp | ft / m | 0.001 to 50 | 2 |  |
| M2 epoch (phase lag kappa) | tide_epoch | deg | 0 to 360 | 90 |  |
| Sea boundary longitude (deg West) | gage_lon | deg | -180 to 180 | 75 |  |
| Start year | year | (none) | 1900 to 2100 | 1988 |  |
| Start month | month | (none) | 1 to 12 | 7 |  |
| Start day | day | (none) | 1 to 31 | 6 |  |
| Start hour | hour | h | 0 to 24 | 0 |  |
| Length of simulation | length_hr | h | 1 to 720 | 30 |  |
| Time step | dt_s | s | 1 to 600 | 60 |  |
| Tabular output interval | out_interval_min | min | 1 to 240 | 15 |  |
| Flood loss coefficient | flood_loss | (none) | 0 to 100 | 4 |  |
| Ebb loss coefficient | ebb_loss | (none) | 0 to 100 | 1 |  |
| Equal-discharge channels per cross-section | n_channels | (none) | 1 to 7 | 4 | each surveyed cross-section is divided into this many channels of equal discharge before the friction is summed (ACES allows up to 7) |
| Manning coefficient C1 | manning_C1 | (none) | 0 to 1 | 0.05 |  |
| Manning coefficient C2 | manning_C2 | (none) | 0 to 1 | 0.0007 |  |
| Bay surface area | bay_area | ft^2 / m^2 | 1 to 1e+15 | 1.8e+09 |  |
| Bay area variation parameter | bay_beta | (none) | 0 to 10 | 0 |  |
| River inflow tabulation interval | river_dt_min | min | 1 to 10000 | 260 |  |
| River / non-inlet inflow series | river | ft^3/s / m^3/s | any | (113.267386, 107.604017, 101.940648, 90.613909, 99.108963, 107.604017, 118.930756, 121.76244, 127.42581) | tabulated discharge (m^3/s) at the river interval; linearly interpolated |
| Inlet cross-sections (bathymetry) | sections | m | any | ((31.6992, 533.4, [0.0, -8.2296, -8.2296, -8.2296, -8.2296, -8.2296, -8.2296, -8.2296, -8.2296, -8.2296, -8.2296, -5.4864, -3.9624, -3.9624, -3.9624, -3.9624, -3.9624, -3.9624, -3.9624, -5.4864, -7.3152, -9.144, -9.7536, -10.3632, -10.3632, -10.3632, -10.3632, -9.7536, -9.7536, -9.7536, -9.7536, -7.3152, -7.3152, -7.3152, -7.3152, -7.62, -7.62, -5.4864, -5.4864, -5.4864, -5.4864, 0.0]), (31.6992, 495.3, [0.0, -9.144, -10.0584, -10.0584, -10.0584, -10.3632, -10.3632, -10.3632, -10.3632, -10.3632, -9.144, -9.144, -6.096, -3.048, 0.0]), (31.6992, 584.3016, [0.0, -3.6576, -5.4864, -6.096, -7.62, -9.144, -10.0584, -10.3632, -10.3632, -10.3632, -10.3632, -10.3632, -10.3632, -9.144, -5.4864, -3.6576, -2.4384, -2.4384, -2.4384, -1.8288, -1.8288, -1.8288, -1.8288, 0.0]), (31.6992, 381.0, [0.0, -5.4864, -11.2776, -11.2776, -15.24, -15.24, -15.24, -10.3632, -10.3632, -10.3632, -10.3632, -7.3152, -5.4864, 0.0]), (31.6992, 0.0, [0.0, -3.3528, -3.3528, -3.3528, -3.6576, -3.6576, -5.1816, -5.1816, -5.1816, -4.572, -4.572, -4.572, -5.4864, -7.62, -7.62, -6.096, -6.096, -6.096, -10.3632, -10.3632, -10.3632, -10.3632, -7.0104, -5.4864, -3.048, -3.048, -3.048, -3.048, -3.048, -3.048, -3.048, -3.048, -3.048, -3.048, -3.048, -3.048, 0.0])) | one row per cross-section: (dX m, along-inlet length dY m, [bed elevations m]); area and width are integrated from the elevation profile relative to datum 0 |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Throat (minimum) cross-section area | throat_area | ft^2 / m^2 | scalar |
| Geometry integral | I_g | ft / m | scalar |
| Bay tidal range | bay_range | ft / m | scalar |
| Peak ebb discharge | max_ebb_Q | ft^3/s / m^3/s | scalar |
| Peak flood discharge | max_flood_Q | ft^3/s / m^3/s | scalar |
| Peak controlling-section velocity | max_vel | ft/s / m/s | scalar |
| Time | t | h | profile |
| Sea elevation | sea_el | ft / m | profile |
| Bay elevation | bay_el | ft / m | profile |
| Inlet discharge | inlet_Q | ft^3/s / m^3/s | profile |
| Controlling-section velocity | control_vel | ft/s / m/s | profile |

*Reference:* Seelig (1977); Seelig, Harris & Herchenroder (1977); Harris & Bodine (1977); Keulegan (1967); Schureman (1971)

*Module:* `backend/applications/chessqc_7_1_inlet_hydraulics.py`

### 7-2 — Wave-Current Interaction in Channels  `[I]`

**Status:** Current.

Originating ACES application: 7-2 "Wave-current Interaction in Channels" (functional area: Inlet Processes; a later ACES addition). Computes how a wave train is modified when it crosses a navigation channel carrying a steady current: the current Doppler- shifts the dispersion relation (changing the wavelength) and changes the wave height through conservation of wave action.

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Wave period | T | s | 0.5 to 30 | 8 | absolute wave period (ground frame) T > 0 |
| Angle of wave orthogonal to current | alpha | deg | 0 to 180 | 0 | 0 = following current, 180 = directly opposing |
| Channel current velocity | V | kt / m/s | -10 to 10 | 1.5 | current speed; combined with the angle (positive magnitude) |
| Channel depth | dT | ft / m | 0.3 to 300 | 10 | still-water channel depth d > 0 |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Dimensionless frequency (omega sqrt(d/g)) | omega_star | (none) | scalar |
| Effective Froude number U/Cg | F | (none) | scalar |
| Wave height factor H/H0 | R_H | (none) | scalar |
| Wavelength factor L/L0 | R_L | (none) | scalar |
| Wavelength (no current) | L0 | ft / m | scalar |
| Wavelength (with current) | L | ft / m | scalar |

*Reference:* Jonsson (1990); Jonsson, Skovgaard & Wang (1970); Peregrine (1976)

*Module:* `backend/applications/chessqc_7_2_wave_current.py`


## Harbor Design

### 8-1 — Properties of Rectangular Basins  `[I]`

**Status:** Current.

Originating ACES grouping: 8-1 "Properties of Rectangular Basins" (functional area: Harbor Design). Computes the natural (resonant) oscillation periods of a rectangular harbor or basin, the seiche modes that a long wave, surge, or tsunami can excite, plus the standing-wave water-particle kinematics at a node.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Basin type | basin_type | (none) | choices: Open (one end open), Closed (both ends), Closed 2-D (rectangular), Helmholtz (basin + channel) | Closed (both ends) | open/closed 1-D, closed 2-D, or Helmholtz pumping |
| Basin length | lB | ft / m | 0.00328084 to 3.28084e+06 | 3000 | length along the resonant (longitudinal) axis |
| Basin width | lC | ft / m | 0.00328084 to 3.28084e+06 | 2000 | 2-D only: transverse dimension |
| Water depth | d | ft / m | 0.00328084 to 32808.4 | 30 | mean basin depth; resonance uses shallow-water c = sqrt(g d) |
| Standing-wave height | H | ft / m | 0 to 3280.84 | 1 | antinode crest-to-trough height; used for the node kinematics |
| Longitudinal mode n | n | (none) | 0 to 50 | 1 | closed/2-D: n>=1; open: n>=0 (n=0 is the fundamental) |
| Transverse mode m | m_mode | (none) | 0 to 50 | 1 | 2-D only: m>=0 (m=0 reduces to the 1-D longitudinal mode) |
| Basin surface area | Ab | ft^2 / m^2 | 0.0107639 to 1.07639e+13 | 6e+06 | Helmholtz only: plan area of the basin |
| Channel cross-section area | Ac | ft^2 / m^2 | 0.0107639 to 1.07639e+11 | 6000 | Helmholtz only: flow area of the entrance channel |
| Channel length | L_ch | ft / m | 0 to 328084 | 500 | Helmholtz only: length of the entrance channel |
| Mouth length correction | L_corr | ft / m | 0 to 328084 | 100 | Helmholtz only: added-mass correction at the mouth |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Shallow-water celerity sqrt(g d) | c_shallow | ft/s / m/s | scalar |
| Resonant period of selected mode | T_mode | s | scalar |
| Modal wavelength (1-D) | L_mode | ft / m | scalar |
| Resonant frequency of selected mode | f_mode | 1/s | scalar |
| Helmholtz (pumping) period | T_helm | s | scalar |
| Max horizontal velocity at node (surface) | Vmax | ft/s / m/s | scalar |
| Horizontal particle semi-excursion at node | excursion | ft / m | scalar |
| Mean horizontal speed at node | Vbar | ft/s / m/s | scalar |
| Mode table: mode number | mode_index | (none) | profile |
| Mode table: resonant period | mode_T | s | profile |

*Reference:* Merian's formula; Helmholtz mode; SPM (1984); Wilson (1972)

*Module:* `backend/applications/chessqc_8_1_rectangular_basins.py`

### 8-2 — Vessel-Generated Waves  `[I]`

**Status:** Current.

Originating ACES application: 8-2 "Vessel-Generated Waves" (functional area: Harbor Design; a later ACES addition). For a design vessel moving at a given speed in a prismatic channel, it returns the depth Froude number, the celerity/period/propagation direction of the generated wave system, and the Schijf (1949) one-dimensional return-current and drawdown (water-level depression alongside the vessel).

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Channel width | b | ft / m | 1 to 10000 | 100 | prismatic channel top width b > 0 |
| Channel depth | d | ft / m | 0.5 to 200 | 6 | still-water channel depth d > 0 |
| Vessel speed | VS | kt / km/h | 0.036 to 108 | 10.8 | vessel speed through the water V_s > 0 |
| Wetted cross-sectional area | Am | ft^2 / m^2 | 0 to 10000 | 25 | submerged midship section area A_m >= 0 |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Depth Froude number | F | (none) | scalar |
| Wave celerity | C | kt / m/s | scalar |
| Wave period | T | s | scalar |
| Wavelength | L | ft / m | scalar |
| Wave crest angle to sailing line | direction | deg | scalar |
| Drawdown (water-level depression) | drawdown | ft / m | scalar |
| Relative drawdown dh/d | D | (none) | scalar |
| Blockage ratio A_m/(b d) | S | (none) | scalar |
| Return current | Vr | kt / m/s | scalar |

*Reference:* Schijf (1949); PIANC (1987); EM 1110-2-1100; Kelvin/Havelock ship-wave theory

*Module:* `backend/applications/chessqc_8_2_vessel_waves.py`

### 8-3 — Surging of a Moored Vessel  `[I]`

**Status:** Current.

Originating ACES application: 8-3 "Surging of a Moored Vessel" (functional area: Harbor Design; a later ACES addition). For a vessel held by a set of mooring lines, it resolves the line stiffnesses onto the surge (longitudinal) axis, adds the hydrodynamic added mass, and returns the natural surge period together with the per-line loading and the forward / reverse / total surge spring constants.

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Vessel mass | m | tonne | 1 to 5e+06 | 5000 | displacement mass m > 0 (stored in kg internally is not used; tonnes here) |
| Surge added-mass coefficient | Ca | (none) | 0 to 1 | 0.1 | virtual mass m_v = m (1 + C_a); surge C_a ~ 0.05-0.25 for ships |
| Mooring lines | lines | (none) | table: Angle to surge axis, Line length, Pretension, Breaking strength, Elongation at break | 4 default rows | one row per line; angle from the +surge (forward) axis to the anchor |
| Safe working load fraction | swl_fraction | (none) | 0.1 to 1 | 0.5 | line-impact flag trips when load exceeds this x breaking |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Virtual (surge) mass | m_v | kg | scalar |
| Forward surge spring constant | k_fwd | N/m | scalar |
| Reverse surge spring constant | k_rev | N/m | scalar |
| Total surge spring constant | k_total | N/m | scalar |
| Natural surge period | T_S | s | scalar |
| Maximum line load (% breaking) | max_load | % | scalar |
| Lines over safe working load | n_overloaded | (none) | scalar |
| Line-impact flag (1 = yes) | impact | (none) | scalar |

*Reference:* EM 1110-2-1100 Part II; PIANC mooring guidelines; ACES manual

*Module:* `backend/applications/chessqc_8_3_moored_vessel_surge.py`


## Storm Surge

### 9-1 — Bathystrophic Storm Surge  `[II]`

**Status:** Screening only. Newer method: ADCIRC (risk assessment).

First Quick Compute tool beyond the original 34 ACES applications (functional area: Storm Surge). Estimates open-coast hurricane surge along a single cross-shelf traverse by the quasi-1D bathystrophic method of Bodine (1971), CERC TM-35 (verification: Pararas-Carayannis 1975, TM-50).

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Traverse bathymetry | bathy | (none) | table: Distance from shore, Depth below SWL | 17 default rows | shelf edge -> shore; one row per station |
| Central pressure | Pc | inHg / hPa | 23.624 to 30.1206 | 27.57 | storm central pressure |
| Peripheral pressure | Pn | inHg / hPa | 28.9394 to 30.4159 | 29.92 | ambient/peripheral pressure |
| Radius of maximum winds | R | nm / km | 1.07991 to 107.991 | 35 | > 0 |
| Forward speed | Vf | kt / km/h | 0 to 77.7538 | 22 | storm translation speed |
| Track offset from traverse | track_offset | nm / km | 0 to 269.978 | 35 | alongshore distance from the traverse to the storm landfall/track |
| Latitude | lat | deg | 0 to 80 | 37 | for the Coriolis parameter |
| Wind model | wind_model | (none) | choices: Bodine TM-35 source curves, Holland (1980), Myers / Bodine (1954) | Bodine TM-35 source curves | Bodine's published source curves, or Holland/Myers parametric sensitivity paths |
| Holland B (peakedness) | B_holland | (none) | 0.5 to 2.5 | 1.5 | Holland shape factor; locked to 1.0 for Myers/Bodine |
| Max wind (optional) | Vmax | mph / km/h | 0 to 268.432 | 0 | if > 0, overrides B via B = rho_a e Vmax^2 / dP |
| Air density | rho_air | kg/m^3 | 1 to 1.3 | 1.2 | ambient air density |
| Bottom friction coefficient | K_bottom | (none) | 0.0001 to 0.02 | 0.0025 | bed friction K (~0.002-0.005) |
| Initial water-level rise | Se | ft / m | -9.84252 to 32.8084 | 0.5 | initial setup at start of computation |
| Astronomical tide | SA | ft / m | -16.4042 to 32.8084 | 2.5 | astronomical tide above MSL datum |
| Time step | dt | hr | 0.0166667 to 2 | 0.5 | integration time step (stored in seconds) |
| Number of time steps | n_steps | (none) | 10 to 400 | 62 | storm is swept past the traverse over these steps |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Peak surge at shore | peak_surge | ft / m | scalar |
|   wind setup (Sx) | S_wind | ft / m | scalar |
|   bathystrophic setup (Sy) | S_bathy | ft / m | scalar |
|   pressure setup | S_press | ft / m | scalar |
| Max wind speed (30 ft) | Vmax_out | mph / km/h | scalar |
| Holland B used | B_used | (none) | scalar |
| Time of peak | t_peak | hr / s | scalar |
| Profile: distance from shore | profile_X | nm / m | profile |
| Profile: setup at peak | profile_eta | ft / m | profile |
| Profile: still-water depth | profile_u | ft / m | profile |
| Profile: total depth at peak | profile_w | ft / m | profile |

*Reference:* Bodine (1971) TM-35; Holland (1980); Myers (1954); TR/CERC

*Module:* `backend/applications/chessqc_9_1_surge_bathystrophic.py`


## Coastal Hazards

### 10-1 — Water Level Detrending  `[II]`

**Status:** Current.

Functional area: Coastal Hazards. Estimates the long-term sea-level trend of a water-level record and removes it, leaving a series that varies about the tidal datum rather than about a moving mean.

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Water-level record | csv | (none) | CSV text / uploaded file | embedded sample | Select a bundled NOAA station or upload your own CSV (column 1 = date, column 2 = water level in m). Header and blank water-level rows are ignored. Hourly records and ready-made monthly means are both accepted. |
| Trend estimator | fit_mode | (none) | choices: Monthly GLS (default), Raw-sample OLS, Specified slope | Monthly GLS (default) | Monthly GLS is NOAA's published CO-OPS 053 procedure: calendar-month means fitted with a seasonal cycle and AR(1) errors, which is the only option here that yields a confidence interval. Raw-sample OLS regresses every sample on time instead (no seasonal term, no serial-correlation correction), for records too short for the monthly fit. Or supply a published rate directly. |
| Level shifts (optional) | level_shifts | (none) | table: Shift year | 0 default rows | Add a row per CONFIRMED step discontinuity in the record - a gauge relocation or datum change - as the decimal year it takes effect (e.g. 1947.5). Each adds an indicator that absorbs the step, so it cannot be soaked up by the trend instead. Leave empty when the record has no documented steps. |
| NTDE start year | ntde_start | yr | 1800 to 2100 | 1983 | National Tidal Datum Epoch start year; with the end year it sets the epoch center the detrended series is referenced to. It plays no part in the trend fit. |
| NTDE end year | ntde_end | yr | 1800 to 2100 | 2001 | National Tidal Datum Epoch end year (inclusive). |
| Specified slope | slope_value | in/yr / mm/yr | -1000 to 1000 | 3 | used only when the estimator is a supplied slope |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Sea-level trend (slope) | slope_per_year | in/yr / mm/yr | scalar |
| Trend 95% CI half-width | ci_halfwidth | in/yr / mm/yr | scalar |
| Residual AR(1) rho | rho_ar1 | (none) | scalar |
| Fitted seasonal range | seasonal_range | ft / m | scalar |
| Complete months in fit | n_months | (none) | scalar |
| Degrees of freedom | dof | (none) | scalar |
| Epoch center (datum) year | epoch_year | yr | scalar |
| Total trend over record | total_trend | ft / m | scalar |
| Record length | record_years | yr | scalar |
| Samples in record | n_samples | (none) | scalar |
| RMS residual about trend | rms_residual | ft / m | scalar |
| Fit converged | converged | (none) | scalar |
| Epoch center | pivot_line | yr | vline |
| Profile: year | profile_year | yr | profile |
| Profile: observed | profile_original | ft / m | profile |
| Profile: linear trend | profile_trend | ft / m | profile |
| Profile: detrended | profile_detrended | ft / m | profile |
| Profile: datum | profile_datum | ft / m | profile |
| handoff | handoff_csv | (none) | data |

*Reference:* Zervas (2009) NOAA CO-OPS 053; Cochrane & Orcutt (1949); NTDE datum convention

*Module:* `backend/applications/chessqc_10_1_water_level_detrending.py`

### 10-2 — Non-Tidal Residual  `[II]`

**Status:** Current.

Functional area: Coastal Hazards. The non-tidal residual is the (detrended) observed water level minus the predicted astronomical tide on the same time grid:

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Water-level record | csv_wl | (none) | CSV text / uploaded file | embedded sample | Observed (ideally detrended, via 10-1) water level: a bundled NOAA station or an uploaded CSV (column 1 = date, column 2 = level in m). |
| Tide prediction | csv_tide | (none) | CSV text / uploaded file | embedded sample | Predicted astronomical tide on the same grid: a bundled NOAA station or an uploaded CSV (column 1 = date, column 2 = tide in m). |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Mean NTR | mean_ntr | ft / m | scalar |
| RMS NTR | rms_ntr | ft / m | scalar |
| Maximum NTR | max_ntr | ft / m | scalar |
| Samples | n_samples | (none) | scalar |
| Profile: year | profile_year | yr | profile |
| Profile: water level | profile_wl | ft / m | profile |
| Profile: tide | profile_tide | ft / m | profile |
| Profile: NTR | profile_ntr | ft / m | profile |
| handoff | handoff_csv | (none) | data |

*Reference:* NOAA CO-OPS tide/residual practice; PyStorm NTR

*Module:* `backend/applications/chessqc_10_2_non_tidal_residual.py`

### 10-3 — Peaks Over Threshold  `[II]`

**Status:** Current.

Functional area: Coastal Hazards. Extracts independent storm peaks from a continuous water-level / non-tidal-residual series: an automatically chosen percentile threshold is raised until the declustered exceedance rate matches a target average number of events per year, then the exceedances are declustered into one peak per independent storm and rank-trimmed to exactly round(target x effective-duration) peaks. These peaks are the input to the Probabilistic Simulation Technique (10-4).

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Water-level / NTR record | csv | (none) | CSV text / uploaded file | embedded sample | Select a bundled NOAA station or upload your own CSV (column 1 = date, column 2 = value). Typically the detrended water level (10-1) or NTR (10-2). |
| Target events per year | target_events_per_year | 1/yr | 0.1 to 365 | 10 | average number of independent peaks per year to retain (matches PST) |
| Inter-event window | interevent_hours | hr | 1 to 2160 | 48 | minimum separation between independent events |
| Declustering method | method | (none) | choices: hydrograph, peak_gap | hydrograph | hydrograph: group + per-group max; peak_gap: sequential gap filter |
| Start percentile | start_percentile | % | 0 to 99.9 | 75 | series percentile where the upward threshold scan begins |
| Percentile step | step_size | % | 0.001 to 5 | 0.05 | percentile increment per scan step (smaller = finer, slower) |
| Convergence tolerance | tolerance | 1/yr | 0.001 to 10 | 0.25 | accept a rate in [target, target + tolerance] |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Threshold | threshold | ft / m | scalar |
| Peaks retained | n_peaks | (none) | scalar |
| Effective rate | events_per_year | 1/yr | scalar |
| Threshold percentile | final_percentile | % | scalar |
| Effective duration | eff_duration | yr | scalar |
| Converged | converged | (none) | scalar |
| Profile: year | profile_year | yr | profile |
| Profile: series | profile_series | ft / m | profile |
| Profile: threshold | profile_threshold | ft / m | profile |
| Peak times | peaks_t | yr | scatter_x |
| Peaks | peaks_v | ft / m | scatter |
| handoff | handoff_csv | (none) | data |

*Reference:* Coles (2001); USACE coastal-hazards practice; PyStorm POT

*Module:* `backend/applications/chessqc_10_3_peaks_over_threshold.py`

### 10-4 — Probabilistic Simulation Technique  `[II]`

**Status:** Current.

Functional area: Coastal Hazards. Turns a Peaks-Over-Threshold sample (10-3) into a hazard curve: response magnitude versus annual exceedance rate (AER), with a bootstrap confidence band. The upper tail is a fitted Generalized Pareto Distribution (GPD); the frequent range is carried empirically; the two are spliced into one continuous curve.

**Inputs** (values in SI units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| POT peak sample | csv | (none) | CSV text / uploaded file | embedded sample | Peaks over threshold (handed off from 10-3) or an uploaded CSV with the peak magnitude in column 2. |
| Events per year | events_per_year | 1/yr | 0.1 to 365 | 10 | POT target rate; sets the auto record length (n_pot / events_per_year) when record length is 0 |
| Record length | record_length_years | yr | 0 to 100000 | 0 | 0 = auto (n_pot / events_per_year) |
| Bootstrap simulations | num_simulations | (none) | 20 to 2000 | 300 | Monte-Carlo realizations for the confidence band |
| Random seed | seed | (none) | 0 to 1e+06 | 12345 | fixed for reproducible confidence bands |
| GPD band: min percentile | threshold_min_percentile | % | 0 to 99 | 50 | lower bound of the GPD-location scan band |
| GPD band: max percentile | threshold_max_percentile | % | 1 to 99.9 | 95 | upper bound of the GPD-location scan band |
| Min exceedances | min_exceedances | (none) | 5 to 100000 | 30 | a candidate location must retain at least this many exceedances |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| POT peaks | n_pot | (none) | scalar |
| GPD location mu | gpd_threshold | ft / m | scalar |
| GPD shape xi | gpd_shape | (none) | scalar |
| GPD scale sigma | gpd_scale | ft / m | scalar |
| Base rate lambda_u | lambda_u | 1/yr | scalar |
| Rate above mu | lambda_mu | 1/yr | scalar |
| 10-yr magnitude | mag_10yr | ft / m | scalar |
| 100-yr magnitude | mag_100yr | ft / m | scalar |
| 500-yr magnitude | mag_500yr | ft / m | scalar |
| Profile: log10 return interval | profile_logmri | log10 yr | profile |
| Profile: best estimate | profile_be | ft / m | profile |
| Profile: 10% confidence | profile_cb10 | ft / m | profile |
| Profile: 90% confidence | profile_cb90 | ft / m | profile |
| Empirical log10 MRI | emp_logmri | log10 yr | scatter_x |
| Empirical peaks | emp_mag | ft / m | scatter |

*Reference:* Coles (2001); Nadal-Caraballo et al. PST; PyStorm PST

*Module:* `backend/applications/chessqc_10_4_probabilistic_simulation.py`


## Miscellaneous Routines

### M-1 — Miscellaneous Breaker and Steepness Routines  `[II]`

**Status:** Current.

Originating ACES grouping: the ACES Technical Reference "Miscellaneous Routines" (a set of shared breaker and steepness utility relations). It is renumbered M-1 here to avoid a clash with the Harbor Design application 8-1 (Properties of Rectangular Basins), and it sorts last. These are the Weggel (1972) breaker index and related limiting-wave criteria that applications 3-1 (Snell) and 5-5 (wave setup) use internally; this app exposes them directly for quick checks.

**Inputs** (values in US units)

| Input | key | units (US/SI) | range | default | notes |
| --- | --- | --- | --- | --- | --- |
| Unrefracted deepwater wave height | H0 | ft / m | 0.000328084 to 3280.84 | 6 |  |
| Wave period | T | s | 0.01 to 1000 | 8 |  |
| Water depth (steepness / flat breaking) | d | ft / m | 0.00328084 to 32808.4 | 15 |  |
| Nearshore slope (tan theta) | m | (none) | 0 to 1 | 0.05 |  |
| Water depth at structure | ds | ft / m | 0.00328084 to 32808.4 | 15 |  |

**Outputs**

| Output | key | units (US/SI) | kind |
| --- | --- | --- | --- |
| Wave length at depth d | L | ft / m | scalar |
| Maximum steepness H/L (Miche) | steepness_max | (none) | scalar |
| Maximum height at steepness limit | Hmax_steep | ft / m | scalar |
| Breaking height, flat slope (McCowan) | Hb_flat | ft / m | scalar |
| Breaking height, finite slope (S&W) | Hb_sloped | ft / m | scalar |
| Breaker depth (Weggel) | db_sloped | ft / m | scalar |
| Weggel breaker index a(m) | a_index | (none) | scalar |
| Weggel breaker index b(m) | b_index | (none) | scalar |
| Breaker height at structure depth | Hb_structure | ft / m | scalar |

*Reference:* Weggel (1972); Miche (1944); McCowan (1894); Singamsetti & Wind (1980)

*Module:* `backend/applications/chessqc_m_1_breaker_routines.py`
