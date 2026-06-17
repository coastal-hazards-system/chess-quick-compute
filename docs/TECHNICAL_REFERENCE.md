# CHESS-QC — Technical Reference
*Coastal Hazards, Engineering, and Structures System (CHESS) — Quick Compute (QC)*

This is the engineering reference for the CHESS-QC (Quick Compute) toolkit. It explains what
each application computes, the physics behind it, the assumptions it rests on, where it is
still the right tool, and where a newer method should be preferred.

It is a from-scratch rewrite. The governing relationships are re-expressed in our own notation
and explained in our own words, drawn from the classical literature (Airy, Stokes,
Korteweg and de Vries, Isobe, Fenton, Hudson, Bodine, Holland, and others) rather than copied
from any single source.

**Built vs Documented.** Each tool's section is tagged one of two ways. A **Built** tool exists
in CHESS-QC and is validated against an independent worked example (the ACES *User's Guide*
problems or the original method paper); it carries a **Validation** note. A **Documented** tool
is one whose method is transcribed from the source and written up here, but whose application
is not yet implemented; it has no Validation note yet, and its inputs/outputs may change when
built. This reference therefore covers the full ACES method set, not only what ships today.

Companion documents, auto-generated from the application contracts and kept in sync by
`common/gen_docs.py`:

- `docs/USER_MANUAL.md` lists how to drive each tool (inputs, units, ranges, defaults, outputs).
- `docs/OUTDATED_APPS.md` collects the apps flagged below as superseded or screening-only.

**Conventions.** `d` is still-water depth; `η` is water-surface elevation measured up from the
still-water line (SWL, `z = 0`); `H` is wave height (crest to trough) and `a = H/2` is
amplitude; `T` is period; `L` is wavelength; `c = L/T` is phase speed (celerity); `k = 2π/L`,
`ω = 2π/T`; `g` is gravity; `ρ` is water density; `z` is measured up from SWL with the bed at
`z = -d`; `s = z + d` is height above the bed; `θ` is wave phase. The tools compute in SI
internally and display US or SI at the user's choice.

---

# Area 1 — Wave Prediction

These tools establish the design wave climate, namely the wind that drives the waves, the
waves a wind generates over a fetch, and the rare large wave a structure must survive, before
any transformation or force calculation. All four methods (1-1 wind/wave growth, 1-2
Beta-Rayleigh height distribution, 1-3 extremal analysis, and 1-4 constituent tide prediction)
are built and validated.

## 1-1 — Windspeed Adjustment and Wave Growth

**What it does.** It takes a single observed wind (speed, the height it was measured at, how it
was observed, the averaging time, the latitude) and does two things: first it converts that
observation into the standard design wind, then it estimates the significant wave height and
peak period the wind builds over a given fetch, in deep or shallow water, for an open coast or
a fetch restricted by surrounding land.

**The physics, briefly.** There are two stages. In the wind-adjustment stage the observation is
corrected for how and where it was taken (a ship-speed bias for shipboard reports), brought to
the standard 10 m height under neutral stability using the constant-stress logarithmic wind
profile (friction velocity and surface roughness solved together), adjusted from its observed
averaging duration to the design duration, and finally converted to the "adjusted wind" used by
the growth laws through a wind-stress (drag) relation. In the wave-growth stage the adjusted
wind feeds the SPM (1984) fetch-limited and duration-limited growth formulas; the smaller of the
two governs, capped at the fully developed sea. Shallow water adds depth limiting. A restricted
fetch is handled by selecting the wave-development direction that maximizes an effective-fetch
versus off-wind-angle trade-off computed from a table of radial fetch lengths.

> **Status and Caveats:** Current for screening, with two important limits.
> (1) Air-sea temperature (stability) defaults to neutral, with an opt-in physical correction.
> The stability function ACES prints (Technical Reference eq 9, unstable branch) is a corrupted
> copy of the standard Businger-Dyer and Paulson (1970) momentum function: it has the `arctan`
> and `π/2` terms sign-flipped and carries spurious extra terms. CHESS-QC therefore offers the
> canonical Businger-Dyer correction as an opt-in option, with the validated neutral path as the
> default. That option is deliberately not validated against ACES Example 3, because the ACES
> example cannot be reconciled with correct physics: it shows a `−3.2 %` change for `ΔT = −3`,
> which is both several times too large and opposite in sign to what standard surface-layer
> theory gives (a small increase for unstable air), and ACES's own equations are internally
> sign-inconsistent. No worked example isolates stability from the observation-type bias
> (Example 1 is shore and neutral, Example 2 is shipboard, Example 3 is overwater with `ΔT ≠ 0`),
> so the correction is shipped as a physically-standard option rather than a validated default.
> (2) The SPM-1984 growth formulas are themselves dated; they were revised after 1989 (for
> example Hurdle and Stive 1989, and the Coastal Engineering Manual). Use 1-1 for quick
> estimates, and prefer a modern spectral hindcast for design.

**Validation.** Reproduces ACES *User's Guide* Example 1-1 (shore-windward, open water, shallow,
with inputs near `25 ft / 45 mph / 26 mi / 13 ft`): equivalent-neutral wind `U_e = 46.4 mph`,
adjusted wind `U_a = 67.9 mph`, `H_mo = 4.2 ft`, `T_p = 4.77 s`. The restricted-fetch geometry
matches Example 3 (effective fetch `26.6 mi`, mean wave direction `93°`), and the restricted
growth reproduces `7.80 ft / 5.74 s` given the manual's adjusted wind.

**References.** SPM (1984) Ch. 3; Resio and Vincent (1977); Smith (1991) for restricted fetch;
Van Dorn (1953) for wind stress.

## 1-2 — Beta-Rayleigh Wave-Height Distribution

**What it does.** Starting from the energy-based significant wave height `H_mo`, the peak period
`T_p`, and the water depth, it estimates the statistical distribution of the individual wave
heights in a sea state and returns the characteristic heights (root-mean-square `H_rms`, median,
and the means of the highest third, tenth, and hundredth: `H_1/3`, `H_1/10`, `H_1/100`). In deep
water it collapses to the classical Rayleigh distribution.

**The physics, briefly.** In deep water the heights of individual waves in a narrow-band sea
follow the Rayleigh distribution (Longuet-Higgins 1952), so every characteristic height is a
fixed multiple of `H_rms`. In shallow water the tallest waves break, which truncates the upper
tail; the distribution becomes the depth-limited Beta-Rayleigh form (Hughes and Borgman 1987), a
Beta-shaped probability density on `0 < H < H_b` where `H_b ≈ 0.9 d` is the breaking height. The
two shape parameters are fixed by matching the sea's root-mean-square and root-mean-quad heights,
which are themselves tied to `H_mo` and the relative depth through empirical fits (Thompson and
Vincent 1985). Because there is no closed form for the characteristic heights, they are obtained
by numerically integrating the fitted density.

> **Status and Caveats:** Current. The method is valid only in the depth-limited regime
> (`d/(gT_p²) ≲ 0.01`); outside it the result reverts to pure Rayleigh. Results depend on the
> chosen breaking-height rule (`0.78 d` versus the ACES-preferred `0.9 d`) and on the empirical
> rms/rmq fits, so treat the tail heights (`H_1/100`) as estimates. One transcription point: the
> relative-depth fit in our equation notes was recorded with the argument inverted; reproducing
> the worked example requires `g T_p² / d`, not `d / (g T_p²)`.

**Validation.** Reproduces ACES *User's Guide* Example 1-2 (`H_mo = 5 ft, T_p = 6.30 s,
d = 10.2 ft`): `H_rms = 3.72 ft`, `H_med = 3.26 ft`, `H_1/3 = 5.18 ft`, `H_1/100 = 7.48 ft`,
plus the Rayleigh-limit ratios in deep water. The manual's `H_1/10 = 6.55 ft` is not reproduced:
the documented Beta-Rayleigh method gives `6.30 ft` (grid-independent and stable across
quadrature schemes), and the manual figure is physically inconsistent with the more strongly
truncated `H_1/3` and `H_1/100`, so it is taken to be a documentation or legacy-code artifact and
the computed value is reported.

**References.** Hughes and Borgman (1987); Longuet-Higgins (1952); Thompson and Vincent (1985);
Hughes and Ebersole (1987); SPM (1984).

## 1-3 — Extremal Significant Wave Height Analysis

**What it does.** Given a sample of storm-peak significant wave heights (one per event), it fits
extreme-value probability distributions and estimates the design wave height for chosen return
periods (for example the 50-year or 100-year wave) together with confidence bounds and the
chance of exceedance over a project lifetime.

**The physics, briefly.** Five candidate distributions are fitted: Fisher-Tippett Type I
(Gumbel) and Weibull with shape `k = 0.75, 1.0, 1.4, 2.0`. Each ranked height is assigned a
non-exceedance probability by Goda's (1988) plotting-position formula, converted to a reduced
variate, and the height-versus-variate relationship is fitted by least squares. The best
distribution is the one with the highest correlation (lowest residual). The return-period height
follows from the fitted line evaluated at the return-period variate, and the confidence band
uses Goda's normalized standard-deviation coefficients.

> **Status and Caveats:** Current. Selection follows the documented least-squares and
> correlation criterion. The original ACES screen also prints a "sum of squares" whose
> normalization is not recoverable from the published equations, so CHESS-QC reports the
> correlation coefficient (which it reproduces exactly) and selects by it. As with any extremal
> fit, results from a small sample have wide confidence bands; treat the band, not just the
> central estimate, as the answer.

**Validation.** Reproduces ACES *User's Guide* Example 1-3 (the EXTDELFT dataset, 15 heights,
`N_T = 20`, record `K = 20 yr`): every return-period height for all five distributions matches
the manual (FT-I 100-year `36.72 ft`, Weibull-2 `34.40 ft`), the correlations match (FT-I
`0.9813`, Weibull-2 `0.9866`), the 90 percent confidence band on the FT-I 100-year value is
`25.2 to 48.2 ft`, and the best fit is Weibull `k = 2.0`.

**References.** Goda (1988); Gringorten (1963); Gumbel (1958); EM 1110-2-1414.

## 1-4 — Constituent Tide Record Generation

**What it does.** Predicts a tide-elevation time series at a coastal site over a chosen period
from the classical harmonic method, given the amplitude and phase (epoch) of each tidal
constituent at that site.

**The physics, briefly.** The astronomical tide is the sum of many cosine harmonics, one per
constituent (the semidiurnal `M_2` and `S_2`, the diurnal `K_1` and `O_1`, and others), each with
its own fixed angular speed (frequency), amplitude, and phase. The predicted height at time `t`
is `h = H_0 + Σ f_n A_n cos(a_n t + (V_0 + u)_n − κ_n)`, summed over the constituents, where the
amplitude `A_n` and local epoch `κ_n` are site inputs, the speeds `a_n` come from the standard
constituent table, and the node factor `f_n` and equilibrium argument `(V_0 + u)_n` are computed
astronomically (Schureman 1971) from the start date and the gage longitude. The semidiurnal lunar
constituent `M_2` usually dominates.

> **Status and Caveats:** Current. It predicts only the astronomical tide; it does not include
> storm surge or meteorological effects (use the surge tools for those). One convention is worth
> recording because it is not stated in the source and was recovered by matching the worked
> example: the equilibrium-argument solar-time term uses `T = 15 H` measured from local midnight,
> and the slow astronomical longitudes are evaluated at `UT = H − longitude_west / 15`. The full
> 37-constituent Table A-5 is supported; the minor constituents outside the worked example use the
> standard Schureman node factors.

**Validation.** Reproduces ACES *User's Guide* Example 1-4 (Buzzards Bay Entrance, MA: 25
constituents, start 1989-01-10 10:00, 120-hour record, gage longitude 70.62° W, datum offset
1.79 ft): the predicted elevations match the published Table 1-4-1 to within 0.04 ft over the
whole record (for example `4.26, 4.35, 4.39 ft` at the first three quarter-hours, and `0.01,
−0.08 ft` at the last two), with no fitted parameters.

**References.** Schureman (1971); Harris (1981); EM 1110-2-1414 Ch. 2.

---

# Area 2 — Wave Theory

These three tools all answer the same question, namely given a wave (height, period, depth, and
a position in the wave) what are its length, speed, energy, water-particle motion, and pressure.
They differ in fidelity and in the range of conditions they cover. The choice among them is
mostly about how nonlinear the wave is, measured by the Ursell number `U_r = H L² / d³` (small
means gentle or deep, large means steep or shallow):

| Tool | Best when | Idea | Fidelity |
| --- | --- | --- | --- |
| 2-1 Linear (Airy) | `U_r` below about 26, gentle waves, any depth | sinusoidal, small amplitude | first order |
| 2-2 Cnoidal | `U_r` above about 26, shallow water, nonlinear | elliptic-function (sharp crest, flat trough) | first or second order |
| 2-3 Fourier (Fenton) | almost anything up to near breaking | numerical steady wave | "exact" (N terms) |

A practical rule of thumb: start with Linear; if the Ursell number is large (the tool warns
you), switch to Cnoidal in shallow water, or use Fourier when you want the most accurate single
answer across the whole range.

A fourth tool, 2-4 Wave Parameters, is a companion to the Linear theory rather than a new
theory. It presents the same first-order solution as the standard dimensionless parameter table
(the functions of relative depth that engineers used to read from printed tables) and adds two
field conveniences: specifying the wave by period or by frequency, and recovering wave height
from a pressure-gauge record.

## 2-1 — Linear (Airy) Wave Theory

**What it does.** It treats the wave as a small, perfectly sinusoidal undulation and returns the
first-order ("Airy") description of its motion: wavelength, celerity, group speed, energy
density and flux, the Ursell number, and, at a chosen depth and phase, the surface elevation,
sub-surface pressure, particle displacements, velocities, and accelerations. It is the workhorse
of coastal engineering and the starting point for almost every other method.

**The physics, briefly.** The water is taken as inviscid and the flow irrotational, so the
velocity field comes from a potential `φ` that satisfies Laplace's equation. Three boundary
conditions close the problem: no flow through the bed, and two conditions at the free surface
(the surface moves with the fluid, and the pressure there is atmospheric). For small waves those
surface conditions are linearized about `z = 0`, and the solution is a single sinusoid.

**Key relationships (our notation).**

- Dispersion, the link between speed and length:
  ```
  ω² = g k tanh(kd)          (equivalently  c² = (g/k) tanh(kd))
  ```
  This is implicit in `k`. The tool solves it with Hunt's (1979) explicit nine-term Pade
  approximation, accurate to better than 0.01 percent at any depth, avoiding iteration:
  ```
  c² = g d · [ y + (1 + Σ_{n=1..9} dₙ yⁿ)⁻¹ ]⁻¹ ,   y = ω² d / g
  ```
- Wavelength `L = cT`; group (energy) speed `C_g = (c/2)[1 + 2kd/sinh(2kd)]`.
- Surface `η = (H/2) cos θ`; energy density `E = ρ g H² / 8`; energy flux `P = E C_g`.
- Pressure and kinematics (displacement, velocity, acceleration) follow as the usual
  `cosh(ks)/sinh(kd)` and `sinh(ks)/sinh(kd)` depth-decay forms evaluated at `s = z + d`.
- Nonlinearity gauge `U_r = H L² / d³`.

**Assumptions and range.** Small steepness and small `H/d`; inviscid, irrotational flow;
constant depth; a single sinusoidal component. Particle velocities below a wave trough are
evaluated by clamping the elevation to the SWL, because the linear profile is undefined above
the trough. Accuracy degrades as the wave steepens or shoals.

> **Status and Caveats:** Current and foundational. Linear theory remains the universal first
> cut and underlies wave transformation, forces, and spectra. However it loses accuracy for
> steep or shallow waves: when `U_r` exceeds about 26 the tool flags the result, and you should
> move to 2-2 Cnoidal in shallow water or 2-3 Fourier. It also says nothing about wave skewness
> or asymmetry (real waves have sharper crests than troughs); use 2-2 or 2-3 if those matter.

**Validation.** Reproduces ACES *User's Guide* Example 2-1 (`H = 6.30 ft, T = 8 s, d = 20 ft`)
to about 0.1 percent: `L = 189.90 ft`, `c = 23.74 ft/s`, `C_g = 20.87 ft/s`, `U_r = 28.4`, and
the sub-surface pressure and kinematics at the example point.

**References.** Airy (1845); Hunt (1979); SPM (1984) Ch. 2.

## 2-2 — Cnoidal Wave Theory

**What it does.** It describes finite-amplitude periodic waves in shallow water, where the real
profile is far from sinusoidal: long flat troughs separated by short peaked crests. It returns
the same family of outputs as 2-1 (length, celerity, energy, pressure, kinematics, profiles) to
first order (Isobe 1985) or second order (Hardy and Kraus 1987), selected by the user.

**The physics, briefly.** "Cnoidal" comes from the Jacobian elliptic cosine `cn`: the surface is
proportional to `cn²`. The theory descends from the Korteweg and de Vries balance between
nonlinear steepening and dispersion, and is written with an auxiliary elliptic modulus `κ` and
the complete elliptic integrals `K(κ)` and `E(κ)`. As `κ` goes to 0 it collapses to linear
theory; as `κ` goes to 1 the crest spacing grows without bound and it approaches a single
solitary wave.

**Key relationships (our notation).** A perturbation parameter `ε = H/d` and the modulus `κ` are
fixed by the first-order dispersion relation:
```
16 κ² K(κ)² / 3 = g H T² / d²
```
(the second-order form multiplies the right side by a small `ε`-correction). With `κ` known, the
surface, celerity, energy, pressure, and velocities are polynomials in `cn²θ` with
elliptic-integral coefficients, for example `η = d (A₀ + A₁ cn²θ)` with `A₀ = ε(λ - μ)`,
`A₁ = ε`, where `λ = κ'²/κ²` and `μ = E/(κ²K)`. The tool solves for `κ` by bisection (the left
side is monotonic) and evaluates `cn, sn, dn` from the modulus.

**Assumptions and range.** Shallow water and long waves (the regime where 2-1 fails); finite but
not breaking amplitude. Best for `U_r` above about 26.

> **Status and Caveats:** Current. This is the correct theory for nonlinear shallow-water waves
> and the natural companion to 2-1. Two things to know.
> (1) Transcription correction (carried in our equation ledger): the original Technical
> Reference's two vertical-acceleration expressions lost the squares on their trigonometric
> factors. CHESS-QC computes those accelerations by analytic differentiation of the correct
> velocity field, which both regenerates the correctly-printed equations and fixes the error,
> confirmed against the manual's `∂w/∂t`.
> (2) Near `κ` equal to 1 the wave is effectively solitary and the period loses meaning; near
> `κ` equal to 0 you are better off in linear theory. Outside shallow water, prefer 2-3 Fourier.

**Validation.** Reproduces ACES *User's Guide* Example 2-2 (first order, `H = 10 ft, T = 15 s,
d = 25 ft`): `L = 455.74 ft`, `c = 30.38 ft/s`, `U_r = 132.9`, the crest and trough elevations,
and the point kinematics and pressure, all to about 0.1 percent.

**References.** Korteweg and de Vries (1895); Isobe (1985); Hardy and Kraus (1987); Abramowitz
and Stegun (1972) for the elliptic functions.

## 2-3 — Fourier Series Wave Theory (Fenton)

**What it does.** It computes a steady progressive wave of permanent form to essentially
arbitrary accuracy by representing it as an N-term Fourier series (up to 25 terms) and solving
the full nonlinear free-surface conditions numerically. It is the most accurate of the three
across the widest range, from deep water to near-breaking shallow water, and it additionally
returns integral properties (impulse, kinetic and potential energy, radiation stress, wave
power, mean flows) that the lower-order theories cannot.

**The physics, briefly.** A stream function is written as a sum of N cosine modes whose vertical
structure satisfies Laplace's equation and the no-flow bed condition automatically. The
remaining unknowns (the modal coefficients, the surface elevation at a set of evenly spaced
points, the wavenumber, and the flow and energy constants) are found by forcing the two
nonlinear free-surface conditions to hold exactly at those points, together with the
wave-height and mean-level constraints. This yields a square system of nonlinear equations
solved by Newton iteration, made robust by wave-height ramping: the height is raised in steps so
each solve starts near the previous answer, avoiding divergence near breaking.

**Key relationships (our notation).** The unknowns are gathered into one vector and the system
`F(z) = 0` (free-surface kinematic and dynamic conditions at each collocation point, plus
height, mean-zero, and the celerity or current condition) is driven to zero. CHESS-QC uses a
compact non-dimensional form (scaling by depth and `√(gd)`) with a finite-difference Jacobian and
`numpy`'s linear solver, so no external solver is needed. A useful internal check: with `N = 1`
it reduces to linear (Airy) theory.

> **Datum note.** The dynamic free-surface (Bernoulli) constant is carried at the still-water
> datum (`r`); the value reported as "Bernoulli constant `R`" adds the `gd` baseline, so
> `R = r + gd`. This matches the original worked example exactly and is documented in the module.

**Assumptions and range.** A single steady wave of permanent form (no spectrum, no wave-current
shear beyond a uniform current); irrotational flow; constant depth. Either of the two classical
celerity definitions (Eulerian mean current, or Stokes mass transport) may be specified.

> **Status and Caveats:** Current and the highest fidelity of the three. Prefer this tool when
> one accurate answer is wanted across the whole depth and steepness range, or when the integral
> wave properties are needed. It is numerical, so it can fail to converge for waves driven past
> the physical breaking limit; keep `H` below the Fenton `H_max` feasibility check the tool
> applies. The very-deep-water branch uses exponential kernels to avoid `cosh` and `sinh`
> overflow.

**Validation.** Reproduces ACES *User's Guide* Example 2-3 (`H = 4.5 ft, T = 9 s, d = 22 ft`,
Eulerian, `N = 16`): `c = 25.62 ft/s`, `L = 230.58 ft`, crest `η = 2.80 ft`, and the full
integral-property table (impulse, energies, radiation stress, wave power) all reproduced; at the
example point `U = 3.15 ft/s`, `a_z = -1.30 ft/s²`, `p = 473.2 psf`.

**References.** Rienecker and Fenton (1981); Fenton (1988a, 1988b, 1990).

## 2-4 — Wave Parameters

**What it does.** Presents linear (Airy) wave theory as the classic dimensionless wave-parameter
table and adds two practical conveniences. The table reports, for the given wave and depth, the
relative depth in several forms, the hyperbolic functions of `kd`, the group-to-phase speed ratio,
the shoaling coefficient, the wave-steepness ratio, the pressure-response factor at the chosen
depth and at the bed, the energy and power, and the second-order mass-transport (drift) velocity.
The two conveniences are: specify the wave by period or by frequency, and recover the wave height
from a measured pressure amplitude (the pressure-transducer use case).

**The physics, briefly.** The engine is identical to 2-1: the same explicit dispersion relation
and the same first-order kinematics. What 2-4 adds is the set of dimensionless ratios that the
Shore Protection Manual tabulated as functions of relative depth. The shoaling coefficient is
`K_s = H/H_0 = sqrt(C_g0 / C_g)`, the ratio of deepwater to local group speed. The
pressure-response factor `K_p(z) = cosh k(z + d) / cosh kd` describes how the wave pressure
signal weakens with depth; the height inversion runs that factor backward, so a measured dynamic
pressure amplitude `p_d` at gauge elevation `z` gives `H = 2 p_d / (ρ g K_p)`. The mass-transport
velocity is the second-order Stokes drift, the small net downwave current the wave leaves behind.

**Assumptions and range.** Same first-order linear theory as 2-1, so the same gentle-wave
assumption applies. The shoaling coefficient assumes refraction is handled separately (it folds
no angle change into `H_0`). The pressure inversion assumes the gauge sits within the water column
and that the measured signal is the wave (dynamic) part only.

> **Status and Caveats:** Current. This is a presentation and utility layer over linear theory,
> identical in fidelity to 2-1, so the same Ursell-number caveat applies: when the wave is strongly
> nonlinear the linear parameters are approximate, and 2-2 or 2-3 should be preferred.

**Validation.** Because 2-4 shares the validated 2-1 core, it reproduces the shared linear-theory
outputs of ACES *User's Guide* Example 2-1 (`H = 6.30 ft, T = 8 s, d = 20 ft, z = -12 ft`):
`L = 189.90 ft`, `C = 23.74 ft/s`, `C_g = 20.87 ft/s`, `p = 767.83 psf`. The height inversion is
checked by a round trip (compute the pressure amplitude from a known height, then recover the
height), and the period and frequency inputs are checked to give identical results.

**References.** Airy (1845); SPM (1984) Appendix C; Hunt (1979).

## 2-5 — Solitary Wave Theory

**What it does.** Returns the kinematics and integral properties of a solitary wave, a single wave
of translation that lies entirely above the still-water level with no trough (the limiting form of
a very long wave, and a useful model for tsunamis and surge bores): the celerity, the surface
profile, particle velocities, the dynamic pressure, the total energy, the McCowan-Munk coefficients
M and N, and the breaking height.

**The physics, briefly.** A solitary wave has infinite length and period, so a single parameter, the
height-to-depth ratio, specifies it. To lowest order the surface is a squared hyperbolic secant
(Boussinesq), the celerity is the square root of g times the sum of depth and height, and the energy
grows with the three-halves power of both height and depth. The McCowan-Munk coefficients M and N,
which ACES reports, follow from a pair of transcendental relations in the height-to-depth ratio.

> **Status and Caveats:** Current (Built), validated analytically only. This application has no
> ACES Technical-Reference chapter and no User's-Guide worked example; it is based on the Shore
> Protection Manual solitary-wave theory, which was not available, so the implementation uses the
> standard, well-established theory from the Coastal Engineering Manual (EM 1110-2-1100, Part II-1)
> together with the classical McCowan-Munk coefficients. Because there is no ACES numeric oracle,
> the values are standard solitary-wave theory and have not been cross-checked against an ACES run.

**Validation.** Analytic, against the closed-form relations and their known limits: the celerity is
the square root of g(d+H) (the first-order McCowan value, confirmed against Zaroodny 1972, "McCowan's
Solitary Wave Expansions," and the OpenFOAM McCowan wave model); the crest elevation equals the wave
height; the McCowan flat-bed breaking height is 0.78 times the depth; the coefficient M approaches the
square root of three times the height-to-depth ratio as that ratio vanishes; and the total energy
scales as height and depth to the three-halves power. The McCowan-Munk M and N relations were
cross-checked against the OpenFOAM McCowan model (whose m equals M divided by the depth).

**References.** McCowan (1891, 1894); Munk (1949); Boussinesq (1872); Zaroodny (1972); SPM (1984);
CEM (EM 1110-2-1100, Part II-1).

---

# Area 3 — Wave Transformation

These tools follow a wave as it moves from deep water toward the shore and changes height and
direction through refraction, shoaling, and eventually breaking.

## 3-1 — Linear Wave Theory with Snell's Law

**What it does.** It transforms a wave of known height, period, and crest angle at one depth both
out to deep water and in to a chosen "subject" depth, and reports the height, crest angle,
length, celerity, group velocity, energy density, and energy flux at all three locations, plus
the Ursell numbers, the deepwater steepness, and the depth and height at which the wave breaks.

**The physics, briefly.** Three classical effects act together over straight, parallel depth
contours. Refraction bends the crest by Snell's law, `sin(α)/c` constant, so the wave turns
toward shore-parallel as it slows in shallower water. Shoaling changes the height by conserving
energy flux, `E C_g cos(α)` constant along a ray, so the height grows as the group speed drops.
Both use the same Hunt (1979) dispersion solver as 2-1. The breaking limit is reported from the
Singamsetti and Wind (1980) breaker-height relation and the Weggel (1972) breaker-depth relation
(beach slope dependent).

**Assumptions and range.** Linear (small-amplitude) transformation; monochromatic wave; straight
and parallel bottom contours, so refraction is purely by Snell's law with no caustics or
crossing rays. The breaker formulas are empirical fits for plane beaches.

> **Status and Caveats:** Current. Suitable for first-order transformation over a regular shelf.
> It does not capture diffraction, reflection, wave-current interaction, irregular bathymetry, or
> directional spreading; for those, use a spectral wave-transformation model. Treat the breaker
> height and depth as screening estimates.

**Validation.** Reproduces ACES *User's Guide* Example 3-1 (`H = 10 ft, T = 7.5 s` at `d = 25 ft`,
crest angle `10°`, to a subject depth `20 ft`): subject height `10.27 ft`, the three-location
length, celerity, and group-velocity table, and the breaker height `12.29 ft` at breaker depth
`15.25 ft`, all to about 0.1 percent.

**References.** O'Brien (1942) for Snell's law in refraction; Weggel (1972) and Singamsetti and
Wind (1980) for breaking; Hunt (1979).

## 3-2 — Irregular Wave Transformation (Goda's Method)

**What it does.** Transforms a whole irregular sea, a spectrum of wave heights and directions
rather than a single wave, across the surf zone, and returns the transformed wave-height
statistics (significant, mean, root-mean-square, `H_1/10`, maximum), the effective shoaling and
refraction coefficients, the wave setup, and the surf-beat amplitude. Where 3-1 transforms one
monochromatic wave, this tool transforms the natural random sea.

**The physics, briefly.** Goda's method (1975, 1984) describes the sea by a
Bretschneider-Mitsuyasu frequency spectrum with a Mitsuyasu directional spread, then computes
spectrum-averaged ("effective") shoaling and refraction coefficients by integrating the
individual-component coefficients weighted by the spectrum. Depth-limited breaking is applied as
a clip on the upper tail of the (Rayleigh) individual-wave-height distribution, with the energy
redistributed; the incipient breaking height follows Goda's empirical relation in beach slope and
relative depth. Near breaking the shoaling is taken as nonlinear (Shuto 1974). Wave setup follows
from the cross-shore radiation-stress gradient, and a surf-beat amplitude is estimated empirically.
The output statistics come from numerically integrating the transformed height distribution.

> **Status and Caveats:** Current (Built), with a documented accuracy limitation. The method
> assumes straight, parallel contours and a narrow-band spectrum, and the ACES procedure restricts
> it to peak period below about 16 s, depth above about 10 ft, and principal incidence within 75° of
> the shore normal. The directional integral uses the Mitsuyasu spreading parameter s_max (default
> 10, wind waves).
>
> The transformation physics is implemented in full, but it does not reproduce the worked example to
> the digit. The shoaling coefficient (0.9133) is exact, and the at-depth significant, mean, and
> root-mean-square heights match to about four percent (significant 17.0 versus 17.7 ft) with surf
> beat to about one percent. Two pieces need source-level detail absent from the public Technical
> Reference: the high quantiles (the highest tenth, highest two percent, and maximum) and a small
> bias on the significant height follow Goda finite-N order-statistic relations (an effective wave
> count, roughly 1200 here) rather than the asymptotic Rayleigh ratios the distribution integration
> yields; and the effective refraction coefficient (computed 0.950 versus the published 0.964)
> depends on the exact directional-integration scheme and the default spreading parameter. The
> deepwater reference column and the high quantiles therefore carry larger residuals (up to about six
> percent). This is the same class of empirical-coefficient gap documented for application 5-4.

**References.** Goda (1975, 1984); Mitsuyasu (1975) for directional spread; Shuto (1974) for
nonlinear shoaling.

## 3-3 — Combined Diffraction and Reflection by a Vertical Wedge

**What it does.** Estimates how the wave height is modified near a fixed structure, such as a
breakwater, jetty, or quay wall represented as a fully reflecting wedge, by the combined action of
diffraction and reflection. It returns the modification factor (modified height divided by
incident height) and the phase at a chosen point; a companion gridded variant evaluates the same
quantities over a map around the structure.

**The physics, briefly.** The linear, monochromatic wave field around a wedge of arbitrary
included angle is solved by an eigenfunction expansion in Bessel functions (Chen 1987), in
cylindrical coordinates centered on the wedge apex. The semi-infinite straight breakwater is the
special case of zero wedge angle (the classical Penny and Price / Wiegel solution). The complex
horizontal potential gives, at any point, the amplitude (which is the modification factor) and the
phase relative to the incident wave; the series is truncated once successive Bessel terms fall
below a small tolerance.

> **Status and Caveats:** Current (Built; the single-point tool). It assumes linear, monochromatic,
> unidirectional waves over constant depth and a perfectly reflecting, vertical structure; real
> structures are partially absorbing and the bathymetry is rarely flat, so treat the result as an
> idealized upper bound on the modification. Fractional-order Bessel functions are evaluated by
> series. The gridded variant (which maps the same quantities) is not yet built.

**Validation.** Reproduces the ACES *User's Guide* Example 1 (a semi-infinite breakwater, wedge
angle zero so nu = 2; incident height 2 ft, 6 s wave in 12 ft of water, wave angle 133 degrees,
field point at (33, -17) ft): wavelength 109.82 ft and a modification factor of 0.58, giving a
modified wave height of 1.16 ft, all matching to better than one percent. The phase comes out
-2.48 rad versus the published -2.58; the roughly 0.1 rad offset is a phase-reference convention in
the original PCDFRAC code and does not affect the modification factor or the modified height.

**References.** Chen (1987); Penny and Price (1952); Wiegel (1962); Stoker (1957); SPM (1984).

## 3-4 — Vertical-Wedge Diffraction/Reflection on a Uniform Grid

**What it does.** The gridded companion to 3-3: it evaluates the same fully-reflecting-wedge
solution at every node of a uniform X-Y grid around the wedge apex, producing maps of the
wave-height modification factor, the modified wave height, and the phase relative to the incident
wave (a contour map of how a structure reshapes the wave field around it).

**The physics, briefly.** Identical to 3-3 (the Chen 1987 Bessel eigenfunction expansion), simply
evaluated over a grid of field points instead of a single one.

> **Status and Caveats:** Current (Built). Shares 3-3's solver, so the same assumptions apply
> (linear, monochromatic, unidirectional waves over constant depth; a perfectly reflecting vertical
> structure). The output is a two-dimensional field; the generic GUI renderer reports scalar
> summaries (the field maxima/minima), with the full grid available from the compute layer.

**Validation.** Reproduces the ACES *User's Guide* Example 3 (semi-infinite breakwater, wedge angle
zero; a 4 ft, 12 s wave in 30 ft of water at 52 degrees) on its X grid from -600 to 200 ft in 200 ft
steps and Y grid from -400 to 200 ft in 100 ft steps: the wavelength is 356.85 ft, and at the grid
point (-600, -400) ft the modification factor is 0.903 and the modified height 3.61 ft, matching the
published Table 3-3-2 to better than one percent.

**References.** Chen (1987); Penny and Price (1952); Wiegel (1962); Stoker (1957); SPM (1984).

---

# Area 4 — Structural Design

These tools size or load common coastal structures: rubble-mound armor, the toe that keeps a
wall from undermining, and the standing-wave forces on a vertical wall.

## 4-1 — Breakwater Design (Hudson)

**What it does.** It sizes the primary armor units of a rubble-mound breakwater or revetment from
the Hudson stability equation, and reports the required individual armor-unit weight, the minimum
crest width, the cover-layer thickness, and the number of armor units placed per unit area.

**The physics, briefly.** The Hudson (1953 to 1961) equation balances the destabilizing wave
force against the submerged weight of an armor unit. The required weight grows with the cube of
the design wave height and falls with the stability coefficient `K_D` (which encodes armor shape,
placement, and whether the structure trunk or head is breaking or non-breaking) and a slope
factor. Crest width, layer thickness, and placement density follow from layer-coefficient and
porosity relations.

> **Status and Caveats:** Current for preliminary sizing, but a newer method is preferred.
> Hudson is a single design-wave-height empirical formula: it ignores wave period, storm
> duration (number of waves), structure permeability, and an explicit damage level. **Van der
> Meer (1988)** stability formulae account for those factors and are preferred for detailed
> design; Hudson remains useful for quick estimates and for consistency with legacy studies.

**Validation.** Reproduces ACES *User's Guide* Example 4-1 (tribar trunk, non-breaking, design
height `11.5 ft`): armor weight `1.59 tons`, crest width `8.21 ft`, cover-layer thickness
`5.47 ft`, and `130.3` units per 1000 ft².

**References.** Hudson (1953 to 1961); SPM (1984) Ch. 7; EM 1110-2-2904; and Van der Meer (1988)
as the preferred successor.

## 4-2 — Toe Protection Design

**What it does.** It designs the toe-apron width and toe-stone weight for a vertical wall,
bulkhead, or revetment, the stone blanket that prevents scour from undermining the structure.

**The physics, briefly.** The apron width is taken as the larger of a geotechnical width (the
Rankine passive wedge that must be confined) and two hydraulic minima tied to the water depth at
the toe. The toe-stone weight uses the Tanimoto, Yagyu, and Goda (1982) stability number in a
Hudson-form sizing equation, evaluated with the local wave conditions from the Hunt (1979)
dispersion solver.

> **Status and Caveats:** Current. Appropriate for the standard vertical-wall and revetment toe
> cases covered by the ACES procedure. As with all rubble sizing, it is empirical; verify
> against site-specific physical or numerical modeling for unusual geometries or extreme waves.

**Validation.** Reproduces ACES *User's Guide* Example 4-2 (case 1, bulkhead toe): apron width
`15.00 ft`, toe-stone weight `12.99 lb`, and water depth at the top of the toe layer `15.50 ft`.

**References.** EM 1110-2-1614; Tanimoto, Yagyu, and Goda (1982); Hunt (1979).

## 4-3 — Nonbreaking Wave Forces at Vertical Walls

**What it does.** It computes the standing-wave (clapotis) forces and overturning moments on a
vertical wall for both the wave crest and the wave trough at the wall, by two classical methods:
Sainflou (1928) and Miche-Rundgren (Miche 1944, Rundgren 1958).

**The physics, briefly.** A non-breaking wave reflecting off a wall forms a standing wave whose
crest pushes harder and whose trough pulls back. The wall surface elevation (height of crest or
trough above the bed) is taken from the Miche-Rundgren second-order result and used by both
methods. The horizontal force and overturning moment are obtained by integrating the Lagrangian
pressure over the wall in 90 increments, weighting each at-rest increment by its stretched
(elevated) thickness, which captures the nonlinear setup of the pressure under the crest.

**Assumptions and range.** Non-breaking, fully reflected standing waves at a vertical,
impermeable wall on a flat bottom. It does not apply to breaking or broken waves, which deliver
much larger short-duration impact (shock) loads.

> **Status and Caveats:** Current for the non-breaking case. For breaking-wave impact forces use
> a dedicated method (for example Goda's formula or Minikin), not this tool; the Technical
> Reference notes that smaller (force) values are recommended when the geometry is near the
> breaking limit.

**Validation.** Reproduces ACES *User's Guide* Example 4-3 (`d = 15 ft, H = 8 ft, T = 10 s`):
crest and trough heights above the bottom `32.95 ft` and `16.95 ft`, and the Miche-Rundgren and
Sainflou crest and trough forces and moments, all to about 0.1 percent.

**References.** Sainflou (1928); Miche (1944); Rundgren (1958); Hunt (1979).

## 4-4 — Rubble-Mound Revetment Design

**What it does.** Sizes a sloping rubble-mound revetment under irregular waves: the median armor
stone weight, the bedding (filter) stone size, the armor and filter layer thicknesses, the stone
gradation (the range of stone sizes), and the expected and conservative wave runup. It evaluates
two stability methods side by side and uses the more demanding of the two.

**The physics, briefly.** The armor weight comes from the Hudson-form stability relation (weight
grows with the cube of the design height and falls with the cube of the stability number), but the
stability number is taken as the larger of two: the CERC zero-damage number (Ahrens 1981), a
simple function of slope, and the Dutch van der Meer (1988) number, which additionally accounts
for the wave period through the surf-similarity (Iribarren) parameter, the storm duration through
the number of waves, the core permeability, and an explicit allowable damage level. Plunging and
surging wave conditions select different van der Meer branches, switching at the surf-similarity
value where they cross. Layer thicknesses and gradation follow standard ratios, and the maximum
runup uses the Ahrens and Heimbaugh (1988) relation in the surf-similarity parameter.

> **Status and Caveats:** Current (Built). This is the revetment counterpart to 4-1 Hudson, and
> unlike 4-1 it already incorporates the modern van der Meer formulas, so it directly addresses the
> limitation flagged there. The van der Meer relations apply to rock armor on a uniform slope with
> little or no overtopping; the ACES procedure applies a fixed 1.2 shallow-water correction and a
> conservative number of waves (7000). Seawater unit weight (64 lb/ft^3) sets the specific gravity.

**Validation.** Reproduces the ACES *User's Guide* worked example (a 1-on-2 slope in 9 ft of water
under a 5 ft, 10 s sea, P = 0.1, S = 2): the CERC stability number governs, giving a median armor
weight of about 2500 lb, a 4.95 ft armor layer and 1.24 ft filter layer, and irregular-wave runup of
10.96 ft (expected) and 13.79 ft (conservative). The runup and the gradation ratios match the
published values exactly; the stone weights match to about half a percent, the residual being the
cubed rounding of the Rayleigh H-one-tenth-to-H-significant factor that the Technical Reference quotes
as 1.27. Two equation-summary errors were corrected against the source: the zero-moment height is the
*smaller* of the depth- and steepness-limited values, and the steepness-limited form uses the peak
period.

**References.** Ahrens (1981); van der Meer (1988a, 1988b); van der Meer and Pilarczyk (1987);
Battjes (1974); Ahrens and Heimbaugh (1988); EM 1110-2-2300.

---

# Area 5 — Wave Runup, Transmission, and Overtopping

This area covers how far waves climb up a beach or structure (runup), how much water passes over
the crest of a low structure (overtopping), how much wave energy gets through a breakwater to the
protected side (transmission), and how much the mean water level itself rises in the surf zone
(setup). All five chapters (5-1 through 5-5) are Built and validated against the ACES worked
examples (5-4's reflection coefficient carries a documented caveat; its transmitted height matches).

## 5-1 — Irregular Wave Runup on Beaches

**What it does.** Estimates how high irregular (real, many-period) waves run up a sloping beach,
reporting several runup statistics: the maximum, the level exceeded by 2 percent of waves, the
average of the highest tenth, the average of the highest third, and the mean.

**The physics, briefly.** Each runup statistic follows a simple power law in the surf-similarity
(Iribarren) number, which combines the beach slope with the deepwater wave steepness. Steeper
slopes and longer (less steep) waves produce higher runup. The power-law coefficients for each
statistic come from the Mase (1989) laboratory fit.

> **Status and Caveats:** Current (Built). The statistic coefficients are taken from Mase (1989)
> directly; they are cited in the ACES source but not tabulated there. Valid for beaches (uniform
> impermeable slopes), not for armored structures (see 5-2).

**Validation.** Reproduces the ACES *User's Guide* worked example (5-1-4) exactly: for a deepwater
significant height of 4.60 ft, a peak period of 9.50 s, and a foreshore cotangent of 13 (Iribarren
number 0.771), the five statistics come out to 8.74, 7.11, 6.50, 5.29, and 3.38 ft (maximum, two
percent, highest tenth, highest third, and mean), matching the published values to 0.01 ft. The
ordering of the statistics and their power-law form are checked as closed-form consistency tests.

**References.** Hunt (1959); Mase (1989); Walton and Ahrens (1989).

## 5-2 — Wave Runup and Overtopping on Impermeable Structures

**What it does.** Computes wave runup on a smooth or rough impermeable slope (such as a seawall or
a revetment face) and the resulting overtopping rate when the runup exceeds the crest, for both
monochromatic (single-period) and irregular waves, with an optional onshore-wind correction.

**The physics, briefly.** Runup on a rough slope uses an empirical relation in the surf-similarity
parameter with coefficients that depend on the armor type; smooth slopes use a separate set of
relations split by surf regime (plunging, transition, nonbreaking). Once runup is known, the
overtopping rate follows the Weggel (1976) and Saville form: it grows with the freeboard ratio (how
far the runup overshoots the crest) and with two empirical coefficients the user supplies from
published figures. The wind correction increases overtopping for onshore winds. For irregular
waves the runup is treated as Rayleigh distributed and the overtopping is averaged over the
distribution.

> **Status and Caveats:** Current (Built). The armor-type runup coefficients (a, b) and the
> overtopping coefficients (alpha, Q*0) are user inputs drawn from the ACES appendix and the SPM
> figures, not computed internally; an optional average alpha = 0.06 - 0.01431 sin(theta) is offered.
> The incident height at the structure toe is deshoaled to deepwater by linear theory for the
> overtopping term.

**Validation.** Reproduces the ACES *User's Guide* worked examples (Examples 1 through 7) on a
1-on-3 structure in 12.5 ft of water under a 7.5 ft, 10 s wave: rough-slope runup 9.421 ft and
smooth-slope runup 21.366 ft (transition regime); monochromatic overtopping with a known runup of
15 ft is 3.565 and of 20 ft is 5.368 ft^3/s/ft; the monochromatic rough-slope runup-plus-overtopping
case gives 0.829 ft^3/s/ft; and the irregular-wave case (Rayleigh runup distribution averaged over
199 quantiles) gives 0.287 ft^3/s/ft. All match the published values to better than 0.01 ft and
0.01 ft^3/s/ft.

**References.** Ahrens and McCartney (1975); Ahrens and Titus (1985); Weggel (1976); Douglass
(1986); SPM (1984).

## 5-3 — Wave Transmission on Impermeable Structures

**What it does.** Estimates the height of waves transmitted past a low impermeable structure by
overtopping, for sloped structures and for vertical or composite (wall-on-mound) structures.

**The physics, briefly.** The transmitted height is the incident height times a transmission
coefficient. For sloped structures the coefficient depends on the freeboard relative to runup and
on the crest width relative to structure height (Seelig 1980). For vertical and composite
structures it follows a Seelig (1976) relation built around the freeboard-to-height ratio with
empirical shape coefficients that depend on the breakwater type and berm depth.

> **Status and Caveats:** Current (Built). The vertical and composite coefficients are valid only
> within a stated depth-to-height range (0.145 to 0.5); outside it a note is issued and the
> transmission coefficient is clamped to its zero or unity limits. Sloped-structure runup reuses the
> 5-2 rough and smooth methods.

**Validation.** Reproduces the ACES *User's Guide* worked examples (Examples 1 through 4): a sloped
structure with a known 15 ft runup transmits 2.275 ft; a vertical wall over a submerged berm
(4.5 s wave, 20 ft depth, 12 ft crest, 6 ft berm) transmits 3.798 ft; a rough riprap slope (runup
9.421 ft) transmits 1.601 ft; and a smooth slope (runup 22.436 ft) transmits 2.652 ft. All match the
published transmitted heights to 0.01 ft.

**References.** Seelig (1976, 1980); Cross and Sollitt (1971); Goda (1969).

## 5-4 — Wave Transmission through Permeable Structures

**What it does.** Estimates the total wave height transmitted past a permeable (rubble) breakwater,
combining the part that passes over the crest by overtopping with the part that seeps through the
porous core.

**The physics, briefly.** The total transmission coefficient is the root-sum-square of an
overtopping coefficient (the same Seelig 1980 form as in 5-3) and a through-transmission
coefficient. The through part uses the Madsen and White (1976) hydraulic model: the layered
trapezoidal structure is first reduced to a hydraulically equivalent rectangle, then long-wave
theory with a linearized (Lorentz) friction term gives the transmitted and reflected amplitudes.
The friction factor itself depends on the flow through the porous medium, so the solution is
iterative.

> **Status and Caveats:** Current (Built), with one documented limitation on the reflection
> coefficient. This is the most involved ACES routine: the multilayer geometry, the
> equivalent-rectangle reduction, the iterative internal-friction model, and a Bessel-function
> seaward-slope solution all have to be reproduced. Complex-argument Bessel functions are
> hand-implemented (numpy plus standard library only). The equivalent-breakwater reference diameter
> is taken as one half the median material diameter (the representative material). Best for periodic,
> relatively long, normally incident, non-breaking waves.
>
> The reflection coefficient is approximate. The transcribed seaward-slope equations give near-total
> reflection for the long-period worked-example wave (the slope barely dissipates a 20 s wave), so
> the model over-predicts the reflection coefficient (about 0.86 versus the published 0.719). The
> additional seaward-slope dissipation needed to match 0.719 is a Madsen and White (1976) empirical
> calibration detail that is not recoverable from the public Technical Reference, and it does not
> affect the transmitted wave height.

**Validation.** Reproduces the *primary* outputs of the ACES *User's Guide* Example 1 (a three-material,
three-layer trapezoidal breakwater: 6.56 ft, 20 s wave in 15.75 ft of water) to better than half a
percent: transmitted height 1.571 ft (published 1.570), total transmission coefficient 0.239,
through-transmission 0.077, and overtopping transmission 0.227. The reflection coefficient is the one
output that does not match (see the caveat above). The hand-coded complex Bessel functions are checked
against their known real-argument values.

**References.** Madsen and White (1976); Seelig (1979, 1980); Ahrens and McCartney (1975).

## 5-5 — Wave Setup

**What it does.** Predicts how the mean water level changes across the surf zone as waves break.
Just seaward of breaking the mean level dips slightly below the still-water line (the set-down);
inside the surf zone it rises steadily (the set-up), lifting the mean water line up and onto the
beach. The tool reports the breaker height and depth, the set-down at breaking, the set-up
gradient, the set-up at the still-water shoreline and its maximum at the displaced waterline, the
surf-zone width, and the horizontal shift of the shoreline.

**The physics, briefly.** Breaking waves transfer momentum to the water column. The cross-shore
gradient of that wave momentum flux (the radiation stress of Longuet-Higgins and Stewart) has to
be balanced by a slope in the mean water surface. Seaward of breaking the balance produces a small
set-down, whose value at the breaker line is `η̄_b = -H_b² k / (8 sinh 2kd)`. Inside the surf zone,
where the broken wave height stays a fixed fraction `γ` of the local water depth, the balance gives
a constant set-up slope, a fixed fraction `β / (1 + β)` of the beach slope with `β = 3γ²/8`.
Marching that slope shoreward from the breaker set-down gives the set-up everywhere, the maximum at
the waterline, and the distance the shoreline is pushed up the beach. The breaker height comes from
the Singamsetti and Wind (1980) relation and the breaker depth and index `γ` from the Weggel (1972)
breaker index.

**Assumptions and range.** A plane beach of constant slope, a single representative wave (not a
spectrum), and a saturated surf zone where height tracks depth. The breaker-height relation is
calibrated for slopes of roughly 1:50 to 1:5. Set-down applies seaward of breaking; the constant
set-up slope applies between the breaker line and the waterline.

> **Status and Caveats:** Current for planning-level estimates of the wave-driven water-level rise
> that adds to storm surge and tide. It does not include wave-by-wave runup (see 5-1 and 5-2),
> infragravity contributions, or alongshore variability. The Weggel breaker index in the source is
> written in US units with gravity folded into the coefficient; CHESS-QC computes in SI using the
> gravity-explicit equivalent, which reproduces the US-unit result exactly.

**Validation.** No ACES *User's Guide* worked example exists for this tool, so the reference is the
Longuet-Higgins and Stewart radiation-stress theory itself. The set-down reproduces the closed form
above and its shallow-water limit `η̄_b → -γ² d_b / 16`; the set-up gradient reproduces the
closed-form surf-zone slope; and the cross-surf-zone table integrates consistently from the breaker
set-down to the maximum set-up. For a 6 ft, 10 s swell on a 1:50 beach the tool gives a breaker
index `γ ≈ 0.89`, a set-down of about `-0.5 ft`, and a maximum set-up of about `2.5 ft` (close to
a quarter of the breaker depth), which is the expected magnitude.

**References.** Longuet-Higgins and Stewart (1962, 1963); Weggel (1972); Singamsetti and Wind
(1980); SPM (1984).

---

# Area 6 — Littoral Processes

This area covers how sand moves along and across the coast: the rate sand is carried along the
shore by waves, how a beach and dune erode during a storm, how to combine and characterize
grain-size samples, and how much extra borrow sand a nourishment project needs. Three of the four
chapters are built and validated (6-1 longshore transport, 6-3 composite grain size, and 6-4 beach
nourishment); the time-dependent beach and dune erosion model (6-2) is Documented but not yet
implemented.

## 6-1 — Longshore Sediment Transport

**What it does.** Estimates the volume rate of sand carried along the shore (the longshore transport
rate) from wave conditions, using the standard CERC energy-flux method, from either breaking-wave
or deepwater inputs.

**The physics, briefly.** Waves approaching the shore at an angle drive a longshore current that
carries sand. The transport rate is proportional to the longshore component of the wave energy flux
at breaking, which scales with the breaker height to the five-halves power and the sine of twice the
breaker angle. An empirical coefficient links that energy flux to the actual sand volume. When only
deepwater conditions are known, refraction and shoaling relations (and Snell's law for the angle)
carry them to the breaker line first.

> **Status and Caveats:** Current. The standard CERC coefficient (about 0.39 for field conditions)
> and the breaker index are long-established but carry real scatter; the method gives a
> planning-level rate, not a precise one. The water density, sediment density, and porosity are
> exposed as inputs, defaulting to seawater, quartz sand, and 40 percent porosity.

**Validation.** With physically standard quartz sand and seawater, the formula reproduces the
well-known CERC literature factor, about 1290 cubic metres per year of transport per watt per metre
of longshore energy flux. The ACES *User's Guide* worked examples (275,234 yd³/yr for the deepwater
case at 1.75 ft and 15°, and 2,662,872 yd³/yr for the breaking case at 3.75 ft and 12°) come out
about 25 percent higher; reproducing them exactly requires an effective sediment density near 2320
kg/m³, which is below quartz, disagrees with the standard CERC factor, and is not stated in the
Technical Reference. CHESS-QC defaults to the physically correct value and documents the example
offset, reproducing the published numbers when the sediment density is set to the lower effective
value.

**References.** SPM (1984); Galvin (1979); Gravens (1988).

## 6-2 — Time-Dependent Beach and Dune Erosion

**What it does.** Simulates how a beach and dune profile erode and recover over the course of a
storm, driven by a time history of water level (surge plus tide) and waves, reporting shoreline
change at several elevation contours. This is the XSHORE cross-shore model.

**The physics, briefly.** The surf-zone profile tends toward an equilibrium shape (depth grows with
the two-thirds power of offshore distance) set by the grain size. A storm surge raises the water
level over that equilibrium beach, putting the profile out of equilibrium; the beach and dune then
recede exponentially toward a new equilibrium. CHESS-QC uses the Kriebel and Dean (1985)
equilibrium-response model in its closed-form limit: the equilibrium (maximum) recession follows from
an equilibrium sand balance, and the time-dependent response is an exponential approach to it with a
characteristic time scale (Kriebel and Dean 1993). This is physically grounded, magnitude-correct,
and has no free parameters. Longshore transport is neglected.

> **Status and Caveats:** Current (Built), as the Kriebel-Dean equilibrium-response model rather than
> the legacy XSHORE finite-difference scheme. It is surge-driven: it represents the dune and beach
> erosion caused by a storm raising the water level over an equilibrium profile (the design use case),
> and for that problem it is magnitude-correct with zero calibration. It does not reproduce the legacy
> ACES XSHORE worked example, which is a no-surge generic-profile readjustment (12 ft) specific to that
> 1970s/80s scheme; with no surge there is no equilibrium shift, so this model returns no recession for
> that case (a deliberate, documented difference, not a deficiency). A faithful XSHORE clone would need
> the ACES source or the Kriebel 1984b EBEACH theory manual, neither available. One-dimensional
> cross-shore only; no longshore effects.

**Validation.** No ACES numeric oracle exists for the surge problem (the only ACES worked example is
the no-surge XSHORE case). Validated analytically against the Kriebel and Dean (1985) paper: for its
Figure-5 case (medium sand so A = 0.118, a 2 m surge, a 4.6 m breaker, a 3 m berm, a 1-on-10 beach)
the equilibrium recession is about 79 m, matching both the independent Bruun-rule sand balance and the
paper's stated tens of metres, and the response time scale is about 19 hours, within the paper's stated
10-to-100-hour storm range. Recession is linear in surge level (the paper's Figure 5), the response is
exponential in time, and with no surge the recession is zero, all as the theory requires.

**References.** Kriebel and Dean (1985, 1993); Dean (1977); Bruun (1954); Moore (1982).

## 6-3 — Composite Grain-Size Distributions

**What it does.** Combines several sediment sieve samples into one composite distribution and reports
standard grain-size statistics (mean, sorting, skewness, kurtosis) by both the Folk graphic method
and the method of moments. This is the application labeled 6-5 in the ACES manual.

**The physics, briefly.** Grain sizes are expressed on the phi scale (the negative base-two logarithm
of the diameter in millimeters), which turns the wide range of sand sizes into convenient numbers.
The graphic statistics read selected percentiles off the cumulative curve; the method of moments
sums over all size classes. Both describe the central size and how well or poorly sorted the sample
is.

> **Status and Caveats:** Current. A descriptive utility rather than a predictive model; results are
> only as good as the sieve data supplied. The current build computes the statistics for one supplied
> distribution (a composite is the same method applied to the averaged samples).

**Validation.** The method of moments is validated analytically on the Panama City CoreSample1 sieve
data shipped with the repo (mean = 2.652 phi, an exact weighted average independently checked). The
ACES *User's Guide* worked example is the composite of two core samples; only CoreSample1 is in the
repo, so that exact composite output cannot be reproduced, but the Folk graphic measures for the
single sample fall in the published composite range (median near 2.6 phi, sorting near 0.5).

**References.** Folk (1974); Krumbein (1934, 1938); SPM (1984).

## 6-4 — Beach Nourishment Overfill Ratio and Volume

**What it does.** Tells a nourishment designer how much borrow sand to place to get a given volume of
usable beach (the overfill ratio) and how much faster or slower the borrow sand will erode compared
to the native sand (the renourishment factor).

**The physics, briefly.** The borrow and native sands are compared through two phi-scale numbers: the
difference in mean grain size and the ratio of sorting (spread). If the borrow sand is finer or more
poorly sorted than the native, more of it washes away, so more must be placed. The James (1975)
formulas turn those two comparison numbers into the overfill ratio (through a normal-distribution
integral, by category of the comparison) and into the renourishment factor.

> **Status and Caveats:** Current. Identical borrow and native sands give a ratio and factor of one
> (no penalty). The method is a screening tool; it does not model the actual profile or wave climate.

**Validation.** Reproduces ACES *User's Guide* Example 6-4 (native phi mean 1.80 and sorting 0.45,
borrow mean 2.25 and sorting 0.76, initial volume 800,000 yd³): overfill ratio 2.00, renourishment
factor 1.08, and design volume about 1.60 million yd³, all within rounding of the published values.

**References.** James (1975); Hobson (1977); SPM (1984).

---

# Area 7/8 — Inlet and Harbor Processes

This area covers tidal inlet hydraulics, a set of shared breaker utilities, and harbor-basin
resonance. A naming note: the ACES technical reference numbers two chapters here, the inlet model
(7-1) and a miscellaneous-routines chapter it labels 8-1, which is not harbor design. The harbor
design applications proper (rectangular basins, vessel waves, moored vessel surge) have no
technical-reference chapter and are sourced from standard theory. Of the sections below, the
inlet model and the breaker utilities are Built and validated; the harbor-basin tool (the manual's
application 8-1, distinct from the technical reference's 8-1 chapter) is also Built and validated.

## 7-1 — Inlet Hydraulics

**What it does.** Models the tidal flow through one or more inlets connecting the sea to a bay,
predicting the discharge through each inlet and the rise and fall of the bay water level over a tidal
cycle, for several sea-inlet-bay configurations.

**The physics, briefly.** Unlike the simpler lumped (Keulegan) inlet model, this is a spatially
integrated one-dimensional model. The momentum equation along the inlet axis, including bottom
friction expressed through Manning's n, is integrated over the inlet length to give an ordinary
differential equation for the discharge. That is coupled with a bay-continuity equation (the bay
level rises as net discharge fills it), and the pair is marched through time with a fourth-order
Runge-Kutta scheme. The inlet is described by a set of surveyed cross-sections: each one's flow
area and top width are integrated by the trapezoidal rule from its bed-elevation profile, the
narrowest section sets the controlling throat, and the sections in series set both a geometry
integral (the inertia coefficient) and the distributed Manning friction. The sea boundary is a
harmonic-constituent tide, synthesized with the same Schureman astronomy as application 1-4.

> **Status and Caveats:** Built (standard). The cross-section area integrator reproduces the published
> flow-net areas exactly (cross-section 1, 100,360 ft2; cross-section 5, 60,112 ft2; throat 40,456
> ft2). The 30-hour time march of User's Guide Example 1 (pure M2 tide) reproduces the bay
> elevation hydrograph to better than 0.02 ft, the first-ebb peak discharge of -207,260 cfs to 0.2
> percent, and the controlling-section velocity (5.05 ft/s) to about 1 percent. The mid-record
> flood and ebb exchange volumes run about 6 percent low: the friction here is computed from each
> section's mean depth rather than the full per-channel flow net, because the User's Guide publishes
> the channel division for only two of the five cross-sections. The headline discharge, velocity,
> and bay-range metrics meet the project bar; the per-channel velocity field (a secondary output)
> would require the full flow-net subdivision.

**References.** Seelig (1977); Seelig, Harris, and Herchenroder (1977); Harris and Bodine (1977);
Keulegan (1967).

## 8-1 — Miscellaneous Routines (Breaker and Steepness Utilities)

**What it does.** Provides the shared breaking-wave and wave-steepness relations that several other
applications call: the maximum stable wave steepness, the breaker height in shallow water, the
breaking depth, and the maximum breaker height in the vicinity of a structure.

**The physics, briefly.** A wave can only get so steep before it breaks (the Miche limit). In shallow
water the breaker height is capped at a fraction of the depth, refined for sloping beds by the
Singamsetti and Wind relation. The Weggel (1972) breaker index gives the breaking depth and the
near-structure breaker height through two slope-dependent coefficients, a and b. Those same Weggel
coefficients are reused by the refraction and wave-setup applications.

> **Status and Caveats:** Current (Built, as application M-1 to distinguish it from the harbor-design
> rectangular-basins app 8-1). Despite the source grouping, this is not harbor design: it is a utilities
> chapter. The harbor-design applications (rectangular basins, vessel waves, moored vessel surge) have
> no technical-reference chapter and are sourced externally. The Weggel a and b coefficients here are
> the canonical source for the breaker index used in 3-1 and in wave setup. One relation is not used
> as transcribed: the elaborate eq-5 closed form for the maximum breaker height in the structure
> vicinity returns non-physical values as written in the public reference (a 79 ft breaker for a 15 ft
> structure depth, growing with depth), and there is no standalone worked example to debug it; the
> depth-consistent inversion of the Weggel breaker-depth relation is used instead, which gives the
> physically expected result.

**Validation.** These are shared utility relations with no standalone ACES worked example, so they are
checked analytically and through the applications that use them (3-1 and 5-5). McCowan flat-slope
breaking is exact (0.78 times the depth); the Miche deepwater steepness limit is 0.142; the Weggel
index has the correct slope limits (a goes to zero and b to 0.78 as the slope vanishes); and the
Weggel breaker depth and height form a mutually consistent pair (inverting the depth relation recovers
the height). For a 1-on-20 structure in 15 ft of water under an 8 s wave the structure breaker height
is 14.2 ft.

**References.** Miche (1944); McCowan (1894); Singamsetti and Wind (1980); Weggel (1972).

## 8-1 (Harbor Design) — Properties of Rectangular Basins

This is the manual's harbor-design application 8-1, which is a different tool from the
miscellaneous-routines chapter 8-1 above; they share a number only by an accident of the ACES
numbering.

**What it does.** Computes the natural oscillation (seiche) periods of a rectangular harbor or
basin, the resonances that a long wave, surge, or tsunami can excite, together with the
water-particle velocity and excursion at a node. It handles a basin open at one end, a fully
enclosed basin in one or two dimensions, and the Helmholtz pumping mode of a basin connected to
the sea by a channel.

**The physics, briefly.** A basin of a given length has a set of natural periods at which water
sloshes back and forth, in the same way a pipe or a guitar string has natural tones. For shallow
water the wave travels at speed `c = sqrt(g d)`, so the periods follow directly from the basin
length and the mode number: a closed basin resonates when a whole number of half-wavelengths fits
its length (Merian's formula), an open basin when an odd number of quarter-wavelengths fits, and a
rectangular basin has a two-dimensional set of modes combining the length and width. A basin joined
to the sea by a short channel also has a Helmholtz mode, the slow in-and-out pumping of the whole
basin through its mouth, whose period grows with the basin area and the channel length and shrinks
with the channel cross-section. At a node (where the surface pivots and the flow is fastest) the
tool reports the peak horizontal velocity, the particle excursion, and the mean speed from the
standing-wave (clapotis) relations.

**Assumptions and range.** A rectangular basin of uniform depth, shallow-water resonance
(wavelength long compared with depth), and an inviscid response with no damping, so the periods are
exact but the amplitudes are idealized. The Helmholtz mode treats the basin as a lumped reservoir
and the channel as a lumped inertia.

> **Status and Caveats:** Current. Real harbors are rarely rectangular or flat-bottomed, so treat
> these periods as first estimates for identifying which incoming wave or surge periods are
> dangerous; a boundary-element or finite-element harbor-resonance model is preferred for a final
> design. Damping (and therefore the actual amplification at resonance) is not modeled.

**Validation.** No ACES *User's Guide* worked example exists for this tool (the ACES manual
confirms the inputs and outputs but prints no numbers), so the reference is the classical seiche
theory. The tool reproduces Merian's closed forms for the open and closed fundamentals
(`T_0 = 4 l_B / sqrt(g d)`, `T_1 = 2 l_B / sqrt(g d)`), the two-dimensional mode reduces exactly to
the one-dimensional mode when the transverse index is zero, a square basin gives the expected
symmetry `T(n,m) = T(m,n)`, and the Helmholtz period matches its closed form. The node kinematics
satisfy the standing-wave identities (excursion equals peak velocity over angular frequency).

**References.** Merian's formula and the seiche literature in SPM (1984); Wilson (1972);
Sorensen (1993); Ippen (1966).

---

# Area 9 — Storm Surge

This area holds tools for the rise of water level at the coast during a storm. It is a new
functional area, added beyond the original eight ACES areas.

## 9-1 — Bathystrophic Storm Surge

**What it does.** It estimates the open-coast hurricane surge along a single cross-shelf traverse
as a storm passes, returning the peak still-water rise at the shore broken into its physical
parts (wind setup, Coriolis or "bathystrophic" setup, atmospheric-pressure setup), plus profiles
of setup and total depth along the traverse.

**The physics, briefly.** The model integrates, from the shelf edge to the shore as the storm
moves by, the simplified depth-integrated equations of motion (Bodine 1971). Onshore wind stress
piles water against the coast (wind setup); the alongshore wind drives an alongshore current
that, turned by the earth's rotation, adds a further rise (the bathystrophic effect); and the low
central pressure lifts the surface (inverse barometer). The wind that drives all this comes from
a parametric hurricane model with two interchangeable choices: Holland (1980) by default, or
Myers (1954), which is exactly the Holland model with shape factor `B = 1`. The wind speed at any
radius follows from gradient-wind balance on the chosen pressure profile, reduced to the surface
and given a forward-speed asymmetry and an inflow angle.

**Key relationships (our notation).**

- Pressure profile: `p(r) = P_c + ΔP exp(-(R/r)^B)`, with `B = 1` recovering Myers.
- Gradient wind: `V(r) = sqrt( (B ΔP / ρ_a)(R/r)^B exp(-(R/r)^B) + (r f / 2)² ) - r f / 2`.
- Setup balance (cross-shore) and transport (alongshore) integrated as Bodine's finite-difference
  analogs, with quadratic bottom friction and a flux limiter; the total shore rise sums wind,
  bathystrophic, pressure, initial-rise, astronomical-tide, and wave-setup parts.

The Holland shape factor `B` (default 1.5) may be entered directly, or back-computed from a
specified maximum wind through `B = ρ_a e V_max² / ΔP`.

> **Status and Caveats:** Screening only. The bathystrophic method is a quasi-one-dimensional
> approximation along a single traverse; it omits the alongshore and two-dimensional dynamics,
> nonlinear advection, and the inlet and bay effects that a full model captures. For structure
> design or risk assessment use **ADCIRC**. Bodine (1971) notes the method can
> be in error by up to a factor of two. CHESS-QC additionally uses a parametric (analytic)
> wind field in place of Bodine's hand-digitized isovel chart, so the surge matches the original
> worked example only in the ballpark.

**Validation.** For the Bodine (1971) Chesapeake Bay Entrance example (central pressure
`27.57 inHg`, peripheral `29.92 inHg`, radius `35 nm`, forward speed `22 kt`, latitude `37°`)
with the Myers wind model, the maximum 30-foot wind is `102.9 mph` (Bodine `V_x = 102`) and the
peak surge is about `16.8 ft` against Bodine's `13.4 ft`. The analytic sub-models are checked
exactly: gradient wind equals `sqrt(B ΔP / (ρ_a e))` at the radius of maximum winds, `B = 1`
collapses Holland onto Myers, the pressure setup is about `1.1 ft` per inch of mercury, and a
specified maximum wind back-computes `B` correctly.

**References.** Bodine (1971), CERC TM-35; Freeman, Baer, and Jung (1957); Holland (1980);
Myers (1954); Van Dorn (1953) for wind stress.

---

*Built areas complete: 1 (partial: 1-1, 1-3), 2 (Wave Theory), 3 (Wave Transformation), 4
(Structural Design), 9 (Storm Surge). The remaining ACES applications and areas will be added to
this reference as they are implemented.*
