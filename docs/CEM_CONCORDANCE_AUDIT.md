# CEM Concordance Audit

Audit source: the Coastal Engineering Manual, Parts I-VI (CHESS-Technical-Guidance, `modules/legacy_guidance/manuals/cem`).  This record distinguishes an implementation of a CEM equation from a retained ACES/SPM or external method.  “Qualified” is not a failure: it means the calculator is deliberately tied to the cited legacy/source method and must not be represented as the newer CEM recommendation.

## Three independent checks

1. **Source check.** The governing equation, method, units, and validity range were compared with the local CEM chapter or, where CEM is not the governing source, the module’s cited source.
2. **Implementation check.** Inputs, outputs, numerical limits, units, method selection, and traceability text were inspected.  This found and corrected two CEM defects in 2-5, a displayed-unit defect in 1-5, an inactive output-control input in 7-1, and the tail-fitting/bootstrap procedure in 10-4.
3. **Numerical check.** Every module’s self-check ran; all applicable ACES User’s Guide examples ran through the independent manual-oracle suite; all 40 modules then built, computed, reset, and converted units through the GUI and web bridge.

## Result by module

| Module | CEM relationship / governing reference | Result after the three checks |
| --- | --- | --- |
| 1-1 Windspeed & wave growth | CEM II-2 supersedes the retained SPM/Resio growth relations. | Qualified legacy method; neutral oracle and physical stability-sign checks pass. |
| 1-2 Beta-Rayleigh | CEM II-1 covers irregular-wave statistics; this is the cited Hughes/Borgman depth-limited distribution. | Qualified standard method; distribution and manual checks pass. |
| 1-3 Extreme Hs | CEM II-8 provides modern extreme-wave context; calculator uses Goda/EM fitting. | Qualified legacy method; all five fits and confidence interval match the ACES example. |
| 1-4 Tide record | CEM water-level guidance is compatible; governing source is Schureman harmonic astronomy. | Verified source method; 120-hour ACES record matches. |
| 1-5 Near-surface wind | CEM II-2 is contextual; ACES TR 1-1 and Silva (2005) eqs. 7.16–7.19 govern the revised PBL resistance law. | Source-completed standard method: Sherlock constants and revision, distinct stability exponents, source-sign μ, ACES KEYPS/PBL iteration, neutral reduction, and exact equation-branch checks pass. |
| 1-6 Holland hurricane | Meteorological parametric model; no CEM replacement equation. | Verified external source method; pressure, radius, and wind-limit identities pass. |
| 2-1 Linear waves | CEM II-1-7 through II-1-60 (Airy/dispersion/kinematics). | Direct CEM-compatible implementation; dispersion, deep/shallow, energy, and ACES checks pass. |
| 2-2 Cnoidal waves | CEM II-1-77 through II-1-82. | Direct CEM-compatible first-order implementation; elliptic identities and ACES example pass. |
| 2-3 Fenton Fourier waves | CEM II-1 stream-function/Fourier discussion. | Direct theory implementation; linear limit and ACES example pass. |
| 2-4 Wave parameters | CEM II-1 linear-wave relations and pressure response. | Direct CEM-compatible implementation; pressure round-trip and 2-1 cross-check pass. |
| 2-5 Solitary waves | CEM II-1-83 through II-1-98. | Corrected direct CEM implementation: CEM Figure II-1-17 M,N values, II-1-92/93 kinematics, II-1-98 sloping breaker, and II-1-97 limit now apply. |
| 3-1 Snell transformation | CEM II-3 nearshore transformation. | Verified compatible method; ACES transformation and breaker example pass. |
| 3-2 Goda transformation | CEM II-3 irregular-wave transformation and complete ACES TR 3-2 equations. | Source-verified standard method; deterministic Goda clipped-distribution integration and explicit finite-wave Hmax definition. |
| 3-3 Wedge diffraction | CEM II-3 context; governing Chen/Penny-Price wedge solution. | Verified legacy source method; manual magnitude/height check passes; phase convention is qualified. |
| 3-4 Wedge grid | Same as 3-3. | Verified legacy source method; ACES grid check passes. |
| 4-1 Hudson breakwater | CEM VI-5 discusses armor design; module retains Hudson/SPM sizing. | Qualified legacy design method; ACES example passes. |
| 4-2 Toe protection | CEM VI-5 context; Tanimoto/Yagyu/Goda and EM rules govern. | Verified source method; ACES example passes. |
| 4-3 Vertical-wall forces | CEM VI-5 context; Sainflou/Miche/Rundgren govern. | Verified source method; ACES example passes. |
| 4-4 Revetment design | CEM VI-5 context; ACES Ahrens/van der Meer method plus selected modern option. | Verified stated methods; ACES runup/design example passes. |
| 5-1 Beach runup | CEM II-4 context; Mase empirical relations govern. | Verified source method; all reported runup statistics match the ACES example. |
| 5-2 Impermeable runup/overtopping | CEM VI-5 context; Ahrens/Weggel/SPM govern. | Verified source method; seven ACES examples pass. |
| 5-3 Impermeable transmission | CEM VI-5-2 includes the selected d’Angremond modern option. | Verified; ACES examples and transmission bounds pass. |
| 5-4 Permeable transmission | CEM VI-5 context; Madsen-White layered method governs. | Source-completed standard method: ACES eq.-64 head-difference closure and Madsen-White Table-2 steep-slope reflection correction are iterated; coefficient identities and bounds pass. |
| 5-5 Wave setup | CEM II-4 radiation-stress setup. | Direct compatible method; sign, limiting, and integration checks pass. |
| 6-1 Longshore transport | CEM III-2 CERC energy-flux transport. | Verified formula; quartz-density default is physical, and ACES density discrepancy is disclosed. |
| 6-2 Dune erosion | CEM III-3 context; Kriebel-Dean governs. | Verified source method; analytical and published-figure checks pass. |
| 6-3 CEDRS transport | CEM III-2 CERC transport over directional climate. | Standard implementation; ACES net transport is within 1 percent and symmetry checks pass. |
| 6-4 Beach nourishment | CEM V-4 context; James TM-60 governs. | Verified source method; ACES example and category identity pass. |
| 6-5 Composite grain size | CEM III-1 sediment-size definitions. | Direct compatible implementation; moment/Folk checks pass. |
| 7-1 Inlet hydraulics | CEM II-6 inlet hydraulics context; complete Seelig/Harris-Bodine source equations govern. | Source-verified standard method within its one-inlet scope. Fixed inactive tabular-output interval without altering integration/extrema. |
| 7-2 Wave-current interaction | CEM II-6 context; Jonsson relations govern. | Verified source method; dispersion/action and current-direction limits pass. |
| 8-1 Rectangular basins | CEM II-7 harbor hydrodynamics context. | Verified closed-form source method; modal/Helmholtz reductions pass. |
| 8-2 Vessel waves | CEM II-7 harbor context; Schijf/ship-wave relations govern. | Verified source method; continuity, energy, and Mach-angle checks pass. |
| 8-3 Moored-vessel surge | CEM II-7 harbor context; mooring dynamics/PIANC govern. | Verified source method; stiffness, mass, symmetry, and loading checks pass. |
| 9-1 Bathystrophic surge | CEM II-5 water-level context; Bodine/Myers/Holland traverse method governs. | Source-completed standard method: TM-35 appendix bathymetry, wind/radius/direction curves, and finite-difference constants are encoded; reference result is reproduced within the published curve precision. Screening-only scope remains. |
| 10-1 Detrending | Coastal-hazards/NOAA method, outside CEM’s governing equations. | Verified external method; slope, pivot, override, gaps, and decimation checks pass. |
| 10-2 Non-tidal residual | NOAA/PyStorm method, outside CEM’s governing equations. | Verified external method; subtraction, alignment, and mean checks pass. |
| 10-3 Peaks over threshold | Coles/USACE/PyStorm method, outside CEM’s governing equations. | Verified external method; threshold, declustering, and ranking checks pass. |
| 10-4 Probabilistic simulation | Coles/Nadal-Caraballo/PyStorm method, outside CEM’s governing equations. | Standard empirical-bulk/GPD-tail implementation: fixed-location maximum-likelihood GPD fits and a nonparametric resampling bootstrap; interval ordering and return-period checks pass. |
| M-1 Breaker routines | CEM II-1, II-3, and II-4 breaker/steepness relations. | Direct/compatible utilities; McCowan, Miche, and Weggel limits are checked. |

## Corrections made during this audit

- **2-5 Solitary Wave Theory:** replaced a lower-order velocity approximation with CEM II-1-92/93, replaced an inconsistent algebraic M estimate with interpolation of CEM Figure II-1-17, enforced the CEM II-1-97 breaking limit, and applied CEM II-1-98 for sloping-bed breaking.
- **1-5 Near-surface wind:** corrected the displayed/input SI wind-speed unit from km/h to m/s. The literature reconciliation also corrected eq. 15's sign, eq. 16's distinct 0.015/0.03 exponents, the sign of `μ=kU*/(fL')`, and the revised neutral limit `A(0)=A_0`; it uses the Sherlock (1996) constants and revision published by Silva (2005).
- **5-4 / 9-1 source completion:** implemented the exact Madsen-White/ACES head-ratio closure, material diameter convention, and coefficient synthesis; and transcribed Bodine TM-35's appendix bathymetry, forcing curves, constants, tide, and original difference scheme, including exact first-step and midpoint-depth logic from the FORTRAN listing.
- **3-2 Goda transformation:** resolved the provisional-source label after verifying every governing relation against the complete ACES TR 3-2 equation set; Hmax is now explicitly documented as a finite-wave statistic.
- **7-1 Inlet Hydraulics:** made the declared tabular-output interval control emitted hydrograph rows while retaining the requested RK4 time step and full-resolution extrema.
- **10-4 Probabilistic Simulation:** replaced method-of-moments fitting and perturbed-value pseudo-bootstrap with fixed-threshold maximum-likelihood GPD fitting and a resampling-with-replacement bootstrap.

## Validation record

The final run completed successfully: all 40 module self-checks; all manual-oracle checks; desktop GUI, table-widget, web-bridge, launcher, and browser smoke checks; Python compilation; and whitespace/diff validation.
