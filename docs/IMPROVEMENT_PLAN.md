# CHESS-QC — Application Improvement Plan
*Coastal Hazards, Engineering, and Structures System (CHESS) — Quick Compute (QC)*

Improvement roadmap for all 35 built applications (the area-10 "Coastal Hazards" placeholders
are excluded). Prepared 2026-06.

Each application is assessed on three axes: **(1) validation shortcomings**, **(2) missing
data**, **(3) a superior / replacement method** that could be adopted. Feasibility tags:
`[drop-in]` = swap a formula or coefficient set in place (same inputs); `[new-calculator]` =
add a better standalone app; `[needs-data]` = blocked on a dataset or source document.

---

## Priorities

**Tier 1 — provisional (class III): fix or supersede first**
- **1-5 Near-Surface Wind** — no numeric oracle at all; stratified resistance law missing its appendix constants.
- **3-2 Goda transformation** — high quantiles off by up to ~6%; relies on an undocumented finite-N largest-wave inflation and spreading scheme.
- **5-4 Permeable transmission** — reflection K_R over-predicted (~0.86 vs 0.719); Madsen-White calibration not recoverable; tuned reference diameter.
- **7-1 Inlet hydraulics** — per-channel velocity field (Table 7-1-2) not reproduced; flow-net subdivision + full bathymetry missing.
- **9-1 Bathystrophic surge** — ballpark only (factor-of-two per Bodine 1971); already flagged superseded by ADCIRC.

**Tier 2 — standard (class II) / known residuals at physical inputs**
- **1-2 Beta-Rayleigh** — H_1/10 6.30 vs manual 6.55 ft, unresolved (suspected source-doc artifact).
- **6-1 Longshore** & **6-3 CEDRS** — only match the ACES examples with a non-physical effective sediment density (~2320 vs 2650); self-chosen binning conventions in 6-3.
- **M-1 breaker routines** — TR eq-5 non-physical as printed; replaced by a self-made substitute; analytic checks only.

**Tier 3 — class I but validated analytically only (no ACES numeric oracle)**
- **1-6, 2-4, 2-5, 5-5, 7-2, 8-1, 8-2, 8-3** (and partly **2-5**, whose SPM source is also absent). Action: give each a published cross-check oracle (literature worked example) so the "exact" label is independently backed, or down-classify.

---

## Cross-cutting themes
1. **Missing ACES appendix constants.** The geostrophic/over-land resistance-law similarity constants (A_0/B_0/B_1) are absent from the public TR, leaving the stratified wind paths in **1-1** and **1-5** provisional. Fix by adopting a self-contained modern bulk scheme (COARE 3.5) rather than chasing the missing constants.
2. **Non-physical effective sediment density.** **6-1** and **6-3** only reproduce the ACES longshore numbers with rho_s ~2320 (not quartz 2650), which the TR never documents. Decide: keep physical density (and accept the ACES examples are ~25% high), or expose the effective density as a documented input.
3. **No numeric oracle (analytic-only).** Ten apps (1-5, 1-6, 2-4, 2-5, 5-5, 7-2, 8-1, 8-2, 8-3, M-1) have no ACES worked example. Several are sound closed forms; the risk is silent regressions. Action: add literature-based cross-check oracles to the test suite.
4. **Missing datasets.** **6-5** (CoreSample2 companion sample), **7-1** (full surveyed bathymetry / flow-net), **9-1** (Bodine's hand-digitized wind+bathymetry), **6-3** (broader CEDRS/WIS station library).
5. **Clear modern supersessions.** Hudson → Van der Meer (4-1); Weggel/Seelig → EurOtop & d'Angremond (5-2, 5-3, 5-4); clapotis → Goda wave pressures (4-3); Ahrens-Heimbaugh → EurOtop runup (4-4); Mase → Stockdon (5-1); CERC → Kamphuis (6-1, 6-3, M-1); bathystrophic → ADCIRC (9-1).

## Quick wins (`[drop-in]`, highest value / lowest effort)
- **4-1** add Van der Meer stability — the kernel already exists in 4-4, so reuse it.
- **5-2 / 5-3 / 5-4** adopt EurOtop (2018) overtopping and d'Angremond/Van der Meer (2005) transmission — documented, well-validated, dimensionless.
- **5-1** add Stockdon et al. (2006) runup; **4-4** swap the runup sub-calc to EurOtop.
- **1-2** add Battjes-Groenendijk (2000) distribution (resolves the H_1/10 ambiguity).
- **6-1 / 6-3 / M-1** offer Kamphuis (1991) as a selectable alternate kernel.

---

## Per-application detail

### Area 1 — Wave Prediction

### 1-1 Windspeed Adjustment and Wave Growth  [class I]
- **Validation shortcomings:** Neutral path reproduces Examples 1-1 (open) and 3 (restricted) exactly. The air-sea stability correction is the weak point: ACES TR eqs 8-9 are internally sign-inconsistent and eq 9 is a corrupted Businger-Dyer/Paulson (1970) psi_m; the shipped canonical correction cannot reproduce ACES Example 3 and is opt-in/unvalidated. The full-PBL geostrophic path (eqs 14-19) needs untranscribed appendix constants and is unvalidated.
- **Missing data:** Appendix similarity constants A_0/B_0/B_1; no worked example decouples stability from observation-type bias.
- **Superior / replacement method:** CEM (EM 1110-2-1100, Part II-2) growth formulas + a COARE 3.5 bulk-flux algorithm (Fairall et al. 2003) for a correctly-signed, validated stability/drag treatment. `[drop-in]` growth coefficients; `[new-calculator]`/`[needs-data]` for the COARE module.

### 1-2 Beta-Rayleigh Distribution  [class II]
- **Validation shortcomings:** Four of five characteristic heights match Example 1-2 to ~1.5%, but **H_1/10 = 6.30 vs manual 6.55 ft** is unreconcilable (grid-independent), argued to be a source-document artifact. Relative-depth fits eq (16)/(19) were transcribed with depth inverted; module silently uses the corrected argument. Single-example validation.
- **Missing data:** None unrecoverable.
- **Superior / replacement method:** Battjes-Groenendijk (2000) composite Weibull depth-limited distribution — modern standard, gives H_1/10/H_2%/H_1% directly, resolves the ambiguity. `[drop-in]`.

### 1-3 Extremal Significant Wave Height  [class I]
- **Validation shortcomings:** None material; reproduces Example 1-3 (correlations, return values, CIs). Methodological limits: 5 fixed candidate distributions, least-squares (not MLE), Goda Table coefficients tied to discrete k.
- **Missing data:** None noted.
- **Superior / replacement method:** L-moments / MLE GEV + GPD peaks-over-threshold (Hosking & Wallis 1997; Coles 2001) for continuous shape estimation and bootstrap CIs. `[new-calculator]` (parallel POT/GEV app).

### 1-4 Constituent Tide Record  [class I]
- **Validation shortcomings:** None material; exact Schureman node factors + equilibrium arguments reproduce Example 1-4 to 0.04 ft over 120 h. Only subtlety is the pinned ACES longitude/time convention.
- **Missing data:** None noted.
- **Superior / replacement method:** UTide / NS_TIDES (Codiga 2011; Foreman et al. 2009) — adds nodal corrections, a far larger constituent set, and harmonic *analysis* of observed records (not just synthesis). `[new-calculator]`.

### 1-5 Near-Surface Wind Speeds  [class III]
- **Validation shortcomings:** Weakest in the set — **no ACES worked example; analytic-only checks** (log-law recovery, C_D band, U_* fraction, stability sign). The stratified resistance law cannot be reproduced (missing constants); stratification enters only via the surface-layer profile, so the modeled response is unverified in magnitude.
- **Missing data:** Resistance-law similarity constants A_0/B_0/B_1; no ACES numeric example.
- **Superior / replacement method:** COARE 3.5 (Fairall et al. 2003) or a Grachev-Fairall stable-PBL scheme — self-contained, data-validated, no missing constants. `[drop-in]` for the flux core.

### 1-6 Holland Hurricane Wind Model  [class I]
- **Validation shortcomings:** Mild; Holland (1980) fully specified, profiles + cyclostrophic V_max validated **analytically** (no ACES numeric example). Fixed air density; axisymmetric, stationary, no translation asymmetry or surface-wind reduction.
- **Missing data:** None noted.
- **Superior / replacement method:** Holland, Belanger & Wright (2010) revised model (radius-varying B) + translation asymmetry & gradient-to-surface reduction (Kepert 2001 / NWS factor) for a realistic forcing field. `[drop-in]` B(r)+factor; `[new-calculator]` for a full asymmetric BL field.

### Area 2 — Wave Theory

### 2-1 Linear Wave Theory  [class I]
- **Validation shortcomings:** None; full Airy theory with Hunt (1979) dispersion reproduces Example 2-1 to the digit. Intrinsic small-amplitude limit reported via Ursell.
- **Missing data:** None noted.
- **Superior / replacement method:** None for the linear regime; optionally auto-recommend the right theory via the Le Méhauté (1976)/CEM validity diagram. `[drop-in]` (advisory).

### 2-2 Cnoidal Wave Theory  [class I]
- **Validation shortcomings:** First-order reproduces Example 2-2; module **corrects a TR transcription artifact** in the dw/dt accelerations (so shipped accelerations deliberately depart from the printed eqs). Second-order branch + corrected accelerations validated analytically, not against an ACES number.
- **Missing data:** None noted.
- **Superior / replacement method:** Fenton (1979/1990) high-order cnoidal, or route to 2-3 for uniformly higher near-breaking accuracy. `[drop-in]` / existing 2-3.

### 2-3 Fenton Fourier (stream-function)  [class I]
- **Validation shortcomings:** None; exact nonlinear steady-wave solver (Rienecker & Fenton 1981) reproduces the example to the digit; manual datum ambiguity resolved correctly. Newton convergence sensitive very near breaking (mitigated by height ramping).
- **Missing data:** None noted.
- **Superior / replacement method:** Already state of the art; optional Clamond & Dutykh (2018) fast/high-precision solver or adaptive-N for extreme steepness. `[drop-in]` (robustness only).

### 2-4 Wave Parameters  [class I]
- **Validation shortcomings:** **No dedicated ACES oracle** — inherits validity from the 2-1 example for shared outputs, plus a pressure-inversion round-trip and analytic limits. SPM App-C dimensionless quantities and the pressure-transducer inversion are validated by construction, not an independent table.
- **Missing data:** None noted.
- **Superior / replacement method:** Replace the linear K_p pressure-transducer inversion with an empirical/nonlinear reconstruction (Bishop & Donelan 1987). `[drop-in]`.

### 2-5 Solitary Wave Theory  [class I]
- **Validation shortcomings:** "Exact" is generous — **no ACES TR chapter, no ACES example, SPM source absent**; analytic-only (celerity, crest, 0.78 breaking, McCowan-Munk M-N). First-order McCowan/Munk values; M,N cross-checked against secondary forms, not the ACES source.
- **Missing data:** SPM solitary-wave chapter and ACES 2-5 example both absent.
- **Superior / replacement method:** Boussinesq/KdV exact soliton or ninth-order Fenton (1972) solitary-wave solution — gives a defensible oracle and accurate kinematics (max H/d ~0.83). `[drop-in]` coefficients; `[new-calculator]` for full kinematics.

### Area 3 — Wave Transformation

### 3-1 Snell's-law Refraction  [class I]
- **Validation shortcomings:** None material; reproduces Example 3-1 to the digit incl. the breaker. Secondary breaker output leans on Singamsetti & Wind (1980)/Weggel (1972).
- **Missing data:** None noted.
- **Superior / replacement method:** For real bathymetry, SWAN (or Boussinesq) handles refraction-diffraction coupling and spreading the straight-contour Snell assumption cannot. `[new-calculator]`.

### 3-2 Goda Random-Wave Transformation/Breaking  [class III]
- **Validation shortcomings:** Shoaling exact and Hs/Hmean/Hrms/surf-beat ~3%, but **high quantiles (H1/10, H1/50, Hmax) up to ~6%** off, using an ACES-specific finite-N largest-wave inflation (N~1200) inconsistent with the manual's own plotted Rayleigh. Effective refraction coefficient depends on an undocumented directional scheme (s_max~13).
- **Missing data:** The finite-N count and the directional-integration scheme / default spreading.
- **Superior / replacement method:** Battjes-Janssen (1978) dissipation model or updated Goda (2010, *Random Seas* 3rd ed.) breaking model — documented, consistent depth-limited statistics. `[drop-in]` core; `[new-calculator]` (SWAN) for full directional transformation.

### 3-3 Wedge Diffraction  [class I]
- **Validation shortcomings:** Minor; exact Chen (1987) eigenfunction solution matches the oracle (Kr=0.58). ~0.1 rad phase offset from a reference convention (does not affect magnitude).
- **Missing data:** None noted.
- **Superior / replacement method:** For multi-structure/partial-reflection fields, Boussinesq (MIKE 21 BW) or SWAN diffraction. `[new-calculator]`.

### 3-4 Wedge Diffraction/Reflection on a Grid  [class I]
- **Validation shortcomings:** None material; same exact solver on a grid reproduces Example 3 to <1%. Generic GUI surfaces only scalar summaries of a 2-D field.
- **Missing data:** None noted.
- **Superior / replacement method:** Phase-resolving Boussinesq / mild-slope for realistic fields. `[new-calculator]`. Also upgrade the GUI to a contour/heatmap renderer `[drop-in]` (web already has a heatmap; desktop heatmap added 2026-06).

### Area 4 — Structural Design

### 4-1 Breakwater Armor (Hudson)  [class I]
- **Validation shortcomings:** Computationally exact, but Hudson ignores period, duration, permeability and damage; accuracy depends entirely on user-selected K_D.
- **Missing data:** None in code (table coefficients are user inputs).
- **Superior / replacement method:** Van der Meer (1988) stability (plunging/surging, permeability P, damage S, storm count N) per the CIRIA Rock Manual (2007) — already named in this app's `superseded_by`, and the kernel exists in 4-4. `[drop-in]`.

### 4-2 Toe Protection  [class I]
- **Validation shortcomings:** None material; Tanimoto/Yagyu/Goda (1982) + EM 1110-2-1614 apron rules reproduce Example 4-2 #1. `cot_phi` collected but unused.
- **Missing data:** None noted.
- **Superior / replacement method:** Van der Meer/Gerding (1995) toe-berm stability (CIRIA Rock Manual 2007; CEM VI-5). `[drop-in]`.

### 4-3 Vertical Wall (Wave) Forces  [class I]
- **Validation shortcomings:** Effectively exact (Sainflou + Miche-Rundgren clapotis, ~0.1%). Scope-limited: nonbreaking/standing only; reflection coefficient is a user input; no impact pressures.
- **Missing data:** None noted.
- **Superior / replacement method:** Goda (1974, rev. 2000/2010) wave-pressure formula for vertical/composite caissons (incl. breaking, λ factors, Takahashi 1994 impulsive coefficient); PROVERBS/Oumeraci (2001) for true impact. `[new-calculator]`.

### 4-4 Revetment Design  [class I]
- **Validation shortcomings:** Exact (W_50 ~0.5%); mixes eras (Hudson 1958, Ahrens 1981, Van der Meer 1988, Ahrens & Heimbaugh 1988). Runup model is dated.
- **Missing data:** None noted.
- **Superior / replacement method:** Keep Van der Meer stability; swap the runup sub-calc to EurOtop (2018) Ru2%. `[drop-in]`.

### Area 5 — Runup, Transmission, Overtopping

### 5-1 Irregular Wave Runup on Beaches  [class I]
- **Validation shortcomings:** None material; exact Mase (1989) power laws reproduce Example 5-1-4 to 0.01 ft. (a,b) cited but not in the TR; published Mase set used and oracle-confirmed. Smooth impermeable slopes only.
- **Missing data:** None unrecoverable.
- **Superior / replacement method:** Stockdon et al. (2006) runup (setup + incident + infragravity swash), the modern natural-beach standard. `[drop-in]` / `[new-calculator]` (returns components).

### 5-2 Runup and Overtopping (impermeable)  [class I]
- **Validation shortcomings:** Exact vs Examples 1-7, but Weggel (1976) Q*0/α coefficients are read by the user from SPM figures (digitized chart values); dimensional Weggel formula largely superseded.
- **Missing data:** Q*0/α not built in (user-supplied from figures).
- **Superior / replacement method:** EurOtop (2018) mean-overtopping formulae (dimensionless, documented coefficients + uncertainty) and EurOtop Ru2%. `[drop-in]`.

### 5-3 Wave Transmission (impermeable)  [class I]
- **Validation shortcomings:** Exact (0.01 ft) but Seelig (1980/1976) regressions have narrow single-dataset validity.
- **Missing data:** None noted.
- **Superior / replacement method:** d'Angremond/Van der Meer/de Jong (1996), refined Van der Meer et al. (2005) (CEM VI-5-2) for low-crested/submerged structures. `[drop-in]`.

### 5-4 Wave Transmission (permeable)  [class III]
- **Validation shortcomings:** Headline K_T <0.5%, but **K_R over-predicted (~0.86 vs 0.719)** — transcribed Madsen & White (1976) seaward-slope eqs give near-total reflection; the needed dissipation calibration is not recoverable. Reference diameter tuned to half the median.
- **Missing data:** Madsen-White seaward-slope calibration.
- **Superior / replacement method:** d'Angremond/Van der Meer et al. (2005) transmission with permeability/Dn50 terms; Zanuttigh & Van der Meer (2008) Kr = tanh(a·ξ^b) for reflection. `[drop-in]`.

### 5-5 Wave Setup  [class I]
- **Validation shortcomings:** Exact radiation-stress theory (Longuet-Higgins & Stewart 1962/63) but **no ACES example/TR chapter** — analytic-only. Breaker height/depth lean on Singamsetti & Wind / Weggel.
- **Missing data:** No ACES numeric oracle.
- **Superior / replacement method:** Drive setup with a Battjes & Stive (1985)/Battjes-Janssen random-wave dissipation model, or fold into a Stockdon (2006) total-water-level framework. `[drop-in]` driver; `[new-calculator]` for a cross-shore setup profile.

### Area 6 — Littoral Processes

### 6-1 Longshore Sediment Transport  [class I]
- **Validation shortcomings:** Matches the literature CERC factor at quartz density, but the ACES examples come out ~25% high and only reproduce with a non-physical effective rho_s ~2320 (not in the TR). The regression test passes only by injecting that density.
- **Missing data:** Unexplained ACES factor / effective density.
- **Superior / replacement method:** Kamphuis (1991) transport (adds period, slope, grain size; better for fine/medium sand). `[drop-in]` (selectable kernel).

### 6-2 Time-Dependent Beach and Dune Erosion (Kriebel-Dean)  [class I]
- **Validation shortcomings:** **No ACES oracle reproduced** — substitutes Kriebel-Dean (1985) equilibrium response for the legacy XSHORE FD scheme; analytic-only vs the source paper. The no-surge ACES Example 2 (12 ft) returns ~0 by construction.
- **Missing data:** XSHORE/EBEACH source and its no-surge worked example unrecoverable.
- **Superior / replacement method:** XBeach (process-based storm dune erosion/overwash) or SBEACH (operational EBEACH successor). `[new-calculator]`, partly `[needs-data]`.

### 6-3 CEDRS Transport  [class II]
- **Validation shortcomings:** Self-chosen bin-midpoint / contributing-fraction conventions; reproduces the net to 0.58% (not to the digit) and needs the same non-physical effective rho_s ~2319.
- **Missing data:** Unstated ACES effective density; only one default station climate (G1033) embedded.
- **Superior / replacement method:** Updated WIS hindcast occurrence tables + Kamphuis (1991) per-bin kernel. `[drop-in]` kernel; `[needs-data]` for broader stations.

### 6-4 Beach Nourishment / Overfill  [class I]
- **Validation shortcomings:** None significant; James (1975) overfill ratio + renourishment factor reproduce Example 6-4 vs a true oracle. Screening estimate; fixed winnowing W=1.
- **Missing data:** None noted.
- **Superior / replacement method:** Equilibrium-profile / closure-depth fill-volume design (Dean overfill + profile translation, Hallermeier closure) and SBEACH/GenCade fill evolution. `[new-calculator]`.

### 6-5 Composite (Native/Borrow) Grain Size  [class I]
- **Validation shortcomings:** The headline ACES example is a **composite of two samples but only one ships** — composite output unreproducible. Moments verified analytically on CoreSample1; Folk measures checked for consistency only.
- **Missing data:** CoreSample2 (Panama City companion) absent.
- **Superior / replacement method:** Full GRADISTAT-style logarithmic+geometric Folk-Ward + moment statistics from laser-diffraction/settling data. `[drop-in]` stats; `[needs-data]` for the second sample.

### Area 7 — Inlet Processes

### 7-1 Spatially Integrated Inlet Hydraulics  [class III]
- **Validation shortcomings:** Peak ebb Q (0.2%), bay hydrograph (<0.02 ft), controlling velocity (~1%) meet bar, but **per-channel velocity field (Table 7-1-2) not reproduced** and mid-record exchange volumes ~6% low — section-mean friction stands in for the flow-net subdivision.
- **Missing data:** Full flow-net channel-subdivision algorithm and cross-section bathymetry (only sections 1 & 5 published).
- **Superior / replacement method:** 2-D depth-averaged inlet model (ADCIRC, Delft3D-FLOW / D-Flow FM). `[new-calculator]`, `[needs-data]`. A like-for-like supplement: Keulegan/DiLorenzo non-linear tidal-prism relations.

### 7-2 Wave-Current Interaction  [class I]
- **Validation shortcomings:** **No ACES oracle** — analytic-only (no-current identity, following/opposing trends, blocking flag). Exact linear theory but restricted to steady, collinear 1-D Doppler interaction.
- **Missing data:** None noted (no published example for this later addition).
- **Superior / replacement method:** Spectral model with wave-current coupling (SWAN/STWAVE with currents) for spreading, shear, partial blocking. `[new-calculator]`.

### Area 8 — Harbor Design

### 8-1 Properties of Rectangular Basins  [class I]
- **Validation shortcomings:** **No TR chapter / no worked example** (externally sourced); oracle is the textbook closed forms only (Merian, 2-D, Helmholtz). Idealized rectangular constant-depth basin; resonant periods only, no amplification/damping.
- **Missing data:** None (no oracle ever existed).
- **Superior / replacement method:** Harbor resonance model on real geometry (CGWAVE or BEM agitation tools) returning amplification + entrance radiation damping. `[new-calculator]`.

### 8-2 Vessel-Generated Waves  [class I]
- **Validation shortcomings:** **No ACES oracle** — analytic-only (drawdown limits, Schijf continuity+energy residual, 35.27° crest angle, Mach angle). Schijf (1949) 1-D prismatic-canal theory; no diverging/transverse wave *height*.
- **Missing data:** None noted.
- **Superior / replacement method:** PIANC WG166 (2014) and EM 1110-2-1100 Part II secondary-wave-height predictors (Sorensen/Weggel) to output wave height. `[drop-in]` (add height output).

### 8-3 Moored Vessel Motions / Surge  [class I]
- **Validation shortcomings:** **No ACES oracle** — analytic-only (single-line spring period, symmetry, virtual mass, stiffness monotonicity). Linear single-DOF surge with user-supplied C_a; no sway/yaw coupling, fender nonlinearity, catenary, or excitation.
- **Missing data:** None noted.
- **Superior / replacement method:** 6-DOF moored-ship dynamics with frequency-dependent added mass/damping and long-wave forcing per PIANC WG115 (2012) (Optimoor/Delft3D-class). `[new-calculator]`; minimally add harbor-oscillation excitation to the spring model.

### Area 9 — Storm Surge

### 9-1 Bathystrophic Storm Surge  [class III]
- **Validation shortcomings:** "Screening only." Analytic sub-models exact, but integrated surge reproduces Bodine TM-35 only to ballpark (~13.4 ft); Bodine (1971) cites a possible factor-of-two. Already flagged superseded by ADCIRC.
- **Missing data:** Bodine's hand-digitized wind-isovel field and bathymetry are not recoverable; app substitutes a parametric (Holland/Myers) wind field.
- **Superior / replacement method:** ADCIRC (named successor) or SLOSH — resolve 2-D shelf geometry/basin response. `[new-calculator]`, `[needs-data]`.

### Miscellaneous

### M-1 Miscellaneous Breaker and Steepness Routines  [class II]
- **Validation shortcomings:** TR eq-5 (max breaker height near a structure) is non-physical as printed (79 ft for a 15 ft depth), so it was replaced by a self-made substitute (depth-consistent inversion of Weggel eq-4). No standalone ACES example; analytic limits only, validated indirectly via 3-1 and 5-5.
- **Missing data:** Correct coefficients/structure for TR eq-5; no worked example.
- **Superior / replacement method:** Kamphuis (1991) breaker index and Battjes-Janssen/Goda random-wave breaking criteria for irregular seas. `[drop-in]` (alternate breaker indices).

---

## Suggested sequencing
1. **Drop-in supersessions for the design apps** (4-1 Van der Meer, 5-2/5-3/5-4 EurOtop & d'Angremond, 4-4 EurOtop runup, 5-1 Stockdon) — highest practical value, low risk, reuse existing kernels.
2. **Resolve the class-C/B residuals** (3-2 breaking model, 5-4 reflection, 1-2 distribution; document or expose the 6-1/6-3 effective density).
3. **Add cross-check oracles** for the ten analytic-only apps so the test suite guards them.
4. **New calculators** where 2-D physics is the real answer (POT/GEV extremes, Goda caisson pressures, SBEACH/XBeach erosion, ADCIRC-class surge/inlet) — larger efforts, schedule individually.
