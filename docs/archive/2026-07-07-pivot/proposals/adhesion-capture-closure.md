> **STATUS: Proposal landed 2026-07-06 via 3dde2bf; implementation Phase B still pending**

# Adhesion-Capture and Resuspension Closure Proposal (D-track, phase B)

Status: PROPOSAL for review. Not yet implemented. Author: physics-literature
research session, 2026-07-06.

Scope: closures for (a) wall attachment/capture and (b) flow-driven
resuspension of ~20 µm near-neutrally-buoyant ORGANIC/SILICONE (PDMS-like),
possibly membrane/shell particles in a water-like carrier, one-way Lagrangian
on a resolved LBM flow. Companion to `docs/DISPERSED_DEPOSITION.md` (D-track,
T18). Every candidate is written to be adoptable under the four-artifact rule
of `.claude/skills/lbmflow-physics-discipline` (Rule 1): (1) citation +
derivation, (2) validity domain, (3) a validation test vs an analytic/reference
result, (4) a PHYSICS.md entry. Draft PHYSICS.md entries are in §7.

**Physics-discipline note up front.** The current example deposits by pure
floor-crossing (a geometric perfect-sink at the plane, capture probability 1
for any particle that reaches z=0). That is defensible ONLY when transport to
the wall is the rate-limiting step AND every arriving particle sticks. For
these particles neither is obviously true: (i) near-neutral buoyancy makes
gravitational transport to the floor ~µm/s, so transport is slow and the
attachment step and the resuspension step both matter to the final pattern;
(ii) "sticks on touch" is an attachment-efficiency = 1 assumption that is only
valid under *favorable* colloid-surface chemistry. This proposal makes both
assumptions explicit, bounds their validity, and gives the minimal validated
way to relax them. No constant below is calibrated to a T18 band; each is a
material property to be supplied or a literature value with a cited range.

---

## 0. Dimensionless framing — which regime are we in?

Before choosing a model we must place the 20 µm particle on the transport map.
Order-of-magnitude with d_p = 20 µm, water (ρ_f = 1000 kg/m³, µ = 1e-3 Pa·s,
ν = 1e-6 m²/s), |Δρ| ≈ 30–50 kg/m³, a representative near-floor flow speed
U ~ 1–10 mm/s and shear rate G ~ 1–100 s⁻¹.

- **Péclet number** Pe = U d_p / D_B, with Brownian diffusivity
  D_B = k_B T / (3π µ d_p) ≈ 2e-14 m²/s at 20 µm, 293 K. Even at U = 1 mm/s,
  Pe ≈ (1e-3·2e-5)/2e-14 ≈ 1e6. **Pe ≫ 1: the particle is NON-Brownian.**
  Diffusion does not carry it to the wall.
- **Stokes settling velocity** v_s = (Δρ) g d_p² / (18 µ)
  ≈ (40·9.81·(2e-5)²)/(18e-3) ≈ 8.7e-6 m/s ≈ **~9 µm/s** (sign follows the
  sign of Δρ; if the particle is lighter than water it creams upward). This is
  the "extremely slow gravity crossing" of the brief, quantified.
- **Gravitational number** N_G = v_s / U ≈ 9e-6 / 1e-3 ≈ 1e-2 to 1e-3.
  Settling is a *weak* transport mechanism relative to advection.
- **Interception number** N_R = d_p / L_c, where L_c is the flow length scale
  (channel half-height / tray depth). For a mm-scale tray, N_R ~ 1e-2.
- **Particle Reynolds** Re_p = v_s d_p / ν ≈ 9e-6·2e-5/1e-6 ≈ 2e-4 ≪ 1
  → **Stokes drag regime** (Schiller-Naumann correction < 0.01%).
- **Stokes number** St = τ_p U / L_c with τ_p = ρ_p d_p²/(18µ) ≈ 2e-5 s → for
  U/L_c ~ 1 s⁻¹, St ~ 2e-5 ≪ 1. **Inertia negligible; particles follow
  streamlines except for the small settling drift.** (Deposition-relevant St
  in DISPERSED_DEPOSITION.md §3.4 is already tracked.)

**Consequence for model choice.** We sit at Pe ≫ 1, St ≪ 1, Re_p ≪ 1,
N_G ≪ 1. This is the *non-Brownian, interception+sedimentation, negligible
inertia* corner. The two implications that drive everything below:

1. **Brownian-diffusion deposition theory (Smoluchowski-Levich, Levich
   convective-diffusion flux, the diffusion term of colloid-filtration
   single-collector efficiency) DOES NOT APPLY** — its flux ∝ Pe^(-2/3) → 0
   here. Adopting it would be a physics error. What survives from colloid
   science is the *attachment-efficiency* concept (DLVO favorable/unfavorable),
   NOT the transport-flux formula.
2. **Transport to the floor is by advection + weak settling + interception**,
   all of which the resolved LBM flow + the existing one-way Lagrangian
   integrator already compute. So the *transport* half is resolved physics; the
   only genuine closures we must add are the **attachment** rule at the wall and
   the **resuspension** threshold. This is the correct minimal surface for new
   closures.

---

## 1. Deposition / attachment at the wall

### 1.1 What the engine already resolves (no closure needed)

Advective transport, the settling drift (Stokes/Schiller-Naumann gravity term),
and interception (a finite-radius particle whose *surface* reaches the wall while
its center is one radius away) are all resolved by the resolved flow + the
Lagrangian integrator. Per Rule 1 row 1 these need no extra artifact. The only
subtlety is the geometric capture criterion (§1.3), which is a definition, not a
closure.

### 1.2 Perfect-sink boundary condition — valid, but state the assumption

**Form.** Concentration/particle perfectly absorbed at contact: the wall is an
irreversible sink, C = 0 at the collector surface (continuum), or "remove on
arrival" (Lagrangian). Attachment efficiency α = 1.

**Source.** Standard in colloid deposition; the "perfect sink" is applied once
the surface-to-surface separation h < ~0.5 nm (primary-minimum contact). See
Elimelech, Gregory, Jia & Williams, *Particle Deposition and Aggregation*
(1995); Adamczyk & van de Ven, *J. Colloid Interface Sci.* 80 (1981) 340.

**Validity.** Correct ONLY under *favorable* deposition (no DLVO energy
barrier: attractive or weakly repulsive double layer) AND when the arriving
particle is not re-entrained. It says nothing about the transport rate — it is a
boundary condition, not a flux law. It is the right default for the D-track
*forward* model IF the user confirms favorable chemistry; otherwise §1.4.

**What is NOT valid here:** pairing the perfect sink with the
Smoluchowski-Levich / Levich diffusion flux to *predict a deposition rate*.
That flux is diffusion-limited (∝ Pe^(-2/3)) and vanishes at our Pe. Our
deposition rate is set by advection+settling+interception delivering particles
to the contact surface, which the trajectory integration already does.

### 1.3 Capture-distance (contact) criterion — the geometric core

**Form.** A particle of radius a is captured when the interpolated trajectory
brings its *center* to within the capture distance h_c of the wall:

    capture when  z_center ≤ a + δ_c

where δ_c is a small contact tolerance (physically the range of the attractive
primary minimum, ~1–10 nm; numerically it can be a fraction of dx). Setting
δ_c = 0 gives "surface touches wall" = pure interception; the current example
uses δ_c effectively = −a (center crosses the plane), which *under*-counts
interception by one radius. **Recommendation: use z_center ≤ a (surface
contact), i.e. capture within one particle radius of the wall.** This is exactly
the working hypothesis in the brief and is a definitional geometric criterion,
not a tunable closure.

**Source.** Trajectory/interception capture: Spielman, *Annu. Rev. Fluid Mech.*
9 (1977) 297; Yao, Habibian & O'Melia, *Environ. Sci. Technol.* 5 (1971) 1105
(interception mechanism η_I). The one-radius contact criterion is standard in
Lagrangian deposition codes.

**Validity.** Any Pe, any St; it is geometry. The physics content is entirely
in what happens *at* contact (α below) and whether the resolved flow near the
wall is accurate (needs adequate near-wall LBM resolution; see §4 test).

### 1.4 DLVO-based attachment efficiency — for unfavorable conditions

If the chemistry is *unfavorable* (like-charged particle and wall, low ionic
strength → a repulsive DLVO barrier), not every contact sticks: α < 1.

**Form (interaction energy).** The DLVO total interaction vs separation h:

    V_T(h) = V_vdW(h) + V_EDL(h)
    V_vdW(h) = − A a / (6 h)                    (sphere–plate, non-retarded)
    V_EDL(h) = 64 π ε a (k_B T / z e)² Γ_p Γ_c exp(−κ h)   (constant-potential, linearized)

with A = Hamaker constant (system, across water), a = particle radius,
ε = medium permittivity, κ = inverse Debye length (κ⁻¹ set by ionic strength),
Γ = tanh(z e ψ / 4 k_B T) reduced surface potentials from zeta potentials ψ.

**Form (attachment efficiency).** Two literature closures map V_T → α:

- **Maxwell/interaction-force-boundary-layer (IFBL)** — Elimelech & O'Melia,
  *Langmuir* 6 (1990) 1153; Spielman & Friedlander, *JCIS* 46 (1974) 22. Gives
  α as a ratio of deposition rate with/without the barrier; α falls
  approximately exponentially with the barrier height V_max/k_B T.
- Empirical correlations of α vs the DLVO barrier and Pe exist but are
  system-specific; do NOT hard-code one.

**Validity.** DLVO itself is quantitative for smooth surfaces, moderate ionic
strength (roughly 1 mM–0.1 M), separations ≳ 1 nm. It is known to
*under-predict* attachment for rough/soft/heterogeneous surfaces (charge
heterogeneity, "patchy" attachment) — PDMS and membrane flakes are soft and
likely rough, so treat α_DLVO as a *lower bound on stickiness* and flag it. At
Pe ≫ 1 the transport is still advective (§0); DLVO only sets α, it does not
resurrect the diffusion flux.

**Material parameters required (from user/experiment):** zeta potential of
particle and of floor, ionic strength / electrolyte composition (→ κ), Hamaker
constant A (or the contact-angle route, §2.4). Without these, α is unknowable
and we must fall back to the α = 1 favorable assumption *and say so*.

**Recommendation for phase B:** implement α as a single scalar attachment
probability with default α = 1 (favorable, perfect sink) and an OPTIONAL
DLVO-computed α when the four electrochemical inputs are supplied. Do not build
the IFBL transport correction (it belongs to the diffusion regime we are not
in). This keeps the added closure minimal and honest.

---

## 2. Resuspension / detachment under flow

Once a particle is on the floor, agitation-driven near-wall shear can detach it.
The literature consensus (Ziskind, Fichman & Gutfinger; Reeks & Hall; Soltani &
Ahmadi) is that **incipient rolling about the downstream contact edge is the
dominant detachment mode** for micro-particles — it beats lift-off and sliding
because the adhesive *moment* arm (the contact radius) is tiny while the
hydrodynamic drag acts with a lever ~particle radius.

### 2.1 Moment-balance rolling criterion (the recommended form)

**Detachment (incipient rolling) when the hydrodynamic moment about the contact
point exceeds the adhesive + gravitational resisting moment:**

    M_drag + M_lift·a  ≥  M_adh + M_grav
    (1.4 a) · F_drag   ≥  a_c · F_adh                        (lift & gravity negligible here)

- **Hydrodynamic drag on a sphere resting on a wall in a wall shear flow**
  (O'Neill 1968; Goldman, Cox & Brenner 1967, *Chem. Eng. Sci.* 22, 637 & 653):

      F_drag = 1.7009 · (6 π µ a) · u_shear ,   u_shear = G·a = (τ_w/µ)·a
             = 1.7009 · 6 π · τ_w · a²  ≈  32.05 · τ_w · a²

  where τ_w = µ G is the wall shear stress (computed by the resolved LBM near-
  wall field), G the local shear rate, a the particle radius. The factor 1.7009
  is O'Neill's exact Stokes-flow correction for a translating sphere touching a
  plane; the accompanying wall torque adds to give an **effective lever ≈ 1.4 a**
  about the contact point (the "1.399 a" of the resuspension literature).
- **Adhesive resisting moment** M_adh = a_c · F_adh, with F_adh the pull-off
  force and a_c the contact radius (lever arm), both from contact mechanics
  (§2.2).

**Resulting critical wall shear stress** (solve M_drag = M_adh):

    τ_w,crit = a_c · F_adh / (1.4 a · 6 π · 1.7009 · a²)  ∝  (a_c F_adh) / a³

Using JKR F_adh ∝ W·a and a_c ∝ (W a²/E\*)^(1/3) ∝ a^(2/3) W^(1/3) E\*^(-1/3):

    τ_w,crit  ∝  W^(4/3) · a^(-4/3) · E\*^(-1/3)
    F_drag,crit ∝ a^(2/3)         (the d^(2/3) scaling reported in the literature)

**Physical reading (behavior anchors for the validation test):**
smaller particles are HARDER to resuspend (τ_crit ↑ as a ↓); stickier surfaces
(larger W) are harder to resuspend; softer particles (smaller E\*, → larger
contact area) are harder to resuspend. All monotone, all falsifiable.

**Source.** Ziskind, Fichman & Gutfinger, "Adhesion moment model for estimating
particle detachment from a surface," *J. Aerosol Sci.* 28 (1997) 623, and their
review *J. Aerosol Sci.* 26 (1995) 613; Reeks & Hall "Rock'n'Roll" model,
*J. Aerosol Sci.* 32 (2001) 1; Soltani & Ahmadi, *J. Adhesion* 44 (1994) 161;
review: Henry & Minier, *Prog. Energy Combust. Sci.* / arXiv:1802.06448 (2018).

**Validity.** Stokes near-wall flow (Re_p ≪ 1 — satisfied here), smooth
single-asperity contact, quasi-static (threshold, not rate). Roughness and
membrane compliance shift the threshold (see §2.3 caveat); treat the smooth-JKR
τ_crit as an *order-of-magnitude threshold*, and expose it as a criterion with a
stated ±factor uncertainty rather than a sharp switch calibrated to data.

### 2.2 Contact mechanics: JKR vs DMT (which one, and why JKR here)

**JKR (Johnson–Kendall–Roberts, 1971,** *Proc. R. Soc. A* 324, 301**):**

    F_pulloff = (3/2) π W R              (pull-off / adhesion force)
    a_c(0)    = (9 π W R² / (2 E\*))^(1/3)   (contact radius at zero external load)

**DMT (Derjaguin–Muller–Toporov, 1975,** *JCIS* 53, 314**):**

    F_pulloff = 2 π W R

with W = work of adhesion (J/m²), R = particle radius (= a; sphere-on-flat
reduced radius = a), E\* = reduced modulus, 1/E\* = (1−ν₁²)/E₁ + (1−ν₂²)/E₂.

**Which model.** The Tabor parameter µ_T = (R W² / (E\*² z₀³))^(1/3) decides:
µ_T ≳ 5 → JKR (compliant, large contact, short-range adhesion); µ_T ≲ 0.1 →
DMT (stiff, small contact). **PDMS is soft (E ≈ 1–3 MPa) and a membrane/shell
is softer still → large µ_T → use JKR.** JKR also predicts a finite contact
radius at zero load and a "neck", which matters for the moment lever a_c. This
is a physics-based selection, not a fit.

### 2.3 Caveat that MUST be reported with any resuspension number

Real thresholds are dominated by (i) nanoscale roughness — asperity contact
reduces the true contact area and F_adh by 1–3 orders of magnitude vs
smooth-sphere JKR, and introduces a *distribution* of thresholds (hence Reeks &
Hall's stochastic "Rock'n'Roll"); (ii) for a soft *membrane flake* the contact
is not Hertzian at all — it is peeling/conformal, governed by a peel energy, not
JKR. So: the JKR τ_crit is a defensible *baseline scalar threshold*; a
membrane's true detachment may be a peel model (Kendall peel, *J. Phys. D* 4
(1971) 1186) which needs the flake bending stiffness. **Do not present the JKR
threshold as exact; present it as the resolved-flow shear compared against a
material-property-derived threshold with an explicit uncertainty band, and gate
only the monotone trends (§4), not an absolute detachment count, until the shape
and roughness are known.**

### 2.4 Material parameters required for resuspension

- **Work of adhesion W** (J/m²): either measured (AFM pull-off, JKR contact
  test) or from the Hamaker constant, W = A / (12 π D₀²), D₀ ≈ 0.165 nm
  (Israelachvili's universal cut-off), A = system Hamaker constant across water.
  For PDMS/organic–water–glass/silica systems, A ≈ (0.3–1.5)×10⁻²⁰ J is the
  defensible literature range (aqueous solid–solid Hamaker constants cluster at
  ~10⁻²⁰ J; polymer latices ~0.3–1×10⁻²⁰ J; silica-in-water ~0.4–0.8×10⁻²⁰ J).
  This yields W ≈ (A/(12π·(0.165e-9)²)) ≈ 1–50 mJ/m² — a wide range, hence W
  must be a *user input with a cited range*, never a hard-coded constant.
- **Reduced modulus E\*** — needs particle Young's modulus (PDMS ~1–3 MPa; a
  membrane much lower effective modulus) and Poisson ratio (~0.5 for elastomer),
  and the floor modulus (glass ~70 GPa → floor compliance negligible, E\* set by
  the particle).
- **Particle radius a** and its distribution.
- (For the peel alternative: flake thickness and bending modulus.)

**Source for W↔A and Hamaker ranges:** Israelachvili, *Intermolecular and
Surface Forces*, 3rd ed. (2011), Ch. 13–17; Hamaker range confirmations from
the aqueous-colloid literature surveyed for this proposal (polymer latex
A ≈ 3.2×10⁻²¹ J; silica-across-water A ~ few×10⁻²¹–10⁻²⁰ J).

---

## 3. Shape corrections for non-spherical (flake) particles

If the particles are thin flakes rather than beads, both the settling drag and
the near-wall detachment lever change. Two regimes.

### 3.1 Terminal-settling drag: Haider–Levenspiel sphericity correlation

**Form** (Haider & Levenspiel, *Powder Technol.* 58 (1989) 63):

    C_D = (24/Re)(1 + A Re^B) + C / (1 + D/Re)
    A = exp(2.3288 − 6.4581 φ + 2.4486 φ²)
    B = 0.0964 + 0.5565 φ
    C = exp(4.905 − 13.8944 φ + 18.4222 φ² − 10.2599 φ³)
    D = exp(1.4681 + 12.2584 φ − 20.7322 φ² + 15.8855 φ³)

with φ = sphericity = (surface area of the volume-equivalent sphere) / (actual
surface area), Re based on the volume-equivalent-sphere diameter.

**Validity.** Re < 2.6×10⁵; correlation fit over φ ∈ [0.026, 1] but most
reliable for φ ≳ 0.67; average error ~24%. **Crucially it gives an
ORIENTATION-AVERAGED drag** — it cannot represent a specific orientation.

**Relevance here.** At our Re_p ≈ 2×10⁻⁴ (deep Stokes), the (1 + A Re^B) term is
essentially 1 and the correlation collapses to a shape-dependent Stokes drag.
For deep-Stokes flakes, the *exact* Stokes results below are preferable to the
empirical Haider–Levenspiel (which was built from higher-Re falling-particle
data). **Recommendation: use Haider–Levenspiel only if Re_p is not ≪ 1; in our
Stokes regime use the analytic disc drag (§3.2).**

### 3.2 Thin-disc drag in Stokes flow and the orientation question

For a thin circular disc of radius R_d in Stokes flow (Re → 0):

    F_broadside (face-on, motion ⟂ disc)  = 16 µ R_d U
    F_edge-on   (motion ∥ disc)           = (32/3) µ R_d U
    ratio  F_broadside / F_edge-on = 3/2

**Source.** Oberbeck / Lamb, classical Stokes-flow disc resistances; Happel &
Brenner, *Low Reynolds Number Hydrodynamics* (1965), §5.

**Orientation, the subtle part.** In *pure* Stokes flow (Re_p → 0) a settling
disc has NO preferred orientation — it is torque-free and settles at whatever
orientation it starts, drifting sideways in general (the mobility tensor is
anisotropic but the torque is zero). At *finite* Re_p, inertia rotates a
settling flake to fall **broadside-on** (max drag, min settling speed). Since we
are at Re_p ≈ 2×10⁻⁴, a free flake keeps its initial orientation while settling
— so the effective settling drag is orientation-dependent and, strictly, a
single scalar drag is an approximation. **Behavior consequence:** a flake
settles 1.0–1.5× *slower* than a volume-equivalent sphere depending on
orientation; this only *reinforces* the "gravity crossing is very slow"
conclusion and makes contact adhesion even more dominant.

**Recommendation.** For phase B, if shape is unconfirmed, model the flake as a
sphere with an *orientation-averaged* Stokes drag correction (a scalar shape
factor K ∈ [1.0, 1.5] applied to Stokes drag, derived from disc theory, NOT
fitted) and DECLARE it an approximation with the validity caveat "torque-free
Stokes orientation not resolved." A full orientation-resolved flake requires
tracking the orientation vector and the anisotropic mobility tensor — defer
unless P3/P4 evidence shows the pattern depends on it. Do not fake orientation
with a case branch.

---

## 4. Candidate summary table + validation-test designs

Each row: form (§ref), source, validity, required material params, validation
test (analytic/reference anchor). Tests follow the two-layer rule (band +
behavior anchor).

| # | Closure | Form | Validity | Params from user | Validation test (analytic anchor) |
|---|---|---|---|---|---|
| C1 | Geometric contact capture (§1.3) | capture when z_center ≤ a+δ_c | any Pe/St; geometry | a (radius) | **Straight-line crossing exactness** (already T18.3c): a particle on a known trajectory is captured exactly at z=a; deposit point = analytic intercept. Behavior anchor: capture count rises when δ_c ↑. |
| C2 | Perfect sink / α=1 (§1.2) | remove on contact | favorable chemistry only | confirm favorable | **Mass conservation**: deposited+suspended = seeded (ledger exact). Behavior: 100% of arriving particles removed. |
| C3 | DLVO attachment efficiency α<1 (§1.4) | α from V_T(h) barrier via IFBL/Maxwell | smooth surfaces, 1mM–0.1M, h≳1nm | ζ_p, ζ_wall, I (ionic str.), A | **α vs barrier**: reproduce published α(V_max) curve (Elimelech & O'Melia 1990) for a chosen ζ,I within band. Behavior: α↓ monotone as ionic strength ↓ (barrier ↑); α→1 as barrier→0. |
| C4 | JKR pull-off + contact radius (§2.2) | F=1.5πWR; a_c=(9πWR²/2E\*)^⅓ | soft, µ_T≳5, smooth | W, E\*, a | **JKR self-consistency**: F_pulloff, a_c match closed form for a test (W,R,E\*); JKR→Hertz as W→0. Behavior: F↑ with W and R. |
| C5 | Rolling resuspension threshold (§2.1) | τ_crit ∝ W^(4/3) a^(−4/3) E\*^(−⅓); O'Neill 1.7009 drag, lever 1.4a | Re_p≪1, single-asperity, threshold | W, E\*, a | **Scaling test**: for a family of (a,W,E\*), computed τ_crit reproduces the analytic exponents (−4/3 in a, +4/3 in W, −1/3 in E\*) AND F_drag,crit ∝ a^(2/3). Behavior: smaller a harder to detach; agitation above τ_crit resuspends, below does not (step response). Reference cross-check: an order-of-magnitude comparison to a published τ_crit(d) curve (Ziskind et al. 1997). |
| C6 | Haider–Levenspiel drag (§3.1) | C_D(Re,φ) 4-coeff form | Re<2.6e5, φ∈[0.026,1], orient-avg | φ (sphericity) | **Sphere limit**: φ=1 recovers standard C_D(Re) (Schiller-Naumann) within band. Behavior: C_D↑ as φ↓ at fixed Re. (Use only if Re_p not ≪1.) |
| C7 | Thin-disc Stokes drag / shape factor (§3.2) | F_bs=16µRU, F_edge=(32/3)µRU | Re_p≪1 | shape=disc, R_d | **Disc resistance**: broadside/edge drag ratio = 3/2 exactly; settling speed 1.0–1.5× slower than volume-equiv sphere. Behavior: flake settles slower than bead. |

**Validation-infrastructure note.** C1, C2, C4, C6, C7 are self-contained
analytic checks (no experimental data needed) and can be implemented
immediately as Rust unit/integration tests in the T18 family. C3 and C5 need a
published reference *curve* to compare against; both are digitizable from the
cited papers, but C5's absolute magnitude is roughness-sensitive so its gate
must be the *scaling exponents and the step-response behavior*, not an absolute
τ_crit value (per §2.3). Every run must emit a visual artifact (deposition
density map / near-wall shear field) per the physics-discipline post-run review.

---

## 5. Recommendation — the minimal defensible stack for phase B

The brief's working hypothesis was: *contact-capture within one particle radius
of the wall, attachment probability 1 under favorable conditions, plus a
JKR-based critical-shear resuspension threshold.* **I confirm this hypothesis as
the correct minimal stack, with three refinements and one explicit dependency.**

**Adopt now (phase B):**

1. **C1 geometric contact capture at z_center ≤ a** (one radius). This corrects
   the current example's crossing test (which under-counts interception by a
   radius) and is pure geometry — no closure debt. *Refinement 1: capture at the
   surface (z=a), not the center-crossing (z=0).*
2. **C2 perfect sink with α = 1** as the DEFAULT attachment rule, explicitly
   labeled "favorable-chemistry assumption." *Refinement 2: make α a named
   scalar parameter (default 1) so C3 can later supply α<1 without a code
   branch keyed to case identity.* This keeps the ban-list clean.
3. **C4 + C5 JKR rolling-resuspension threshold** as the detachment criterion,
   using the resolved near-wall τ_w from the LBM field, O'Neill's 1.7009 drag,
   the 1.4a lever, and JKR (justified by PDMS/membrane softness, §2.2).
   *Refinement 3: gate it on the SCALING and STEP-RESPONSE behavior (C5 test),
   not on an absolute detachment count, and report τ_crit with a stated ±factor
   from roughness/shape uncertainty (§2.3).* Resuspension only activates under
   agitation (τ_w > τ_crit); quiescent settle never resuspends.

**Explicitly DEFER (do NOT build speculatively):**

- **Smoluchowski-Levich / Levich diffusion flux** — WRONG regime (Pe≫1); would
  be a physics error. Do not implement.
- **IFBL transport correction** — same diffusion-regime reason.
- **DLVO α<1 (C3)** — build only when the user supplies ζ-potentials + ionic
  strength; until then α=1 is the honest default and C3 is a plug-in.
- **Orientation-resolved flake dynamics** — defer to a P3/P4 trigger (build only
  if evidence shows the pattern depends on flake orientation); until then a
  scalar disc-derived shape factor K∈[1,1.5] (C7), declared approximate.
- **Peel/membrane detachment model** — defer until shape is confirmed as
  membrane; JKR is the baseline, with the §2.3 caveat attached.

**Why this is defensible under the discipline.** Transport (advection, settling,
interception) is resolved physics. The only genuine closures are α (default 1,
a stated assumption, upgradeable to validated DLVO) and the JKR rolling
threshold (a literature model with all four artifacts). No constant is
calibrated to a T18 band; every material number is a user input with a cited
range. No case-identity branch (α and W are parameters, not per-sample
switches). No transport-absorbing clamp (capture removes the particle and logs
it — accounted, not clamped). If, at gate time, a T18 band cannot be met without
tuning W or α to hit it, that is a STOP-RULE event (spec revision), not a fit.

**Known open risk to surface now:** if the particles are actually lighter than
water (Δρ < 0, they cream), the floor deposition the tool targets may be
physically minor compared to *ceiling*/free-surface accumulation, and the whole
"deposit on floor" framing needs the user's confirmation of the sign of Δρ
before phase B is worth building. Flag this to the PM/user.

---

## 6. Measurements to request from the user

Ranked by leverage on the model:

1. **Sign and magnitude of Δρ (particle vs carrier density).** Decides whether
   deposition is on the floor at all (§5 open risk). Highest priority.
2. **Particle shape** (bead / shell / membrane flake) and size distribution.
   Selects §3 treatment and JKR-vs-peel (§2.3).
3. **Work of adhesion W** OR the inputs to compute it: **Hamaker constant** (or
   contact angles of carrier on particle and on floor → surface energies →
   W via Dupré/Young), for C4/C5. Without this the resuspension threshold is
   only a scaling law, not a number.
4. **Particle Young's modulus E** (PDMS 1–3 MPa? membrane much lower) and
   Poisson ratio, for E\* in JKR.
5. **Zeta potential of particle and floor + ionic strength / electrolyte** of
   the carrier — only needed to move off α=1 (C3). If the process runs in
   near-DI water with like charges, α<1 matters; in salty/attractive conditions
   α=1 is fine.
6. **Any observed macroscopic timescale** — an observed settling/clearing time
   or an observed resuspension onset (agitation speed at which deposit lifts).
   This is a *validation cross-check* (compare predicted v_s and τ_crit to
   observation), NOT a calibration target — used to confirm the model is in the
   right order of magnitude, per the two-layer discipline.

---

## 7. Draft PHYSICS.md entries (to be added on implementation, not now)

```markdown
### <date> Contact-capture attachment (particles.rs / deposition step)
- Form: capture particle when z_center ≤ a + δ_c (surface contact, one radius);
  attachment probability α applied on contact (default α = 1, favorable sink).
- Source: Spielman Annu. Rev. Fluid Mech. 9 (1977) 297; Elimelech et al.,
  Particle Deposition & Aggregation (1995) — perfect-sink primary-minimum.
- Validity domain: Pe ≫ 1 (non-Brownian; transport advective+settling, NOT
  diffusive — Smoluchowski-Levich explicitly inapplicable); α=1 valid only for
  favorable DLVO chemistry.
- Validation: t18_x_capture.rs — straight-line crossing exact (existing T18.3c);
  mass ledger exact; α sweep changes deposit count monotonically.
- Replaces / interacts with: replaces center-crossing test (under-counted
  interception by one radius); α is a named parameter, not a case branch.

### <date> JKR rolling-resuspension threshold (particles.rs)
- Form: detach when 1.4a·F_drag ≥ a_c·F_adh, F_drag = 1.7009·6π·τ_w·a²
  (O'Neill 1968), F_adh = 1.5πWa, a_c = (9πWa²/2E*)^(1/3) (JKR 1971);
  τ_crit ∝ W^(4/3) a^(−4/3) E*^(−1/3).
- Source: Ziskind, Fichman & Gutfinger, J. Aerosol Sci. 28 (1997) 623;
  Reeks & Hall, J. Aerosol Sci. 32 (2001) 1; O'Neill, Chem. Eng. Sci. 23
  (1968) 1293; Johnson, Kendall & Roberts, Proc. R. Soc. A 324 (1971) 301.
- Validity domain: Re_p ≪ 1 (Stokes near-wall), single-asperity smooth contact,
  soft particle (Tabor µ_T ≳ 5 → JKR); threshold is order-of-magnitude
  (roughness/membrane compliance shift it by up to orders — reported ± factor).
- Validation: t18_x_resuspension.rs — scaling exponents (−4/3 in a, +4/3 in W,
  −1/3 in E*), F_drag,crit ∝ a^(2/3); step response (τ_w>τ_crit resuspends,
  below does not); JKR→Hertz as W→0.
- Replaces / interacts with: adds detachment to the deposit-on-contact model;
  requires resolved near-wall τ_w from the LBM field.
```

---

## 8. Sources

- Elimelech, Gregory, Jia & Williams, *Particle Deposition and Aggregation:
  Measurement, Modelling and Simulation* (Butterworth-Heinemann, 1995).
- Elimelech & O'Melia, "Kinetics of deposition of colloidal particles in porous
  media," *Environ. Sci. Technol.* 24 (1990) 1528; *Langmuir* 6 (1990) 1153.
- Yao, Habibian & O'Melia, "Water and waste water filtration: concepts and
  applications," *Environ. Sci. Technol.* 5 (1971) 1105 (interception η_I).
- Spielman, "Particle capture from low-speed laminar flows," *Annu. Rev. Fluid
  Mech.* 9 (1977) 297.
- Tufenkji & Elimelech, "Correlation equation for predicting single-collector
  efficiency," *Environ. Sci. Technol.* 38 (2004) 529 (shows η_diffusion → 0 at
  high Pe; interception/gravity dominate for large particles).
- Ziskind, Fichman & Gutfinger, "Adhesion moment model for estimating particle
  detachment from a surface," *J. Aerosol Sci.* 28 (1997) 623; review
  *J. Aerosol Sci.* 26 (1995) 613.
  https://www.sciencedirect.com/science/article/abs/pii/S0021850296004600
- Reeks & Hall, "Kinetic models for particle resuspension in turbulent flows:
  theory and measurement," *J. Aerosol Sci.* 32 (2001) 1.
  https://www.sciencedirect.com/science/article/abs/pii/S002185020000063X
- Soltani & Ahmadi, "On particle adhesion and removal mechanisms in turbulent
  flows," *J. Adhesion Sci. Technol.* 8 (1994) 763.
- Henry & Minier, "Progress in particle resuspension from rough surfaces by
  turbulent flows" / "Particle resuspension from complex surfaces: current
  knowledge and limitations," arXiv:1802.06448 (2018).
  https://arxiv.org/pdf/1802.06448
- O'Neill, "A sphere in contact with a plane wall in a slow linear shear flow,"
  *Chem. Eng. Sci.* 23 (1968) 1293; Goldman, Cox & Brenner, *Chem. Eng. Sci.*
  22 (1967) 637 & 653 (near-wall drag/lift/torque).
  https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/lift-and-drag-forces-acting-on-a-particle-moving-with-zero-slip-in-a-linear-shear-flow-near-a-wall/83AABC94313E41AC0C190AB88264BC86
- Johnson, Kendall & Roberts, "Surface energy and the contact of elastic
  solids," *Proc. R. Soc. A* 324 (1971) 301 (JKR).
- Derjaguin, Muller & Toporov, "Effect of contact deformations on the adhesion
  of particles," *J. Colloid Interface Sci.* 53 (1975) 314 (DMT).
- Kendall, "Thin-film peeling — the elastic term," *J. Phys. D* 4 (1971) 1186
  (membrane peel alternative).
- Israelachvili, *Intermolecular and Surface Forces*, 3rd ed. (Academic Press,
  2011) — W = A/(12πD₀²), Hamaker constants across water, Tabor parameter.
  https://www.sciencedirect.com/topics/engineering/hamaker-constant
- Haider & Levenspiel, "Drag coefficient and terminal velocity of spherical and
  nonspherical particles," *Powder Technol.* 58 (1989) 63.
  https://www.semanticscholar.org/paper/7f86e05b6c74884e08491f84943584d48331fb9c
- Happel & Brenner, *Low Reynolds Number Hydrodynamics* (Prentice-Hall, 1965)
  — thin-disc Stokes resistances, anisotropic mobility.
- Review — Drag Coefficients of Non-spherical and Irregularly-Shaped Particles,
  ResearchGate 368865234 (Haider–Levenspiel accuracy ~24%, Re<2.6e5 range).
