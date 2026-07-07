# D3Q19 cumulant / central-moment Galilean-invariance finding: derivation and fix decision

Status: PROPOSAL (theorist analysis, no code changed by this document)
Author role: LBM theorist
Scope: `CollisionKind::CentralMoment` on **D3Q19 only** (D3Q27 passes all holdouts)
Related: `docs/PHYSICS.md` 2026-07-06 / 2026-07-07 entries, ANOM-P4-008;
`crates/lbm-core/tests/cumulant_holdout.rs`; `crates/lbm-core/tests/cumulant_acceptance.rs`.

---

## 0. What the live code actually does (audited against source, 2026-07-07)

The task brief quotes `omega_eff = omega_shear * (1 + 0.0025 - 0.16 |u|^2)`. **That form
is stale.** ANOM-P4-008 already removed the `+0.0025` offset from all three emission sites.
The live term, identical on CPU-scalar, CPU-SIMD and generated GPU WGSL, is:

- `crates/lbm-core/src/kernels.rs:397-404`
- `crates/lbm-core/src/backend_simd.rs:634-641`
- `crates/lbm-core/src/gpu/wgsl.rs:641-649`

```
os_base = omega_shear                       (or per-cell WALE omega)
velocity_correction = -0.16 * |u|^2         (|u|^2 = ux^2 + uy^2 + uz^2)   [unless ablation flag]
os = min(2.0, os_base * (1.0 + velocity_correction))
```

`os` is the relaxation rate applied to **all** second-order deviatoric central moments
(the diagonal deviator via `kernels.rs:436-441`, the off-diagonals `xy`,`xz`,`yz` via the
`order == 2 => os` arm at `kernels.rs:419-424`). The bulk trace relaxes at rate 1.0.
This is emitted for **both** D3Q19 and D3Q27; the term is not gated by lattice in the
current code — but the finding, calibration, and this analysis are about D3Q19, because
D3Q27 has no third-order defect to correct (see §1) and passes the holdouts either way.

The effective shear viscosity that this rate produces is the standard LBM relation
applied to the corrected rate:

```
nu_eff(u) = ( 1/os(u) - 1/2 ) * cs^2 ,   cs^2 = 1/3
```

Because `velocity_correction < 0`, `os(u) < omega_shear`, so `nu_eff(u) > nu_target`, and the
excess **grows as `|u|^2`**. The correction was calibrated so that, in the *rest-frame*
advected-TGV Galilean-defect metric of `cumulant_acceptance.rs`
(`galilean_defect = |rate(u=0.05) - rate(u=0)| / rate(u=0)`), the D3Q19 CentralMoment defect
(`9.996e-4`) sits below half the BGK baseline (`2.539e-3`). That is a **training metric**:
`-0.16` was chosen to pass it. The holdout `advected_tgv3d_decay_rate_is_frame_independent`
is a different, adversarial probe and it FAILS.

**Provenance verdict (physics-discipline decision table):** `-0.16` is a constant calibrated
to pass an acceptance band. It has no recorded derivation, no validity domain, and no
independent validation test — it is a **banned ad-hoc constant**, on exactly the same footing
as the `+0.0025` that ANOM-P4-008 already removed. The only reason it survived that sweep is
that it is "compile-visible ablatable" and normal builds keep it "so the Galilean-defect
acceptance stays covered" — i.e. it is retained *because* it passes its own calibration
target, which is precisely the prohibited pattern.

---

## 1. WHY: derivation of the frame-dependent viscous error

### 1.1 The exact object: the third-order equilibrium moment tensor on D3Q19

The Chapman–Enskog viscous stress in LBM is set by the relaxation of the second-order
non-equilibrium moment `Pi^(1)_{ab}`. Its spatial-gradient closure is fed by the
**third-order equilibrium moment** through the streaming operator `c_g d_g`:

```
d_t Pi^(0)_{ab} + d_g Q^(0)_{abg} = -omega ( Pi^(1)_{ab} ) / dt + ...
```

where `Q^(0)_{abg} = sum_i c_ia c_ib c_ig f_i^(eq)`. For the second-order Hermite (or discrete)
equilibrium with mean velocity `u`, the *continuum* target third moment is the Galilean-covariant

```
Q^(0),cont_{abg} = rho cs^2 ( u_a delta_bg + u_b delta_ag + u_g delta_ab )        (Eq. 1)
```

which is fully symmetric and isotropic in `(a,b,g)`. Substituting Eq. 1 into the C–E closure
gives the exact Navier–Stokes stress with `nu = cs^2 (1/omega - 1/2)`, **independent of `u`**.
That is Galilean invariance of the viscosity.

### 1.2 The D3Q19 defect: the diagonal cubic moment `Q_{aaa}` is not representable

The lattice can only reproduce Eq. 1 for the moments its velocity set spans isotropically.
On a lattice, `sum_i w_i c_ia c_ib c_ig c_id` (the 4th-order weight tensor that generates the
`u`-linear part of `Q^(0)`) must equal the isotropic `cs^4 (delta_ab delta_gd + delta_ag delta_bd
+ delta_ad delta_bg)` for Eq. 1 to hold. For **D3Q27** it does. For **D3Q19** it does **not**,
because the 8 body-diagonal velocities `(+-1,+-1,+-1)` are absent (`central_basis` drops exactly
these: `lattice.rs` omits `ax>0 && ay>0 && az>0`; the basis correspondingly drops the `xyz`,
`xxy`, ... corner moments). The concrete failure is in the **pure-diagonal** components:

```
D3Q27:  sum_i w_i c_ix^4 = 1            (isotropic value cs^4*3 = 1/3 ... consistent set)
D3Q19:  sum_i w_i c_ix^4 = 1  BUT  the cross weight sum_i w_i c_ix^2 c_iy^2  is short,
        so the 4th-order tensor is NOT isotropic: it carries a residual
        diagonal excess  Delta_{aaaa} = sum_i w_i c_ia^4 - 3 (sum_i w_i c_ia^2 c_ib^2)  != 0.
```

Equivalently, the *third-order equilibrium moment* on D3Q19 acquires a lattice error whose
`u`-linear leading term is **not** the isotropic Eq. 1 but a form with an anisotropic diagonal
correction. Writing the defect as `Q^(0),D3Q19_{abg} = Q^(0),cont_{abg} + dQ_{abg}`, the leading
`dQ` is diagonal-dominant: it modifies `Q_{aaa}` (the `c_a^3` moment) differently from the
mixed `Q_{aab}` moments. This is the **classical D3Q19 cubic defect** (Geier et al. 2015,
cumulant paper, appendix on the "defective" third-order moments of D3Q19; also White & Chong,
Dellar; the standard remedy is the cubic equilibrium correction of the diagonal
`c_a^3`-moments).

### 1.3 Consequence: the error in `nu_eff` is a TENSOR in `u_a u_b`, not a scalar in `|u|^2`

Propagating `dQ_{abg}` through the C–E closure, the spurious contribution to the deviatoric
stress `Pi^(1)_{ab}` is second order in velocity and **inherits the diagonal anisotropy of `dQ`**.
For a single Fourier shear mode `u ~ exp(i k.x)` superposed on a uniform frame velocity `U`,
the frame velocity enters the closure through `d_g Q^(0)_{abg}` where the `u`-linear part of `Q`
is evaluated at the *total* velocity `U + u'`. The frame-dependent piece of the effective
viscosity tensor is therefore

```
delta nu_eff_{ab}(U) = -cs^2 * C_lat * ( alpha * U_a U_b  +  beta * |U|^2 delta_ab )    (Eq. 2)
```

with `C_lat` a lattice constant fixed by `Delta_{aaaa}` (D3Q19: nonzero; D3Q27: exactly zero,
which is why D3Q27 passes). The **directional part `alpha U_a U_b`** is the essence: the spurious
viscosity is larger along the frame-motion axis than transverse to it. **A scalar
`-0.16 |u|^2` correction has only the `beta |U|^2 delta_ab` structure — it is isotropic and
therefore cannot cancel the `alpha U_a U_b` anisotropy.** It can null the *trace* of the defect
at one calibration amplitude, which is exactly why the rest-frame Galilean-defect training
metric (which averages a mode over all orientations) can be driven to ~1e-3, while the holdout
(a fixed x-directed frame `U = (u_frame,0,0)` acting on a fixed-orientation TGV shear) sees the
uncancelled anisotropic residual and the error grows `~ u_frame^2`. The measured holdout errors
`{1.83e-3, 2.91e-3, 6.04e-3}` at `u_frame = {0, 0.05, 0.1}` fit `err = err_0 + kappa*u_frame^2`
with `err_0 ≈ 1.83e-3` and `kappa ≈ 0.42`, i.e. clean quadratic-in-frame growth — the signature
of an **uncorrected `U_a U_b` term**, not of a mis-tuned scalar coefficient.

**Sign/magnitude of the "right" scalar:** even the best scalar can only match `beta`. Fitting the
isotropic part of Eq. 2 to the D3Q19 weight defect gives a coefficient of the correct sign
(positive `delta nu`, hence negative `velocity_correction`) and O(0.1) magnitude in `omega`-space
— which is why `-0.16` "works" for the trace and is not crazy — but the coefficient is
**not uniquely derivable as a scalar**, because the true object is a rank-2 tensor. So the
current term has the right *sign* and *rough magnitude* for the trace, but the **wrong structure**.

### 1.4 Why this is a collision fix, not a relaxation-rate fix, in principle

The defect lives in the **equilibrium third moment** `Q^(0)`. The physically-correct remedy is to
**repair `Q^(0)` itself** (add the missing diagonal cubic term to the equilibrium / third-order
central-moment target), not to bend the *second-order relaxation rate* `omega` to compensate its
downstream effect. Bending `omega` (what `-0.16|u|^2` does) is a lumped, amplitude- and
orientation-specific approximation to a gradient closure error — structurally it can never be
exact, only calibrated. This is the root cause of the finding.

---

## 2. Candidate fixes

### (a) Derivation-backed tensorial relaxation correction

Replace the scalar `os` with a per-shear-moment rate carrying the `U_a U_b` structure of Eq. 2:
`os_{ab} = omega_shear * (1 - C_lat(alpha u_a u_b/|u|... + beta |u|^2))`, applied component-wise
to the five deviatoric second-order moments (`xx-yy`, `xx-zz`, `xy`, `xz`, `yz`).
**Derivability:** the tensor structure IS derivable (Eq. 2); `C_lat`, `alpha`, `beta` are fixed
by `Delta_{aaaa}` of D3Q19. **But** the mapping from "gradient-closure stress defect" to
"per-moment relaxation-rate multiplier" is itself only leading-order: it assumes the neq moment
is quasi-static and the mode is a single wavevector. Off that assumption (multi-mode, boundaries,
non-diffusive scaling) the compensation drifts. So (a) is *more* derived than the scalar but still
a **relaxation-rate proxy for an equilibrium defect** — it narrows the residual to higher order in
`(k dx)` and `Ma`, not to zero, and it introduces per-moment rate machinery (breaks the single-`os`
structure that `backend_simd_equiv` and the WGSL emitter assume; touches the diagonal-deviator
split at `kernels.rs:436-441` and the GPU codegen at `wgsl.rs:651-662`). High implementation cost,
still not exact.

### (b) Best achievable scalar `os` correction

Any scalar can only cancel `beta |U|^2 delta_ab`, leaving the full anisotropic residual
`-cs^2 C_lat alpha U_a U_b`. Predicted residual frame error for an x-directed frame acting on the
fixed-orientation TGV: `~ alpha/(alpha+3 beta)` of the *current* uncorrected defect — i.e. the
holdout would still show quadratic-in-`u_frame` growth with the same `kappa` (the anisotropic part
is untouched by construction). **Predicted holdout spread stays O(4e-3), still above the
1.157e-3 band.** No scalar retune (including re-deriving `beta` from the weight defect instead of
calibrating it) passes the holdout. This is a dead end for the holdout by construction, and any
nonzero scalar remains a banned calibrated constant unless `beta` is derived — and even the derived
`beta` fails the gate. **Reject.**

### (c) Drop the empirical term on D3Q19; route Galilean-sensitive cases to D3Q27

Set `velocity_correction = 0` for D3Q19 (equivalently make the existing ablation flag the default
for D3Q19). What is lost at rest-frame TGV: the `cumulant_acceptance` Galilean-defect metric for
D3Q19 CentralMoment rises from `9.996e-4` back toward the **BGK-equivalent** level, because with no
correction the D3Q19 central-moment operator has the *same* uncorrected third-order defect as BGK
in its trace (the central-moment machinery does not by itself repair `Q^(0)`). From the recorded
numbers, expect the D3Q19 CentralMoment rest-frame defect to land near the BGK baseline
`~2.5e-3` (the ablation E1 run will pin the exact value). That would **fail** the current
`cumulant_improves...galilean_invariance` "< 0.5 * BGK" claim for D3Q19 — so this path is coupled to
narrowing that claim (D3Q27 still passes it; D3Q19's central-moment advantage over BGK for
Galilean invariance is simply **not real** without a corrected equilibrium). Nothing is lost for the
*rest-frame decay-rate accuracy* (T15 class): the `cumulant_acceptance` viscosity smoke is already
**re-frozen to the UNCORRECTED N=32 value** (`nu_eff = 2.0455e-2`) per ANOM-P4-008, so removing the
term does not regress the viscosity gate at all. This is the honest state.

### (d) Corrected equilibrium (add the missing third-order moment where representable)

Repair `Q^(0)` directly: add the cubic diagonal correction to the equilibrium so that the
D3Q19 third moment matches Eq. 1 in its representable subspace. For D3Q19 the standard result
(Geier 2015 cumulant appendix; the "3rd-order cumulant" correction, and the older cubic-defect
fixes of Dellar / White-Chong) is to add `-omega_3 * (Delta_{aaaa} term)` acting on the
**third-order central moments** `m_{aab}` — i.e. relax the *third-order* cumulants toward the
Galilean-corrected target rather than toward the plain discrete-Hermite value the code currently
uses. Crucially the D3Q19 basis **does** retain the mixed third moments `xxy, xxz, yyx, yyz, zzx,
zzy` (only the `xyz` corner is dropped), and the cubic Galilean defect Eq. 2 is generated precisely
by those representable mixed third moments. So the correction **is** expressible on D3Q19's own 19
moments — it is a change to the third-order relaxation targets (currently hardwired to rate 1.0
toward the discrete equilibrium at `kernels.rs:419-424`, `order >= 3 => rate 1.0`), not new storage.
This is the **only** candidate that attacks the root cause (§1.4) and is representable on D3Q19.
Its coefficients are fully derived from the lattice weight tensor (no calibration). Cost: it changes
the third-order target/rate, so it touches the same three emission sites plus the target
construction, and must re-clear `backend_simd_equiv` (bit/threshold) and T13. It is the correct
long-term fix but is a genuine collision-operator extension, not a one-line change.

---

## 3. RECOMMENDATION — exactly one path

**Recommended: (c) now — remove the empirical `-0.16 |u|^2` term on D3Q19 and narrow the claim —
with (d) recorded as the derivation-backed follow-up for anyone who needs Galilean-invariant
viscosity on D3Q19.**

Rationale, decided not surveyed:
- (a) and (b) keep a relaxation-rate proxy for an equilibrium defect; (b) provably fails the holdout;
  (a) does not close cleanly (residual at higher order, high cost, still a proxy).
- (d) closes it at the derivation level but is a real collision-operator change; shipping it *now*
  under time pressure would itself risk an under-tested constant. It is the right *next* order, not
  this document's emergency action.
- The physics-discipline prime directive forces the immediate action regardless of (d)'s timeline:
  `-0.16` is a banned calibrated constant (§0). It must be removed, exactly as `+0.0025` was, and it
  cannot be "retained because it passes its calibration metric." The `-0.16` term does **not**
  establish Galilean-invariant viscosity (holdout FAILS) and does **not** improve rest-frame
  viscosity accuracy (that gate is frozen to the uncorrected value). It buys only a
  calibrated-to-pass number on its own training metric. Remove it.

### 3.1 Implementation plan for (c)

Files (D3Q19 path only; leave D3Q27 mathematically unaffected — see note):
1. `crates/lbm-core/src/params.rs:16` — this is the cleanest lever: the ablation flag already
   zeroes the term. Make the **default behavior** for D3Q19 be "no velocity correction." Because the
   term is emitted lattice-agnostically today, the correct change is to gate the correction on
   `L::Q == 27` (D3Q27 has `Delta_{aaaa}=0`, so the term is identically unjustified there too and
   should also go — see note) OR simply set the term to `0.0` unconditionally and delete the flag.
   Decision: **set `velocity_correction = 0.0` unconditionally and remove the term + flag**, because
   for D3Q27 the term corrects a defect that is exactly zero (so it is pure calibration there as well)
   and for D3Q19 it fails the holdout. This makes the central-moment operator carry **no** calibrated
   velocity constant on any lattice.
2. `crates/lbm-core/src/kernels.rs:398-404` — replace the `velocity_correction` block with
   `let os = os_base.min(2.0);`.
3. `crates/lbm-core/src/backend_simd.rs:635-641` — same replacement.
4. `crates/lbm-core/src/gpu/wgsl.rs:641-649` — emit `let os = min(2.0f, os_base);`.
5. `crates/lbm-core/src/params.rs` — delete `CENTRAL_MOMENT_DISABLE_VELOCITY_CORRECTION_FOR_ABLATION`
   and its references (it exists only to toggle the removed term).

Note (D3Q27): removing the term on D3Q27 changes its numbers slightly; D3Q27 currently PASSES the
holdouts *with* the term, so confirm D3Q27 still passes after removal (it must, since the term was
correcting a null defect — expect D3Q27 holdout spread to move by O(1e-4), staying inside band).
If for any reason a reviewer wants to minimize churn, gating the removal to `L::Q == 19` is an
acceptable fallback, but the principled action is total removal.

Claim narrowing (the `cumulant_acceptance` Galilean claim must change, since D3Q19 will no longer
beat 0.5*BGK):
- `cumulant_improves_advected_tgv3d_galilean_invariance` must be re-scoped to assert the
  improvement **only for D3Q27** (which genuinely has no cubic defect). For D3Q19 the assertion
  becomes: the CentralMoment defect is **not worse** than BGK (it should be ~equal), NOT
  "< 0.5 * BGK". This is the honest, derivation-consistent claim.

Gates (must pass WITHOUT retuning any band):
- `cumulant_holdout.rs` — all three holdouts at their **already-derived** bands. The advected holdout
  must now pass its `Ma^2 (k dx)^2 = 1.157e-3` band (prediction in §4). Off-Re and D3Q19-vs-D3Q27
  already pass and must stay passing.
- `cumulant_acceptance.rs` — viscosity smoke is frozen to the **uncorrected** value already, so it
  passes unchanged. `assert_rest_exact` / `assert_uniform_exact` unaffected (no velocity gradient,
  term was ~0 there anyway). The Galilean-improvement test must be updated per the narrowing above.
- `backend_simd_equiv.rs` and T13 — removing a term that was identical across CPU/SIMD/GPU keeps them
  bit/threshold-consistent; still run them.
- Rest-frame T15 TGV (`t15_3d.rs` T15.4) — must not regress. It uses `u=0` frame so the term was
  inactive there; expect **no change**.

### 3.2 PHYSICS.md entry text (to be added by the implementer, not by this document)

```
### 2026-07-07 — ANOM-P4-009: D3Q19 cumulant velocity term (-0.16|u|^2) removed

Form removed: omega_eff = omega_shear * (1 - 0.16 |u|^2), clamped <= 2, applied
to all 2nd-order deviatoric central moments on CPU scalar/SIMD + generated GPU WGSL.
Provenance: -0.16 was calibrated to pass the rest-frame advected-TGV Galilean-defect
band in cumulant_acceptance.rs. It is a banned calibrated constant (same class as the
+0.0025 removed in ANOM-P4-008): no derivation, no validity domain, and it FAILS the
adversarial holdout advected_tgv3d_decay_rate_is_frame_independent (frame spread
4.195e-3 > derived band 1.157e-3, error growing ~u_frame^2).
Root cause (derivation): D3Q19 lacks the 8 body-diagonal velocities, so its 4th-order
weight tensor is anisotropic (Delta_{aaaa} != 0), making the third-order equilibrium
moment Q^(0)_{abg} deviate from the Galilean-covariant Eq. rho cs^2 (u_a d_bg + ...).
The resulting frame-dependent viscosity error is a TENSOR ~ U_a U_b (Eq. 2), not a
scalar |u|^2; a scalar omega correction can null only its trace at one calibration
amplitude, leaving the anisotropic U_a U_b residual that the holdout measures.
D3Q27 has Delta_{aaaa} = 0 (no defect) and is unaffected; the term was pure calibration
there too and was removed on all lattices.
Claim narrowed: cumulant Galilean-invariance improvement over BGK is asserted for D3Q27
only. On D3Q19 the central-moment operator does NOT establish Galilean-invariant
viscosity without a corrected equilibrium; its defect is ~ the BGK level.
Follow-up (derivation-backed, not yet implemented): repair Q^(0) via the D3Q19 cubic
third-order-moment correction (Geier et al. 2015 cumulant appendix; Dellar; White-Chong),
expressible on the retained mixed 3rd moments xxy,xxz,yyx,yyz,zzx,zzy (only xyz is
dropped). Coefficients derived from the lattice weight tensor, zero calibration.
Validity domain: no D3Q19 velocity-dependent viscosity closure is live. tau = 3 nu + 0.5
is the only viscosity relation. Galilean-sensitive advected-flow use cases should select
D3Q27 until the corrected-equilibrium follow-up lands.
```

---

## 4. Predicted holdout outcomes for the recommended path (falsifiable, NOT tuned-to-pass)

Prediction basis: with the calibrated term removed, the D3Q19 central-moment operator relaxes the
second-order moments at the clean `omega_shear` with the discrete-Hermite equilibrium — identical in
structure to the operator whose *rest-frame* uncorrected decay was measured in ANOM-P4-008 at
`nu_eff = 2.0455e-2` (rel err `2.27e-2`). Removing a `-0.16|u|^2` term that at the holdout amplitudes
(`u0 = 0.012`, frames up to 0.1) contributes at most `-0.16 * (0.1^2) ≈ -1.6e-3` relative to `omega`
means the frame-*dependent* part of the rate is no longer being amplified by an anisotropic residual;
the remaining frame dependence is the genuine second-order spatial truncation of the single Fourier
mode, `O((k dx)^2)`, exactly the band derivation in the holdout.

Quantitative predictions (the implementation order should treat these as the falsifiable target):

- **Advected holdout `advected_tgv3d_decay_rate_is_frame_independent`:**
  - The three per-frame rates converge toward a common value near the uncorrected rest-frame rate.
    Expect all three within `+-0.3%` of each other.
  - **Predicted frame spread `(max-min)/mean` ≈ 3e-4 to 9e-4**, i.e. **below** the derived band
    `Ma^2 (k dx)^2 = 1.156594266e-3`. (Center estimate ~6e-4; the point is it drops from 4.195e-3 to
    sub-band.) **PASS is the prediction.** If the removed term were actually correcting real physics,
    removal would *increase* the spread above 4.195e-3 — that is the falsifier.
  - Rest-frame (`u_frame=0`) rate stays essentially unchanged from `4.6349e-3` (the term is ~0 at
    zero mean velocity and small `u0`); the *advected* rates move DOWN toward it instead of drifting
    up, because the spurious `|u|^2`-amplified viscosity bump is gone.
- **Off-calibration-Re holdout:** unchanged pass; the term contributed negligibly at `u_frame=0`
  (rate stays `~9.223e-3`, rel err `~3.17e-3 << 2e-2`).
- **D3Q19-vs-D3Q27 holdout:** unchanged pass; both at `u_frame=0`, D3Q19 rel err `~3.2e-3` remains a
  non-outlier vs D3Q27 `~5.1e-3` + `2e-2` band.
- **`cumulant_acceptance` Galilean-improvement test:** D3Q19 CentralMoment defect rises from
  `9.996e-4` toward the BGK level `~2.5e-3` — this is why the claim must be narrowed to D3Q27-only
  (§3.1). D3Q27 CentralMoment defect stays `~1.16e-3 < 0.5*2.473e-3`, still passing.

If the advected holdout spread after removal does **not** fall below `1.157e-3`, the tensorial-defect
diagnosis in §1 is wrong and the finding must be re-opened (this is the built-in falsification test).

---
---

# Round 2 (post-falsification)

Status: PROPOSAL — REVISION. The Round-1 recommendation (path (c), remove the `-0.16 |u|^2` term)
was **experimentally falsified**. This section supersedes §3–§4 above; §0–§2 remain the audited
record of the live code and the qualitative defect picture. The falsifier fired: removing the term
made the advected-TGV3D frame spread **worse**, not better. The built-in falsification test at the
end of §4 therefore mandates re-opening the finding, which this section does.

The falsification data (recorded verbatim in `docs/PHYSICS.md` 2026-07-07 falsification record):

| quantity | value |
|---|---|
| frame spread **WITH** term `(max-min)/mean` | `4.195075506e-3` |
| frame spread **WITHOUT** term | `1.051034711e-2` (2.5x WORSE) |
| per-frame rel. decay-rate error WITH term, `u_frame={0,0.05,0.1}` | `{1.834009427e-3, 2.908688268e-3, 6.044168839e-3}` |
| derived band | `1.156594266e-3` |
| case | `N=32, nu=0.02, u0=0.012, 160 steps, D3Q19 CentralMoment` |

These are the **four measured numbers** any surviving model must reproduce before it may predict a
post-fix spread.

## R2.0 The observable the model must target (read from `cumulant_holdout.rs`)

`measure_decay` computes, over `steps=160`:
`rate = -ln(E1/E0)/steps`, `analytic_rate = 6*nu*k^2`, `k = 2*pi/N`,
`rel_err = |rate - analytic_rate| / analytic_rate`. Since the fluctuation KE of a single
diffusive Fourier mode decays as `exp(-2 * nu_eff * k^2 * t)` and the code's `analytic_rate`
uses `6*nu` (the TGV3D factor for the summed x/y modes), the measured `rate` is proportional to
the **effective viscosity** the operator applies to that mode:

```
nu_eff / nu = rate / analytic_rate  =>  rel_err = | nu_eff/nu - 1 |     (signed: d := nu_eff/nu - 1)
```

So the entire holdout is a measurement of **d(u_frame) = nu_eff/nu - 1** at three frames. The
model below is written in `d`, exactly this observable. `spread = (max rate - min rate)/mean rate`
= `(max d - min d)` to first order (the mean is ~1). The observable is the **frame dependence of
`d`**, nothing else.

## R2.1 WHY the Round-1 prediction failed — the corrected error model

Round-1's error was a **sign/magnitude misattribution**, not (only) a tensoriality claim. It
assumed the scalar term merely nulls a trace it can partially reach, so removing it would expose a
*smaller* isotropic residual and drop below band. The data say the opposite: the term is
**over-compensating a real, negative, `|u|^2`-scaling native defect**, and removing it un-masks the
full native defect, which is *larger* than the net residual the term leaves.

### R2.1.1 Linearized transfer from the `os` factor to `nu_eff`

The term is `os = omega_shear * (1 - C |u|^2)` with `C = 0.16`, and
`nu_eff = (1/os - 1/2) * cs^2`. At the calibration point `omega_shear = 1/(3*0.02+0.5) = 1.7857143`
(`tau = 0.56`), the sensitivity of `nu_eff` to a *fractional* change `epsilon` in `os`
(`os -> omega(1+epsilon)`) is, to first order,

```
d(nu_eff)/nu_eff = -epsilon * (1/os) / (1/os - 1/2)
                 = -epsilon * A ,      A = (1/omega)/(1/omega - 1/2) = 0.56/0.06 = 9.3333.   (R2-1)
```

The amplification `A = 9.333` is the crux Round-1 missed: because `tau = 0.56` is close to `1/2`,
a *small* fractional `os` change is amplified ~9x into `nu_eff`. The term's `epsilon = -C |u|^2`
therefore shifts `d = nu_eff/nu - 1` by `+C |u|^2 * A` (positive: lowering `os` raises `nu_eff`).

### R2.1.2 The spatially-averaged `|u|^2` the term actually sees

The term reads the **local** `|u|^2 = |u_frame + u_TGV(x)|^2` at every cell. For the holdout field
(`u_TGV,x = u0 sinX cosY cosZ`, `u_TGV,y = -u0 cosX sinY cosZ`, x-directed frame `u_frame`):

```
|u|^2 = u_frame^2 + 2 u_frame u_TGV,x(x) + |u_TGV(x)|^2 .
```

- Spatial average `<|u_TGV|^2> = u0^2/4 = 3.6e-5` (each of `<ux^2>,<uy^2> = u0^2/8`).
- **Linear cross term `2 u_frame u_TGV,x`**: its plain average is zero (`sinX` odd). The brief asks
  whether the *dissipation-weighted* correlation of this asymmetric term is what breaks the scalar.
  It is **not**, for this field: the decay rate is weighted by the mode's local energy `~ |u_TGV|^2`,
  which is **even** in X, while `u_TGV,x ~ sinX` is **odd** in X — so
  `<u_TGV,x * |u_TGV|^2> = 0` by parity. The cross term contributes **no linear-in-`u_frame`**
  asymmetry to the measured rate; the frame dependence is purely **even** (`~ u_frame^2`), exactly
  as the monotone measured sequence `{1.83, 2.91, 6.04}e-3` shows. **So the cross-term-asymmetry
  hypothesis in the brief is falsified for this flow; the residual is the genuine even
  `u_frame^2` defect, and tensoriality/orientation is the only remaining structural question (R2.2).**

The effective `|u|^2` the term applies to the mode is thus `u_frame^2 + u0^2/4`.

### R2.1.3 Reconstructing the WITHOUT-term errors — the model reproduces the falsification

Define the **native defect** `d_native(u_frame)` = the frame dependence with the term OFF. By (R2-1)
with the averaged `|u|^2`:

```
d_native(f) = d_with(f) - C * (f^2 + u0^2/4) * A ,   C=0.16, A=9.3333.               (R2-2)
```

Plugging the three measured `d_with = {1.834, 2.909, 6.044}e-3`:

```
f=0.00:  d_native = 1.834e-3 - 0.16*(3.6e-5)*9.333       =  1.780e-3
f=0.05:  d_native = 2.909e-3 - 0.16*(2.536e-3)*9.333     = -0.878e-3
f=0.10:  d_native = 6.044e-3 - 0.16*(1.0036e-2)*9.333    = -8.943e-3
```

Reconstructed rates `= (1 + d_native)*analytic_rate` give a predicted WITHOUT-term
**spread = 1.075e-2**, against the **measured `1.051034711e-2`** — agreement to **2.3%** with a model
carrying *no free parameter* (C, A, and `u0^2/4` are all fixed a priori). **This is the decisive
result: the corrected model reproduces the falsifier.** The native defect frame-growth (fit to the
three `d_native`) is

```
kappa_native = -1.0718  per u_frame^2      (negative: the D3Q19 cubic defect LOWERS nu_eff at frame).
```

The term supplies `+C*A = +0.16*9.3333 = +1.4933` per `u_frame^2`. The net WITH-term frame growth is
`kappa_with = kappa_native + C*A = -1.0718 + 1.4933 = +0.4215` per `u_frame^2`, which reproduces the
measured WITH-term frame residual `d_with - d_with[0] = {0, 1.075e-3, 4.210e-3}` **exactly**
(`0.4215 * {0, 0.0025, 0.01}`). So:

- Round-1's *quantitative* claim `err ~ err_0 + kappa*u_frame^2, kappa ~ 0.42` was numerically right
  (that is `kappa_with`), but it **mis-attributed** the sign of the native defect. The native defect
  is `-1.07 u_frame^2`; the term over-corrects to `+1.49`; the net `+0.42` is *smaller in magnitude
  than the native `-1.07`* — which is precisely why removal (exposing `-1.07`) is **2.5x worse**
  (`1.07/0.42 = 2.55`, matching the measured `1.051e-2/4.195e-3 = 2.51`). Round-1 assumed removal
  would shrink the residual; the arithmetic shows removal enlarges it.

### R2.1.4 What an isotropic scalar CAN and CANNOT do — resolved by the data

The Round-1 tensoriality argument (§1.3) claimed a scalar can only null the **trace** `beta |U|^2`
and never the **anisotropic** `alpha U_a U_b`. The falsification data now let us test this. Model
Eq. 2: `delta nu_eff_{ab}(U) = -cs^2 C_lat (alpha U_a U_b + beta |U|^2 delta_ab)`. For the **single
geometry the holdout uses** (x-directed frame `U=(u_frame,0,0)` acting on a fixed-orientation xy-TGV
shear), the measured decay rate is one linear functional of `(alpha, beta)`:

```
kappa_native = -cs^2 C_lat (alpha * P + 3 beta) ,                                    (R2-3)
```

where `P` is the fixed projection of `U_aU_b = u_frame^2 delta_{ax}delta_{bx}` onto the deviatoric
moments the xy-TGV excites. **A single orientation cannot separate `alpha` from `beta`** — the data
constrain only the combination `(alpha P + 3 beta)`. Consequently:

- A scalar `-C|u|^2` term supplies `+C A` per `u_frame^2` and can be tuned so that
  `C A = -kappa_native`, i.e. it **fully nulls the frame growth for THIS geometry** regardless of
  whether the underlying defect is trace or anisotropic — because at fixed orientation both project
  to the same `u_frame^2`. **This directly refutes Round-1's claim that no scalar can pass the
  holdout.** The holdout, being single-orientation, cannot see the anisotropy Round-1 invoked.
- What a scalar **cannot** do is null the frame growth *simultaneously across orientations* unless
  `alpha = 0` (pure-trace defect). That is the real, still-open tensoriality question — but it is
  **not** what the current holdout measures, and Round-1 conflated "fails a multi-orientation test"
  (plausible, untested) with "fails this holdout" (false).

## R2.2 Re-evaluating the candidate fixes against the corrected model

### (i) Retuned / re-derived scalar coefficient

Setting `C A = -kappa_native` gives the frame-nulling coefficient

```
C_null = -kappa_native / A = 1.0718 / 9.3333 = 0.11484 .                              (R2-4)
```

Predicted holdout with `C = 0.11484`: the three `d(f)` collapse to `{1.819, 1.840, 1.814}e-3`
(the `u0^2/4` rest shift is `+3.9e-5`, negligible), giving **spread = 2.6e-5**, far below the
`1.157e-3` band. **This scalar passes the holdout by construction and reproduces all four numbers**
(with-term spread from `C=0.16`; without-term spread from `C=0` via R2-2; the three per-frame errors;
and predicts the retuned spread).

Is `C_null = 0.1148` **derivable** rather than calibrated? Partially, and this is the honesty crux:

- The *structure* — a negative `|u|^2` correction of `O(0.1)` in `omega`-space with the amplification
  `A` — is derived. `kappa_native = -1.07` is a **measured** frame-growth of the native operator, not
  itself derived from the D3Q19 weight tensor in this document; I attempted the standard 4th-order
  weight-defect route and found `Delta_4 = sum_i w_i c_ix^4 - 3 sum_i w_i c_ix^2 c_iy^2 = 1/3 - 3(1/9)
  = 0` for D3Q19 — i.e. the defect is **not** in the 4th-order weight tensor (Round-1 §1.2's stated
  mechanism is imprecise). The genuine D3Q19 cubic defect lives in the **third-order moment closure**:
  the central basis (`kernels.rs:161-182`, `pow_upto2`) caps every axis exponent at 2, so the
  pure-cubic moment `Q_{aaa}=<c_a^3 f^eq>` is **absent from the moment space entirely**, and the
  retained mixed third moments `xxy,xxz,yyx,yyz,zzx,zzy` (only `xyz` is dropped) carry an equilibrium
  short-fall because the corner velocities are missing. That short-fall, propagated through
  `d_g Q^(0)_{abg}`, is the `O(Ma^2)` frame-dependent `nu_eff` error.
- **Verdict on (i):** `C_null = 0.1148` is derivable *in structure* but its *value* rests on the
  measured `kappa_native`. Under the physics-discipline decision table it is a **calibrated
  constant unless `kappa_native` is derived from the third-moment closure**, AND — critically — it is
  a scalar tuned on a **single orientation**, so shipping it would risk exactly a second falsified
  prediction the moment a multi-orientation holdout is added (if `alpha != 0`). A scalar retune is
  therefore **not** safe to recommend as a physics fix; it is at best a documented interim with a
  hard validity domain "x-aligned / weakly-tensorial frames".

### (ii) Cubic third-moment equilibrium correction (Geier 2015 / Dellar) — CONCRETE terms

Repair `Q^(0)` at the source. The defect is in the equilibrium of the **retained mixed third central
moments**. In the code these are the basis vectors with `order == 3` and max axis-power 2:
`m_{xxy} (2,1,0)`, `m_{xxz} (2,0,1)`, `m_{yyx} (1,2,0)`, `m_{yyz} (0,2,1)`, `m_{zzx} (1,0,2)`,
`m_{zzy} (0,1,2)`. Currently (`kernels.rs:372-377`, `backend_simd.rs` mirror) these relax at
`rate = 1.0` toward the discrete-Hermite equilibrium `eq[m]`. The Galilean-covariant continuum
target for the *raw* third moment is `Q^(0)_{aab} = rho cs^2 u_b` (from Eq. 1 with `a=a, g=b`);
in **central-moment** coordinates the equilibrium mixed third central moment should be **zero**
(a Maxwellian has vanishing central third moments), but the D3Q19 discrete equilibrium leaves a
non-zero residual `kappa`-generating value because the corner populations that would cancel it are
absent. The correction is the standard cumulant/Dellar cubic term: add to the mixed-third
central-moment equilibrium the compensation

```
Delta_eq[m_{aab}] = - kappa_3 * rho * u_a^2 * u_b ,     for each retained mixed triple (a,a,b),   (R2-5)
```

with the coefficient fixed by the D3Q19 velocity set (the missing-corner deficit),

```
kappa_3 = 1 - 3 * ( sum_i w_i c_ia^2 c_ib^2 ) / ( sum_i w_i c_ia^2 )^2   evaluated on D3Q19,
        = the same third-order isotropy deficit that generates kappa_native (its value must be
          computed from the 19-velocity tensor and cross-checked to reproduce kappa_native = -1.0718).
```

Equivalently (Geier 2015, Eq. for the D3Q19 3rd-order cumulant correction), relax the mixed third
cumulants toward the Galilean-corrected target `C_{aab} = 0` at a rate `omega_3` with the added
equilibrium contribution above, rather than toward the plain discrete-Hermite value. Because this
repairs the *equilibrium*, it removes the `alpha U_aU_b` AND `beta |U|^2` parts of Eq. 2 together —
it is orientation-general, unlike (i).

- **Expected residual after correction:** the leading `O(Ma^2)` frame defect is cancelled; the
  remainder is `O(Ma^2 (k dx)^2)` from the finite-mode gradient closure plus `O(Ma^4)`. Estimate
  `|kappa_native| * (k dx)^2 = 1.07 * (2 pi/32)^2 = 1.07 * 0.03855 = 0.0413` per `u_frame^2`, giving a
  predicted spread `~ 0.0413 * 0.01 = 4.1e-4 < 1.157e-3` band. Passes with margin.
- **Implementation sites** (all three must stay bit/threshold-consistent under
  `backend_simd_equiv.rs`): the mixed-third `eq[m]` target and its relaxation rate in
  `kernels.rs:372-377` (the `order` match arm; add the `Delta_eq` of R2-5 to `eq[m]` for the six
  mixed triples, or switch their rate from `1.0` to `omega_3` toward the corrected target),
  the SIMD mirror in `backend_simd.rs:644-652`, and the GPU emitter
  `gpu/wgsl.rs::emit_central_moment_collide` (the `order` arm that currently emits rate `1.0` for
  order>=3). The `solve_moment_system` back-transform (`kernels.rs:204`) is unchanged — no new storage,
  the moments already exist in the D3Q19 basis. D3Q27 has the corner velocities, so its `kappa_3 = 0`
  and it is untouched (verify unchanged).

### (iii) Honest retention with the narrowed claim

Keep `-0.16 |u|^2` as-is; keep the holdout `#[ignore = FINDING]` at its derived band; keep the
narrowed claim (valid for non-advected / weakly-framed decay, Galilean invariance at finite frame
NOT established on D3Q19). This is the current standing state per the falsification record. The
corrected model now *explains* why the term helps at all (it over-compensates a real defect and the
net residual `+0.42 u_frame^2` is smaller than the native `-1.07 u_frame^2`) — so retention is no
longer "a mystery calibration"; it is a documented partial compensation with a known,
model-reproduced residual. Its liability is unchanged: `0.16` is still a calibrated constant and the
holdout still fails its band.

## R2.3 RECOMMENDATION — path (ii), the cubic third-moment equilibrium correction

Decision: **(ii)**. Rationale, decided not surveyed:

- (i) a retuned scalar `C_null = 0.1148` *would* pass this holdout and reproduces all four numbers,
  but it is single-orientation-calibrated and would risk a **second falsified prediction** the moment
  a differently-oriented frame is tested if the defect has any `alpha U_aU_b` content. The brief
  explicitly warns that a second falsified prediction is worse than an honest open finding. Reject as
  a *fix* (acceptable only as an explicitly-bounded interim, which we do not need since retention
  already covers the interim).
- (iii) retention is the safe null, but it leaves a banned calibrated constant live and the holdout
  red. It is the fallback **only if the model below fails its reproduction gate**.
- (ii) is the sole orientation-general, derivation-backed fix. The corrected model gives it a
  concrete coefficient target (`kappa_3` must reproduce `kappa_native = -1.0718`) and a falsifiable
  residual prediction. It attacks the root cause (§1.4) at the equilibrium, not via a relaxation-rate
  proxy.

### R2.3.1 NEW falsifiable prediction (must reproduce ALL FOUR measurements first)

The implementation order's model is accepted **only if**, before predicting any post-fix number, a
one-off diagnostic (term ON with `C=0.16`; term OFF; and the derived `kappa_3`) reproduces:

1. WITH-term (`C=0.16`) frame spread `4.195075506e-3` — tolerance **+-3%** (this is the live code; it
   must match to measurement noise). 
2. WITHOUT-term (`C=0`) frame spread `1.051034711e-2` — tolerance **+-5%** (model R2-2 predicts
   `1.075e-2`, a `+2.3%` deviation, inside tolerance). 
3. The three per-frame WITH-term errors `{1.834, 2.909, 6.044}e-3` — each tolerance **+-3%** (these are
   the model's inputs; reproduction is the sanity check that `A=9.333`, `u0^2/4`, and the fit are
   wired correctly). 
4. The frame-growth decomposition `kappa_native = -1.0718 +-5%` extracted from the term-OFF run, and
   `kappa_native + 0.16*A` reproducing `kappa_with = +0.4215 +-5%`.

**Only after all four pass** does the post-fix prediction stand:

> With the cubic third-moment equilibrium correction (R2-5) applied on the six retained D3Q19 mixed
> third moments, with `kappa_3` derived from the 19-velocity tensor and verified to reproduce
> `kappa_native = -1.0718` (independent check, not fit to the holdout), the advected-TGV3D frame
> spread drops to **`3e-4` to `5e-4`** (center `~4.1e-4`), **below** the `1.156594266e-3` band —
> **PASS**. Rest-frame (`u_frame=0`) rate is essentially unchanged; the off-Re and D3Q19-vs-D3Q27
> holdouts stay passing (both at `u_frame=0`, where the third-moment correction is `~0` for the small
> `u0`). D3Q27's `kappa_3 = 0`, so its holdouts are unchanged.

**Falsifier:** if the derived `kappa_3` does **not** independently reproduce `kappa_native = -1.0718`
within `+-5%`, or if the post-correction spread does **not** fall below `1.157e-3`, then the
third-moment-closure diagnosis is wrong. In that event **do NOT ship a retuned scalar** (that would be
the second falsified prediction). Fall back to **(iii) honest retention** and record the open finding.
If even reproduction gates 1–4 cannot be met, the corrected model itself is wrong: report it as an
open finding and recommend (iii). A second falsified prediction is worse than an honest open finding.

### R2.3.2 Implementation-order skeleton (recommended path (ii))

```
ORDER: D3Q19 cumulant cubic third-moment Galilean correction (ANOM-P4-009 follow-up)
Scope: crates/lbm-core central-moment collision, D3Q19 only (D3Q27 kappa_3=0, must stay unchanged).
Physics-discipline: embed lbmflow-physics-discipline Step 1.5 clauses. kappa_3 MUST be derived from
  the D3Q19 19-velocity tensor with the derivation recorded in the order output; NO value fitted to
  the holdout band. This order adds a physics term => PHYSICS.md entry mandatory.

STEP 0 (model-reproduction gate — BLOCKING, do this FIRST):
  Add a diagnostic run (not a shipped test) that measures, at N=32,nu=0.02,u0=0.012,160 steps,
  frames {0,0.05,0.1}, the D3Q19 CentralMoment decay rates for (a) C=0.16 [live], (b) C=0
  [ablation flag]. Confirm ALL FOUR against docs/proposals R2.3.1:
    (1) WITH spread = 4.195e-3 +-3%   (2) WITHOUT spread = 1.051e-2 +-5%
    (3) per-frame errs {1.834,2.909,6.044}e-3 each +-3%
    (4) kappa_native = -1.0718 +-5% and kappa_native+0.16*A = +0.4215 +-5% (A=9.3333).
  If any gate 1-4 fails: STOP. Report the model as wrong; recommend path (iii) retention. Do not
  proceed to STEP 1. Do not ship a scalar retune.

STEP 1 (derive kappa_3): compute the D3Q19 third-order isotropy deficit from the 19-velocity set
  and the mixed-third central-moment equilibrium; show it reproduces kappa_native=-1.0718 within
  +-5% as an INDEPENDENT check (not a fit). Record the derivation + validity domain in the order log.
  If it does NOT reproduce kappa_native: STOP, recommend (iii), do not ship.

STEP 2 (implement R2-5): add Delta_eq to the six retained mixed-third-moment equilibria (xxy,xxz,
  yyx,yyz,zzx,zzy) OR relax them toward the corrected target at omega_3. Sites, all three consistent:
    - crates/lbm-core/src/kernels.rs:372-377 (order>=3 arm + eq[m] target)
    - crates/lbm-core/src/backend_simd.rs:644-652 (mirror)
    - crates/lbm-core/src/gpu/wgsl.rs emit_central_moment_collide (order>=3 emission)
  Gate the correction on L::Q == 19 (D3Q27 kappa_3=0 -> no-op; assert D3Q27 numbers unchanged).
  Remove the -0.16|u|^2 scalar term and the CENTRAL_MOMENT_DISABLE_VELOCITY_CORRECTION_FOR_ABLATION
  flag ONLY after STEP 3's advected holdout passes (until then keep it, so a failed fix falls back
  cleanly to the documented interim retention).

GATES (must pass WITHOUT retuning any band):
  - cumulant_holdout.rs advected_tgv3d_decay_rate_is_frame_independent: un-ignore ONLY if spread
    <= 1.156594266e-3 (predicted 3e-4..5e-4). If it does not pass, RE-IGNORE, keep the scalar term,
    record failure, recommend (iii). Off-Re and D3Q19-vs-D3Q27 holdouts stay green.
  - cumulant_acceptance.rs: viscosity smoke frozen to uncorrected N=32 value -> the third-moment
    correction is ~0 at u_frame=0,small u0 so expect no smoke regression; Galilean-improvement claim
    re-scoped per §3.1 (D3Q27 genuinely; D3Q19 now via the derived correction if it passes).
  - backend_simd_equiv.rs (bit/threshold) + T13 (partition invariance): MUST stay consistent across
    CPU scalar/SIMD/GPU. This is the primary structural risk of touching the order>=3 arm.
  - t15_3d.rs T15.4 rest-frame TGV: u_frame=0 => correction ~0 => expect no change; confirm.
  - Full lbmflow-build-verify (cargo test --workspace --release, then --include-ignored for the
    cumulant tests) green before reporting done. A codex order finishing is NOT evidence of green.

STOP-RULE (physics-discipline template):
  STOP and report to PM (do not merge, do not fake a band) if ANY of:
    - STEP 0 model-reproduction gate fails (model wrong -> recommend (iii)).
    - STEP 1 kappa_3 does not reproduce kappa_native within +-5% (diagnosis wrong -> recommend (iii)).
    - advected holdout spread does not fall below 1.157e-3 after correction (fix insufficient ->
      re-ignore holdout, retain scalar as documented interim, recommend (iii), open finding).
    - backend_simd_equiv or T13 goes red (structural regression -> revert, re-scope).
    - the only way to pass any band is to adjust kappa_3 away from its derived value (that is a
      banned calibrated constant -> STOP, the spec is revised, the physics is not faked).
  A second falsified prediction is worse than an honest open finding: when in doubt, fall back to
  (iii) and report.

PHYSICS.md entry (implementer writes, not this document): record ANOM-P4-009 as the derived cubic
  third-moment correction superseding the -0.16 scalar; the derivation of kappa_3; the four
  reproduced falsification numbers; the post-fix spread; validity domain (D3Q19 mixed-third moments;
  D3Q27 unaffected); and the narrowed-then-restored Galilean claim.
```

## R2.4 Summary of the corrected model (for the PHYSICS.md cross-reference)

- Observable: `d(u_frame) = nu_eff/nu - 1`, measured by `cumulant_holdout.rs` as `rate/analytic - 1`.
- Amplification: `A = (1/omega)/(1/omega - 1/2) = 9.3333` at `tau=0.56` — the factor Round-1 omitted.
- Term-seen `|u|^2 = u_frame^2 + u0^2/4` (linear cross term vanishes by TGV parity — the brief's
  asymmetry hypothesis is falsified for this field).
- Native defect frame growth `kappa_native = -1.0718` per `u_frame^2` (real, negative `|u|^2` error).
- Term supplies `+0.16 A = +1.4933`; net WITH-term `+0.4215` (reproduces measured residual exactly).
- Model reproduces the without-term spread `1.075e-2` vs measured `1.051e-2` (2.3%) with **no free
  parameter** — the falsification is quantitatively explained.
- Frame-nulling scalar would be `C_null = 0.1148` (single-orientation, rejected as a fix).
- Recommended fix: cubic third-moment equilibrium correction (R2-5), orientation-general, residual
  `~4.1e-4 < 1.157e-3` band, gated on independently reproducing all four falsification numbers.
