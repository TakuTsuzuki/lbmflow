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
