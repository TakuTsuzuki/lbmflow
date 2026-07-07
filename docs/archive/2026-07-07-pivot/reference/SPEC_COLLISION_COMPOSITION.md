# SPEC — Collision-composition trait boundary (B-1 / B-8 enabler, MF-α)

> Research-agent deliverable for the R-Phase 2 B-1/B-8 trait design (the MF-α enabler:
> D3Q27 + central-moment/cumulant collision). Ready-to-dispatch spec section, same depth /
> format as SPEC_UNIT_CONVERTER.md. Derived from firsthand reads of the three live collision
> sites — `kernels.rs::collide_row`, `backend_simd.rs::{collide_span_flat,collide_span_blocked}`,
> `gpu/wgsl.rs::generate` — plus OpenLB `Tuple<MOMENTA,EQ,COLLISION,FORCING>` and Palabos
> comprehensive moment-space templates (both re-derived from published methods; no GPL/AGPL code).
>
> **The design goal in one sentence:** define the collision operator **once** so that adding MRT,
> cumulant, regularized, LES-effective-omega, or non-Newtonian collision is one new zero-sized
> type that CPU-scalar, CPU-SIMD (flat + blocked), and the WGSL generator all pick up — with **no
> per-operator GPU shader hand-authoring** and **no regression** in the backend-equivalence gates.

---

## 1. Scope & the problem being solved

The TRT+Guo collision math exists today in **three hand-synchronised copies**, kept
operand-for-operand identical by the equivalence gates:

| Site | Form | File / fn |
|---|---|---|
| Scalar reference | per-cell, per-q TRT pair | `kernels.rs::collide_row` (L127–212) |
| Fused CPU-SIMD (D2Q9) | flat, pair-shared `base/r3/r45` regrouping | `backend_simd.rs::collide_span_flat` (L263–358) |
| Fused CPU-SIMD (D3Q19) | block-staged (`[T;BLOCK]`) for vectorization | `backend_simd.rs::collide_span_blocked` (L365–473) |
| GPU | generated WGSL SSA (`e{q}`, `s{q}`, `fc{a}`…) | `gpu/wgsl.rs::generate` step section (L316–362) |

These are held equal by `tests/backend_simd_equiv.rs` (CpuScalar↔CpuSimd, ≤1e-11 f64 / ≤1e-6 f32)
and T14 (CPU↔GPU, ≤1e-5). The SIMD form is a deliberate algebraic regrouping of `collide_row`
that "differs from it by last-ulp rounding only" (backend_simd.rs L196–199, L60–66).

**The problem:** MF-α needs MRT + central-moment + cumulant; MF-β needs LES; the roadmap needs
non-Newtonian and per-cell omega. Naïvely, each new operator would be written **3–4 more times**
(scalar, flat, blocked, WGSL) and the GPU one would be **hand-authored per operator** — exactly the
fork OpenLB (CSE codegen) and Palabos (separate `acceleratedLattice`) suffered, which our
COMPETITOR_ANALYSIS §D said NOT to adopt. This spec defines the boundary that collapses the copies
to **one operator definition + N interpreters**, so operator count and backend count multiply
instead of the operator×backend matrix.

---

## 2. The composition model (MOMENTA × EQ × COLLISION × FORCING)

Mirror OpenLB's four-slot `Tuple`, as zero-sized types (ZSTs) composed at the type level and
monomorphised — no runtime dispatch in the hot loop. The four slots and who owns what:

| Slot | Owns | Today = | Future variants |
|---|---|---|---|
| **MOMENTA** | compute `(rho, u)` from `f` for streaming/output, incl. the Guo **F/2 velocity correction** | `moments_row` / `moments_span_*` | unchanged for hydrodynamics; cumulant/CM still expose `(rho,u)`, higher moments stay internal to COLLISION |
| **EQUILIBRIUM** | `feq_q` (deviation form `feq−w`) | `equilibrium` (kernels.rs L101–119) | incompressible-`feq`, cumulant moment-space `feq`, complete-`feq` |
| **COLLISION** | the relaxation: read `f`, `feq`, `source`, per-cell `omega`, write post-collide `f` | inline TRT pair split (`collide_row` L189–210) | BGK (single ω), TRT (pair ω±), MRT (per-moment `s_i`), CM/cumulant (moment transform → relax → back), regularized (project onto Π then relax) |
| **FORCING** | `source_q` + the velocity shift MOMENTA applies | Guo (`src[q]`, `cp/cm` prefactors) | Shan-Chen shift, He, plain-Guo, Liang family (phase-field) |

Composition rule: **COLLISION is the outer driver**; it calls `EQUILIBRIUM::feq` and
`FORCING::source` as sub-operations. A concrete dynamics is a type alias, OpenLB-style:

```rust
// today's default, expressed in the new model:
type TrtGuo = Collide<Trt, SecondOrderEq, GuoForcing>;
// MF-α targets:
type CumulantGuo = Collide<Cumulant, MomentSpaceEq, GuoForcing>;   // needs D3Q27
type MrtGuo      = Collide<Mrt,      SecondOrderEq, GuoForcing>;
// MF-β / roadmap, as *wrappers* (decorators) around an inner collision:
type SmagorinskyTrt = Collide<EffectiveOmega<Smagorinsky, Trt>, SecondOrderEq, GuoForcing>;
type CarreauTrt     = Collide<EffectiveOmega<PowerLaw,    Trt>, SecondOrderEq, GuoForcing>;
```

`EffectiveOmega<Model, Inner>` is the OpenLB `SmagorinskyEffectiveOmega<COLLISION>` / COMPETITOR
§B5 trick: compute `omega_eff` from the non-equilibrium stress `Π_neq = Σ c c (f − feq)` **already
available inside the collision** (zero extra field/pass for plain Smagorinsky; non-Newtonian reads
the same `Π_neq` → strain rate → `μ(γ̇)` → omega), then delegate to `Inner` with `omega_eff`.

---

## 3. The abstract-arithmetic boundary (the one mechanism that prevents the GPU fork)

**Key decision.** Write each operator **once** against an abstract numeric sink, `trait Arith`, and
provide three interpreters. The operator never sees `f32`/`f64`/WGSL directly — it composes values
through `Arith`, so the **operand order is defined exactly once, in the operator**, and every
interpreter reproduces it by construction. This is what keeps CPU and GPU bit-close *and* lets the
WGSL generator emit any operator without hand-authoring.

```rust
/// Abstract arithmetic over a value handle `V`. Lane-parametric: one impl does
/// scalar CPU, one does block-staged SIMD, one emits WGSL SSA text.
trait Arith {
    type V: Copy;                       // CPU: [T; LANES];  WGSL: an SSA id
    fn lit(&mut self, x: f64) -> Self::V;      // numeric constant (WGSL: `lit()` f32 literal)
    fn add(&mut self, a: Self::V, b: Self::V) -> Self::V;
    fn sub(&mut self, a: Self::V, b: Self::V) -> Self::V;
    fn mul(&mut self, a: Self::V, b: Self::V) -> Self::V;
    fn recip(&mut self, a: Self::V) -> Self::V;         // moments' 1/rho; cumulant normalise
    // context accessors (interpreter decides load vs recompute):
    fn pop(&mut self, q: usize) -> Self::V;             // f_q in deviation form
    fn set_pop(&mut self, q: usize, v: Self::V);        // post-collide write
    fn rho(&mut self) -> Self::V;
    fn u(&mut self, d: usize) -> Self::V;
    fn force(&mut self, d: usize) -> Self::V;
    fn omega_p(&mut self) -> Self::V;                   // may be per-cell (B-6)
    fn omega_m(&mut self) -> Self::V;
    fn w(&self, q: usize) -> f64;                       // lattice weight (compile-time)
    fn c(&self, q: usize, d: usize) -> i8;              // lattice velocity (compile-time)
}
```

Three interpreters — the collision matrix collapses to these three, forever:

| Interpreter | `V` | Purpose | Replaces |
|---|---|---|---|
| `ScalarArith<T, 1>` | `[T;1]` | CPU scalar reference + D2Q9 flat (LLVM vectorises the outer cell loop) | `collide_row`, `collide_span_flat` inline math |
| `BlockArith<T, BLOCK>` | `[T;BLOCK]` | D3Q19 block-staged; lane-wise ops keep the per-pair sweeps vectorising | `collide_span_blocked` inline math |
| `WgslArith` | SSA id (`Val(u32)`) | append `let vN = a + b;` lines; `pop(q)`→`f{q}`, `set_pop`→`f_out`/register | `wgsl.rs` hand-written `step` collision block (L316–362) |

**Why lane-parametric `V` matters (performance, not just cleanliness):** the D3Q19 blocked kernel
(`collide_span_blocked`) stages `base/r3/r45` across `[T;BLOCK]` and sweeps once per TRT pair
specifically because a flat per-cell form left D3Q19 "essentially scalar: 18 vs 285 vector
instructions and 22 vs 43 MLUPS at 128³" (backend_simd.rs L170–178). If the operator were written
against a scalar `V`, re-deriving the block form would re-fork it. With `V = [T; BLOCK]`, every
`add`/`mul` is a lane-wise loop the compiler vectorises, and shared-subexpression hoisting inside
the operator (compute `base/r3/r45` before the pair loop) stages them across the block **exactly
like the current code** — one operator, both the flat and blocked performance profiles preserved.

**Why this is the anti-fork mechanism for GPU:** after the refactor, `wgsl::generate` becomes
`generate::<L, C: Collision>()`; its collision section is produced by running
`C::relax::<WgslArith>(...)`. The streaming push, edge-stash, BC passes, and `emit_cell_prologue`
(wgsl.rs L185–225) are operator-independent and stay hand-written. **Adding MRT/cumulant/LES = one
`Collision` impl; the shader for it is generated, never authored.** The existing
`generated_wgsl_parses_and_validates_with_naga` test (wgsl.rs L556) auto-covers each new operator's
shader; T14 covers its numerics.

---

## 4. Trait definitions (Rust sketch — labelled as sketch, not final)

```rust
/// EQUILIBRIUM slot: feq_q in deviation form (feq − w), operand order = kernels::equilibrium.
trait Equilibrium { fn feq<A: Arith>(a: &mut A, q: usize) -> A::V; }

/// FORCING slot: source_q added after relaxation; also names the velocity shift MOMENTA applies.
trait Forcing {
    fn source<A: Arith>(a: &mut A, q: usize) -> A::V;       // Guo: w_q(3(cf−uf)+9 cu cf)
    const HALF_FORCE_VELOCITY_SHIFT: bool;                   // Guo/He = true → MOMENTA adds F/2
}

/// COLLISION slot: the relaxation. Reads pop/feq/source/omega via `a`, writes set_pop.
trait Collision {
    type Eq: Equilibrium;
    type Force: Forcing;
    fn relax<A: Arith, L: Lattice>(a: &mut A);              // the whole per-cell (or per-block) update
}
```

How each concern composes through this boundary:

- **BGK / TRT** — `relax` loops `L::PAIRS` (compile-time table already used by all three sites),
  builds `feq` from `Self::Eq`, `source` from `Self::Force`, splits into symmetric/antisymmetric
  parts and relaxes with `omega_p`/`omega_m`. BGK = `Trt` with `omega_m == omega_p` (as today,
  params.rs `omegas()` L30–41). **Transcribe the existing operand order** into `Arith` calls so the
  gates stay green (see §6).
- **TRT magic (COMPETITOR Palabos B6 / OpenLB tunable)** — already a field of
  `CollisionKind::Trt{magic}` (params.rs L18–23) feeding `omegas()`. No trait change; just surface
  `magic` in the scenario JSON (default `MAGIC_STD = 3/16`, params.rs L27). The operator only ever
  sees `omega_p/omega_m`.
- **Guo forcing** — `GuoForcing: Forcing`. `source` = the `src[q]` expression (kernels.rs L184);
  `HALF_FORCE_VELOCITY_SHIFT = true` so MOMENTA keeps the `(m + F/2)/rho` correction (kernels.rs
  L348, invariant in CLAUDE.md). The `cp = 1−ω_p/2`, `cm = 1−ω_m/2` prefactors (params.rs L127–128)
  stay as the FORCING↔COLLISION contract.
- **Per-cell omega (B-6 / LES / non-Newtonian)** — `a.omega_p()` returns a **value**, not a global.
  For uniform viscosity it's a splat of the global (today). For LES/non-Newtonian, the
  `EffectiveOmega<Model, Inner>` wrapper computes `omega_eff` from `Π_neq = Σ_q c_q c_q (f_q −
  feq_q)` (the operator already has `f` and `feq` in hand — zero extra traffic, COMPETITOR §B5),
  then calls `Inner::relax` with `omega_eff` substituted. Storage for an explicit per-cell omega
  field is optional and reuses the existing per-cell `force_field` storage pattern
  (`SoaFields`, backend_simd.rs L1344 `force_field`).
- **Future cumulant slot (MF-α)** — `Cumulant: Collision` with `Eq = MomentSpaceEq`. `relax`
  transforms `f` to the cumulant basis (Chimera transform, D3Q27), relaxes each cumulant toward its
  target, transforms back — all through `Arith` ops (`add/sub/mul/recip` suffice; enumerate any
  extra primitive when the transform is written). Requires the planned `Q_MAX = 27` bump; MOMENTA
  still returns `(rho,u)`. This is precisely the OpenLB comprehensive-template / Palabos
  `comprehensiveModelsTemplates` payoff: once the transform machinery exists behind `Arith`, RM / HM
  / CM / CHM / K / RR are near-free additional `Collision` impls.

---

## 5. Integration into the existing code (what each site becomes)

- `kernels::collide_row` → thin caller of `TrtGuo::relax::<ScalarArith<T,1>, L>` per cell. Keeps the
  scalar reference semantics; **removes one hand-synced copy**.
- `backend_simd::collide_span_flat` → caller of `relax::<ScalarArith<T,1>, L>` over the span (LLVM
  vectorises the cell loop, as it does today for D2Q9).
- `backend_simd::collide_span_blocked` → caller of `relax::<BlockArith<T,BLOCK>, L>` over each block
  (lane-wise ops preserve the per-pair vectorised sweeps). `use_blocked::<L>()` (L183) still selects
  flat vs blocked — now just selects the interpreter, not a re-implemented kernel.
- `gpu/wgsl::generate` → `generate::<L, C: Collision>()`; the collision block (L316–362) is emitted
  by `C::relax::<WgslArith>`; prologue/push/stash/BC unchanged.
- `StepParams`/`KParams` (params.rs) gain the operator selection (a `CollisionKind`-style enum at the
  boundary) and, when present, the per-cell omega source. The `moments_*` functions become the
  MOMENTA slot (they already are, essentially).

**Net:** the collision×backend matrix (currently 4 sites for 1 operator) becomes **1 operator
definition + 3 fixed interpreters**. Operator N and backend M stop multiplying.

---

## 6. Invariants & bit-equivalence guardrails (non-negotiable)

- **This touches the hottest code and the strictest gates.** The migration MUST be a *mechanical,
  operand-preserving* transcription first: port **only** the existing TRT+Guo into the `Arith` form,
  transcribing expressions in their current order, and prove **all three gates stay green** BEFORE
  adding any new operator:
  - `tests/backend_simd_equiv.rs` (CpuScalar↔CpuSimd ≤1e-11 f64 / ≤1e-6 f32),
  - T14 backend equivalence (CPU↔GPU ≤1e-5),
  - `collide_feq_matches_equilibrium_bitwise_*` (kernels.rs L770–780 — the equilibrium fixed-point
    bit test) and the SIMD `v1_match`-class reassociation tolerance (backend_simd.rs L60–66).
- The flat/blocked/WGSL forms are already *different algebraic regroupings* accepted at last-ulp;
  the `Arith` boundary does not tighten or loosen that — each interpreter keeps its current
  regrouping (flat = pair-shared `base/r3/r45`; blocked = block-staged; WGSL = the generated SSA).
  The operator defines the *sequence*; the interpreter defines the *rounding*. Keep the SIMD
  pair-shared form as the canonical operator body (it is the faster one — 219 vs 182 MLUPS,
  backend_simd.rs L61–64 — and CpuScalar already tolerates the drift).
- Deviation storage (`f−w`), q-major SoA `f[q*np+i]`, one-step pass order (collide→halo→stream→open
  BC→moments), and the D2Q9 direction ordering (CLAUDE.md invariants) are untouched — the boundary
  is strictly *inside* the collide phase.
- Determinism (R4): all three interpreters are pure functions of their inputs; no interpreter may
  introduce data-dependent reordering (WGSL probe atomics remain the only tolerated nondeterminism,
  unchanged, wgsl.rs L300).

---

## 7. Operator cost (after boundary)

Each new operator = one ZST, live on CPU-scalar, CPU-SIMD (flat+blocked), and GPU
without further backend porting. Concrete targets: MRT (S/M), cumulant D3Q27 (M),
Smagorinsky/WALE LES (S/S–M), non-Newtonian Carreau/Herschel-Bulkley (S), RLB (S/M).

---

## 8. Acceptance / adversarial test matrix (for codex, written from §3/§6)

1. **Operand-preserving port (regression):** after TRT+Guo is expressed via `Arith`, `backend_simd_equiv`,
   T14, and the equilibrium-fixed-point bit test pass with **unchanged tolerances** — no gate is
   relaxed to accommodate the refactor.
2. **Interpreter agreement:** for random valid `(rho,u,force,ω±)`, `ScalarArith<_,1>` and
   `BlockArith<_,BLOCK>` produce results within the existing SIMD tolerance; `WgslArith` output,
   run on GPU, agrees within T14.
3. **Cross-lane invariance:** `BlockArith` with `BLOCK` vs `BLOCK/2` gives bitwise-identical
   per-cell results (lane count must not change rounding).
4. **New-operator drop-in (MRT):** adding `MrtGuo` requires no change to scalar/flat/blocked/WGSL
   plumbing; its generated WGSL passes `naga` validation (wgsl.rs L556 test) and a TGV/Poiseuille
   MRT run matches its analytic reference in the existing validation bands.
5. **EffectiveOmega correctness:** `EffectiveOmega<Const,Trt>` (a no-op model returning the global
   omega) is **bit-identical** to plain `TrtGuo` — proving the wrapper adds no arithmetic when the
   model is trivial.
6. **Per-cell omega path:** a scenario with a uniform per-cell omega field equals the global-omega
   run to the SIMD tolerance (the field path and the scalar path agree).
7. **feq/source single-source:** `Equilibrium::feq` and `Forcing::source` used by COLLISION are the
   same expressions the MOMENTA/GPU-prologue reads — extend the existing
   `collide_feq_matches_equilibrium_bitwise` property to the trait form.
8. **No-forcing fast path:** `Forcing = NoForce` elides all source terms (matches the `force_on ==
   false` branch in collide_span_flat L350–355 and wgsl push), verified equal to today's forceless
   run bitwise.
