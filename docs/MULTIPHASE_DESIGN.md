# Multiphase Flow Module Design (Phase 4)

**Status (2026-07-07)**: **Landed**: Shan-Chen SCMP, MCMP, and contact-angle validation exist on the 2D compat facade path.
**Scope**: scenario/CLI multiphase builds `compat::Simulation<T>` plus optional `compat::multiphase::ShanChen<T>`; WASM exposes f32 SCMP through the same facade.
**Superseded**: claims that multiphase kernels run on D3Q27/CpuSimd/Wgpu; GPU scenario dispatch rejects multiphase and D3Q27/cumulant are core primitives, not the landed multiphase path.
**Current intent**: Allen-Cahn W-VOF / MF-gamma remains pending outside this Phase-4 Shan-Chen design.

## Approach

(landed 2026-07-07 — Shan-Chen is implemented in `crates/lbm-core/src/compat/multiphase.rs`; product path is 2D compat/D2Q9 CPU, with WASM f32 SCMP)

Adopt the Shan-Chen pseudopotential method. Rationale:
- Implementation is **orthogonal** to collision/streaming (it injects a force
  field), so it does not perturb the single-phase core.
- Interface tracking is unnecessary (diffuse interface); strong for demos.
- Naturally extends to both single-component multiphase (droplets, bubbles,
  condensation) and two-component multiphase (RT instability, two-phase flow).

Constraints (documented weaknesses):
- Density ratio ~50 with the classic ψ, ~100-1000 with CS-EOS.
- Spurious currents near the interface at O(1e-2).
- Surface tension and density ratio couple through G (independent control
  requires a multi-range ψ — future work).

## Phase 4a: Single-component multiphase (SCMP) — landed

(landed 2026-07-07 — `ShanChen<T>` updates a compat per-cell force field before `Simulation::step`)

### Engine changes (lbm-core)

(landed with changed scope 2026-07-07 — public compat API is `force_field_mut`/`clear_force_field`; V2 core stores a 3-component `force_field`)

1. **Per-cell force field**: `force_field: Option<Vec<[T; 2]>>` on
   `Simulation`. `F_local = force + force_field[i]` in collide's Guo term and
   in `update_moments`' F/2 correction. Public API:
   `sim.set_force_field(Some(vec))` / `sim.force_field_mut() -> &mut [...]`
   (allocation reuse for callers that rewrite every step).
2. The uniform `force` stays as-is (for gravity).

### multiphase module

(landed with API drift 2026-07-07 — actual `Psi::Exponential` also carries `psi0`, and `ShanChen` includes `wall_rho`/`with_wall_rho`)

```rust
pub struct ShanChen<T> {
    pub g: T,                 // fluid-fluid interaction strength (negative = attractive)
    pub g_wall: T,            // wall adhesion strength (contact-angle control)
    pub psi: Psi,             // choice of potential function
}
pub enum Psi { Classic /* 1 - exp(-rho) */, Exponential { rho0: f64 } }

impl ShanChen<T> {
    /// Computes the SC force field from the current rho field and sets it on sim.
    /// Usage: loop { sc.update_force(&mut sim); sim.step(); }
    pub fn update_force(&self, sim: &mut Simulation<T>);
}
```

- F(x) = −G ψ(x) Σ_q w_q ψ(x+c_q) c_q (fluid neighbors)
- Wall: F_ads(x) = −G_w ψ(x) Σ_q w_q s(x+c_q) c_q (s=1 if solid)
- Periodic boundaries wrap; open-boundary cells use zero-gradient
  extrapolation. SC combined with open boundaries is unsupported.

### Validation (T11)

(landed with frozen bands 2026-07-07 — T11/T11b/T11c live in `validation_multiphase.rs` and `validation_contact_angle.rs`; old contact-angle target range below is historical)

- Flat-interface coexistence density: G=−5.0, ψ=Classic, 128×64 periodic;
  steady-state density within ±3% of the Maxwell construction.
- Laplace's law: droplets R ∈ {12, 16, 20, 24}, linear fit Δp = σ/R with R² ≥ 0.99.
- Spurious currents: max|u| ≤ 0.05 (G=−5, τ=1).
- Contact angle: sweeping G_w hits the full θ ∈ {~30°, 60°, 90°, 120°, ~150°}
  range within ±10° (spherical-cap fit).

## Phase 4b: Two-component multiphase (MCMP) — landed

(landed 2026-07-07 — `MultiComponent<T>` exists in compat and T12 validates separation plus Rayleigh-Taylor growth)

- `MultiComponent<T>`: two sets of distribution functions with shared
  solid/geometry.
- Interaction: F_σ = −G_AB ψ_σ(x) Σ w_q ψ_σ̄(x+c_q) c_q (σ̄ = counterpart).
- Composite velocity: standard Shan-Chen common-u for the collision step.
- Validation (T12): RT instability linear growth rate vs theory (Atwood 0.5,
  ±20%); mass conservation across droplet separation/coalescence.

## MF-alpha stage (D3Q27 + cumulant)

(superseded 2026-07-07 — D3Q27 and `CollisionKind::Cumulant` landed in the V2 core, but Shan-Chen multiphase remains on the 2D compat/D2Q9 path and GPU dispatch rejects multiphase)

The multiphase kernels now run on D3Q27 with the cumulant central-moment
collision on CpuSimd and Wgpu (commit `20d0e10`, cx/cumulant-s3). This lifts
the higher-order isotropy of the lattice into SC/MCMP and reduces the
velocity-dependent truncation error of BGK at strong density gradients.
(Correction 2026-07-07: "removes the Galilean-invariance error" is NOT
established — the D3Q19 cumulant correction has an open finite-frame
Galilean-invariance holdout finding; see the PHYSICS.md 2026-07-07 entry.)
Acceptance lives under the MF-alpha rows in `docs/paper/claims-ledger.md`.

## Still pending

(current intent 2026-07-07 — not implemented in this Shan-Chen module)

- **W-VOF proper** (Allen-Cahn free surface, MF-γ): gated on W-GRAV proper
  (already landed). Owned by the QA-sweep session; see
  `docs/HANDOFF-PM-2026-07-07.md` §4.
