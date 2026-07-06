# Multiphase Flow Module Design (Phase 4)

Status: Phase-4 SCMP + MCMP + full-range contact angles all landed on trunk.
MF-alpha (D3Q27 + cumulant central-moment collision) landed 2026-07-06
(commit `20d0e10`, cx/cumulant-s3 — see docs/HANDOFF-PM-2026-07-07.md §3),
lifting the multiphase workload onto the higher-isotropy lattice and the
Galilean-invariant collision. W-VOF proper (Allen-Cahn free surface, MF-γ) is
still pending — owned by the QA-sweep session per handoff §4.

## Approach

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

### Engine changes (lbm-core)

1. **Per-cell force field**: `force_field: Option<Vec<[T; 2]>>` on
   `Simulation`. `F_local = force + force_field[i]` in collide's Guo term and
   in `update_moments`' F/2 correction. Public API:
   `sim.set_force_field(Some(vec))` / `sim.force_field_mut() -> &mut [...]`
   (allocation reuse for callers that rewrite every step).
2. The uniform `force` stays as-is (for gravity).

### multiphase module

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

- Flat-interface coexistence density: G=−5.0, ψ=Classic, 128×64 periodic;
  steady-state density within ±3% of the Maxwell construction.
- Laplace's law: droplets R ∈ {12, 16, 20, 24}, linear fit Δp = σ/R with R² ≥ 0.99.
- Spurious currents: max|u| ≤ 0.05 (G=−5, τ=1).
- Contact angle: sweeping G_w hits the full θ ∈ {~30°, 60°, 90°, 120°, ~150°}
  range within ±10° (spherical-cap fit).

## Phase 4b: Two-component multiphase (MCMP) — landed

- `MultiComponent<T>`: two sets of distribution functions with shared
  solid/geometry.
- Interaction: F_σ = −G_AB ψ_σ(x) Σ w_q ψ_σ̄(x+c_q) c_q (σ̄ = counterpart).
- Composite velocity: standard Shan-Chen common-u for the collision step.
- Validation (T12): RT instability linear growth rate vs theory (Atwood 0.5,
  ±20%); mass conservation across droplet separation/coalescence.

## MF-alpha stage (D3Q27 + cumulant)

The multiphase kernels now run on D3Q27 with the cumulant central-moment
collision on CpuSimd and Wgpu (commit `20d0e10`, cx/cumulant-s3). This lifts
the higher-order isotropy of the lattice into SC/MCMP and removes the
Galilean-invariance error of BGK at strong density gradients. Acceptance
lives under the MF-alpha rows in `docs/paper/claims-ledger.md`.

## Still pending

- **W-VOF proper** (Allen-Cahn free surface, MF-γ): gated on W-GRAV proper
  (already landed). Owned by the QA-sweep session; see
  `docs/HANDOFF-PM-2026-07-07.md` §4.
