# Multiphase Flow Module Design (Phase 4)

## Approach

Adopt the Shan-Chen pseudopotential method. Rationale:
- The implementation is **orthogonal** to collision/streaming (it merely injects a force field),
  so it does not break the already-validated single-phase core
- Interface tracking is unnecessary (diffuse interface), which makes it strong for
  beginner-oriented demos
- It naturally extends to both single-component multiphase (droplets, bubbles, condensation)
  and two-component multiphase (RT instability, two-phase flow)

Constraints (known weaknesses; must be documented explicitly):
- The density ratio is ~50 with the classic ψ, ~100-1000 once CS-EOS is introduced (Phase 4 starts from the former)
- Spurious currents appear near the interface at O(1e-2)
- Surface tension and density ratio are coupled through G (independent control requires a
  multi-range ψ → future work)

## Phase 4a: Single-component multiphase (SCMP)

### Engine changes (lbm-core)

1. **Per-cell force field**: add `force_field: Option<Vec<[T; 2]>>` to `Simulation`.
   - `F_local = force + force_field[i]` in collide's Guo term / update_moments' F/2 correction
   - Public API: `sim.set_force_field(Some(vec))` / `sim.force_field_mut() -> &mut [...]`
     (allocation reuse for callers who rewrite it every step)
2. The existing uniform `force` stays as-is (for gravity).

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
- Periodic boundaries wrap; open-boundary cells use zero-gradient extrapolation (in Phase 4,
  it is acceptable for SC combined with open boundaries to be unsupported — specify this explicitly)

### Validation (T11 concretization)

- Flat-interface coexistence density: G=−5.0, ψ=Classic, 128×64 periodic, top/bottom halves
  initialized with ρ_l/ρ_v → steady-state density within ±3% of the theoretical Maxwell
  construction (reference value computed separately via numerical integration)
- Laplace's law: droplets with R ∈ {12, 16, 20, 24}, linear fit of Δp = σ/R with R² ≥ 0.99
- Spurious currents: max|u| ≤ 0.05 (G=−5, τ=1)
- Contact angle: sweeping G_w to hit θ ∈ {~60°, 90°, ~120°} within ±10° (measured via
  spherical-cap fit from droplet height/radius)

## Phase 4b: Two-component multiphase (MCMP)

- `MultiComponent<T>`: two sets of distribution functions (to be decided at implementation time —
  either holding two internal Simulation instances, or a dedicated struct holding f as [2][N*Q].
  Shared solid/geometry is mandatory)
- Interaction: F_σ = −G_AB ψ_σ(x) Σ w_q ψ_σ̄(x+c_q) c_q (σ̄ is the counterpart component)
- Common velocity u' = (Σ_σ m_σ ω_σ + ...)/... (collision uses the standard Shan-Chen composite velocity)
- Validation (T12): RT instability linear growth rate vs theory (Atwood number 0.5, ±20%),
  mass conservation for droplet separation/coalescence

## Implementation order

1. Force-field API (+ unit test: uniform force_field ≡ uniform force agreement)
2. SCMP + flat interface/Laplace (codex validation commissioned)
3. Contact angle
4. MCMP + RT (codex validation commissioned)
