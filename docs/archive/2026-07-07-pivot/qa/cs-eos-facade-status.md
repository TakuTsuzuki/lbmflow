# CS-EOS Compat Facade Status

Lifecycle: snapshot, created 2026-07-07. Scope: radar row 37 and
code-to-spec-diff §B rows A16/A17.

## Implementation Status

No Carnahan-Starling (CS-EOS) helper implementation was found in
`crates/lbm-core/src` or `crates/lbm-core/tests`.

Search terms checked: `Carnahan`, `Starling`, `CS-EOS`, `cs_eos`, `cseos`,
`Peng`, `Robinson`, `EOS`, `equation of state`, `pseudopotential`, and `psi`.
The only implemented single-component multiphase pressure helper found in the
compat path is:

```rust
pub fn pressure(&self, rho: T) -> T {
    let cs2 = T::r(CS2);
    let psi = T::r(self.psi.eval(rho.as_f64()));
    cs2 * rho + T::r(0.5) * self.g * cs2 * psi * psi
}
```

This is the Shan-Chen pseudopotential bulk pressure,
`p = cs^2 rho + (G cs^2 / 2) psi(rho)^2`, not a Carnahan-Starling EOS.

The compat `Psi` enum has exactly two variants:

```rust
pub enum Psi {
    Classic,
    Exponential { psi0: f64, rho0: f64 },
}
```

`Psi::Exponential` is the Shan-Chen 1994 pseudopotential form. It is a
non-default Shan-Chen pseudopotential, but it is not a named CS-EOS helper and
does not expose a separate tunable thermodynamic EOS.

The native V2 `Solver` has:

```rust
pub fn update_shan_chen_force(&mut self, g: T, psi: impl Fn(T) -> T)
pub fn update_shan_chen_force_with_walls(
    &mut self,
    g: T,
    g_wall: T,
    psi_wall: T,
    psi: impl Fn(T) -> T,
)
```

That native API accepts a caller-provided Shan-Chen pseudopotential closure, but
there is no native named CS-EOS helper to re-export through the compat facade.

## Facade Exposure Status

`crates/lbm-core/src/compat/mod.rs` exposes `pub mod multiphase`, so users can
reach `lbm_core::compat::multiphase::{Psi, ShanChen, MultiComponent}`.

The convenient compat prelude is inline in `compat/mod.rs`; there is no
`crates/lbm-core/src/compat/prelude.rs` file in this tree. The prelude exports
only domain/config types, `Real`, and `Simulation`:

```rust
pub mod prelude {
    pub use super::domain::{Collision, ConfigError, Edge, EdgeBC, Edges, SimConfig, MAX_SPEED};
    pub use super::real::Real;
    pub use super::sim::Simulation;
}
```

Therefore:

- CS-EOS is not exposed through `compat::multiphase`.
- CS-EOS is not exposed through `compat::prelude`.
- No native CS-EOS helper exists elsewhere in `lbm-core` for a small compat
  re-export.
- The only non-default EOS-adjacent compat knob is `Psi::Exponential`, which
  remains within the Shan-Chen pseudopotential model.

## Gap Classification

Classification: missing capability plus stale documentation claim.

This is not the A16/A17 "native helper exists but compat forgot to re-export it"
case. The search found no native Carnahan-Starling implementation at all. The
gap is that PHYSICS.md currently claims "Carnahan-Starling (CS) EOS helpers"
and "Density-ratio ceiling from CS-EOS", while the code exposes only classic
and exponential Shan-Chen pseudopotentials.

Radar row 37 is directionally correct that the compat facade lacks CS-EOS, but
the sharper classification is "missing capability", not merely "missing facade
re-export".

## Recommended Action

Recommended near-term action: remove or condition the CS-EOS claim in
`docs/PHYSICS.md` and update row 37 / code-to-spec A16/A17 wording to say that
no CS-EOS helper exists anywhere in `lbm-core` as of 2026-07-07.

Recommended implementation action, if CS-EOS remains desired: add a real
capability with derivation and validation before advertising it. The minimum
acceptable work is:

1. Define the Carnahan-Starling pressure relation and the derived
   pseudopotential mapping used by the Shan-Chen force, including parameter
   validity domains.
2. Add a typed API, for example a new `Psi::CarnahanStarling { ... }` or a
   separate EOS helper type, only after the derivation is recorded in
   `docs/PHYSICS.md`.
3. Add T11-style validation covering coexistence, pressure equilibrium, Laplace
   behavior, stability envelope, and at least one holdout that is not just a
   regression pin.
4. Then expose the helper consistently through the native path and
   `compat::multiphase`; decide separately whether multiphase belongs in the
   convenient prelude.

No facade-only diff is proposed. A small re-export cannot fix this because
there is no implementation to re-export.
