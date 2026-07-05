//! V1 API facade over the V2 core — a supported public API, not a shim to
//! be removed: the scenario / CLI / wasm 2D paths and the inherited V1 test
//! suite all run through it.
//!
//! Module tree, type names and behaviour mirror the retired V1 engine
//! (`crates/lbm-core` until 2026-07-05; V1 client code ported with the
//! textual substitution `lbm_core::` → `lbm_core::compat::`). The physics
//! is the V2 `Solver<D2Q9, T, CpuScalar, LocalPeriodic>`, which reproduced
//! V1 trajectories bit-exactly — proven by `tests/v1_match.rs`, whose final
//! frozen measurements live in its header in branch history (deleted with
//! V1).
//!
//! R5 (COMPETITIVE_SPEC.md): the pre-existing 2D suite must stay green,
//! unmodified, on top of the new core. This module is that guarantee.

pub mod domain;
pub mod lattice;
pub mod multiphase;
pub mod real;
pub mod sim;

/// Convenient glob import: `use lbm_core::compat::prelude::*;`.
pub mod prelude {
    pub use super::domain::{Collision, ConfigError, Edge, EdgeBC, Edges, SimConfig, MAX_SPEED};
    pub use super::real::Real;
    pub use super::sim::Simulation;
}
