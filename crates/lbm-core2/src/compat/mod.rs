//! V1 API facade (`lbm_core` drop-in) over the V2 core.
//!
//! Module tree, type names and behaviour mirror V1 exactly, so client code
//! ports with the textual substitution `lbm_core::` → `lbm_core2::compat::`
//! (`scripts/sync-tests.sh` applies it to the V1 test suite). The physics is
//! the V2 `Solver<D2Q9, T, CpuScalar, LocalPeriodic>`, which reproduces V1
//! trajectories bit-exactly (see `tests/v1_match.rs`).
//!
//! R5 (COMPETITIVE_SPEC.md): the pre-existing 2D suite must stay green,
//! unmodified, on top of the new core. This module is that guarantee.

pub mod domain;
pub mod lattice;
pub mod multiphase;
pub mod real;
pub mod sim;

/// Convenient glob import: `use lbm_core2::compat::prelude::*;`.
pub mod prelude {
    pub use super::domain::{Collision, ConfigError, Edge, EdgeBC, Edges, SimConfig, MAX_SPEED};
    pub use super::real::Real;
    pub use super::sim::Simulation;
}
