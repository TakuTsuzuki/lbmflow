//! # lbm-core2
//!
//! Lattice Boltzmann core V2: the dimension / lattice / precision / backend /
//! decomposition axes are orthogonal (docs/ARCHITECTURE_V2.md). The physics
//! kernels are written once, generically over a [`lattice::Lattice`] and a
//! [`real::Real`], and specialise at compile time.
//!
//! Layer map (docs/ARCHITECTURE_V2.md §1):
//!
//! - [`lattice`] — compile-time velocity sets (D2Q9, D3Q19) with derived
//!   tables (TRT pairs, per-face unknown sets).
//! - [`fields`] — q-major SoA deviation storage over halo-padded local boxes.
//! - [`kernels`] — the physics (collide/stream/moments/BCs), written once,
//!   generic over lattice and precision; V1-faithful arithmetic.
//! - [`backend`] — the compute-target trait and the `CpuScalar` reference.
//! - [`subdomain`] / [`halo`] — decomposition and halo exchange.
//! - [`solver`] — the orchestrator (V1 step sequence over parts).
//!
//! V1 (`crates/lbm-core`) is frozen as the reference implementation; the
//! equivalence test suite in `tests/` compares V2 against it field-by-field.

pub mod backend;
pub mod fields;
pub mod halo;
mod kernels;
pub mod lattice;
pub mod params;
pub mod real;
pub mod solver;
pub mod subdomain;

/// Convenient glob import for the V2 API.
pub mod prelude {
    pub use crate::backend::{Backend, CellRange, CpuScalar, HostMoments, PARALLEL_MIN_CELLS};
    pub use crate::fields::{LocalGeom, SoaFields};
    pub use crate::halo::{HaloExchange, LocalPeriodic};
    pub use crate::lattice::{Face, Lattice, D2Q9, D3Q19};
    pub use crate::params::{CollisionKind, FaceBC, Reduction, StepParams};
    pub use crate::real::Real;
    pub use crate::solver::{build_wall_rims, partition, GlobalSpec, Solver, WallSpec};
    pub use crate::subdomain::Subdomain;
}
