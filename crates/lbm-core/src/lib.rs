//! # lbm-core
//!
//! Lattice Boltzmann core (the V2 architecture, sole engine since the V1
//! retirement on 2026-07-05): the dimension / lattice / precision / backend /
//! decomposition axes are orthogonal (docs/ARCHITECTURE_V2.md). The physics
//! kernels are written once, generically over a [`lattice::Lattice`] and a
//! [`real::Real`], and specialise at compile time.
//!
//! Layer map (docs/ARCHITECTURE_V2.md §1):
//!
//! - [`lattice`] — compile-time velocity sets (D2Q9, D3Q19, D3Q27) with derived
//!   tables (TRT pairs, per-face unknown sets).
//! - [`fields`] — q-major SoA deviation storage over halo-padded local boxes.
//! - `kernels` (private) — the physics (collide/stream/moments/BCs), written
//!   once, generic over lattice and precision; V1-faithful arithmetic.
//! - [`backend`] — the compute-target trait and the `CpuScalar` reference.
//! - [`subdomain`] / [`halo`] — decomposition and halo exchange.
//! - [`solver`] — the orchestrator (V1 step sequence over parts).
//! - [`compat`] — the retired V1 engine's public API as a supported facade
//!   (scenario/CLI/wasm run through it); V1↔V2 equivalence was proven by
//!   `tests/v1_match.rs`, frozen and removed with V1 (see branch history).

pub mod backend;
pub mod backend_simd;
pub mod bouzidi;
pub mod compat;
#[cfg(feature = "mpi")]
pub mod dist;
pub mod fields;
#[cfg(feature = "gpu")]
pub mod gpu;
pub mod halo;
mod kernels;
pub mod lattice;
pub mod les;
pub mod params;
pub mod particles;
pub mod real;
pub mod rotating_ibm;
pub mod solver;
pub mod subdomain;

/// Convenient glob import for the V2 API.
pub mod prelude {
    pub use crate::backend::{Backend, CellRange, CpuScalar, HostMoments, PARALLEL_MIN_CELLS};
    pub use crate::backend_simd::CpuSimd;
    pub use crate::bouzidi::{BouzidiLink, BouzidiLinks};
    pub use crate::fields::{LocalGeom, SoaFields};
    pub use crate::halo::{HaloExchange, InProcess, LocalPeriodic};
    pub use crate::lattice::{Face, Lattice, D2Q9, D3Q19, D3Q27};
    pub use crate::les::{WaleLes, WALE_CW};
    pub use crate::params::{
        CollisionKind, FaceBC, FacePatch, Reduction, SourceKind, SourceRegion, StepParams,
        VolumeSource,
    };
    pub use crate::real::Real;
    pub use crate::rotating_ibm::{DirectForcingConfig, IbmDiagnostics, IbmMarker, RotatingBody};
    pub use crate::solver::{
        build_wall_rims, partition, CheckpointError, Diverged, GlobalSpec, Solver, SpecError,
        WallSpec,
    };
    pub use crate::subdomain::Subdomain;

    #[cfg(feature = "gpu")]
    #[allow(deprecated)]
    pub use crate::gpu::{GpuContext, GpuSolver, WgpuBackend};

    #[cfg(feature = "mpi")]
    pub use crate::dist::{MpiExchange, MpiSolver};
}
