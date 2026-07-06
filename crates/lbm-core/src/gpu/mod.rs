//! wgpu compute backend (feature `gpu`): D2Q9, f32, device-resident fields.
//!
//! Layer map:
//!
//! - [`wgsl`] *(private)* — kernel generation from the [`crate::lattice`]
//!   tables (fused push collide+stream, face BCs, moments, probe clear).
//! - [`backend`] — [`GpuContext`], [`WgpuBackend`] (the [`crate::backend::
//!   Backend`] implementation with deferred submission) and [`GpuFields`]
//!   (ping-pong device buffers).
//! - [`solver`] — [`GpuSolver`], the batch-first driver with a host mirror
//!   for setup and explicit-readback accessors.
//!
//! Capability envelope (asserted, not silently degraded): 2D D2Q9, `f32`,
//! monolithic decomposition. TRT/BGK + Guo forcing (uniform and per-cell),
//! half-way bounce-back walls (still and moving) with force probes, Zou–He
//! velocity/pressure faces (with inlet profiles), Outflow and
//! ConvectiveOutflow. Equivalence against `CpuScalar` is frozen as T14
//! (`tests/t14_backend_equiv.rs`).

mod backend;
mod solver;
mod wgsl;

pub use backend::{GpuContext, GpuError, GpuFields, GpuInitError, GpuStorage, KernelCfg, WgpuBackend};
#[allow(deprecated)]
pub use solver::GpuSolver;
