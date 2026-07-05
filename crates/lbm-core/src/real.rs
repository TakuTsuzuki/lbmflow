//! Scalar abstraction so the engine runs in either `f32` or `f64`.
//!
//! Identical contract to the retired V1 engine's `real::Real` (the V2 core
//! was forbidden from build-depending on V1, hence the duplicate; V1 itself
//! was deleted 2026-07-05 — `compat::real` re-exports this trait).

use num_traits::{Float, FromPrimitive};

/// Floating-point scalar used by the simulation (`f32` or `f64`).
///
/// The accuracy/speed trade-off axis #1: `f64` for validation-grade accuracy,
/// `f32` for roughly 2x memory-bandwidth-bound speed. The deviation storage
/// (populations hold `f - w`) keeps `f32` runs at validation grade for the
/// standard benchmarks (docs/PHYSICS.md).
pub trait Real:
    Float
    + FromPrimitive
    + Default
    + std::fmt::Debug
    + std::fmt::Display
    + std::iter::Sum
    + Send
    + Sync
    + 'static
{
    /// Convert an `f64` constant into `Self` (lossy for `f32`).
    #[inline(always)]
    fn r(v: f64) -> Self {
        Self::from_f64(v).expect("finite f64 constant")
    }

    /// Convert into `f64` (exact for both supported types).
    #[inline(always)]
    fn as_f64(self) -> f64 {
        self.to_f64().expect("real is representable as f64")
    }
}

impl Real for f32 {}
impl Real for f64 {}
