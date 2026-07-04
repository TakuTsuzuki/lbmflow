//! Scalar abstraction so the engine runs in either `f32` or `f64`.

use num_traits::{Float, FromPrimitive};

/// Floating-point scalar used by the simulation (`f32` or `f64`).
///
/// The accuracy/speed trade-off axis #1: `Simulation<f64>` for validation-grade
/// accuracy, `Simulation<f32>` for roughly 2x memory-bandwidth-bound speed.
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
