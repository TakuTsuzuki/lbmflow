//! # lbm-core
//!
//! Lattice Boltzmann method (D2Q9) fluid simulation engine.
//!
//! - Collision: BGK or TRT (magic 3/16 ⇒ exact half-way walls for parabolic
//!   flows), with 2nd-order Guo forcing.
//! - Boundaries: periodic, half-way bounce-back (stationary/moving walls,
//!   realised as one-cell solid rims for domain edges), Zou–He velocity /
//!   pressure open edges, zero-gradient outflow, interior solid obstacles.
//! - Precision: generic over `f32` / `f64` (see [`real::Real`]).
//! - Parallelism: rayon row-parallel loops behind the `parallel` feature
//!   (enabled by default; disable for WASM).
//!
//! See `docs/VALIDATION.md` at the repository root for the validation matrix
//! this engine is required to pass.
//!
//! ```
//! use lbm_core::prelude::*;
//!
//! // Lid-driven cavity, Re = U*L/nu = 0.1*62/0.02 = 310
//! let mut sim: Simulation<f64> = SimConfig {
//!     nx: 64,
//!     ny: 64,
//!     nu: 0.02,
//!     edges: Edges {
//!         left: EdgeBC::BounceBack,
//!         right: EdgeBC::BounceBack,
//!         bottom: EdgeBC::BounceBack,
//!         top: EdgeBC::MovingWall { u: [0.1, 0.0] },
//!     },
//!     ..Default::default()
//! }
//! .build()
//! .unwrap();
//! sim.run(100);
//! assert!(sim.ux(32, 60) != 0.0); // fluid is being dragged by the lid
//! ```

pub mod domain;
pub mod lattice;
pub mod real;
pub mod sim;

/// Convenient glob import: `use lbm_core::prelude::*;`.
pub mod prelude {
    pub use crate::domain::{Collision, ConfigError, Edge, EdgeBC, Edges, SimConfig, MAX_SPEED};
    pub use crate::real::Real;
    pub use crate::sim::Simulation;
}
