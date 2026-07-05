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
//!
//! V1 (`crates/lbm-core`) is frozen as the reference implementation; the
//! equivalence test suite in `tests/` compares V2 against it field-by-field.

pub mod fields;
pub mod lattice;
pub mod real;
