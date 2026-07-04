//! Shan–Chen pseudopotential multiphase model (single component).
//!
//! The model injects an interaction force computed from the current density
//! field before each step:
//!
//! ```text
//! F(x)     = -G  psi(x) Σ_q w_q psi(x + c_q) c_q       (fluid–fluid cohesion)
//! F_ads(x) = -Gw psi(x) Σ_q w_q  s(x + c_q) c_q        (wall adhesion, s = solid)
//! ```
//!
//! With `G` sufficiently negative the fluid separates into a dense (liquid)
//! and a light (vapour) phase. The equation of state is
//! `p = cs² rho + (G cs² / 2) psi(rho)²` — use [`ShanChen::pressure`] when
//! comparing pressures (e.g. the Laplace law), NOT the bare `cs² rho`.
//!
//! Usage:
//! ```
//! use lbm_core::prelude::*;
//! use lbm_core::multiphase::{Psi, ShanChen};
//!
//! let mut sim: Simulation<f64> = SimConfig {
//!     nx: 64, ny: 64, nu: 1.0 / 6.0, ..Default::default()
//! }.build().unwrap();
//! // dense blob in light vapour
//! sim.init_with(|x, y| {
//!     let inside = (x as f64 - 32.0).powi(2) + (y as f64 - 32.0).powi(2) < 12.0f64.powi(2);
//!     (if inside { 2.0 } else { 0.15 }, 0.0, 0.0)
//! });
//! let sc = ShanChen::new(-5.0);
//! for _ in 0..100 {
//!     sc.update_force(&mut sim);
//!     sim.step();
//! }
//! assert!(sim.total_mass().is_finite());
//! ```

use crate::lattice::{CS2, CX, CY, Q, W};
use crate::real::Real;
use crate::sim::Simulation;

/// Pseudopotential ψ(ρ).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Psi {
    /// `psi = 1 - exp(-rho)` — the classic Shan–Chen form. Critical point at
    /// `G = -4` (ρ_c = ln 2); use `G ≈ -5` for a clear two-phase state with
    /// density ratio ~20–60.
    Classic,
    /// `psi = psi0 * exp(-rho0 / rho)` — Shan–Chen 1994 form.
    Exponential {
        /// Amplitude ψ0.
        psi0: f64,
        /// Reference density ρ0.
        rho0: f64,
    },
}

impl Psi {
    #[inline]
    fn eval(&self, rho: f64) -> f64 {
        match *self {
            Psi::Classic => 1.0 - (-rho).exp(),
            Psi::Exponential { psi0, rho0 } => psi0 * (-rho0 / rho).exp(),
        }
    }
}

/// Single-component Shan–Chen model driver.
///
/// Owns no simulation state; call [`ShanChen::update_force`] before every
/// [`Simulation::step`].
#[derive(Clone, Debug)]
pub struct ShanChen<T: Real> {
    /// Fluid–fluid interaction strength (negative = cohesion).
    pub g: T,
    /// Wall adhesion strength: negative values wet the wall (small contact
    /// angle), positive values de-wet it. Zero = neutral (~90°).
    pub g_wall: T,
    /// Pseudopotential form.
    pub psi: Psi,
    /// Scratch buffer for ψ(ρ), reused between calls.
    psi_buf: std::cell::RefCell<Vec<T>>,
}

impl<T: Real> ShanChen<T> {
    /// New model with the given cohesion strength, neutral walls, classic ψ.
    pub fn new(g: f64) -> Self {
        Self {
            g: T::r(g),
            g_wall: T::zero(),
            psi: Psi::Classic,
            psi_buf: std::cell::RefCell::new(Vec::new()),
        }
    }

    /// Set the wall adhesion strength (contact-angle control).
    pub fn with_wall(mut self, g_wall: f64) -> Self {
        self.g_wall = T::r(g_wall);
        self
    }

    /// Set the pseudopotential form.
    pub fn with_psi(mut self, psi: Psi) -> Self {
        self.psi = psi;
        self
    }

    /// Shan–Chen equation of state: `p = cs² rho + (G cs²/2) psi(rho)²`.
    pub fn pressure(&self, rho: T) -> T {
        let cs2 = T::r(CS2);
        let psi = T::r(self.psi.eval(rho.as_f64()));
        cs2 * rho + T::r(0.5) * self.g * cs2 * psi * psi
    }

    /// Compute the interaction + adhesion force from the current density
    /// field and store it into the simulation's per-cell force field.
    ///
    /// Out-of-domain neighbours on non-periodic edges contribute nothing
    /// (zero-gradient approximation); solid neighbours contribute to the
    /// adhesion term only.
    pub fn update_force(&self, sim: &mut Simulation<T>) {
        let (nx, ny) = (sim.nx(), sim.ny());
        let per_x = sim.is_periodic_x();
        let per_y = sim.is_periodic_y();
        let n = nx * ny;

        // psi field (0 on solids so cohesion skips walls cleanly)
        let mut psi_buf = self.psi_buf.borrow_mut();
        psi_buf.resize(n, T::zero());
        for i in 0..n {
            psi_buf[i] = if sim.solid_field()[i] {
                T::zero()
            } else {
                T::r(self.psi.eval(sim.rho_field()[i].as_f64()))
            };
        }

        // borrow rules: compute into a local, then write
        let mut wrap = |x: isize, y: isize| -> Option<usize> {
            let mut x = x;
            let mut y = y;
            if x < 0 || x >= nx as isize {
                if per_x {
                    x = (x + nx as isize) % nx as isize;
                } else {
                    return None;
                }
            }
            if y < 0 || y >= ny as isize {
                if per_y {
                    y = (y + ny as isize) % ny as isize;
                } else {
                    return None;
                }
            }
            Some(y as usize * nx + x as usize)
        };

        let mut forces = vec![[T::zero(); 2]; n];
        for y in 0..ny {
            for x in 0..nx {
                let i = y * nx + x;
                if sim.solid_field()[i] {
                    continue;
                }
                let psi_i = psi_buf[i];
                if psi_i == T::zero() {
                    continue;
                }
                let mut sx = T::zero();
                let mut sy = T::zero();
                let mut ax = T::zero();
                let mut ay = T::zero();
                for q in 1..Q {
                    let Some(j) = wrap(x as isize + CX[q] as isize, y as isize + CY[q] as isize)
                    else {
                        continue;
                    };
                    let w = T::r(W[q]);
                    let (cx, cy) = (T::r(CX[q] as f64), T::r(CY[q] as f64));
                    if sim.solid_field()[j] {
                        ax = ax + w * cx;
                        ay = ay + w * cy;
                    } else {
                        let pj = psi_buf[j];
                        sx = sx + w * pj * cx;
                        sy = sy + w * pj * cy;
                    }
                }
                forces[i] = [
                    -psi_i * (self.g * sx + self.g_wall * ax),
                    -psi_i * (self.g * sy + self.g_wall * ay),
                ];
            }
        }
        sim.force_field_mut().copy_from_slice(&forces);
    }
}
