//! Shan–Chen pseudopotential multiphase model (V1 `lbm_core::multiphase`
//! facade; identical numerics, driven through the facade `Simulation`).
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
//! use lbm_core::compat::prelude::*;
//! use lbm_core::compat::multiphase::{Psi, ShanChen};
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

use super::lattice::{CS2, CX, CY, Q, W};
use super::real::Real;
use super::sim::Simulation;

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

/// Two-component (immiscible) Shan–Chen driver.
///
/// Components A and B are two [`Simulation`]s sharing the same grid and
/// solid geometry. A cross-repulsion force `F_A(x) = -G_ab rho_A(x)
/// Σ_q w_q rho_B(x+c_q) c_q` (and symmetrically for B) separates the fluids;
/// the force pairs are equal-and-opposite per link, so total momentum is
/// conserved. Optional per-component gravity `F = rho_sigma(x) * g_sigma`
/// enables buoyancy-driven flows (Rayleigh–Taylor).
///
/// Usage (both sims must be stepped together):
/// ```no_run
/// use lbm_core::compat::prelude::*;
/// use lbm_core::compat::multiphase::MultiComponent;
///
/// let mut a: Simulation<f64> = SimConfig::default().build().unwrap();
/// let mut b: Simulation<f64> = SimConfig::default().build().unwrap();
/// let mc = MultiComponent::new(2.6).with_gravity([0.0, -5e-5], [0.0, 0.0]);
/// for _ in 0..1000 {
///     mc.update_forces(&mut a, &mut b);
///     a.step();
///     b.step();
/// }
/// ```
#[derive(Clone, Debug)]
pub struct MultiComponent<T: Real> {
    /// Cross-component repulsion strength (positive separates the fluids;
    /// ~1.0–2.0 with psi = rho and background densities O(1)).
    pub g_ab: T,
    /// Gravity applied to component A as `rho_A(x) * g_a` per cell.
    pub g_a: [T; 2],
    /// Gravity applied to component B as `rho_B(x) * g_b` per cell.
    pub g_b: [T; 2],
    /// Wall affinity for A (negative attracts A to walls).
    pub g_wall_a: T,
    /// Wall affinity for B.
    pub g_wall_b: T,
}

impl<T: Real> MultiComponent<T> {
    /// New driver with the given cross-repulsion, no gravity, neutral walls.
    pub fn new(g_ab: f64) -> Self {
        Self {
            g_ab: T::r(g_ab),
            g_a: [T::zero(); 2],
            g_b: [T::zero(); 2],
            g_wall_a: T::zero(),
            g_wall_b: T::zero(),
        }
    }

    /// Set per-component gravity vectors.
    pub fn with_gravity(mut self, g_a: [f64; 2], g_b: [f64; 2]) -> Self {
        self.g_a = [T::r(g_a[0]), T::r(g_a[1])];
        self.g_b = [T::r(g_b[0]), T::r(g_b[1])];
        self
    }

    /// Set wall affinities (negative = wetting for that component).
    pub fn with_walls(mut self, g_wall_a: f64, g_wall_b: f64) -> Self {
        self.g_wall_a = T::r(g_wall_a);
        self.g_wall_b = T::r(g_wall_b);
        self
    }

    /// Compute cross-interaction + gravity forces into both simulations'
    /// force fields. Panics if the two grids differ; in debug builds also
    /// verifies that the solid geometry and the axis periodicities agree
    /// (the stencil below evaluates both components over component A's
    /// mask and wrap rules, so a divergence would silently misplace forces —
    /// A-10e).
    pub fn update_forces(&self, a: &mut Simulation<T>, b: &mut Simulation<T>) {
        assert_eq!(
            (a.nx(), a.ny()),
            (b.nx(), b.ny()),
            "components must share the grid"
        );
        debug_assert_eq!(
            a.solid_field(),
            b.solid_field(),
            "components must share the solid geometry (forces use A's mask)"
        );
        debug_assert_eq!(
            (a.is_periodic_x(), a.is_periodic_y()),
            (b.is_periodic_x(), b.is_periodic_y()),
            "components must share axis periodicity (forces use A's wrap)"
        );
        let (nx, ny) = (a.nx(), a.ny());
        let per_x = a.is_periodic_x();
        let per_y = a.is_periodic_y();
        let n = nx * ny;

        let wrap = |x: isize, y: isize| -> Option<usize> {
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

        // psi = rho (0 on solids), captured before mutating force fields
        let solid: Vec<bool> = a.solid_field().to_vec();
        let psi_a: Vec<T> = a
            .rho_field()
            .iter()
            .zip(&solid)
            .map(|(&r, &s)| if s { T::zero() } else { r })
            .collect();
        let psi_b: Vec<T> = b
            .rho_field()
            .iter()
            .zip(&solid)
            .map(|(&r, &s)| if s { T::zero() } else { r })
            .collect();

        let mut fa = vec![[T::zero(); 2]; n];
        let mut fb = vec![[T::zero(); 2]; n];
        for y in 0..ny {
            for x in 0..nx {
                let i = y * nx + x;
                if solid[i] {
                    continue;
                }
                let mut sum_b = [T::zero(); 2];
                let mut sum_a = [T::zero(); 2];
                let mut adh = [T::zero(); 2];
                for q in 1..Q {
                    let Some(j) = wrap(x as isize + CX[q] as isize, y as isize + CY[q] as isize)
                    else {
                        continue;
                    };
                    let w = T::r(W[q]);
                    let (cxq, cyq) = (T::r(CX[q] as f64), T::r(CY[q] as f64));
                    if solid[j] {
                        adh[0] = adh[0] + w * cxq;
                        adh[1] = adh[1] + w * cyq;
                    } else {
                        sum_b[0] = sum_b[0] + w * psi_b[j] * cxq;
                        sum_b[1] = sum_b[1] + w * psi_b[j] * cyq;
                        sum_a[0] = sum_a[0] + w * psi_a[j] * cxq;
                        sum_a[1] = sum_a[1] + w * psi_a[j] * cyq;
                    }
                }
                fa[i] = [
                    -psi_a[i] * (self.g_ab * sum_b[0] + self.g_wall_a * adh[0])
                        + psi_a[i] * self.g_a[0],
                    -psi_a[i] * (self.g_ab * sum_b[1] + self.g_wall_a * adh[1])
                        + psi_a[i] * self.g_a[1],
                ];
                fb[i] = [
                    -psi_b[i] * (self.g_ab * sum_a[0] + self.g_wall_b * adh[0])
                        + psi_b[i] * self.g_b[0],
                    -psi_b[i] * (self.g_ab * sum_a[1] + self.g_wall_b * adh[1])
                        + psi_b[i] * self.g_b[1],
                ];
            }
        }
        a.force_field_mut().copy_from_slice(&fa);
        b.force_field_mut().copy_from_slice(&fb);
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
    /// Virtual wall density: when set, solid neighbours contribute
    /// `psi(wall_rho)` to the *cohesion* sum (instead of 0). This gives full
    /// contact-angle control: `wall_rho` near the liquid density wets the
    /// wall (θ → 0°), near the vapour density de-wets it (θ → 180°),
    /// intermediate values interpolate through 90°. Preferred over `g_wall`.
    pub wall_rho: Option<f64>,
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
            wall_rho: None,
            psi: Psi::Classic,
            psi_buf: std::cell::RefCell::new(Vec::new()),
        }
    }

    /// Set the wall adhesion strength (contact-angle control).
    pub fn with_wall(mut self, g_wall: f64) -> Self {
        self.g_wall = T::r(g_wall);
        self
    }

    /// Set the virtual wall density (full-range contact-angle control).
    pub fn with_wall_rho(mut self, wall_rho: f64) -> Self {
        self.wall_rho = Some(wall_rho);
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
        let wrap = |x: isize, y: isize| -> Option<usize> {
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
                let psi_wall = T::r(self.wall_rho.map_or(0.0, |r| self.psi.eval(r)));
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
                        // virtual wall density feeds the cohesion sum;
                        // g_wall adds the legacy adhesion term on top
                        sx = sx + w * psi_wall * cx;
                        sy = sy + w * psi_wall * cy;
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
