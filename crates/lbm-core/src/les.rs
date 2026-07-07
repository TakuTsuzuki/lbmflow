//! Large-eddy simulation closures.
//!
//! The MF-beta subset implements WALE (Nicoud & Ducros, 1999) as a
//! solver-level relaxation-field driver. The driver computes an effective
//! `omega_plus = 1 / (3 * (nu0 + nu_t) + 0.5)` field from the current velocity
//! gradient and applies it to the next collision, so the model has a one-step
//! lag.

use crate::backend::Backend;
use crate::halo::HaloExchange;
use crate::lattice::Lattice;
use crate::real::Real;
use crate::solver::Solver;

/// WALE model coefficient recommended by Nicoud & Ducros (1999).
pub const WALE_CW: f64 = 0.325;

/// Diagnostics from the most recent WALE relaxation-field update.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaleLesDiagnostics<T: Real> {
    /// Cells where raw WALE `nu_t` exceeded the configured active bound.
    pub clipped_cells: usize,
    /// `clipped_cells / total_cells` for the update.
    pub clipped_fraction: T,
    /// Maximum raw WALE `nu_t` before applying the optional upper bound.
    pub max_nu_t_before_clipping: T,
    /// Configured upper bound on `tau_eff`, if clipping is active.
    pub tau_eff_max: Option<T>,
    /// Equivalent active upper bound on `nu_t`, if clipping is active.
    pub nu_t_max: Option<T>,
}

impl<T: Real> Default for WaleLesDiagnostics<T> {
    fn default() -> Self {
        Self {
            clipped_cells: 0,
            clipped_fraction: T::zero(),
            max_nu_t_before_clipping: T::zero(),
            tau_eff_max: None,
            nu_t_max: None,
        }
    }
}

/// WALE SGS-viscosity driver.
#[derive(Clone, Debug)]
pub struct WaleLes<T: Real> {
    cw: T,
    delta: T,
    tau_eff_max: Option<T>,
    omega: Vec<T>,
    nu_t: Vec<T>,
    diagnostics: WaleLesDiagnostics<T>,
}

impl<T: Real> Default for WaleLes<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Real> WaleLes<T> {
    /// Build WALE with `Cw = 0.325` and lattice filter width `Delta = 1`.
    pub fn new() -> Self {
        Self {
            cw: T::r(WALE_CW),
            delta: T::one(),
            tau_eff_max: None,
            omega: Vec::new(),
            nu_t: Vec::new(),
            diagnostics: WaleLesDiagnostics::default(),
        }
    }

    /// Build WALE with an explicit upper bound on
    /// `tau_eff = 0.5 + 3 * (nu_0 + nu_t)`.
    pub fn with_tau_eff_max(mut self, tau_eff_max: T) -> Self {
        self.set_tau_eff_max(Some(tau_eff_max));
        self
    }

    /// Set or clear the explicit upper bound on `tau_eff`.
    ///
    /// `None` keeps the historical unclipped WALE path. A configured bound is
    /// validated again during [`update`](Self::update), where the laminar
    /// viscosity is known.
    pub fn set_tau_eff_max(&mut self, tau_eff_max: Option<T>) {
        if let Some(value) = tau_eff_max {
            assert!(
                value.is_finite() && value > T::r(0.5),
                "tau_eff_max must be finite and greater than 0.5"
            );
        }
        self.tau_eff_max = tau_eff_max;
    }

    /// Configured upper bound on `tau_eff`, if clipping is active.
    pub fn tau_eff_max(&self) -> Option<T> {
        self.tau_eff_max
    }

    /// Last computed eddy viscosity in global compact order.
    pub fn nu_t(&self) -> &[T] {
        &self.nu_t
    }

    /// Diagnostics from the most recent update.
    pub fn diagnostics(&self) -> WaleLesDiagnostics<T> {
        self.diagnostics
    }

    /// Recompute `nu_t` and install the next-step `omega_plus` field.
    ///
    /// Standard WALE definition:
    /// `S_ij = (g_ij + g_ji)/2`,
    /// `S^d_ij = (g_ik g_kj + g_jk g_ki)/2 - delta_ij tr(g^2)/3`, and
    /// `nu_t = (Cw Delta)^2 (S^d:S^d)^(3/2) /
    /// ((S:S)^(5/2) + (S^d:S^d)^(5/4))`.
    /// The zero-gradient `0/0` limit is defined as `nu_t = 0`.
    pub fn update<L, B, H>(&mut self, solver: &mut Solver<L, T, B, H>)
    where
        L: Lattice,
        B: Backend<L, T, Fields = crate::fields::SoaFields<T>>,
        H: HaloExchange<T>,
    {
        let grad = solver.gather_velocity_gradient();
        let n = grad.len();
        if self.omega.len() != n {
            self.omega.resize(n, T::zero());
            self.nu_t.resize(n, T::zero());
        }
        let cw_delta_sq = (self.cw * self.delta) * (self.cw * self.delta);
        let nu0 = T::r(solver.nu());
        let tau0 = T::r(0.5) + T::r(3.0) * nu0;
        let nu_t_max = self.tau_eff_max.map(|tau_eff_max| {
            assert!(
                tau_eff_max >= tau0,
                "tau_eff_max must be at least the laminar tau for this solver"
            );
            (tau_eff_max - tau0) / T::r(3.0)
        });
        let mut clipped_cells = 0usize;
        let mut max_nu_t_before_clipping = T::zero();
        for (i, g) in grad.iter().enumerate() {
            let mut s = [[T::zero(); 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    s[a][b] = T::r(0.5) * (g[a][b] + g[b][a]);
                }
            }
            let mut g2 = [[T::zero(); 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    for k in 0..3 {
                        g2[a][b] = g2[a][b] + g[a][k] * g[k][b];
                    }
                }
            }
            let tr_g2 = g2[0][0] + g2[1][1] + g2[2][2];
            let mut sd = [[T::zero(); 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    sd[a][b] = T::r(0.5) * (g2[a][b] + g2[b][a]);
                    if a == b {
                        sd[a][b] = sd[a][b] - tr_g2 / T::r(3.0);
                    }
                }
            }
            let mut ss = T::zero();
            let mut sdsd = T::zero();
            for a in 0..3 {
                for b in 0..3 {
                    ss = ss + s[a][b] * s[a][b];
                    sdsd = sdsd + sd[a][b] * sd[a][b];
                }
            }
            let denom = ss.powf(T::r(2.5)) + sdsd.powf(T::r(1.25));
            let nut = if denom > T::zero() {
                cw_delta_sq * sdsd.powf(T::r(1.5)) / denom
            } else {
                T::zero()
            };
            max_nu_t_before_clipping = max_nu_t_before_clipping.max(nut);
            let nut = if let Some(max) = nu_t_max {
                if nut > max {
                    clipped_cells += 1;
                    max
                } else {
                    nut
                }
            } else {
                nut
            };
            self.nu_t[i] = nut;
            self.omega[i] = T::one() / (T::r(3.0) * (nu0 + nut) + T::r(0.5));
        }
        self.diagnostics = WaleLesDiagnostics {
            clipped_cells,
            clipped_fraction: if n == 0 {
                T::zero()
            } else {
                T::from_usize(clipped_cells).expect("cell count fits Real")
                    / T::from_usize(n).expect("cell count fits Real")
            },
            max_nu_t_before_clipping,
            tau_eff_max: self.tau_eff_max,
            nu_t_max,
        };
        solver.set_omega_field(Some(&self.omega));
    }

    /// Remove the installed relaxation field.
    pub fn clear<L, B, H>(&mut self, solver: &mut Solver<L, T, B, H>)
    where
        L: Lattice,
        B: Backend<L, T, Fields = crate::fields::SoaFields<T>>,
        H: HaloExchange<T>,
    {
        solver.set_omega_field(None);
    }
}
