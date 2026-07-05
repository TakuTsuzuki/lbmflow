//! Rotating-impeller volume penalization for the 2D compat API.
//!
//! [`Rotor::update_force`] must be called once per step before
//! [`Simulation::step`](crate::compat::sim::Simulation::step). It adds its
//! contribution into [`Simulation::force_field_mut`], so the caller owns force
//! composition and any per-step zeroing order.
//!
//! With Guo forcing, the physical velocity is
//! `u_phys = u_star + F_total / (2 rho)`, where `u_star` is the bare first
//! moment velocity. This module computes
//! `F = 2 rho chi (u_target - u_star)`. With no other forces,
//! `u_phys = u_star + chi (u_target - u_star)` exactly after the force is
//! applied: `chi = 1` pins the blade cells to solid-body rotation with no
//! overshoot by construction, and `0 < chi < 1` gives a monotone convex blend.
//! `Simulation::ux()` / `uy()` are not used for the feedback velocity because
//! they already include the half-force correction from the force field present
//! when moments were last updated; the implementation reads the bare first
//! moment from the populations through a crate-private compat helper.

use super::real::Real;
use super::sim::Simulation;

/// Force formula shared by 2D compat and future 3D callers.
pub fn penalization_force<T: Real>(rho: T, chi: T, u_target: [T; 2], u_star: [T; 2]) -> [T; 2] {
    let two = T::r(2.0);
    [
        two * rho * chi * (u_target[0] - u_star[0]),
        two * rho * chi * (u_target[1] - u_star[1]),
    ]
}

/// Rotating blade volume-penalization source.
#[derive(Clone, Debug)]
pub struct Rotor<T: Real> {
    cx: T,
    cy: T,
    n_blades: usize,
    r_hub: T,
    r_blade: T,
    blade_thickness: T,
    omega: T,
    chi: T,
    omega_ramp_steps: u64,
    theta0: T,
    accumulated_angle: T,
    last_torque: T,
    torque_integral: T,
}

impl<T: Real> Rotor<T> {
    /// Build a rotor with conservative defaults from the 2026-07-06
    /// stability-envelope experiment.
    pub fn new(cx: T, cy: T) -> Self {
        Self {
            cx,
            cy,
            n_blades: 4,
            r_hub: T::r(4.0),
            r_blade: T::r(40.0),
            blade_thickness: T::r(3.0),
            omega: T::zero(),
            chi: T::one(),
            omega_ramp_steps: 200,
            theta0: T::zero(),
            accumulated_angle: T::zero(),
            last_torque: T::zero(),
            torque_integral: T::zero(),
        }
    }

    pub fn center(mut self, cx: T, cy: T) -> Self {
        self.cx = cx;
        self.cy = cy;
        self
    }

    pub fn n_blades(mut self, n_blades: usize) -> Self {
        assert!(n_blades > 0, "n_blades must be positive");
        self.n_blades = n_blades;
        self
    }

    pub fn r_hub(mut self, r_hub: T) -> Self {
        assert!(r_hub >= T::zero(), "r_hub must be non-negative");
        self.r_hub = r_hub;
        self
    }

    pub fn r_blade(mut self, r_blade: T) -> Self {
        assert!(r_blade > T::zero(), "r_blade must be positive");
        self.r_blade = r_blade;
        self
    }

    pub fn blade_thickness(mut self, blade_thickness: T) -> Self {
        assert!(
            blade_thickness > T::zero(),
            "blade_thickness must be positive"
        );
        self.blade_thickness = blade_thickness;
        self
    }

    pub fn omega(mut self, omega: T) -> Self {
        assert!(omega.is_finite(), "omega must be finite");
        self.omega = omega;
        self
    }

    pub fn chi(mut self, chi: T) -> Self {
        assert!(chi > T::zero() && chi <= T::one(), "chi must be in (0, 1]");
        self.chi = chi;
        self
    }

    pub fn omega_ramp_steps(mut self, omega_ramp_steps: u64) -> Self {
        self.omega_ramp_steps = omega_ramp_steps;
        self
    }

    pub fn theta0(mut self, theta0: T) -> Self {
        assert!(theta0.is_finite(), "theta0 must be finite");
        self.theta0 = theta0;
        self
    }

    /// Last-step reaction torque on the rotor, `sum r x (-F)`.
    ///
    /// Positive `omega` gives positive torque on the fluid during spin-up and
    /// therefore a negative reaction torque on the rotor.
    pub fn torque(&self) -> T {
        self.last_torque
    }

    /// Time integral of [`Rotor::torque`] over all calls to
    /// [`Rotor::update_force`] with `dt = 1`.
    pub fn torque_integral(&self) -> T {
        self.torque_integral
    }

    pub fn accumulated_angle(&self) -> T {
        self.accumulated_angle
    }

    fn omega_eff(&self, time: u64) -> T {
        if self.omega_ramp_steps == 0 {
            return self.omega;
        }
        let s = time.min(self.omega_ramp_steps);
        self.omega * T::r(s as f64 / self.omega_ramp_steps as f64)
    }

    fn blade_chi(&self, x: T, y: T) -> T {
        let dx = x - self.cx;
        let dy = y - self.cy;
        let r2 = dx * dx + dy * dy;
        if r2 < self.r_hub * self.r_hub || r2 > self.r_blade * self.r_blade {
            return T::zero();
        }
        let half = T::r(0.5) * self.blade_thickness;
        let two_pi = T::r(std::f64::consts::TAU);
        for b in 0..self.n_blades {
            let theta = self.theta0
                + self.accumulated_angle
                + two_pi * T::r(b as f64 / self.n_blades as f64);
            let nx = -theta.sin();
            let ny = theta.cos();
            let perp = (dx * nx + dy * ny).abs();
            if perp <= half {
                return self.chi;
            }
        }
        T::zero()
    }

    /// Add this step's penalization force into `sim.force_field_mut()`.
    pub fn update_force(&mut self, sim: &mut Simulation<T>) {
        let omega = self.omega_eff(sim.time());
        let (nx, ny) = (sim.nx(), sim.ny());
        let mut additions = vec![[T::zero(); 2]; nx * ny];
        let mut reaction_torque = T::zero();
        let mut fluid_torque = T::zero();

        for y in 0..ny {
            for x in 0..nx {
                if sim.is_solid(x, y) {
                    continue;
                }
                let xf = T::r(x as f64);
                let yf = T::r(y as f64);
                let chi = self.blade_chi(xf, yf);
                if chi == T::zero() {
                    continue;
                }
                let rx = xf - self.cx;
                let ry = yf - self.cy;
                let u_target = [-omega * ry, omega * rx];
                let u_star = sim.bare_velocity(x, y);
                let force = penalization_force(sim.rho(x, y), chi, u_target, u_star);
                let i = y * nx + x;
                additions[i] = force;
                let tq_fluid = rx * force[1] - ry * force[0];
                fluid_torque = fluid_torque + tq_fluid;
                reaction_torque = reaction_torque - tq_fluid;
            }
        }

        let field = sim.force_field_mut();
        for (dst, add) in field.iter_mut().zip(additions.iter()) {
            dst[0] = dst[0] + add[0];
            dst[1] = dst[1] + add[1];
        }
        self.last_torque = reaction_torque;
        self.torque_integral = self.torque_integral + reaction_torque;
        let _ = fluid_torque;
        self.accumulated_angle = self.accumulated_angle + omega;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compat::prelude::*;

    fn closed_tank(nx: usize, ny: usize) -> Simulation<f64> {
        SimConfig {
            nx,
            ny,
            nu: 0.02,
            collision: Collision::Trt {
                magic: Collision::MAGIC_STD,
            },
            edges: Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        }
        .build()
        .unwrap()
    }

    #[test]
    fn penalization_formula_has_no_overshoot() {
        let rho: f64 = 1.2;
        let u_star: [f64; 2] = [0.03, -0.02];
        let u_target: [f64; 2] = [0.08, 0.04];
        for chi in [0.25, 0.5, 1.0] {
            let f = penalization_force(rho, chi, u_target, u_star);
            let u_phys = [
                u_star[0] + f[0] / (2.0 * rho),
                u_star[1] + f[1] / (2.0 * rho),
            ];
            let expected = [
                u_star[0] + chi * (u_target[0] - u_star[0]),
                u_star[1] + chi * (u_target[1] - u_star[1]),
            ];
            assert!((u_phys[0] - expected[0]).abs() < 1e-15);
            assert!((u_phys[1] - expected[1]).abs() < 1e-15);
            assert!((u_phys[0] - u_star[0]).abs() <= (u_target[0] - u_star[0]).abs());
            assert!((u_phys[1] - u_star[1]).abs() <= (u_target[1] - u_star[1]).abs());
        }
    }

    #[test]
    fn accumulated_angle_uses_ramped_sum() {
        let mut rotor = Rotor::new(8.0, 8.0)
            .omega(0.02)
            .omega_ramp_steps(4)
            .r_hub(1.0)
            .r_blade(5.0)
            .blade_thickness(2.0);
        let mut sim = closed_tank(16, 16);
        for _ in 0..6 {
            rotor.update_force(&mut sim);
            sim.step();
        }
        let expected = 0.02 * (0.0 / 4.0 + 1.0 / 4.0 + 2.0 / 4.0 + 3.0 / 4.0 + 1.0 + 1.0);
        assert!((rotor.accumulated_angle() - expected).abs() < 1e-15);
    }

    #[test]
    fn spinup_fluid_torque_sign_is_positive() {
        let mut rotor = Rotor::new(32.0, 32.0)
            .omega(0.0025)
            .chi(1.0)
            .r_hub(3.0)
            .r_blade(16.0)
            .blade_thickness(3.0)
            .omega_ramp_steps(0);
        let mut sim = closed_tank(64, 64);
        rotor.update_force(&mut sim);
        assert!(rotor.torque() < 0.0);
    }

    #[test]
    fn solid_body_tracking_inside_blades() {
        let mut rotor = Rotor::new(32.0, 32.0)
            .omega(0.0025)
            .chi(1.0)
            .r_hub(3.0)
            .r_blade(16.0)
            .blade_thickness(5.0)
            .omega_ramp_steps(200);
        let mut sim = closed_tank(64, 64);
        for _ in 0..5000 {
            sim.force_field_mut().fill([0.0; 2]);
            rotor.update_force(&mut sim);
            sim.step();
        }
        let u_tip = 0.0025 * 16.0;
        let mut max_rel: f64 = 0.0;
        for y in 1..63 {
            for x in 1..63 {
                let xf = x as f64;
                let yf = y as f64;
                let dx = xf - 32.0;
                let dy = yf - 32.0;
                let r = (dx * dx + dy * dy).sqrt();
                if !(5.0..=15.0).contains(&r) {
                    continue;
                }
                let mut strictly_inside = false;
                for b in 0..4 {
                    let theta = rotor.accumulated_angle() + std::f64::consts::TAU * b as f64 / 4.0;
                    let perp = (dx * -theta.sin() + dy * theta.cos()).abs();
                    if perp <= 1.5 {
                        strictly_inside = true;
                    }
                }
                if !strictly_inside {
                    continue;
                }
                let u_target = [-0.0025 * dy, 0.0025 * dx];
                let err = ((sim.ux(x, y) - u_target[0]).powi(2)
                    + (sim.uy(x, y) - u_target[1]).powi(2))
                .sqrt()
                    / u_tip;
                max_rel = max_rel.max(err);
            }
        }
        assert!(max_rel < 0.01, "max_rel={max_rel}");
    }
}
