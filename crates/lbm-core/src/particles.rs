//! One-way Lagrangian particles for resolved-flow callers.
//!
//! This module is deliberately engine-agnostic: particles are advanced from a
//! caller-supplied sampler closure and never hold a reference to a solver,
//! backend, or simulation type. It implements only the one-way FR-PART-01
//! subset: particles feel the sampled fluid velocity, buoyancy-reduced gravity,
//! and Schiller-Naumann drag, but they do not apply any reaction force back to
//! the fluid. Two-way/four-way coupling, Saffman/Basset/Faxen forces, collision
//! models, and stochastic LES dispersion are not implemented here.
//!
//! No stochastic terms are used. In particular, this module intentionally does
//! not include the uniform random kick pseudo-turbulence anti-pattern called out
//! by FR-PART-03. Exposure accumulation is deterministic and caller-supplied:
//! pass a shear-rate or other resolved-only exposure field through the
//! `exposure_rate` closure.
//!
//! Solid handling uses a simple staircase-wall model. If a proposed end point
//! is solid, the step segment is subdivided so no sub-step spans more than one
//! cell in any coordinate, preventing tunneling through 1-cell walls. The first
//! solid sub-step defines the blocked axis as the coordinate whose one-axis
//! advance first enters solid, falling back to the largest component of the
//! attempted displacement. The blocked velocity component is reflected with the
//! configured restitution; tangential motion is preserved when the
//! axis-corrected point is fluid. With restitution 0, particles can rest in
//! contact under gravity while still responding to later tangential drag.

/// State of one Lagrangian particle in lattice units.
#[derive(Clone, Debug, PartialEq)]
pub struct Particle {
    pub pos: [f64; 3],
    pub vel: [f64; 3],
    pub d: f64,
    pub rho_p: f64,
    pub exposure: f64,
}

/// Fluid sample at a particle position.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sample {
    pub u: [f64; 3],
    pub solid: bool,
}

/// One-way particle container and fluid parameters in lattice units.
#[derive(Clone, Debug, PartialEq)]
pub struct ParticleSet {
    pub particles: Vec<Particle>,
    pub rho_f: f64,
    pub nu: f64,
    pub g: [f64; 3],
    pub restitution: f64,
}

impl ParticleSet {
    /// Builds a particle set with zero restitution.
    pub fn new(particles: Vec<Particle>, rho_f: f64, nu: f64, g: [f64; 3]) -> Self {
        Self {
            particles,
            rho_f,
            nu,
            g,
            restitution: 0.0,
        }
    }

    /// Sets the wall restitution coefficient used for blocked components.
    pub fn with_restitution(mut self, restitution: f64) -> Self {
        self.restitution = restitution;
        self
    }

    /// Advances all particles by one lattice time step.
    ///
    /// `sample` supplies the fluid velocity and solid mask at arbitrary
    /// positions. `exposure_rate`, when present, supplies a deterministic local
    /// rate accumulated as `exposure += rate(pos) * dt` with `dt = 1`.
    pub fn step<F, E>(&mut self, sample: F, exposure_rate: Option<E>)
    where
        F: Fn([f64; 3]) -> Sample,
        E: Fn([f64; 3]) -> f64,
    {
        assert!(self.rho_f > 0.0, "fluid density must be positive");
        assert!(self.nu > 0.0, "kinematic viscosity must be positive");
        assert!(
            (0.0..=1.0).contains(&self.restitution),
            "restitution must be in [0, 1]"
        );

        for p in &mut self.particles {
            assert!(p.d > 0.0, "particle diameter must be positive");
            assert!(p.rho_p > 0.0, "particle density must be positive");

            let s = sample(p.pos);
            if let Some(rate) = &exposure_rate {
                p.exposure += rate(p.pos);
            }

            let v_new = particle_velocity(p.vel, s.u, p.d, p.rho_p, self.rho_f, self.nu, self.g);
            let pos_new = add(p.pos, v_new);

            if sample(pos_new).solid {
                let (pos, vel) =
                    resolve_solid_contact(p.pos, pos_new, v_new, self.restitution, &sample);
                p.pos = pos;
                p.vel = vel;
            } else {
                p.pos = pos_new;
                p.vel = v_new;
            }
        }
    }
}

impl Default for ParticleSet {
    fn default() -> Self {
        Self {
            particles: Vec::new(),
            rho_f: 1.0,
            nu: 1.0 / 6.0,
            g: [0.0; 3],
            restitution: 0.0,
        }
    }
}

fn particle_velocity(
    v: [f64; 3],
    u: [f64; 3],
    d: f64,
    rho_p: f64,
    rho_f: f64,
    nu: f64,
    g: [f64; 3],
) -> [f64; 3] {
    let slip = sub(u, v);
    let re = norm(slip) * d / nu;
    debug_assert!(
        re < 800.0,
        "Schiller-Naumann drag is used outside its Re_p < 800 range: {re}"
    );
    let re = re.min(800.0);
    let drag_correction = 1.0 + 0.15 * re.powf(0.687);
    let tau_p = rho_p * d * d / (18.0 * rho_f * nu * drag_correction);
    let g_eff = scale(g, 1.0 - rho_f / rho_p);

    let mut out = [0.0; 3];
    for a in 0..3 {
        out[a] = (tau_p * v[a] + u[a] + tau_p * g_eff[a]) / (tau_p + 1.0);
    }
    out
}

fn resolve_solid_contact<F>(
    start: [f64; 3],
    end: [f64; 3],
    vel: [f64; 3],
    restitution: f64,
    sample: &F,
) -> ([f64; 3], [f64; 3])
where
    F: Fn([f64; 3]) -> Sample,
{
    let delta = sub(end, start);
    let n_sub = max_abs(delta).ceil().max(1.0) as usize;
    let mut prev = start;

    for i in 1..=n_sub {
        let t = i as f64 / n_sub as f64;
        let cand = add(start, scale(delta, t));
        if sample(cand).solid {
            let axis = blocked_axis(prev, cand, vel, sample);
            let mut corrected = cand;
            corrected[axis] = prev[axis];
            if sample(corrected).solid {
                corrected = prev;
            }

            let mut reflected = vel;
            reflected[axis] = -restitution * reflected[axis];
            return (corrected, reflected);
        }
        prev = cand;
    }

    let axis = largest_abs_axis(delta);
    let mut reflected = vel;
    reflected[axis] = -restitution * reflected[axis];
    (prev, reflected)
}

fn blocked_axis<F>(prev: [f64; 3], cand: [f64; 3], vel: [f64; 3], sample: &F) -> usize
where
    F: Fn([f64; 3]) -> Sample,
{
    let mut best = None;
    for axis in 0..3 {
        if (cand[axis] - prev[axis]).abs() == 0.0 {
            continue;
        }
        let mut one_axis = prev;
        one_axis[axis] = cand[axis];
        if sample(one_axis).solid {
            let mag = vel[axis].abs();
            if best.map_or(true, |(_, best_mag)| mag > best_mag) {
                best = Some((axis, mag));
            }
        }
    }
    best.map(|(axis, _)| axis)
        .unwrap_or_else(|| largest_abs_axis(sub(cand, prev)))
}

/// Samples a regular grid with trilinear interpolation.
///
/// `dims` is `[nx, ny, nz]`; positions are lattice-node coordinates and are
/// clamped to the grid bounds. When `nz == 1`, interpolation is bilinear in
/// `x,y`. The accessor is engine-free and returns a node velocity plus a solid
/// flag. Solid neighbors contribute `u = 0` with their normal interpolation
/// weight, matching the half-way wall convention used by callers. The returned
/// `solid` flag is the flag of the containing clamped lower node, intended for
/// contact tests rather than volume-fraction interpolation.
pub fn sample_grid<F>(pos: [f64; 3], dims: [usize; 3], accessor: F) -> Sample
where
    F: Fn(usize, usize, usize) -> ([f64; 3], bool),
{
    assert!(
        dims[0] > 0 && dims[1] > 0 && dims[2] > 0,
        "grid dimensions must be nonzero"
    );

    let (x0, x1, tx) = bracket(pos[0], dims[0]);
    let (y0, y1, ty) = bracket(pos[1], dims[1]);
    let (z0, z1, tz) = if dims[2] == 1 {
        (0, 0, 0.0)
    } else {
        bracket(pos[2], dims[2])
    };

    let mut u = [0.0; 3];
    for (ix, wx) in [(x0, 1.0 - tx), (x1, tx)] {
        for (iy, wy) in [(y0, 1.0 - ty), (y1, ty)] {
            for (iz, wz) in [(z0, 1.0 - tz), (z1, tz)] {
                let w = wx * wy * wz;
                if w == 0.0 {
                    continue;
                }
                let (node_u, solid) = accessor(ix, iy, iz);
                if !solid {
                    for a in 0..3 {
                        u[a] += w * node_u[a];
                    }
                }
            }
        }
    }

    let (_, solid) = accessor(x0, y0, z0);
    Sample { u, solid }
}

fn bracket(x: f64, n: usize) -> (usize, usize, f64) {
    if n == 1 {
        return (0, 0, 0.0);
    }
    let x = x.clamp(0.0, (n - 1) as f64);
    let lo = x.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    (lo, hi, x - lo as f64)
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn max_abs(a: [f64; 3]) -> f64 {
    a[0].abs().max(a[1].abs()).max(a[2].abs())
}

fn largest_abs_axis(a: [f64; 3]) -> usize {
    if a[1].abs() > a[0].abs() && a[1].abs() >= a[2].abs() {
        1
    } else if a[2].abs() > a[0].abs() && a[2].abs() > a[1].abs() {
        2
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fluid(u: [f64; 3]) -> impl Fn([f64; 3]) -> Sample {
        move |_| Sample { u, solid: false }
    }

    #[test]
    fn quiescent_settling_matches_implicit_terminal_velocity() {
        let rho_f = 1.0;
        let rho_p = 1.05;
        let d = 4.0;
        let nu = 0.1;
        let g = [0.0, -1e-5, 0.0];
        let mut set = ParticleSet::new(
            vec![Particle {
                pos: [0.0, 100.0, 0.0],
                vel: [0.0; 3],
                d,
                rho_p,
                exposure: 0.0,
            }],
            rho_f,
            nu,
            g,
        );

        for _ in 0..200_000 {
            set.step(fluid([0.0; 3]), None::<fn([f64; 3]) -> f64>);
        }

        let g_eff = (1.0_f64 - rho_f / rho_p) * g[1].abs();
        let mut vt = 0.0_f64;
        for _ in 0..10_000 {
            let re = vt * d / nu;
            let f = 1.0 + 0.15 * re.min(800.0).powf(0.687);
            let tau_p = rho_p * d * d / (18.0 * rho_f * nu * f);
            let next = tau_p * g_eff;
            if (next - vt).abs() < 1e-15 {
                vt = next;
                break;
            }
            vt = next;
        }

        let got = -set.particles[0].vel[1];
        assert!(
            (got - vt).abs() / vt < 1e-6,
            "terminal velocity got {got:e}, expected {vt:e}"
        );
    }

    #[test]
    fn tracer_limit_follows_step_change_in_one_step() {
        let mut set = ParticleSet::new(
            vec![Particle {
                pos: [0.0; 3],
                vel: [0.0; 3],
                d: 1e-8,
                rho_p: 1.0,
                exposure: 0.0,
            }],
            1.0,
            0.1,
            [0.0; 3],
        );

        let u = [0.2, -0.1, 0.05];
        set.step(fluid(u), None::<fn([f64; 3]) -> f64>);
        for (got, want) in set.particles[0].vel.iter().zip(u) {
            assert!((*got - want).abs() < 1e-12, "got {got:e}, want {want:e}");
        }
    }

    #[test]
    fn stiff_response_stays_bounded() {
        let rho_f = 1.0;
        let rho_p = 1.0;
        let nu = 1.0 / 18.0;
        let d = 1e-3;
        let tau_p = rho_p * d * d / (18.0 * rho_f * nu);
        assert!((tau_p - 1e-6_f64).abs() < 1e-18_f64);

        let u = [0.4, 0.0, 0.0];
        let g = [0.0, -0.2, 0.0];
        let mut set = ParticleSet::new(
            vec![Particle {
                pos: [0.0; 3],
                vel: [0.0; 3],
                d,
                rho_p,
                exposure: 0.0,
            }],
            rho_f,
            nu,
            g,
        );
        set.step(fluid(u), None::<fn([f64; 3]) -> f64>);

        let speed = norm(set.particles[0].vel);
        let bound = norm(u) + norm(g) * tau_p + 1e-12;
        assert!(speed <= bound, "speed {speed:e} > bound {bound:e}");
    }

    #[test]
    fn restitution_zero_rests_on_floor_and_later_moves_by_drag() {
        let floor = |p: [f64; 3], u: [f64; 3]| Sample {
            u,
            solid: p[1] < 0.0,
        };
        let mut set = ParticleSet::new(
            vec![Particle {
                pos: [0.0, 0.0, 0.0],
                vel: [0.0; 3],
                d: 1.0,
                rho_p: 2.0,
                exposure: 0.0,
            }],
            1.0,
            0.1,
            [0.0, -0.01, 0.0],
        );

        for _ in 0..1000 {
            set.step(|p| floor(p, [0.0; 3]), None::<fn([f64; 3]) -> f64>);
            assert!(
                set.particles[0].pos[1] >= 0.0,
                "particle tunneled through floor"
            );
            assert!(
                set.particles[0].vel[1].abs() <= 1e-15,
                "normal velocity did not rest"
            );
        }
        let y = set.particles[0].pos[1];
        let x = set.particles[0].pos[0];
        for _ in 0..20 {
            set.step(|p| floor(p, [0.1, 0.0, 0.0]), None::<fn([f64; 3]) -> f64>);
        }
        assert!((set.particles[0].pos[1] - y).abs() < 1e-12);
        assert!(
            set.particles[0].pos[0] > x,
            "horizontal drag did not move rested particle"
        );
    }

    #[test]
    fn interpolation_reproduces_linear_fields() {
        let linear = |x: usize, y: usize, z: usize| {
            let x = x as f64;
            let y = y as f64;
            let z = z as f64;
            (
                [
                    1.0 + 2.0 * x - 3.0 * y + 0.5 * z,
                    -2.0 + x + 4.0 * z,
                    7.0 - y + z,
                ],
                false,
            )
        };
        let p = [1.25, 2.5, 3.75];
        let s = sample_grid(p, [5, 6, 7], linear);
        let want = [
            1.0 + 2.0 * p[0] - 3.0 * p[1] + 0.5 * p[2],
            -2.0 + p[0] + 4.0 * p[2],
            7.0 - p[1] + p[2],
        ];
        for (got, want) in s.u.iter().zip(want) {
            assert!(
                (*got - want).abs() < 1e-14,
                "trilinear got {got:e}, want {want:e}"
            );
        }

        let p = [2.25, 1.5, 0.0];
        let s = sample_grid(p, [5, 6, 1], linear);
        let want = [1.0 + 2.0 * p[0] - 3.0 * p[1], -2.0 + p[0], 7.0 - p[1]];
        for (got, want) in s.u.iter().zip(want) {
            assert!(
                (*got - want).abs() < 1e-14,
                "bilinear got {got:e}, want {want:e}"
            );
        }
    }

    #[test]
    fn fast_particle_does_not_tunnel_through_one_cell_wall() {
        let wall = |p: [f64; 3]| Sample {
            u: [0.0; 3],
            solid: (1.0..2.0).contains(&p[0]),
        };
        let mut set = ParticleSet::new(
            vec![Particle {
                pos: [0.0, 0.0, 0.0],
                vel: [2.0, 0.0, 0.0],
                d: 1000.0,
                rho_p: 1.0,
                exposure: 0.0,
            }],
            1.0,
            0.1,
            [0.0; 3],
        );
        set.step(wall, None::<fn([f64; 3]) -> f64>);
        assert!(
            set.particles[0].pos[0] < 1.0,
            "particle tunneled to x={}",
            set.particles[0].pos[0]
        );
    }

    #[test]
    fn exposure_constant_rate_is_exact() {
        let mut set = ParticleSet::new(
            vec![Particle {
                pos: [0.0; 3],
                vel: [0.0; 3],
                d: 1.0,
                rho_p: 1.0,
                exposure: 0.0,
            }],
            1.0,
            0.1,
            [0.0; 3],
        );
        for _ in 0..37 {
            set.step(fluid([0.0; 3]), Some(|_| 2.5));
        }
        assert_eq!(set.particles[0].exposure, 92.5);
    }
}
