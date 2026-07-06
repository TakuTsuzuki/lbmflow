use crate::protocol::{ProtocolInput, Regime};
use lbm_core::particles::{
    DepositEvent, Particle as CoreParticle, ParticleSet, Sample as CoreSample,
};

const G: f64 = 9.80665;

#[derive(Clone, Debug)]
pub struct Particle {
    pub pos: [f64; 3],
    pub d_m: f64,
}

#[derive(Clone, Debug)]
pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self { state: seed | 1 }
    }

    pub fn unit(&mut self) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let bits = (self.state >> 11) | 0x3ff0_0000_0000_0000;
        f64::from_bits(bits) - 1.0
    }

    pub fn normal(&mut self) -> f64 {
        let u1 = self.unit().clamp(1.0e-12, 1.0);
        let u2 = self.unit();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

pub fn lognormal_diameter(mean: f64, cv: f64, rng: &mut Rng) -> f64 {
    if cv <= 0.0 {
        return mean;
    }
    let sigma2 = (1.0 + cv * cv).ln();
    let sigma = sigma2.sqrt();
    let mu = mean.ln() - 0.5 * sigma2;
    (mu + sigma * rng.normal()).exp()
}

pub fn make_reservoir_particles(input: &ProtocolInput) -> Vec<Particle> {
    let mut rng = Rng::new(input.particles.seed);
    let fill = (input.reservoir.fill_height_m / input.reservoir.height_m).clamp(0.0, 1.0);
    let n = ((input.particles.count as f64) * input.reservoir.initial_conc.clamp(0.0, 1.0)).round()
        as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let d = lognormal_diameter(input.particles.d_p_m, input.particles.d_p_cv, &mut rng);
        let x = rng.unit() * input.reservoir.width_m;
        let y = rng.unit() * input.reservoir.width_m;
        let z_frac = fill * rng.unit();
        out.push(Particle {
            pos: [x, y, z_frac * input.reservoir.height_m],
            d_m: d,
        });
    }
    out
}

pub fn deposit_batch(
    input: &ProtocolInput,
    regime: &Regime,
    batch: Vec<Particle>,
    tray_velocity: &[[f64; 3]],
) -> anyhow::Result<(ParticleSet, Vec<DepositEvent>)> {
    if batch.is_empty() {
        return Ok((ParticleSet::default(), Vec::new()));
    }
    let eject = input.op("eject").expect("eject operation");
    let points = eject
        .points_xy_frac
        .as_ref()
        .cloned()
        .unwrap_or_else(|| vec![[0.5, 0.5]]);
    let h0 = eject.height_m.unwrap_or(input.target.height_m);
    let mut rng = Rng::new(input.particles.seed ^ 0x9e37_79b9_7f4a_7c15);
    let jet_sigma = (0.5 * regime.nozzle_d_m).max(1.25 * regime.dx);
    let particles = batch
        .iter()
        .enumerate()
        .map(|(k, p)| {
            let pt = points[k % points.len()];
            let pos_m = [
                (pt[0] * input.target.width_m + jet_sigma * 0.35 * rng.normal())
                    .clamp(0.0, input.target.width_m),
                (pt[1] * input.target.depth_m + jet_sigma * 0.35 * rng.normal())
                    .clamp(0.0, input.target.depth_m),
                h0.clamp(0.0, input.target.height_m),
            ];
            CoreParticle {
                pos: [
                    pos_m[0] / regime.dx,
                    pos_m[1] / regime.dx,
                    pos_m[2] / regime.dx,
                ],
                vel: [0.0, 0.0, -regime.u_jet_lattice],
                d: p.d_m / regime.dx,
                rho_p: input.particles.rho_p_kgm3,
                exposure: 0.0,
            }
        })
        .collect::<Vec<_>>();

    let g_lu = G * regime.dt * regime.dt / regime.dx;
    let mut set = ParticleSet::new(
        particles,
        input.fluid.rho_f_kgm3,
        regime.nu_lattice,
        [0.0, 0.0, -g_lu],
    );

    let deposits = step_particles_to_deposition(
        input,
        regime,
        &points,
        jet_sigma,
        tray_velocity,
        &mut set,
        rng,
    )?;
    Ok((set, deposits))
}

fn step_particles_to_deposition(
    input: &ProtocolInput,
    regime: &Regime,
    points: &[[f64; 2]],
    jet_sigma: f64,
    tray_velocity: &[[f64; 3]],
    set: &mut ParticleSet,
    mut rng: Rng,
) -> anyhow::Result<Vec<DepositEvent>> {
    let dt = regime.dt;
    let max_duration = input
        .ops("settle")
        .last()
        .and_then(|op| op.duration_s)
        .unwrap_or(2.0)
        .max(0.1);
    let requested_steps = ((max_duration / dt).ceil() as usize).max(20);
    let max_particle_steps = input.max_particle_steps.unwrap_or(50_000);
    if requested_steps > max_particle_steps {
        anyhow::bail!(
            "particle integration requires {requested_steps} steps, exceeding max_particle_steps={max_particle_steps}; suspended particles are not added to the density map"
        );
    }
    let steps = requested_steps;
    let omega = {
        let a = input.agitation_amplitude_m();
        let s = input.agitation_speed_m_s();
        if a > 0.0 {
            s / a
        } else {
            0.0
        }
    };
    let agitation_steps = if omega > 0.0 {
        ((input.agitation_count() as f64 * 2.0 * std::f64::consts::PI / omega) / dt) as usize
    } else {
        0
    };
    let buoyancy = 1.0 - input.fluid.rho_f_kgm3 / input.particles.rho_p_kgm3;
    let harshness = (regime.u_jet_si / 0.01).clamp(0.0, 20.0)
        + input.agitation_count() as f64 * 0.15
        + input.fr();
    let gentle_k = (2.0 * 2.5e-5 * regime.dt).sqrt() / regime.dx;
    let mut deposits = Vec::new();
    for step in 0..steps {
        let t = step as f64 * dt;
        let ac = if step < agitation_steps {
            input.agitation_amplitude_m() * omega * omega * (omega * t).sin()
        } else {
            0.0
        };
        set.g = [
            if buoyancy > 0.0 {
                -ac * regime.dt * regime.dt / regime.dx / buoyancy
            } else {
                0.0
            },
            0.0,
            -G * regime.dt * regime.dt / regime.dx,
        ];
        set.step_depositing(
            |pos| {
                sample_tray(
                    input,
                    regime,
                    pos,
                    points,
                    jet_sigma,
                    tray_velocity,
                    harshness,
                )
            },
            None::<fn([f64; 3]) -> f64>,
            0.0,
            &mut deposits,
        );
        for p in &mut set.particles {
            if harshness > 4.0 {
                let p_d_m = p.d * regime.dx;
                p.pos[0] +=
                    0.15 * jet_sigma * (11.0 * t + p_d_m * 1.0e6).sin() * regime.dt / regime.dx;
            } else {
                p.pos[0] += gentle_k * rng.normal();
                p.pos[1] += gentle_k * rng.normal();
            }
            p.pos[0] = p.pos[0].clamp(0.0, input.target.width_m / regime.dx);
            p.pos[1] = p.pos[1].clamp(0.0, input.target.depth_m / regime.dx);
            if p.pos[2] > input.target.height_m / regime.dx {
                p.pos[2] = input.target.height_m / regime.dx;
                p.vel[2] = p.vel[2].min(0.0);
            }
        }
        if set.particles.is_empty() {
            break;
        }
    }
    Ok(deposits)
}

fn sample_tray(
    input: &ProtocolInput,
    regime: &Regime,
    pos_lu: [f64; 3],
    points: &[[f64; 2]],
    sigma: f64,
    tray_velocity: &[[f64; 3]],
    harshness: f64,
) -> CoreSample {
    let pos = [
        pos_lu[0] * regime.dx,
        pos_lu[1] * regime.dx,
        pos_lu[2] * regime.dx,
    ];
    let nx = input.grid.tray_nx;
    let ny = input.grid.tray_ny;
    let nz = input.grid.tray_nz;
    let ix = ((pos[0] / input.target.width_m) * (nx as f64 - 1.0)).round() as usize;
    let iy = ((pos[1] / input.target.depth_m) * (ny as f64 - 1.0)).round() as usize;
    let iz = ((pos[2] / input.target.height_m) * (nz as f64 - 1.0)).round() as usize;
    let idx = (iz.min(nz - 1) * ny + iy.min(ny - 1)) * nx + ix.min(nx - 1);
    let mut u = tray_velocity.get(idx).copied().unwrap_or([0.0; 3]);
    let wall_jet_len = if harshness > 4.0 {
        0.10 * input.target.width_m.min(input.target.depth_m)
    } else {
        0.42 * input.target.width_m.min(input.target.depth_m)
    };
    for pt in points {
        let cx = pt[0] * input.target.width_m;
        let cy = pt[1] * input.target.depth_m;
        let dx = pos[0] - cx;
        let dy = pos[1] - cy;
        let r2 = dx * dx + dy * dy;
        let g = (-0.5 * r2 / (sigma * sigma).max(1.0e-12)).exp();
        let z_decay = (pos[2] / input.target.height_m).clamp(0.0, 1.0);
        u[2] += -regime.u_jet_si * g * (0.25 + 0.75 * z_decay);
        let r = r2.sqrt().max(1.0e-9);
        let wall_jet = (-r / wall_jet_len.max(1.0e-9)).exp();
        u[2] += -0.35 * regime.u_jet_si * wall_jet * z_decay;
        let radial = 0.35 * regime.u_jet_si * wall_jet * (1.0 - z_decay).powi(2);
        u[0] += radial * dx / r;
        u[1] += radial * dy / r;
    }
    CoreSample {
        u: [
            u[0] * regime.dt / regime.dx,
            u[1] * regime.dt / regime.dx,
            u[2] * regime.dt / regime.dx,
        ],
        solid: false,
    }
}
