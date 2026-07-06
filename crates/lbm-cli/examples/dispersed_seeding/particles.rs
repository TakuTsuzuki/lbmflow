use crate::protocol::{ProtocolInput, Regime};
use crate::Sim3;
use lbm_core::fields::SoaFields;
use lbm_core::particles::{
    sample_grid, DepositEvent, Particle as CoreParticle, ParticleSet, Sample as CoreSample,
};

const G: f64 = 9.80665;

#[derive(Clone, Debug)]
pub struct Particle {
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
    let n = ((input.particles.count as f64) * input.reservoir.initial_conc).round() as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let d = lognormal_diameter(input.particles.d_p_m, input.particles.d_p_cv, &mut rng);
        out.push(Particle { d_m: d });
    }
    out
}

pub fn deposit_batch(
    input: &ProtocolInput,
    regime: &Regime,
    batch: Vec<Particle>,
    tray_sim: &mut Sim3,
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
    let nozzle_radius = 0.5 * regime.nozzle_d_m;
    let particles = batch
        .iter()
        .enumerate()
        .map(|(k, p)| {
            let pt = points[k % points.len()];
            let r = nozzle_radius * rng.unit().sqrt();
            let theta = 2.0 * std::f64::consts::PI * rng.unit();
            // Frozen §3.2 defines each eject point as a nozzle disk. Sampling
            // r = R*sqrt(U), theta = 2πV gives a uniform area distribution
            // over that disk; no extra jitter or case-dependent spread is used.
            let pos_m = [
                pt[0] * input.target.width_m + r * theta.cos(),
                pt[1] * input.target.depth_m + r * theta.sin(),
                h0,
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

    let deposits = step_particles_to_deposition(input, regime, tray_sim, &mut set)?;
    Ok((set, deposits))
}

fn step_particles_to_deposition(
    input: &ProtocolInput,
    regime: &Regime,
    tray_sim: &mut Sim3,
    set: &mut ParticleSet,
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
    let mut deposits = Vec::new();
    let settling_lu = regime.settling_m_s * regime.dt / regime.dx;
    let quiet_threshold = 0.1 * settling_lu;
    let mut flow_frozen = false;
    for step in 0..steps {
        let t = step as f64 * dt;
        let ac = if step < agitation_steps {
            input.agitation_amplitude_m() * omega * omega * (omega * t).sin()
        } else {
            0.0
        };
        let accel_lu = ac * regime.dt * regime.dt / regime.dx;
        tray_sim.set_gravity([-accel_lu, 0.0, 0.0]);
        if !flow_frozen {
            tray_sim.step();
        }
        let fields = tray_sim.backend_fields(0);
        if step >= agitation_steps && step % 128 == 0 {
            let max_u = max_tray_speed(fields);
            // Numerical guard only: once the resolved tray velocity is below
            // one tenth of the Stokes settling speed, further fluid evolution
            // is negligible for particle transport over this example's late
            // quiescent tail. Particles still continue settling under gravity.
            if max_u < quiet_threshold {
                flow_frozen = true;
            }
        }
        // In the translating non-inertial frame, the fluid receives the
        // uniform pseudo-acceleration -a_frame. ParticleSet applies the same
        // body acceleration through the gravity path, which multiplies by
        // (1-rho_f/rho_p); dividing by that buoyancy factor here yields the
        // required density-weighted pseudo-force on the particles.
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
            |pos| sample_tray(input, fields, pos),
            None::<fn([f64; 3]) -> f64>,
            0.0,
            &mut deposits,
        );
        if set.particles.is_empty() {
            break;
        }
    }
    Ok(deposits)
}

fn sample_tray(input: &ProtocolInput, fields: &SoaFields<f64>, pos_lu: [f64; 3]) -> CoreSample {
    let nx = input.grid.tray_nx;
    let ny = input.grid.tray_ny;
    let nz = input.grid.tray_nz;
    sample_grid(pos_lu, [nx, ny, nz], |x, y, z| {
        let c = fields.geom.cidx(x, y, z);
        let p = fields.geom.pidx(x, y, z);
        ([fields.ux[c], fields.uy[c], fields.uz[c]], fields.solid[p])
    })
}

fn max_tray_speed(fields: &SoaFields<f64>) -> f64 {
    (0..fields.ux.len())
        .map(|i| {
            (fields.ux[i] * fields.ux[i]
                + fields.uy[i] * fields.uy[i]
                + fields.uz[i] * fields.uz[i])
                .sqrt()
        })
        .fold(0.0, f64::max)
}
