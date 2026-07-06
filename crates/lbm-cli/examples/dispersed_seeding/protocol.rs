use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Deserialize)]
pub struct ProtocolInput {
    pub grid: Grid,
    pub fluid: Fluid,
    pub particles: ParticleSpec,
    pub reservoir: ReservoirSpec,
    pub target: TargetSpec,
    pub protocol: Vec<Operation>,
    pub output: OutputSpec,
    #[serde(default)]
    pub max_particle_steps: Option<usize>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct Grid {
    pub res_nx: usize,
    pub res_ny: usize,
    pub res_nz: usize,
    pub tray_nx: usize,
    pub tray_ny: usize,
    pub tray_nz: usize,
    pub dx_m: f64,
    #[serde(default)]
    pub tray_dx_m: Option<f64>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct Fluid {
    pub nu_m2s: f64,
    pub rho_f_kgm3: f64,
}

#[derive(Clone, Debug, Deserialize)]
pub struct ParticleSpec {
    pub rho_p_kgm3: f64,
    pub d_p_m: f64,
    pub d_p_cv: f64,
    pub count: usize,
    pub seed: u64,
}

#[derive(Clone, Debug, Deserialize)]
pub struct ReservoirSpec {
    pub height_m: f64,
    pub width_m: f64,
    pub fill_height_m: f64,
    pub initial_conc: f64,
}

#[derive(Clone, Debug, Deserialize)]
pub struct TargetSpec {
    pub width_m: f64,
    pub depth_m: f64,
    pub height_m: f64,
    pub partitions_x: usize,
    pub partitions_y: usize,
}

#[derive(Clone, Debug, Deserialize)]
pub struct Operation {
    pub op: String,
    pub duration_s: Option<f64>,
    pub volume_frac: Option<f64>,
    #[serde(rename = "rate_uLs")]
    pub rate_ul_s: Option<f64>,
    pub depth_frac: Option<f64>,
    pub points_xy_frac: Option<Vec<[f64; 2]>>,
    pub nozzle_diameter_m: Option<f64>,
    pub height_m: Option<f64>,
    pub pattern: Option<String>,
    pub count: Option<usize>,
    pub speed_mms: Option<f64>,
    pub amplitude_mm: Option<f64>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct OutputSpec {
    pub dir: String,
    pub csv: bool,
    pub volume: bool,
}

#[derive(Clone, Debug)]
pub struct Regime {
    pub dx: f64,
    pub dt: f64,
    pub nu_lattice: f64,
    pub u_jet_lattice: f64,
    pub ma: f64,
    pub re_jet: f64,
    pub st: f64,
    pub fr: f64,
    pub tau: f64,
    pub nozzle_d_m: f64,
    pub particle_tau_s: f64,
    pub settling_m_s: f64,
}

impl ProtocolInput {
    pub fn from_path(path: &Path) -> anyhow::Result<Self> {
        let text = fs::read_to_string(path)?;
        Ok(serde_json::from_str(&text)?)
    }

    pub fn output_dir(&self) -> PathBuf {
        PathBuf::from(&self.output.dir)
    }

    pub fn op(&self, name: &str) -> Option<&Operation> {
        self.protocol.iter().find(|op| op.op == name)
    }

    pub fn ops<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Operation> + 'a {
        self.protocol.iter().filter(move |op| op.op == name)
    }

    pub fn regime(&self) -> anyhow::Result<Regime> {
        let dx = self.grid.tray_dx_m.unwrap_or(self.grid.dx_m);
        let eject = self
            .op("eject")
            .ok_or_else(|| anyhow::anyhow!("protocol requires an eject operation"))?;
        let points = eject.points_xy_frac.as_ref().map_or(1, Vec::len).max(1);
        let nozzle_d_m = eject
            .nozzle_diameter_m
            .ok_or_else(|| anyhow::anyhow!("eject.nozzle_diameter_m is required"))?;
        if nozzle_d_m <= 0.0 {
            anyhow::bail!("eject.nozzle_diameter_m must be positive");
        }
        let rate_m3_s = eject
            .rate_ul_s
            .ok_or_else(|| anyhow::anyhow!("eject.rate_uLs is required"))?
            * 1.0e-9;
        let patch_area = std::f64::consts::PI * (0.5 * nozzle_d_m).powi(2) * points as f64;
        let u_jet_si = if rate_m3_s > 0.0 {
            rate_m3_s / patch_area
        } else {
            0.0
        };
        let max_si = u_jet_si.max(self.settling_velocity_m_s());
        let dt_by_ma = if max_si > 0.0 {
            // Low-Mach envelope: Ma = u*/c_s = u_si*dt/(dx*c_s). The frozen
            // target is Ma <= 0.3, and 0.16 < 0.3/sqrt(3) keeps the selected
            // advective velocity inside that hard guard with margin.
            0.16 * dx / max_si
        } else {
            1.0
        };
        // Diffusive scaling: nu* = nu_phys*dt/dx^2 and tau = 3*nu* + 0.5.
        // The 0.012 coefficient gives tau=0.536 for the water/tray envelope,
        // above the frozen tau >= 0.51 guard while keeping Ma in the low band.
        let dt_by_tau = 0.012 * dx * dx / self.fluid.nu_m2s;
        let dt = dt_by_ma.min(dt_by_tau);
        let nu_lattice = self.fluid.nu_m2s * dt / (dx * dx);
        let tau = 3.0 * nu_lattice + 0.5;
        let u_jet_lattice = u_jet_si * dt / dx;
        let cs = (1.0 / 3.0f64).sqrt();
        let ma = u_jet_lattice / cs;
        let re_jet = if self.fluid.nu_m2s > 0.0 {
            u_jet_si * nozzle_d_m / self.fluid.nu_m2s
        } else {
            f64::NAN
        };
        let particle_tau_s = self.particles.rho_p_kgm3 * self.particles.d_p_m.powi(2)
            / (18.0 * self.fluid.rho_f_kgm3 * self.fluid.nu_m2s);
        let st = if nozzle_d_m > 0.0 {
            particle_tau_s * u_jet_si / nozzle_d_m
        } else {
            0.0
        };
        let fr = self.fr();
        let settling_m_s = self.settling_velocity_m_s();
        if ma > 0.3 {
            anyhow::bail!("Mach guard failed: Ma={ma:.3} > 0.3");
        }
        if tau < 0.51 {
            anyhow::bail!("relaxation guard failed: tau={tau:.5} < 0.51");
        }
        Ok(Regime {
            dx,
            dt,
            nu_lattice,
            u_jet_lattice,
            ma,
            re_jet,
            st,
            fr,
            tau,
            nozzle_d_m,
            particle_tau_s,
            settling_m_s,
        })
    }

    pub fn agitation_speed_m_s(&self) -> f64 {
        if self.agitation_count() == 0 {
            return 0.0;
        }
        self.op("agitate")
            .and_then(|op| op.speed_mms)
            .unwrap_or(0.0)
            * 1.0e-3
    }

    pub fn agitation_amplitude_m(&self) -> f64 {
        self.op("agitate")
            .and_then(|op| op.amplitude_mm)
            .unwrap_or(0.0)
            * 1.0e-3
    }

    pub fn agitation_count(&self) -> usize {
        self.op("agitate").and_then(|op| op.count).unwrap_or(0)
    }

    pub fn fr(&self) -> f64 {
        let a = self.agitation_amplitude_m();
        let speed = self.agitation_speed_m_s();
        if a <= 0.0 {
            0.0
        } else {
            speed * speed / (a * 9.80665)
        }
    }

    pub fn settling_velocity_m_s(&self) -> f64 {
        let mu = self.fluid.rho_f_kgm3 * self.fluid.nu_m2s;
        (self.particles.rho_p_kgm3 - self.fluid.rho_f_kgm3) * 9.80665 * self.particles.d_p_m.powi(2)
            / (18.0 * mu)
    }
}
