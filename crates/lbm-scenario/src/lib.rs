//! Scenario JSON (v0): the single execution contract shared by the CLI,
//! the MCP server and (in spirit) the GUI presets.
//!
//! See `docs/AGENT_MODE_DESIGN.md` for the schema rationale. Field names are
//! camelCase in JSON.

use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use serde::{Deserialize, Serialize};

mod units;
pub use units::{
    report as unit_report, resolve, ConversionFactors, DimensionlessNumbers, FlowParams,
    LatticeUnits, UnitConstructor, UnitDiagnostic, UnitInputsEcho, UnitReport, UnitSuggestion,
    UnitVerdict, GRID_RE_WARN_THRESHOLD, LATTICE_SPEED_WARN_THRESHOLD, TAU_HIGH_WARN_THRESHOLD,
    TAU_LOW_WARN_THRESHOLD,
};

// ---------------------------------------------------------------- schema

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct Scenario {
    #[serde(default)]
    pub version: u32,
    pub name: String,
    pub grid: Grid,
    pub physics: Physics,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub units: Option<FlowParams>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub compute: Option<ComputeSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall: Option<WallModel>,
    pub edges: EdgesSpec,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inlet_profile: Option<InletProfile>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub obstacles: Vec<Obstacle>,
    #[serde(default)]
    pub init: InitSpec,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub multiphase: Option<MultiphaseSpec>,
    /// Rotating-impeller volume penalization (MF-delta interim; 2D only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rotor: Option<RotorSpec>,
    /// One-way Lagrangian particles (FR-PART-01 subset; 2D only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub particles: Option<ParticlesSpec>,
    pub run: RunSpec,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub probes: Vec<ProbeSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub outputs: Vec<OutputSpec>,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Grid {
    pub nx: usize,
    pub ny: usize,
    /// Cells along z. Omitted or 1 = 2D (D2Q9); > 1 = 3D (D3Q19, runs on
    /// the V2 core). Not serialised for 2D scenarios, so existing files
    /// round-trip byte-identically.
    #[serde(default = "default_nz", skip_serializing_if = "is_default_nz")]
    pub nz: usize,
}

fn default_nz() -> usize {
    1
}

#[allow(clippy::trivially_copy_pass_by_ref)]
fn is_default_nz(nz: &usize) -> bool {
    *nz == 1
}

impl Scenario {
    /// Whether this scenario runs on the 3D (D3Q19) engine.
    pub fn is_3d(&self) -> bool {
        self.grid.nz > 1
    }
}

/// Compute-target selection (ARCHITECTURE_V2 §3). All fields optional.
/// Explicit backend requests must be honored or rejected; only `auto` may
/// silently choose an available fallback.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct ComputeSpec {
    #[serde(default)]
    pub backend: BackendSpec,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum BackendSpec {
    #[default]
    Auto,
    Cpu,
    Gpu,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum WallModel {
    Bouzidi,
}

const GPU_FEATURE_UNAVAILABLE: &str = "requested backend \"gpu\" is unavailable: this binary was built without GPU scenario dispatch (--features gpu). Use compute.backend \"cpu\" or rebuild with --features gpu.";
const GPU_F64_UNSUPPORTED: &str = "requested backend \"gpu\" is unsupported for physics.precision f64; GPU scenario dispatch currently supports f32 storage/compute only. Use physics.precision \"f32\" or compute.backend \"cpu\".";
const GPU_3D_UNSUPPORTED: &str = "requested backend \"gpu\" is unsupported for this 3D scenario: unsupported 3D GPU scenario combinations are f64 precision, multiphase, rotor, particles, non-rest init, force probes, and non-f32 storage. Use compute.backend \"cpu\" or simplify the scenario.";

/// Backend selected after applying explicit/auto policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "camelCase")]
pub enum BackendChoice {
    Cpu,
    Gpu,
}

pub fn requested_backend(sc: &Scenario) -> BackendSpec {
    sc.compute.map(|c| c.backend).unwrap_or_default()
}

pub fn auto_gpu_threshold(sc: &Scenario) -> bool {
    if sc.is_3d() {
        sc.grid.nx.saturating_mul(sc.grid.ny).saturating_mul(sc.grid.nz) >= 64usize.pow(3)
    } else {
        sc.grid.nx.saturating_mul(sc.grid.ny) >= 256usize.pow(2)
    }
}

pub fn selected_backend(sc: &Scenario) -> BackendChoice {
    match requested_backend(sc) {
        BackendSpec::Cpu => BackendChoice::Cpu,
        BackendSpec::Gpu => BackendChoice::Gpu,
        BackendSpec::Auto => {
            if auto_gpu_threshold(sc) && gpu_capability_error(sc).is_none() {
                BackendChoice::Gpu
            } else {
                BackendChoice::Cpu
            }
        }
    }
}

pub fn gpu_capability_error(sc: &Scenario) -> Option<&'static str> {
    if sc.physics.precision == Precision::F64 {
        return Some(GPU_F64_UNSUPPORTED);
    }
    if sc.is_3d()
        && (sc.multiphase.is_some()
            || sc.rotor.is_some()
            || sc.particles.is_some()
            || !matches!(sc.init, InitSpec::Rest)
            || sc
                .probes
                .iter()
                .any(|p| matches!(p, ProbeSpec::Force { .. })))
    {
        return Some(GPU_3D_UNSUPPORTED);
    }
    None
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct Physics {
    /// Kinematic viscosity (lattice units); tau = 3 nu + 0.5.
    pub nu: f64,
    #[serde(default)]
    pub collision: CollisionSpec,
    #[serde(default)]
    pub force: [f64; 2],
    /// Per-mass gravity g (lattice units): lowered to the solver's single
    /// Guo-force composition point as `rho(x)*g`, additive with `force` and
    /// any caller-owned per-cell force. W-VOF will substitute `rho(phi)*g`
    /// at that same solver point; this schema field is intentionally just
    /// the acceleration vector.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gravity: Option<[f64; 3]>,
    #[serde(default)]
    pub precision: Precision,
}

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
pub enum CollisionSpec {
    Bgk,
    #[default]
    #[serde(rename_all = "camelCase")]
    Trt,
}

impl CollisionSpec {
    pub fn to_core(self) -> Collision {
        match self {
            CollisionSpec::Bgk => Collision::Bgk,
            CollisionSpec::Trt => Collision::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum Precision {
    F32,
    #[default]
    F64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
pub enum EdgeSpec {
    Periodic,
    BounceBack,
    MovingWall {
        u: [f64; 2],
    },
    VelocityInlet {
        u: [f64; 2],
    },
    PressureOutlet {
        rho: f64,
    },
    Outflow,
    /// Convective (radiation) outflow: far less pressure-reflective than
    /// `Outflow`. `uConv` is the expected mean outflow speed, in (0, 1].
    #[serde(rename_all = "camelCase")]
    ConvectiveOutflow {
        u_conv: f64,
    },
}

impl EdgeSpec {
    fn to_core<T: Real>(self) -> EdgeBC<T> {
        match self {
            EdgeSpec::Periodic => EdgeBC::Periodic,
            EdgeSpec::BounceBack => EdgeBC::BounceBack,
            EdgeSpec::MovingWall { u } => EdgeBC::MovingWall {
                u: [T::r(u[0]), T::r(u[1])],
            },
            EdgeSpec::VelocityInlet { u } => EdgeBC::VelocityInlet {
                u: [T::r(u[0]), T::r(u[1])],
            },
            EdgeSpec::PressureOutlet { rho } => EdgeBC::PressureOutlet { rho: T::r(rho) },
            EdgeSpec::Outflow => EdgeBC::Outflow,
            EdgeSpec::ConvectiveOutflow { u_conv } => EdgeBC::ConvectiveOutflow {
                u_conv: T::r(u_conv),
            },
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdgesSpec {
    pub left: EdgeSpec,
    pub right: EdgeSpec,
    pub bottom: EdgeSpec,
    pub top: EdgeSpec,
    /// z = 0 face (3D only; ignored in 2D). Omitted = periodic.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub front: Option<EdgeSpec>,
    /// z = nz - 1 face (3D only; ignored in 2D). Omitted = periodic.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub back: Option<EdgeSpec>,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct InletProfile {
    pub edge: EdgeName,
    pub kind: ProfileKind,
    pub umax: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "camelCase")]
pub enum EdgeName {
    Left,
    Right,
    Bottom,
    Top,
}

impl EdgeName {
    pub fn to_core(self) -> Edge {
        match self {
            EdgeName::Left => Edge::Left,
            EdgeName::Right => Edge::Right,
            EdgeName::Bottom => Edge::Bottom,
            EdgeName::Top => Edge::Top,
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "camelCase")]
pub enum ProfileKind {
    /// Poiseuille parabola with the given peak velocity along the edge normal.
    Parabolic,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(tag = "shape", rename_all = "camelCase")]
pub enum Obstacle {
    /// 2D: a disk. 3D: extruded along z (a cylinder through the domain).
    Circle { cx: f64, cy: f64, r: f64 },
    /// 2D: a rectangle. 3D: extruded along z (a box through the domain).
    Rect {
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
    },
    /// 3D only: a solid ball (staircase approximation).
    Sphere { cx: f64, cy: f64, cz: f64, r: f64 },
}

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "camelCase")]
pub enum InitSpec {
    #[default]
    Rest,
    /// Dense liquid disk in vapour (pairs with `multiphase`).
    #[serde(rename_all = "camelCase")]
    Droplet {
        cx: f64,
        cy: f64,
        r: f64,
        rho_liquid: f64,
        rho_vapor: f64,
    },
    /// Liquid layer at the bottom (pairs with `multiphase` + gravity force).
    #[serde(rename_all = "camelCase")]
    Pool {
        height_frac: f64,
        rho_liquid: f64,
        rho_vapor: f64,
    },
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct MultiphaseSpec {
    /// Shan-Chen cohesion strength (negative; -5.0 is the validated default).
    pub g: f64,
    #[serde(default)]
    pub g_wall: f64,
    /// Virtual wall density for full-range contact-angle control (preferred
    /// over `gWall`): values near the liquid density wet the wall (θ → 0°),
    /// near the vapour density de-wet it (θ → 180°). See VALIDATION.md T11c.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_rho: Option<f64>,
}

/// Rotating-impeller volume penalization (MF-delta interim). The runner
/// constructs `lbm_core::compat::rotor::Rotor` from this each run; defaults
/// mirror the measured stability envelope in docs/PHYSICS.md.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct RotorSpec {
    pub cx: f64,
    pub cy: f64,
    pub n_blades: usize,
    #[serde(default)]
    pub r_hub: f64,
    pub r_blade: f64,
    pub thickness: f64,
    /// Angular velocity in rad/step. Tip speed omega*rBlade obeys the same
    /// low-Mach limits as walls (warn > 0.15, reject > 0.3).
    pub omega: f64,
    /// Penalization strength in (0, 1]; 1 pins blade cells to solid-body
    /// rotation exactly.
    #[serde(default = "default_chi")]
    pub chi: f64,
    /// Linear motor spin-up ramp on omega, in steps (0 = instant start).
    #[serde(default)]
    pub ramp_steps: usize,
    #[serde(default)]
    pub theta0: f64,
}

fn default_chi() -> f64 {
    1.0
}

/// One-way Lagrangian particles: deterministic grid seeding inside `seed`,
/// Schiller-Naumann drag, buoyancy-reduced gravity from `physics.gravity`.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct ParticlesSpec {
    pub count: usize,
    /// Particle diameter (lattice units).
    pub d: f64,
    /// Particle density relative to the fluid (rho_f = 1).
    pub rho_p: f64,
    #[serde(default)]
    pub restitution: f64,
    pub seed: SeedRegion,
    /// Write particles_<step>.csv every N steps (0 = end of run only).
    #[serde(default)]
    pub output_every: usize,
}

/// Axis-aligned seeding region [x0, x1] x [y0, y1] (lattice units).
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct SeedRegion {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct RunSpec {
    pub steps: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stop_when_steady: Option<SteadySpec>,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct SteadySpec {
    pub epsilon: f64,
    #[serde(default = "default_check_every")]
    pub check_every: usize,
}

fn default_check_every() -> usize {
    500
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
pub enum ProbeSpec {
    /// Momentum-exchange force on all obstacle cells.
    #[serde(rename_all = "camelCase")]
    Force { every: usize },
    /// Point time series of (ux, uy, rho); 3D also logs uz. `z` is 3D-only
    /// (omitted = mid-plane nz/2).
    #[serde(rename_all = "camelCase")]
    Point {
        x: usize,
        y: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        z: Option<usize>,
        every: usize,
    },
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "camelCase")]
pub enum FieldKind {
    Speed,
    Ux,
    Uy,
    Rho,
    Vorticity,
    /// Strain-rate invariant `gamma_dot = sqrt(2 S:S)` from the native
    /// non-equilibrium stress (exact — not a finite-difference reconstruction).
    /// Consumed by shear-threshold analyses; paired with `DissipationRate`.
    ShearRate,
    /// Viscous/turbulent dissipation rate `epsilon = nu * gamma_dot^2`
    /// (= `2 nu S:S`), from the same native gather as `ShearRate`. Consumed by
    /// the SCALEUP S-Fingerprint layer for `<eps>_vol` and the Kolmogorov length
    /// `eta = (nu^3 / eps)^(1/4)`.
    DissipationRate,
    /// Vorticity magnitude `|omega|` (central-difference velocity gradient).
    VorticityMag,
    /// Q-criterion `Q = 0.5 (|Omega|^2 - |S|^2)` (central-difference gradient).
    QCriterion,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct OutputSpec {
    pub field: FieldKind,
    pub format: OutputFormat,
    /// "end" or a step number is expressed via `every`/`at_end`; v0 keeps it
    /// simple: snapshots every N steps (0 = only at the end).
    #[serde(default)]
    pub every: usize,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "camelCase")]
pub enum OutputFormat {
    Png,
    Csv,
    /// VTK legacy structured points (ASCII), openable in ParaView etc.
    Vtk,
}

// ---------------------------------------------------------------- validation

/// A non-fatal advisory produced by [`validate`].
#[derive(Clone, Debug, Serialize)]
pub struct Warning {
    pub field: String,
    pub message: String,
}

/// Validate scenario semantics beyond what serde enforces. Returns warnings;
/// hard errors come from `SimConfig::build` when the scenario is applied.
pub fn validate(sc: &Scenario) -> Vec<Warning> {
    let mut warnings = Vec::new();
    let mut warn = |field: &str, message: String| {
        warnings.push(Warning {
            field: field.to_string(),
            message,
        });
    };
    let tau = 3.0 * sc.physics.nu + 0.5;
    if tau < TAU_LOW_WARN_THRESHOLD {
        warn(
            "physics.nu",
            format!("tau = {tau:.3} is close to the stability limit (below {TAU_LOW_WARN_THRESHOLD}). Increase the viscosity or the resolution"),
        );
    }
    let max_edge_speed = edge_speeds(&sc.edges).into_iter().fold(0.0, f64::max);
    if max_edge_speed > LATTICE_SPEED_WARN_THRESHOLD {
        warn(
            "edges",
            format!("inlet/wall velocity {max_edge_speed:.3} is at a level where compressibility error is noticeable (above {LATTICE_SPEED_WARN_THRESHOLD})"),
        );
    }
    if max_edge_speed > 0.0 && sc.physics.nu > 0.0 {
        let grid_re = max_edge_speed / sc.physics.nu;
        if grid_re > GRID_RE_WARN_THRESHOLD {
            warn(
                "physics",
                format!(
                    "grid Reynolds number U/ν = {grid_re:.1} > {GRID_RE_WARN_THRESHOLD}: risk of divergence (see PHYSICS.md)"
                ),
            );
        }
    }
    if let Some(r) = &sc.rotor {
        let tip = (r.omega * r.r_blade).abs();
        if tip > 0.3 {
            warn(
                "rotor.omega",
                format!("tip speed {tip:.3} exceeds the hard low-Mach limit 0.3 (the runner will reject it)"),
            );
        } else if tip > 0.15 {
            warn(
                "rotor.omega",
                format!("tip speed {tip:.3} is at a level where compressibility error is noticeable (above 0.15)"),
            );
        }
        if !(r.chi > 0.0 && r.chi <= 1.0) {
            warn("rotor.chi", format!("chi = {} is outside (0, 1]", r.chi));
        }
        if r.r_blade as usize * 2 >= sc.grid.nx.min(sc.grid.ny) {
            warn(
                "rotor.rBlade",
                "impeller diameter reaches the domain edge (leave clearance to the walls)"
                    .to_string(),
            );
        }
    }
    if let Some(p) = &sc.particles {
        if p.count == 0 || p.d <= 0.0 || p.rho_p <= 0.0 {
            warn(
                "particles",
                "count, d and rhoP must be positive".to_string(),
            );
        }
        if sc.physics.gravity.is_none() {
            warn(
                "particles",
                "no physics.gravity set: particles will neither settle nor feel buoyancy"
                    .to_string(),
            );
        }
        let (w, h) = (sc.grid.nx as f64, sc.grid.ny as f64);
        if p.seed.x0 < 0.0
            || p.seed.y0 < 0.0
            || p.seed.x1 >= w
            || p.seed.y1 >= h
            || p.seed.x0 > p.seed.x1
            || p.seed.y0 > p.seed.y1
        {
            warn(
                "particles.seed",
                "seed region is empty or outside the grid".to_string(),
            );
        }
    }
    if sc.multiphase.is_some() && sc.physics.precision == Precision::F32 {
        warn(
            "physics.precision",
            "f64 is recommended for multiphase flow (headroom for the steep interface gradient)"
                .to_string(),
        );
    }
    if let Some(mp) = &sc.multiphase {
        if mp.g > -4.0 {
            warn(
                "multiphase.g",
                format!(
                    "G = {} is weaker than the critical value -4 and will not phase-separate (recommended -5.0)",
                    mp.g
                ),
            );
        }
    }
    if selected_backend(sc) == BackendChoice::Gpu {
        #[cfg(not(feature = "gpu"))]
        if requested_backend(sc) == BackendSpec::Gpu {
            warn("compute.backend", GPU_FEATURE_UNAVAILABLE.to_string());
        }
        #[cfg(feature = "gpu")]
        if let Some(e) = gpu_capability_error(sc) {
            warn("compute.backend", e.to_string());
        } else {
            warn(
                "compute.backend",
                format!(
                    "GPU backend selected by {:?} policy for {} cells",
                    requested_backend(sc),
                    sc.grid.nx * sc.grid.ny * sc.grid.nz
                ),
            );
        }
    }
    if sc.is_3d() {
        if sc.multiphase.is_some() {
            warn(
                "multiphase",
                "3D (nz > 1) does not support multiphase flow (will error at build time)"
                    .to_string(),
            );
        }
    } else {
        if sc.edges.front.is_some() || sc.edges.back.is_some() {
            warn(
                "edges",
                "front/back are for 3D (nz > 1) only and are ignored in 2D".to_string(),
            );
        }
    }
    let mut named_edges = vec![
        ("edges.left", sc.edges.left),
        ("edges.right", sc.edges.right),
        ("edges.bottom", sc.edges.bottom),
        ("edges.top", sc.edges.top),
    ];
    if let Some(front) = sc.edges.front {
        named_edges.push(("edges.front", front));
    }
    if let Some(back) = sc.edges.back {
        named_edges.push(("edges.back", back));
    }
    for (name, spec) in named_edges {
        if let EdgeSpec::ConvectiveOutflow { u_conv } = spec {
            if !(u_conv > 0.0 && u_conv <= 1.0) {
                warn(
                    name,
                    format!(
                        "uConv = {u_conv} is outside (0,1] and will error at build time. \
                         Specify the expected mean outflow velocity (e.g. 0.05-0.15, comparable to the inlet velocity)"
                    ),
                );
            }
        }
    }
    warnings
}

pub fn unit_report_for(sc: &Scenario) -> Result<Option<UnitReport>, String> {
    sc.units.as_ref().map(unit_report).transpose()
}

fn edge_speeds(e: &EdgesSpec) -> [f64; 6] {
    [
        e.left,
        e.right,
        e.bottom,
        e.top,
        e.front.unwrap_or(EdgeSpec::Periodic),
        e.back.unwrap_or(EdgeSpec::Periodic),
    ]
    .map(|s| match s {
        EdgeSpec::MovingWall { u } | EdgeSpec::VelocityInlet { u } => {
            (u[0] * u[0] + u[1] * u[1]).sqrt()
        }
        _ => 0.0,
    })
}

fn named_edge_spec(edges: &EdgesSpec, edge: EdgeName) -> EdgeSpec {
    match edge {
        EdgeName::Left => edges.left,
        EdgeName::Right => edges.right,
        EdgeName::Bottom => edges.bottom,
        EdgeName::Top => edges.top,
    }
}

fn validate_inlet_profile_request(sc: &Scenario) -> Result<(), ConfigError> {
    let Some(p) = &sc.inlet_profile else {
        return Ok(());
    };
    if !matches!(
        named_edge_spec(&sc.edges, p.edge),
        EdgeSpec::VelocityInlet { .. }
    ) {
        return Err(ConfigError::InvalidParameter {
            what: "inletProfile edge must reference a velocityInlet boundary",
            value: 0.0,
        });
    }
    let speed = p.umax.abs();
    if !(speed <= MAX_SPEED) {
        return Err(ConfigError::VelocityTooHigh { speed });
    }
    Ok(())
}

// ---------------------------------------------------------------- build

/// A built simulation, precision-erased for the runner.
pub enum SimHandle {
    F32(Simulation<f32>, Option<ShanChen<f32>>),
    F64(Simulation<f64>, Option<ShanChen<f64>>),
}

#[cfg(feature = "gpu")]
pub type GpuSim2 = lbm_core::gpu::GpuSolver<lbm_core::lattice::D2Q9>;

/// Build error for 2D scenarios: either a core configuration error or a
/// requested scenario capability the 2D compat execution path cannot provide.
#[derive(Debug)]
pub enum BuildError {
    /// Invalid physical/boundary configuration.
    Core(ConfigError),
    /// Invalid SI unit conversion at the scenario boundary.
    Units(String),
    /// Feature not available on the 2D scenario path (message is user-facing).
    Unsupported(&'static str),
    /// Explicit or auto-selected backend that cannot be honored.
    BackendUnavailable(&'static str),
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BuildError::Core(e) => write!(f, "{e}"),
            BuildError::Units(e) => write!(f, "{e}"),
            BuildError::Unsupported(what) => write!(f, "{what}"),
            BuildError::BackendUnavailable(what) => write!(f, "{what}"),
        }
    }
}

impl std::error::Error for BuildError {}

impl From<ConfigError> for BuildError {
    fn from(e: ConfigError) -> Self {
        BuildError::Core(e)
    }
}

// ---------------------------------------------------------------- build (3D)

/// The 3D engine type behind a scenario: V2 core, D3Q19, CPU backend,
/// monolithic decomposition (ARCHITECTURE_V2; `compute.backend: "cpu"`).
pub type Solver3<T> = lbm_core::solver::Solver<
    lbm_core::lattice::D3Q19,
    T,
    lbm_core::backend::CpuScalar,
    lbm_core::halo::LocalPeriodic,
>;

/// A built 3D simulation, precision-erased for the runner.
pub enum Sim3Handle {
    F32(Solver3<f32>),
    F64(Solver3<f64>),
}

/// Build error for 3D scenarios: either a core configuration error (same
/// semantics as the 2D `SimConfig::build`) or a scenario feature the 3D
/// engine does not support yet.
#[derive(Debug)]
pub enum Build3Error {
    /// Invalid physical/boundary configuration (compat-facade error kind, kept
    /// for the scenario-level checks that have no core equivalent — e.g.
    /// unpaired periodic faces).
    Core(ConfigError),
    /// Invalid native `GlobalSpec`, as reported by `GlobalSpec::validate`
    /// (A-4): uncovered face, ν ≤ 0, periodic × open, out-of-range BC, …
    Spec(lbm_core::solver::SpecError),
    /// Invalid SI unit conversion at the scenario boundary.
    Units(String),
    /// Explicit backend request that this scenario path cannot honor.
    BackendUnavailable(&'static str),
    /// Feature not available on the 3D engine (message is user-facing).
    Unsupported(&'static str),
}

impl std::fmt::Display for Build3Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Build3Error::Core(e) => write!(f, "{e}"),
            Build3Error::Spec(e) => write!(f, "{e}"),
            Build3Error::Units(e) => write!(f, "{e}"),
            Build3Error::BackendUnavailable(what) => write!(f, "{what}"),
            Build3Error::Unsupported(what) => write!(f, "unsupported in 3D (nz > 1): {what}"),
        }
    }
}

impl std::error::Error for Build3Error {}

impl From<ConfigError> for Build3Error {
    fn from(e: ConfigError) -> Self {
        Build3Error::Core(e)
    }
}

impl From<lbm_core::solver::SpecError> for Build3Error {
    fn from(e: lbm_core::solver::SpecError) -> Self {
        Build3Error::Spec(e)
    }
}

/// Dimension-dispatching build check for validators (CLI `validate`, MCP
/// `validate_scenario`): construct the simulation the same way `run` would
/// (2D or 3D) and report the error text, discarding the handle.
pub fn build_check(sc: &Scenario) -> Result<(), String> {
    let resolved = match resolve(sc) {
        Ok(Some(r)) => r.scenario,
        Ok(None) => sc.clone(),
        Err(e) => return Err(e),
    };
    let sc = &resolved;
    if selected_backend(sc) == BackendChoice::Gpu {
        #[cfg(not(feature = "gpu"))]
        if requested_backend(sc) == BackendSpec::Gpu {
            return Err(GPU_FEATURE_UNAVAILABLE.to_string());
        }
        #[cfg(feature = "gpu")]
        if let Some(e) = gpu_capability_error(sc) {
            return Err(e.to_string());
        }
    }
    if sc.is_3d() {
        build3d(sc).map(|_| ()).map_err(|e| e.to_string())
    } else {
        build(sc).map(|_| ()).map_err(|e| e.to_string())
    }
}

/// The six face BCs of a 3D scenario in `Face::index()` order
/// (left, right, bottom, top, front, back); omitted z faces are periodic.
fn face_specs(e: &EdgesSpec) -> [EdgeSpec; 6] {
    [
        e.left,
        e.right,
        e.bottom,
        e.top,
        e.front.unwrap_or(EdgeSpec::Periodic),
        e.back.unwrap_or(EdgeSpec::Periodic),
    ]
}

/// Build a 3D (D3Q19) simulation from a scenario with `grid.nz > 1`.
///
/// Feature scope (minimal wiring, COMPETITIVE_SPEC M-C): single phase,
/// `init: rest`, CPU backend. Boundary semantics mirror the 2D contract:
/// walls are one-cell solid rims (half-way bounce-back), periodic faces must
/// pair, open faces (Zou–He / outflow / convective) must all lie on one axis.
pub fn build3d(sc: &Scenario) -> Result<Sim3Handle, Build3Error> {
    let resolved;
    let sc = match resolve(sc) {
        Ok(Some(r)) => {
            resolved = r.scenario;
            &resolved
        }
        Ok(None) => sc,
        Err(e) => return Err(Build3Error::Units(e)),
    };
    validate_inlet_profile_request(sc)?;
    Ok(match sc.physics.precision {
        Precision::F32 => Sim3Handle::F32(build3d_t::<f32>(sc)?),
        Precision::F64 => Sim3Handle::F64(build3d_t::<f64>(sc)?),
    })
}

fn build3d_t<T: lbm_core::real::Real>(sc: &Scenario) -> Result<Solver3<T>, Build3Error> {
    use lbm_core::prelude::{
        build_wall_rims, CollisionKind, CpuScalar, Face, FaceBC, GlobalSpec, LocalPeriodic, Solver,
        WallSpec,
    };

    assert!(sc.is_3d(), "build3d requires grid.nz > 1");
    if sc.multiphase.is_some() {
        return Err(Build3Error::Unsupported("multiphase (multiphase flow)"));
    }
    if sc.rotor.is_some() {
        return Err(Build3Error::Unsupported(
            "rotor (2D only in this increment; 3D rotor lands with the z-extruded evaluator)",
        ));
    }
    if sc.particles.is_some() {
        return Err(Build3Error::Unsupported(
            "particles (2D only in this increment)",
        ));
    }
    if !matches!(sc.init, InitSpec::Rest) {
        return Err(Build3Error::Unsupported("init must be rest only"));
    }
    if selected_backend(sc) == BackendChoice::Gpu {
        #[cfg(not(feature = "gpu"))]
        if requested_backend(sc) == BackendSpec::Gpu {
            return Err(Build3Error::BackendUnavailable(GPU_FEATURE_UNAVAILABLE));
        }
        #[cfg(feature = "gpu")]
        if let Some(e) = gpu_capability_error(sc) {
            return Err(Build3Error::BackendUnavailable(e));
        }
    }
    let dims = [sc.grid.nx, sc.grid.ny, sc.grid.nz];
    let specs = face_specs(&sc.edges);
    // Periodic *pairing* per axis is a scenario-level concern (two separate
    // EdgeSpecs collapse into one `GlobalSpec::periodic` bool), so it stays
    // here; the extents / open-axis / BC-range / coverage checks are delegated
    // to `GlobalSpec::validate` below (A-4, no duplication).
    let mut periodic = [false; 3];
    for (axis, name) in [(0usize, "x"), (1, "y"), (2, "z")] {
        let lo = matches!(specs[2 * axis], EdgeSpec::Periodic);
        let hi = matches!(specs[2 * axis + 1], EdgeSpec::Periodic);
        if lo != hi {
            return Err(ConfigError::UnpairedPeriodic { axis: name }.into());
        }
        periodic[axis] = lo && hi;
    }
    // Open faces must not share a domain edge (V1's corner rule, lifted to 3D:
    // all open faces on one axis). Kept here so callers keep seeing the
    // `AdjacentOpenEdges` kind; `GlobalSpec::validate` re-checks it as
    // `OpenFacesOnMultipleAxes`.
    let is_open = |s: &EdgeSpec| {
        matches!(
            s,
            EdgeSpec::VelocityInlet { .. }
                | EdgeSpec::PressureOutlet { .. }
                | EdgeSpec::Outflow
                | EdgeSpec::ConvectiveOutflow { .. }
        )
    };
    if (0..3)
        .filter(|a| is_open(&specs[2 * a]) || is_open(&specs[2 * a + 1]))
        .count()
        > 1
    {
        return Err(ConfigError::AdjacentOpenEdges.into());
    }
    // A wall's velocity lives in `WallSpec`, not in `GlobalSpec::faces`, so its
    // low-Mach limit is checked here (validate only covers inlet faces).
    for s in &specs {
        if let EdgeSpec::MovingWall { u } = *s {
            let sp = (u[0] * u[0] + u[1] * u[1]).sqrt();
            if sp > MAX_SPEED {
                return Err(ConfigError::VelocityTooHigh { speed: sp }.into());
            }
        }
    }

    // Walls and open-face BCs. The scenario's 2D velocity vectors embed as
    // (ux, uy, 0) — z-face inlets/lids thus carry in-plane velocity only.
    let mut walls = WallSpec::<T>::default();
    let mut faces = [FaceBC::<T>::Closed; 6];
    for (i, s) in specs.iter().enumerate() {
        match *s {
            EdgeSpec::Periodic => {}
            EdgeSpec::BounceBack => walls.is_wall[i] = true,
            EdgeSpec::MovingWall { u } => {
                walls.is_wall[i] = true;
                walls.u[i] = [T::r(u[0]), T::r(u[1]), T::zero()];
            }
            EdgeSpec::VelocityInlet { u } => {
                faces[i] = FaceBC::Velocity {
                    u: [T::r(u[0]), T::r(u[1]), T::zero()],
                }
            }
            EdgeSpec::PressureOutlet { rho } => faces[i] = FaceBC::Pressure { rho: T::r(rho) },
            EdgeSpec::Outflow => faces[i] = FaceBC::Outflow,
            EdgeSpec::ConvectiveOutflow { u_conv } => {
                faces[i] = FaceBC::Convective {
                    u_conv: T::r(u_conv),
                }
            }
        }
    }
    let spec = GlobalSpec::<T> {
        dims,
        nu: sc.physics.nu,
        collision: match sc.physics.collision {
            CollisionSpec::Bgk => CollisionKind::Bgk,
            CollisionSpec::Trt => CollisionKind::Trt {
                magic: CollisionKind::MAGIC_STD,
            },
        },
        periodic,
        faces,
        force: [
            T::r(sc.physics.force[0]),
            T::r(sc.physics.force[1]),
            T::zero(),
        ],
        sources: Vec::new(),
        face_patches: Vec::new(),
    };
    let (solid, wall_u) = build_wall_rims::<T>(3, dims, &walls);
    // A-4: the single native config gate (extents, ν, periodic × open,
    // open-axis, uncovered faces, inlet/pressure/convective ranges, force[2]).
    // Surfaces a typed error here; `Solver::new` re-checks as a last-line
    // panic guard.
    spec.validate(3, &solid)?;
    let mut s: Solver3<T> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );

    if let Some(g) = sc.physics.gravity {
        validate_gravity_vector(g, true)?;
        s.set_gravity([T::r(g[0]), T::r(g[1]), T::r(g[2])]);
    }

    // Obstacles: 2D shapes extrude along z; spheres are native 3D.
    let mut any_obstacle = false;
    for ob in &sc.obstacles {
        let mut set_region = |pred: &dyn Fn(usize, usize, usize) -> bool| {
            for z in 0..dims[2] {
                for y in 0..dims[1] {
                    for x in 0..dims[0] {
                        if pred(x, y, z) {
                            s.set_solid(x, y, z);
                            any_obstacle = true;
                        }
                    }
                }
            }
        };
        match *ob {
            Obstacle::Circle { cx, cy, r } => {
                let r2 = r * r;
                set_region(&move |x, y, _| {
                    let (dx, dy) = (x as f64 - cx, y as f64 - cy);
                    dx * dx + dy * dy <= r2
                });
            }
            Obstacle::Rect { x0, y0, x1, y1 } => {
                set_region(&move |x, y, _| x >= x0 && x <= x1 && y >= y0 && y <= y1);
            }
            Obstacle::Sphere { cx, cy, cz, r } => {
                let r2 = r * r;
                set_region(&move |x, y, z| {
                    let (dx, dy, dz) = (x as f64 - cx, y as f64 - cy, z as f64 - cz);
                    dx * dx + dy * dy + dz * dz <= r2
                });
            }
        }
    }

    // Parabolic inlet profile: duct-type product profile
    // u(t1, t2) = umax f(t1) f(t2) along the inward normal, where f is the
    // half-way-wall parabola on a walled tangent axis and 1 on a periodic
    // one (so a z-periodic slab degenerates to the 2D parabola exactly).
    if let Some(p) = &sc.inlet_profile {
        let face = match p.edge {
            EdgeName::Left => Face::XNeg,
            EdgeName::Right => Face::XPos,
            EdgeName::Bottom => Face::YNeg,
            EdgeName::Top => Face::YPos,
        };
        if !matches!(faces[face.index()], FaceBC::Velocity { .. }) {
            return Err(Build3Error::Unsupported(
                "inletProfile can only be set on a velocityInlet edge",
            ));
        }
        let (t1, t2) = face.tangents();
        let n = face.n_in();
        let normal = [n[0] as f64, n[1] as f64, n[2] as f64];
        let umax = p.umax;
        let factor = move |axis: usize, c: usize| -> f64 {
            if periodic[axis] {
                return 1.0;
            }
            let h = (dims[axis] - 2) as f64;
            if c == 0 || c as f64 >= h + 1.0 {
                return 0.0;
            }
            let w = c as f64 - 0.5;
            4.0 * w * (h - w) / (h * h)
        };
        s.set_inlet_profile_with(face, move |c1, c2| {
            let mag = umax * factor(t1, c1) * factor(t2, c2);
            [
                T::r(mag * normal[0]),
                T::r(mag * normal[1]),
                T::r(mag * normal[2]),
            ]
        });
    }

    if sc
        .probes
        .iter()
        .any(|p| matches!(p, ProbeSpec::Force { .. }))
    {
        if !any_obstacle {
            return Err(Build3Error::Unsupported(
                "the force probe requires obstacles",
            ));
        }
        // Probe all obstacle solids (cells strictly inside the domain box,
        // rims excluded) — 2D convention lifted to 3D.
        let solid: Vec<bool> = (0..dims[2])
            .flat_map(|z| (0..dims[1]).flat_map(move |y| (0..dims[0]).map(move |x| (x, y, z))))
            .map(|(x, y, z)| s.is_solid(x, y, z))
            .collect();
        let (nx, ny) = (dims[0], dims[1]);
        let rim = move |c: usize, n: usize| c == 0 || c == n - 1;
        s.set_force_probe(move |x, y, z| {
            !rim(x, dims[0]) && !rim(y, dims[1]) && !rim(z, dims[2]) && solid[(z * ny + y) * nx + x]
        });
    }

    Ok(s)
}

/// Build the 2D simulation (+ optional multiphase driver) from a scenario.
/// Scenarios with `grid.nz > 1` must go through [`build3d`] instead.
pub fn build(sc: &Scenario) -> Result<SimHandle, BuildError> {
    let resolved;
    let sc = match resolve(sc) {
        Ok(Some(r)) => {
            resolved = r.scenario;
            &resolved
        }
        Ok(None) => sc,
        Err(e) => return Err(BuildError::Units(e)),
    };
    validate_inlet_profile_request(sc)?;
    if sc.is_3d() {
        return Err(ConfigError::InvalidParameter {
            what: "grid.nz (2D build requires nz == 1; the runner dispatches 3D to build3d)",
            value: sc.grid.nz as f64,
        }
        .into());
    }
    if selected_backend(sc) == BackendChoice::Gpu {
        #[cfg(not(feature = "gpu"))]
        if requested_backend(sc) == BackendSpec::Gpu {
            return Err(BuildError::BackendUnavailable(GPU_FEATURE_UNAVAILABLE));
        }
        #[cfg(feature = "gpu")]
        if let Some(e) = gpu_capability_error(sc) {
            return Err(BuildError::BackendUnavailable(e));
        }
    }
    Ok(match sc.physics.precision {
        Precision::F32 => {
            let (sim, mp) = build_t::<f32>(sc)?;
            SimHandle::F32(sim, mp)
        }
        Precision::F64 => {
            let (sim, mp) = build_t::<f64>(sc)?;
            SimHandle::F64(sim, mp)
        }
    })
}

fn build_t<T: Real>(sc: &Scenario) -> Result<(Simulation<T>, Option<ShanChen<T>>), ConfigError> {
    let mut sim: Simulation<T> = SimConfig {
        nx: sc.grid.nx,
        ny: sc.grid.ny,
        nu: sc.physics.nu,
        collision: sc.physics.collision.to_core(),
        edges: Edges {
            left: sc.edges.left.to_core(),
            right: sc.edges.right.to_core(),
            bottom: sc.edges.bottom.to_core(),
            top: sc.edges.top.to_core(),
        },
        force: [T::r(sc.physics.force[0]), T::r(sc.physics.force[1])],
    }
    .build()?;

    if let Some(g) = sc.physics.gravity {
        validate_gravity_vector(g, false)?;
        sim.set_gravity([T::r(g[0]), T::r(g[1])]);
    }

    for ob in &sc.obstacles {
        match *ob {
            Obstacle::Circle { cx, cy, r } => {
                let r2 = r * r;
                sim.set_solid_region(|x, y| {
                    let dx = x as f64 - cx;
                    let dy = y as f64 - cy;
                    dx * dx + dy * dy <= r2
                });
                if sc.wall == Some(WallModel::Bouzidi) {
                    sim.set_bouzidi_circle(cx, cy, r);
                }
            }
            Obstacle::Rect { x0, y0, x1, y1 } => {
                sim.set_solid_region(|x, y| x >= x0 && x <= x1 && y >= y0 && y <= y1);
            }
            Obstacle::Sphere { r, .. } => {
                return Err(ConfigError::InvalidParameter {
                    what: "obstacles: sphere requires a 3D grid (nz > 1)",
                    value: r,
                });
            }
        }
    }

    if let Some(p) = &sc.inlet_profile {
        let edge = p.edge.to_core();
        let (nx, ny) = (sim.nx(), sim.ny());
        let len = match edge {
            Edge::Left | Edge::Right => ny,
            Edge::Bottom | Edge::Top => nx,
        };
        let h = (len - 2) as f64;
        let umax = p.umax;
        let normal_sign: [f64; 2] = match edge {
            Edge::Left => [1.0, 0.0],
            Edge::Right => [-1.0, 0.0],
            Edge::Bottom => [0.0, 1.0],
            Edge::Top => [0.0, -1.0],
        };
        sim.set_inlet_profile(edge, move |c| {
            if c == 0 || c as f64 >= h + 1.0 {
                return [T::zero(); 2];
            }
            let yw = c as f64 - 0.5;
            let mag = 4.0 * umax * yw * (h - yw) / (h * h);
            [T::r(mag * normal_sign[0]), T::r(mag * normal_sign[1])]
        });
    }

    match sc.init {
        InitSpec::Rest => {}
        InitSpec::Droplet {
            cx,
            cy,
            r,
            rho_liquid,
            rho_vapor,
        } => {
            let r2 = r * r;
            sim.init_with(|x, y| {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let rho = if dx * dx + dy * dy <= r2 {
                    rho_liquid
                } else {
                    rho_vapor
                };
                (T::r(rho), T::zero(), T::zero())
            });
        }
        InitSpec::Pool {
            height_frac,
            rho_liquid,
            rho_vapor,
        } => {
            let ny = sim.ny();
            let cut = (height_frac * ny as f64) as usize;
            sim.init_with(|_, y| {
                let rho = if y < cut { rho_liquid } else { rho_vapor };
                (T::r(rho), T::zero(), T::zero())
            });
        }
    }

    if sc
        .probes
        .iter()
        .any(|p| matches!(p, ProbeSpec::Force { .. }))
    {
        // probe all obstacle solids (rims excluded: only cells strictly inside)
        let (nx, ny) = (sim.nx(), sim.ny());
        let solid: Vec<bool> = sim.solid_field().to_vec();
        sim.set_force_probe(move |x, y| {
            x > 0 && y > 0 && x < nx - 1 && y < ny - 1 && solid[y * nx + x]
        });
    }

    let mp = sc.multiphase.as_ref().map(|m| {
        let mut model = ShanChen::<T>::new(m.g).with_wall(m.g_wall);
        if let Some(rho_w) = m.wall_rho {
            model = model.with_wall_rho(rho_w);
        }
        model
    });
    Ok((sim, mp))
}

fn validate_gravity_vector(g: [f64; 3], is_3d: bool) -> Result<(), ConfigError> {
    for (i, v) in g.iter().enumerate() {
        if !v.is_finite() {
            return Err(ConfigError::NonFiniteParameter {
                what: match i {
                    0 => "physics.gravity[0]",
                    1 => "physics.gravity[1]",
                    _ => "physics.gravity[2]",
                },
            });
        }
    }
    if !is_3d && g[2] != 0.0 {
        return Err(ConfigError::InvalidParameter {
            what: "physics.gravity[2] (2D scenarios require gz == 0)",
            value: g[2],
        });
    }
    Ok(())
}

#[cfg(feature = "gpu")]
pub fn build_gpu2d(sc: &Scenario) -> Result<GpuSim2, BuildError> {
    use lbm_core::lattice::{D2Q9, Face};
    use lbm_core::prelude::{build_wall_rims, CollisionKind, FaceBC, GlobalSpec, WallSpec};
    use std::sync::Arc;

    let resolved;
    let sc = match resolve(sc) {
        Ok(Some(r)) => {
            resolved = r.scenario;
            &resolved
        }
        Ok(None) => sc,
        Err(e) => return Err(BuildError::Units(e)),
    };
    validate_inlet_profile_request(sc)?;
    if sc.is_3d() {
        return Err(BuildError::Unsupported("build_gpu2d requires nz == 1"));
    }
    if let Some(e) = gpu_capability_error(sc) {
        return Err(BuildError::BackendUnavailable(e));
    }
    if sc.multiphase.is_some() {
        return Err(BuildError::Unsupported("GPU scenario dispatch does not support multiphase"));
    }
    if sc.rotor.is_some() {
        return Err(BuildError::Unsupported("GPU scenario dispatch does not support rotor"));
    }
    if sc.particles.is_some() {
        return Err(BuildError::Unsupported("GPU scenario dispatch does not support particles"));
    }
    if !matches!(sc.init, InitSpec::Rest) {
        return Err(BuildError::Unsupported("GPU scenario dispatch supports init rest only"));
    }
    let dims = [sc.grid.nx, sc.grid.ny, 1];
    let specs = [sc.edges.left, sc.edges.right, sc.edges.bottom, sc.edges.top];
    let mut periodic = [false; 3];
    for (axis, name) in [(0usize, "x"), (1, "y")] {
        let lo = matches!(specs[2 * axis], EdgeSpec::Periodic);
        let hi = matches!(specs[2 * axis + 1], EdgeSpec::Periodic);
        if lo != hi {
            return Err(ConfigError::UnpairedPeriodic { axis: name }.into());
        }
        periodic[axis] = lo && hi;
    }
    periodic[2] = true;
    let mut walls = WallSpec::<f32>::default();
    let mut faces = [FaceBC::<f32>::Closed; 6];
    for (i, s) in specs.iter().enumerate() {
        match *s {
            EdgeSpec::Periodic => {}
            EdgeSpec::BounceBack => walls.is_wall[i] = true,
            EdgeSpec::MovingWall { u } => {
                walls.is_wall[i] = true;
                walls.u[i] = [u[0] as f32, u[1] as f32, 0.0];
            }
            EdgeSpec::VelocityInlet { u } => {
                faces[i] = FaceBC::Velocity {
                    u: [u[0] as f32, u[1] as f32, 0.0],
                };
            }
            EdgeSpec::PressureOutlet { rho } => faces[i] = FaceBC::Pressure { rho: rho as f32 },
            EdgeSpec::Outflow => faces[i] = FaceBC::Outflow,
            EdgeSpec::ConvectiveOutflow { u_conv } => {
                faces[i] = FaceBC::Convective {
                    u_conv: u_conv as f32,
                }
            }
        }
    }
    let spec = GlobalSpec::<f32> {
        dims,
        nu: sc.physics.nu,
        collision: match sc.physics.collision {
            CollisionSpec::Bgk => CollisionKind::Bgk,
            CollisionSpec::Trt => CollisionKind::Trt {
                magic: CollisionKind::MAGIC_STD,
            },
        },
        periodic,
        faces,
        force: [sc.physics.force[0] as f32, sc.physics.force[1] as f32, 0.0],
        ..Default::default()
    };
    let (mut solid, wall_u) = build_wall_rims::<f32>(2, dims, &walls);
    let idx = |x: usize, y: usize| y * dims[0] + x;
    for ob in &sc.obstacles {
        match *ob {
            Obstacle::Circle { cx, cy, r } => {
                let r2 = r * r;
                for y in 0..dims[1] {
                    for x in 0..dims[0] {
                        let dx = x as f64 - cx;
                        let dy = y as f64 - cy;
                        if dx * dx + dy * dy <= r2 {
                            solid[idx(x, y)] = true;
                        }
                    }
                }
            }
            Obstacle::Rect { x0, y0, x1, y1 } => {
                for y in y0..=y1.min(dims[1] - 1) {
                    for x in x0..=x1.min(dims[0] - 1) {
                        solid[idx(x, y)] = true;
                    }
                }
            }
            Obstacle::Sphere { r, .. } => {
                return Err(ConfigError::InvalidParameter {
                    what: "obstacles: sphere requires a 3D grid (nz > 1)",
                    value: r,
                }
                .into());
            }
        }
    }
    spec.validate(2, &solid)
        .map_err(|e| BuildError::Unsupported(Box::leak(e.to_string().into_boxed_str())))?;
    let ctx = lbm_core::gpu::GpuContext::new()
        .map_err(|_| BuildError::BackendUnavailable("requested backend \"gpu\" is unavailable: no usable GPU adapter was found"))?;
    let mut s = lbm_core::gpu::GpuSolver::<D2Q9>::new(&spec, &solid, &wall_u, Arc::clone(&ctx));
    if let Some(p) = &sc.inlet_profile {
        let face = match p.edge {
            EdgeName::Left => Face::XNeg,
            EdgeName::Right => Face::XPos,
            EdgeName::Bottom => Face::YNeg,
            EdgeName::Top => Face::YPos,
        };
        let (t1, _t2) = face.tangents();
        let len = dims[t1];
        let h = (len - 2) as f64;
        let n = face.n_in();
        let normal = [n[0] as f64, n[1] as f64, 0.0];
        let values = (0..len)
            .map(|c| {
                if c == 0 || c as f64 >= h + 1.0 {
                    [0.0, 0.0, 0.0]
                } else {
                    let w = c as f64 - 0.5;
                    let mag = 4.0 * p.umax * w * (h - w) / (h * h);
                    [
                        (mag * normal[0]) as f32,
                        (mag * normal[1]) as f32,
                        0.0,
                    ]
                }
            })
            .collect::<Vec<_>>();
        s.set_inlet_profile(face, &values);
    }
    if sc
        .probes
        .iter()
        .any(|p| matches!(p, ProbeSpec::Force { .. }))
    {
        let solid_probe = solid.clone();
        let (nx, ny) = (dims[0], dims[1]);
        s.set_force_probe(move |x, y, _| {
            x > 0 && x + 1 < nx && y > 0 && y + 1 < ny && solid_probe[y * nx + x]
        });
    }
    Ok(s)
}

// ---------------------------------------------------------------- presets

/// Built-in presets: (name, description, scenario JSON factory).
pub fn presets() -> Vec<(&'static str, &'static str, Scenario)> {
    let cavity = Scenario {
        version: 0,
        name: "cavity".into(),
        grid: Grid {
            nx: 128,
            ny: 128,
            nz: 1,
        },
        physics: Physics {
            nu: 0.02,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            gravity: None,
            precision: Precision::F64,
        },
        units: None,
        compute: None,
        wall: None,
        edges: EdgesSpec {
            left: EdgeSpec::BounceBack,
            right: EdgeSpec::BounceBack,
            bottom: EdgeSpec::BounceBack,
            top: EdgeSpec::MovingWall { u: [0.1, 0.0] },
            front: None,
            back: None,
        },
        inlet_profile: None,
        obstacles: vec![],
        init: InitSpec::Rest,
        multiphase: None,
        rotor: None,
        particles: None,
        run: RunSpec {
            steps: 20_000,
            stop_when_steady: Some(SteadySpec {
                epsilon: 1e-8,
                check_every: 500,
            }),
        },
        probes: vec![],
        outputs: vec![OutputSpec {
            field: FieldKind::Speed,
            format: OutputFormat::Png,
            every: 0,
        }],
    };
    let cylinder = Scenario {
        version: 0,
        name: "cylinder-karman".into(),
        grid: Grid {
            nx: 440,
            ny: 164,
            nz: 1,
        },
        physics: Physics {
            nu: 0.04,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            gravity: None,
            precision: Precision::F64,
        },
        units: None,
        compute: None,
        wall: None,
        edges: EdgesSpec {
            left: EdgeSpec::VelocityInlet { u: [0.1, 0.0] },
            right: EdgeSpec::PressureOutlet { rho: 1.0 },
            bottom: EdgeSpec::BounceBack,
            top: EdgeSpec::BounceBack,
            front: None,
            back: None,
        },
        inlet_profile: Some(InletProfile {
            edge: EdgeName::Left,
            kind: ProfileKind::Parabolic,
            umax: 0.15,
        }),
        obstacles: vec![Obstacle::Circle {
            cx: 80.0,
            cy: 80.0,
            r: 20.0,
        }],
        init: InitSpec::Rest,
        multiphase: None,
        rotor: None,
        particles: None,
        run: RunSpec {
            steps: 40_000,
            stop_when_steady: None,
        },
        probes: vec![ProbeSpec::Force { every: 10 }],
        outputs: vec![
            OutputSpec {
                field: FieldKind::Vorticity,
                format: OutputFormat::Png,
                every: 10_000,
            },
            OutputSpec {
                field: FieldKind::Speed,
                format: OutputFormat::Png,
                every: 0,
            },
        ],
    };
    let droplet = Scenario {
        version: 0,
        name: "two-phase-droplet".into(),
        grid: Grid {
            nx: 128,
            ny: 128,
            nz: 1,
        },
        physics: Physics {
            nu: 1.0 / 6.0,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            gravity: None,
            precision: Precision::F64,
        },
        units: None,
        compute: None,
        wall: None,
        edges: EdgesSpec {
            left: EdgeSpec::Periodic,
            right: EdgeSpec::Periodic,
            bottom: EdgeSpec::Periodic,
            top: EdgeSpec::Periodic,
            front: None,
            back: None,
        },
        inlet_profile: None,
        obstacles: vec![],
        init: InitSpec::Droplet {
            cx: 64.0,
            cy: 64.0,
            r: 20.0,
            rho_liquid: 2.0,
            rho_vapor: 0.15,
        },
        multiphase: Some(MultiphaseSpec {
            g: -5.0,
            g_wall: 0.0,
            wall_rho: None,
        }),
        rotor: None,
        particles: None,
        run: RunSpec {
            steps: 20_000,
            stop_when_steady: None,
        },
        probes: vec![],
        outputs: vec![OutputSpec {
            field: FieldKind::Rho,
            format: OutputFormat::Png,
            every: 0,
        }],
    };
    // T11c geometry: half-disk on the bottom wall, virtual wall density 1.0
    // relaxes to a spherical cap with contact angle ~63 deg.
    let droplet_on_wall = Scenario {
        version: 0,
        name: "droplet-on-wall".into(),
        grid: Grid {
            nx: 160,
            ny: 100,
            nz: 1,
        },
        physics: Physics {
            nu: 1.0 / 6.0,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            gravity: None,
            precision: Precision::F64,
        },
        units: None,
        compute: None,
        wall: None,
        edges: EdgesSpec {
            left: EdgeSpec::Periodic,
            right: EdgeSpec::Periodic,
            bottom: EdgeSpec::BounceBack,
            top: EdgeSpec::BounceBack,
            front: None,
            back: None,
        },
        inlet_profile: None,
        obstacles: vec![],
        init: InitSpec::Droplet {
            cx: 80.0,
            cy: 1.0,
            r: 22.0,
            rho_liquid: 2.0,
            rho_vapor: 0.15,
        },
        multiphase: Some(MultiphaseSpec {
            g: -5.0,
            g_wall: 0.0,
            wall_rho: Some(1.0),
        }),
        rotor: None,
        particles: None,
        run: RunSpec {
            steps: 30_000,
            stop_when_steady: None,
        },
        probes: vec![],
        outputs: vec![
            OutputSpec {
                field: FieldKind::Rho,
                format: OutputFormat::Png,
                every: 0,
            },
            OutputSpec {
                field: FieldKind::Rho,
                format: OutputFormat::Vtk,
                every: 0,
            },
        ],
    };
    vec![
        (
            "cavity",
            "lid-driven cavity (with steady-state detection)",
            cavity,
        ),
        (
            "cylinder-karman",
            "Kármán vortex street around a cylinder + drag probe",
            cylinder,
        ),
        (
            "two-phase-droplet",
            "equilibration of a Shan-Chen two-phase droplet",
            droplet,
        ),
        (
            "droplet-on-wall",
            "contact-angle demo of a droplet on a wall (virtual wall density wallRho=1.0 → θ≈63°)",
            droplet_on_wall,
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_roundtrip_and_build() {
        for (name, _, sc) in presets() {
            let json = serde_json::to_string_pretty(&sc).unwrap();
            let back: Scenario = serde_json::from_str(&json).unwrap();
            assert_eq!(back.name, sc.name, "{name} roundtrip");
            build(&back).unwrap_or_else(|e| panic!("{name}: {e}"));
        }
    }

    /// Backward compatibility of the 3D-era schema: 2D scenarios neither
    /// require nor emit the new fields (`grid.nz`, `edges.front/back`,
    /// `compute`), so pre-existing JSON files and their serialised forms are
    /// unchanged.
    #[test]
    fn schema_2d_backward_compat() {
        // Old-style JSON (no new fields) parses, defaults to 2D.
        let sc: Scenario = serde_json::from_str(
            r#"{
                "name": "legacy",
                "grid": { "nx": 16, "ny": 12 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "periodic" }, "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" }
                },
                "run": { "steps": 1 }
            }"#,
        )
        .unwrap();
        assert_eq!(sc.grid.nz, 1);
        assert!(!sc.is_3d());
        // New fields stay invisible on serialisation of 2D scenarios.
        let json = serde_json::to_string(&sc).unwrap();
        for key in ["\"nz\"", "\"front\"", "\"back\"", "\"compute\"", "\"z\""] {
            assert!(
                !json.contains(key),
                "2D JSON must not contain {key}: {json}"
            );
        }
        // deny_unknown_fields still rejects typos.
        assert!(serde_json::from_str::<Scenario>(
            r#"{ "name": "x", "grid": { "nx": 3, "ny": 3, "nw": 4 },
                 "physics": { "nu": 0.05 },
                 "edges": { "left": {"type":"periodic"}, "right": {"type":"periodic"},
                            "bottom": {"type":"periodic"}, "top": {"type":"periodic"} },
                 "run": { "steps": 1 } }"#
        )
        .is_err());
    }

    #[test]
    fn web_export_shape_roundtrips_stably() {
        let exported = r#"{
            "version": 0,
            "name": "web-export-smoke",
            "grid": { "nx": 24, "ny": 18 },
            "physics": {
                "nu": 0.04,
                "collision": { "type": "bgk" },
                "force": [0.000001, 0.0],
                "precision": "f64"
            },
            "edges": {
                "left": { "type": "periodic" },
                "right": { "type": "periodic" },
                "bottom": { "type": "bounceBack" },
                "top": { "type": "movingWall", "u": [0.05, 0.0] }
            },
            "obstacles": [
                { "shape": "rect", "x0": 8, "y0": 7, "x1": 10, "y1": 9 }
            ],
            "init": { "kind": "rest" },
            "run": { "steps": 20000 },
            "outputs": [
                { "field": "speed", "format": "png", "every": 0 }
            ]
        }"#;
        let sc: Scenario = serde_json::from_str(exported).unwrap();
        build(&sc).unwrap();

        let first = serde_json::to_string(&sc).unwrap();
        let back: Scenario = serde_json::from_str(&first).unwrap();
        let second = serde_json::to_string(&back).unwrap();
        assert_eq!(first, second);
        for key in [
            "\"nz\"",
            "\"compute\"",
            "\"front\"",
            "\"back\"",
            "\"probes\"",
        ] {
            assert!(!second.contains(key), "web export gained {key}: {second}");
        }
    }

    fn duct3d() -> Scenario {
        serde_json::from_str(
            r#"{
                "name": "duct3d",
                "grid": { "nx": 12, "ny": 10, "nz": 10 },
                "physics": { "nu": 0.1, "force": [1e-6, 0.0] },
                "compute": { "backend": "cpu" },
                "edges": {
                    "left": { "type": "periodic" }, "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" },
                    "front": { "type": "bounceBack" }, "back": { "type": "bounceBack" }
                },
                "run": { "steps": 10 }
            }"#,
        )
        .unwrap()
    }

    #[test]
    fn build3d_runs_and_guards() {
        let sc = duct3d();
        assert!(sc.is_3d());
        // Builds and steps on the V2 core.
        match build3d(&sc).unwrap() {
            Sim3Handle::F64(mut s) => {
                s.run(3);
                let u = s.u(6, 5, 5);
                assert!(u[0].is_finite());
            }
            _ => panic!("expected f64"),
        }
        // 2D build refuses 3D scenarios; the dispatching check accepts them.
        assert!(build(&sc).is_err());
        assert!(build_check(&sc).is_ok());
        // gpu backend is rejected until the wgpu backend lands.
        let mut gpu = duct3d();
        gpu.compute = Some(ComputeSpec {
            backend: BackendSpec::Gpu,
        });
        assert!(matches!(
            build3d(&gpu),
            Err(Build3Error::BackendUnavailable(_))
        ));
        assert!(build_check(&gpu).is_err());
        assert!(validate(&gpu).iter().any(|w| w.field == "compute.backend"));
        // multiphase is 2D-only for now.
        let mut mp = duct3d();
        mp.multiphase = Some(MultiphaseSpec {
            g: -5.0,
            g_wall: 0.0,
            wall_rho: None,
        });
        assert!(matches!(build3d(&mp), Err(Build3Error::Unsupported(_))));
        // Unpaired z periodicity is a config error.
        let mut unpaired = duct3d();
        unpaired.edges.back = Some(EdgeSpec::Periodic);
        assert!(matches!(
            build3d(&unpaired),
            Err(Build3Error::Core(ConfigError::UnpairedPeriodic {
                axis: "z"
            }))
        ));
        // Open faces on two axes violate the corner rule.
        let mut cross = duct3d();
        cross.edges.left = EdgeSpec::VelocityInlet { u: [0.05, 0.0] };
        cross.edges.right = EdgeSpec::PressureOutlet { rho: 1.0 };
        cross.edges.front = Some(EdgeSpec::Outflow);
        cross.edges.back = Some(EdgeSpec::Outflow);
        assert!(matches!(
            build3d(&cross),
            Err(Build3Error::Core(ConfigError::AdjacentOpenEdges))
        ));
        // Spheres require a 3D grid.
        let mut sphere2d = duct3d();
        sphere2d.grid.nz = 1;
        sphere2d.obstacles = vec![Obstacle::Sphere {
            cx: 6.0,
            cy: 5.0,
            cz: 5.0,
            r: 2.0,
        }];
        assert!(build(&sphere2d).is_err());
    }

    #[test]
    fn explicit_gpu_backend_is_rejected_for_2d() {
        let mut sc = preset("cavity");
        sc.compute = Some(ComputeSpec {
            backend: BackendSpec::Gpu,
        });

        let err = match build(&sc) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("explicit GPU backend should fail for 2D scenarios"),
        };
        #[cfg(not(feature = "gpu"))]
        assert!(
            err.contains("requested backend \"gpu\" is unavailable"),
            "{err}"
        );
        #[cfg(feature = "gpu")]
        assert!(
            err.contains("requested backend \"gpu\" is unsupported for physics.precision f64"),
            "{err}"
        );
        #[cfg(not(feature = "gpu"))]
        assert!(err.contains("--features gpu"), "{err}");

        let check = build_check(&sc).unwrap_err();
        assert_eq!(check, err);
        assert!(validate(&sc).iter().any(|w| {
            w.field == "compute.backend"
                && (w.message
                    .contains("requested backend \"gpu\" is unavailable")
                    || w.message
                        .contains("requested backend \"gpu\" is unsupported"))
        }));
    }

    #[test]
    fn build_check_rejects_hard_scenario_physics_errors() {
        let mut bad_nu = preset("cavity");
        bad_nu.physics.nu = 0.0;
        let err = build_check(&bad_nu).unwrap_err();
        assert!(err.contains("kinematic viscosity must be > 0"), "{err}");
        assert!(err.contains("tau = 3*nu + 0.5"), "{err}");

        let mut too_fast_wall = preset("cavity");
        too_fast_wall.edges.top = EdgeSpec::MovingWall { u: [0.31, 0.0] };
        let err = build_check(&too_fast_wall).unwrap_err();
        assert!(err.contains("low-Mach limit"), "{err}");

        let mut adjacent_open = preset("cavity");
        adjacent_open.edges.left = EdgeSpec::VelocityInlet { u: [0.05, 0.0] };
        adjacent_open.edges.right = EdgeSpec::PressureOutlet { rho: 1.0 };
        adjacent_open.edges.bottom = EdgeSpec::Outflow;
        let err = build_check(&adjacent_open).unwrap_err();
        assert!(
            err.contains("orthogonal") || err.contains("adjacent") || err.contains("perpendicular"),
            "{err}"
        );

        let mut profile_not_on_inlet = preset("cavity");
        profile_not_on_inlet.inlet_profile = Some(InletProfile {
            edge: EdgeName::Top,
            kind: ProfileKind::Parabolic,
            umax: 0.05,
        });
        let err = build_check(&profile_not_on_inlet).unwrap_err();
        assert!(err.contains("inletProfile edge"), "{err}");

        let mut too_fast_profile = preset("cylinder-karman");
        too_fast_profile.inlet_profile = Some(InletProfile {
            edge: EdgeName::Left,
            kind: ProfileKind::Parabolic,
            umax: 0.31,
        });
        let err = build_check(&too_fast_profile).unwrap_err();
        assert!(err.contains("low-Mach limit"), "{err}");
    }

    #[test]
    fn unsupported_source_sink_schema_fields_are_rejected() {
        let with_sources = serde_json::json!({
            "name": "unsupported-sources",
            "grid": { "nx": 16, "ny": 16 },
            "physics": { "nu": 0.05 },
            "edges": {
                "left": { "type": "periodic" },
                "right": { "type": "periodic" },
                "bottom": { "type": "bounceBack" },
                "top": { "type": "bounceBack" }
            },
            "sources": [],
            "run": { "steps": 1 }
        });
        let err = serde_json::from_value::<Scenario>(with_sources)
            .unwrap_err()
            .to_string();
        assert!(err.contains("unknown field `sources`"), "{err}");

        let with_sinks = serde_json::json!({
            "name": "unsupported-sinks",
            "grid": { "nx": 16, "ny": 16 },
            "physics": { "nu": 0.05 },
            "edges": {
                "left": { "type": "periodic" },
                "right": { "type": "periodic" },
                "bottom": { "type": "bounceBack" },
                "top": { "type": "bounceBack" }
            },
            "sinks": [],
            "run": { "steps": 1 }
        });
        let err = serde_json::from_value::<Scenario>(with_sinks)
            .unwrap_err()
            .to_string();
        assert!(err.contains("unknown field `sinks`"), "{err}");
    }

    #[test]
    fn auto_backend_still_builds_and_runs_for_2d() {
        let mut sc = preset("cavity");
        sc.compute = Some(ComputeSpec {
            backend: BackendSpec::Auto,
        });

        match build(&sc).unwrap() {
            SimHandle::F64(mut sim, None) => {
                sim.run(2);
                assert_eq!(sim.time(), 2);
            }
            _ => panic!("expected an f64 single-phase CPU compat build"),
        }
        assert!(build_check(&sc).is_ok());
    }

    #[test]
    fn scenario_gravity_validation_rejects_bad_vectors_with_reason() {
        let wrong_dim = serde_json::from_str::<Scenario>(
            r#"{
                "name": "bad-gravity-dim",
                "grid": { "nx": 8, "ny": 8 },
                "physics": { "nu": 0.1, "gravity": [0.0, -1e-6] },
                "edges": {
                    "left": { "type": "periodic" }, "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" }
                },
                "run": { "steps": 1 }
            }"#,
        )
        .unwrap_err()
        .to_string();
        assert!(
            wrong_dim.contains("invalid length") || wrong_dim.contains("expected an array"),
            "{wrong_dim}"
        );

        let mut nan = preset("cavity");
        nan.physics.gravity = Some([f64::NAN, 0.0, 0.0]);
        let err = match build(&nan) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("NaN gravity should fail"),
        };
        assert!(err.contains("physics.gravity[0]"), "{err}");
        assert!(err.contains("finite"), "{err}");

        let mut bad_z = preset("cavity");
        bad_z.physics.gravity = Some([0.0, 0.0, 1.0e-6]);
        let err = match build(&bad_z) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("2D nonzero gz should fail"),
        };
        assert!(err.contains("physics.gravity[2]"), "{err}");
        assert!(err.contains("2D scenarios require gz == 0"), "{err}");
    }

    /// The z-periodic 3D parabolic inlet degenerates to the 2D parabola:
    /// the built profile must drive the same inlet-node velocities.
    #[test]
    fn inlet_profile_3d_product_form() {
        let sc: Scenario = serde_json::from_str(
            r#"{
                "name": "duct-inlet",
                "grid": { "nx": 12, "ny": 10, "nz": 10 },
                "physics": { "nu": 0.1 },
                "edges": {
                    "left": { "type": "velocityInlet", "u": [0.0, 0.0] },
                    "right": { "type": "pressureOutlet", "rho": 1.0 },
                    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" },
                    "front": { "type": "bounceBack" }, "back": { "type": "bounceBack" }
                },
                "inletProfile": { "edge": "left", "kind": "parabolic", "umax": 0.1 },
                "run": { "steps": 2 }
            }"#,
        )
        .unwrap();
        match build3d(&sc).unwrap() {
            Sim3Handle::F64(mut s) => {
                s.run(2);
                // Duct-type product profile: node (y, z) carries
                // umax f(y) f(z), enforced exactly by the Zou-He face.
                let ny = 10usize;
                let h = (ny - 2) as f64;
                let fac = |c: usize| {
                    let w = c as f64 - 0.5;
                    4.0 * w * (h - w) / (h * h)
                };
                for (y, z) in [(4, 5), (1, 1), (5, 2)] {
                    let expect = 0.1 * fac(y) * fac(z);
                    let got = s.u(0, y, z)[0];
                    assert!(
                        (got - expect).abs() < 1e-13,
                        "inlet ({y},{z}): got {got}, expect {expect}"
                    );
                }
            }
            _ => panic!("expected f64"),
        }
    }

    #[test]
    fn validate_flags_dangerous_settings() {
        let (_, _, mut sc) = presets().remove(0);
        sc.physics.nu = 0.005;
        sc.edges.top = EdgeSpec::MovingWall { u: [0.2, 0.0] };
        let warnings = validate(&sc);
        assert!(
            warnings.iter().any(|w| w.field == "physics"),
            "{warnings:?}"
        );
        assert!(warnings.iter().any(|w| w.field == "edges"), "{warnings:?}");
    }

    fn preset(name: &str) -> Scenario {
        presets()
            .into_iter()
            .find(|(n, _, _)| *n == name)
            .unwrap_or_else(|| panic!("preset {name} not found"))
            .2
    }

    #[test]
    fn convective_outflow_roundtrip_and_hints() {
        // camelCase JSON tag/field
        let spec: EdgeSpec =
            serde_json::from_str(r#"{ "type": "convectiveOutflow", "uConv": 0.08 }"#).unwrap();
        assert!(matches!(spec, EdgeSpec::ConvectiveOutflow { u_conv } if u_conv == 0.08));
        let text = serde_json::to_string(&spec).unwrap();
        assert!(text.contains("\"uConv\":0.08"), "{text}");

        // valid uConv: builds, no edge warnings
        let mut sc = preset("cylinder-karman");
        sc.edges.right = EdgeSpec::ConvectiveOutflow { u_conv: 0.1 };
        build(&sc).unwrap();
        assert!(
            validate(&sc).iter().all(|w| !w.field.starts_with("edges.")),
            "{:?}",
            validate(&sc)
        );

        // uConv out of (0,1]: validate warns with a hint, core build rejects
        for bad in [0.0, -0.1, 1.5] {
            sc.edges.right = EdgeSpec::ConvectiveOutflow { u_conv: bad };
            let warnings = validate(&sc);
            assert!(
                warnings.iter().any(|w| w.field == "edges.right"),
                "uConv={bad}: {warnings:?}"
            );
            assert!(build(&sc).is_err(), "uConv={bad} should fail to build");
        }
    }

    #[test]
    fn wall_rho_wires_into_shan_chen() {
        let sc = preset("droplet-on-wall");
        match build(&sc).unwrap() {
            SimHandle::F64(_, Some(mp)) => assert_eq!(mp.wall_rho, Some(1.0)),
            _ => panic!("expected an f64 multiphase build"),
        }
        // omitted wallRho stays None (legacy scenarios unchanged)
        let sc = preset("two-phase-droplet");
        match build(&sc).unwrap() {
            SimHandle::F64(_, Some(mp)) => assert_eq!(mp.wall_rho, None),
            _ => panic!("expected an f64 multiphase build"),
        }
    }
}
