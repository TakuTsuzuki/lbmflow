//! Scenario JSON (v0): the single execution contract shared by the CLI,
//! the MCP server and (in spirit) the GUI presets.
//!
//! See `docs/AGENT_MODE_DESIGN.md` for the schema rationale. Field names are
//! camelCase in JSON.

use lbm_core::multiphase::ShanChen;
use lbm_core::prelude::*;
use serde::{Deserialize, Serialize};

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
    pub compute: Option<ComputeSpec>,
    pub edges: EdgesSpec,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inlet_profile: Option<InletProfile>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub obstacles: Vec<Obstacle>,
    #[serde(default)]
    pub init: InitSpec,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub multiphase: Option<MultiphaseSpec>,
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

/// Compute-target selection (ARCHITECTURE_V2 §3). All fields optional; the
/// 3D engine currently runs on the CPU backend only ("gpu" is rejected at
/// build time until the wgpu backend lands).
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

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct Physics {
    /// Kinematic viscosity (lattice units); tau = 3 nu + 0.5.
    pub nu: f64,
    #[serde(default)]
    pub collision: CollisionSpec,
    #[serde(default)]
    pub force: [f64; 2],
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
    MovingWall { u: [f64; 2] },
    VelocityInlet { u: [f64; 2] },
    PressureOutlet { rho: f64 },
    Outflow,
    /// Convective (radiation) outflow: far less pressure-reflective than
    /// `Outflow`. `uConv` is the expected mean outflow speed, in (0, 1].
    #[serde(rename_all = "camelCase")]
    ConvectiveOutflow { u_conv: f64 },
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
    Circle {
        cx: f64,
        cy: f64,
        r: f64,
    },
    /// 2D: a rectangle. 3D: extruded along z (a box through the domain).
    Rect {
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
    },
    /// 3D only: a solid ball (staircase approximation).
    Sphere {
        cx: f64,
        cy: f64,
        cz: f64,
        r: f64,
    },
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
    if tau < 0.55 {
        warn(
            "physics.nu",
            format!("tau = {tau:.3} は安定限界に近い（0.55 未満）。粘性を上げるか解像度を上げてください"),
        );
    }
    let max_edge_speed = edge_speeds(&sc.edges).into_iter().fold(0.0, f64::max);
    if max_edge_speed > 0.15 {
        warn(
            "edges",
            format!("流入/壁速度 {max_edge_speed:.3} は圧縮性誤差が目立つ水準（0.15 超）"),
        );
    }
    if max_edge_speed > 0.0 && sc.physics.nu > 0.0 {
        let grid_re = max_edge_speed / sc.physics.nu;
        if grid_re > 15.0 {
            warn(
                "physics",
                format!(
                    "グリッドレイノルズ数 U/ν = {grid_re:.1} > 15: 発散の恐れ（PHYSICS.md 参照）"
                ),
            );
        }
    }
    if sc.multiphase.is_some() && sc.physics.precision == Precision::F32 {
        warn(
            "physics.precision",
            "多相流は f64 を推奨（界面の急勾配に対する余裕）".to_string(),
        );
    }
    if let Some(mp) = &sc.multiphase {
        if mp.g > -4.0 {
            warn(
                "multiphase.g",
                format!(
                    "G = {} は臨界値 -4 より弱く、相分離しません（推奨 -5.0）",
                    mp.g
                ),
            );
        }
    }
    if sc.is_3d() {
        if sc.multiphase.is_some() {
            warn(
                "multiphase",
                "3D (nz > 1) は多相流未対応です（構築時エラーになります）".to_string(),
            );
        }
        if matches!(
            sc.compute,
            Some(ComputeSpec {
                backend: BackendSpec::Gpu
            })
        ) {
            warn(
                "compute.backend",
                "gpu バックエンドは未提供です（cpu / auto を指定。構築時エラーになります）"
                    .to_string(),
            );
        }
    } else if sc.edges.front.is_some() || sc.edges.back.is_some() {
        warn(
            "edges",
            "front/back は 3D (nz > 1) 専用で、2D では無視されます".to_string(),
        );
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
                        "uConv = {u_conv} は (0,1] の範囲外で、構築時にエラーになります。\
                         期待される平均流出速度（例: 流入速度と同程度の 0.05〜0.15）を指定してください"
                    ),
                );
            }
        }
    }
    warnings
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

// ---------------------------------------------------------------- build

/// A built simulation, precision-erased for the runner.
pub enum SimHandle {
    F32(Simulation<f32>, Option<ShanChen<f32>>),
    F64(Simulation<f64>, Option<ShanChen<f64>>),
}

// ---------------------------------------------------------------- build (3D)

/// The 3D engine type behind a scenario: V2 core, D3Q19, CPU backend,
/// monolithic decomposition (ARCHITECTURE_V2; `compute.backend: "cpu"`).
pub type Solver3<T> = lbm_core2::solver::Solver<
    lbm_core2::lattice::D3Q19,
    T,
    lbm_core2::backend::CpuScalar,
    lbm_core2::halo::LocalPeriodic,
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
    /// Invalid physical/boundary configuration.
    Core(ConfigError),
    /// Feature not available on the 3D engine (message is user-facing).
    Unsupported(&'static str),
}

impl std::fmt::Display for Build3Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Build3Error::Core(e) => write!(f, "{e}"),
            Build3Error::Unsupported(what) => write!(f, "3D (nz > 1) では未対応: {what}"),
        }
    }
}

impl std::error::Error for Build3Error {}

impl From<ConfigError> for Build3Error {
    fn from(e: ConfigError) -> Self {
        Build3Error::Core(e)
    }
}

/// Dimension-dispatching build check for validators (CLI `validate`, MCP
/// `validate_scenario`): construct the simulation the same way `run` would
/// (2D or 3D) and report the error text, discarding the handle.
pub fn build_check(sc: &Scenario) -> Result<(), String> {
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
    Ok(match sc.physics.precision {
        Precision::F32 => Sim3Handle::F32(build3d_t::<f32>(sc)?),
        Precision::F64 => Sim3Handle::F64(build3d_t::<f64>(sc)?),
    })
}

fn build3d_t<T: lbm_core2::real::Real>(sc: &Scenario) -> Result<Solver3<T>, Build3Error> {
    use lbm_core2::prelude::{
        build_wall_rims, CollisionKind, CpuScalar, Face, FaceBC, GlobalSpec, LocalPeriodic,
        Solver, WallSpec,
    };

    assert!(sc.is_3d(), "build3d requires grid.nz > 1");
    if sc.multiphase.is_some() {
        return Err(Build3Error::Unsupported("multiphase（多相流）"));
    }
    if !matches!(sc.init, InitSpec::Rest) {
        return Err(Build3Error::Unsupported("init は rest のみ"));
    }
    if let Some(c) = &sc.compute {
        if c.backend == BackendSpec::Gpu {
            return Err(Build3Error::Unsupported(
                "compute.backend \"gpu\"（cpu / auto を指定してください）",
            ));
        }
    }
    let dims = [sc.grid.nx, sc.grid.ny, sc.grid.nz];
    if dims[0] < 3 || dims[1] < 3 {
        return Err(ConfigError::DomainTooSmall {
            nx: dims[0],
            ny: dims[1],
        }
        .into());
    }
    if dims[2] < 3 {
        return Err(ConfigError::InvalidParameter {
            what: "grid.nz (3D requires nz >= 3)",
            value: dims[2] as f64,
        }
        .into());
    }
    if sc.physics.nu <= 0.0 {
        return Err(ConfigError::NonPositiveViscosity { nu: sc.physics.nu }.into());
    }

    let specs = face_specs(&sc.edges);
    // Periodic pairing per axis.
    let mut periodic = [false; 3];
    for (axis, name) in [(0usize, "x"), (1, "y"), (2, "z")] {
        let lo = matches!(specs[2 * axis], EdgeSpec::Periodic);
        let hi = matches!(specs[2 * axis + 1], EdgeSpec::Periodic);
        if lo != hi {
            return Err(ConfigError::UnpairedPeriodic { axis: name }.into());
        }
        periodic[axis] = lo && hi;
    }
    // Open faces must not share a domain edge (V1's corner rule, lifted to
    // 3D): all open faces on one axis only.
    let is_open = |s: &EdgeSpec| {
        matches!(
            s,
            EdgeSpec::VelocityInlet { .. }
                | EdgeSpec::PressureOutlet { .. }
                | EdgeSpec::Outflow
                | EdgeSpec::ConvectiveOutflow { .. }
        )
    };
    let open_axes: Vec<usize> = (0..3)
        .filter(|a| is_open(&specs[2 * a]) || is_open(&specs[2 * a + 1]))
        .collect();
    if open_axes.len() > 1 {
        return Err(ConfigError::AdjacentOpenEdges.into());
    }
    // Speed / density / parameter limits (2D `SimConfig::build` semantics).
    let speed_of = |u: [f64; 2]| (u[0] * u[0] + u[1] * u[1]).sqrt();
    for s in &specs {
        match *s {
            EdgeSpec::MovingWall { u } | EdgeSpec::VelocityInlet { u } => {
                let sp = speed_of(u);
                if sp > MAX_SPEED {
                    return Err(ConfigError::VelocityTooHigh { speed: sp }.into());
                }
            }
            EdgeSpec::PressureOutlet { rho } => {
                if rho <= 0.0 {
                    return Err(ConfigError::NonPositiveDensity { rho }.into());
                }
            }
            EdgeSpec::ConvectiveOutflow { u_conv } => {
                if !(u_conv > 0.0 && u_conv <= 1.0) {
                    return Err(ConfigError::InvalidParameter {
                        what: "u_conv",
                        value: u_conv,
                    }
                    .into());
                }
            }
            _ => {}
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
    };
    let (solid, wall_u) = build_wall_rims::<T>(3, dims, &walls);
    let mut s: Solver3<T> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );

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
                "inletProfile は velocityInlet の辺にのみ指定できます",
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
                "force プローブには obstacles が必要です",
            ));
        }
        // Probe all obstacle solids (cells strictly inside the domain box,
        // rims excluded) — 2D convention lifted to 3D.
        let solid: Vec<bool> = (0..dims[2])
            .flat_map(|z| {
                (0..dims[1]).flat_map(move |y| (0..dims[0]).map(move |x| (x, y, z)))
            })
            .map(|(x, y, z)| s.is_solid(x, y, z))
            .collect();
        let (nx, ny) = (dims[0], dims[1]);
        let rim = move |c: usize, n: usize| c == 0 || c == n - 1;
        s.set_force_probe(move |x, y, z| {
            !rim(x, dims[0]) && !rim(y, dims[1]) && !rim(z, dims[2])
                && solid[(z * ny + y) * nx + x]
        });
    }

    Ok(s)
}

/// Build the 2D simulation (+ optional multiphase driver) from a scenario.
/// Scenarios with `grid.nz > 1` must go through [`build3d`] instead.
pub fn build(sc: &Scenario) -> Result<SimHandle, ConfigError> {
    if sc.is_3d() {
        return Err(ConfigError::InvalidParameter {
            what: "grid.nz (2D build requires nz == 1; the runner dispatches 3D to build3d)",
            value: sc.grid.nz as f64,
        });
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

    for ob in &sc.obstacles {
        match *ob {
            Obstacle::Circle { cx, cy, r } => {
                let r2 = r * r;
                sim.set_solid_region(|x, y| {
                    let dx = x as f64 - cx;
                    let dy = y as f64 - cy;
                    dx * dx + dy * dy <= r2
                });
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

// ---------------------------------------------------------------- presets

/// Built-in presets: (name, description, scenario JSON factory).
pub fn presets() -> Vec<(&'static str, &'static str, Scenario)> {
    let cavity = Scenario {
        version: 0,
        name: "cavity".into(),
        grid: Grid { nx: 128, ny: 128, nz: 1 },
        physics: Physics {
            nu: 0.02,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            precision: Precision::F64,
        },
        compute: None,
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
        grid: Grid { nx: 440, ny: 164, nz: 1 },
        physics: Physics {
            nu: 0.04,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            precision: Precision::F64,
        },
        compute: None,
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
        grid: Grid { nx: 128, ny: 128, nz: 1 },
        physics: Physics {
            nu: 1.0 / 6.0,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            precision: Precision::F64,
        },
        compute: None,
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
        grid: Grid { nx: 160, ny: 100, nz: 1 },
        physics: Physics {
            nu: 1.0 / 6.0,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            precision: Precision::F64,
        },
        compute: None,
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
        ("cavity", "リッド駆動キャビティ（定常判定つき）", cavity),
        (
            "cylinder-karman",
            "円柱まわりのカルマン渦列 + 抗力プローブ",
            cylinder,
        ),
        ("two-phase-droplet", "Shan-Chen 二相液滴の平衡化", droplet),
        (
            "droplet-on-wall",
            "壁上液滴の接触角デモ（仮想壁密度 wallRho=1.0 → θ≈63°）",
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
            assert!(!json.contains(key), "2D JSON must not contain {key}: {json}");
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
        assert!(matches!(build3d(&gpu), Err(Build3Error::Unsupported(_))));
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
            Err(Build3Error::Core(ConfigError::UnpairedPeriodic { axis: "z" }))
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
