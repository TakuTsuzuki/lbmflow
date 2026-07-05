//! Scenario JSON (v0): the single execution contract shared by the CLI,
//! the MCP server and (in spirit) the GUI presets.
//!
//! See `docs/AGENT_MODE_DESIGN.md` for the schema rationale. Field names are
//! camelCase in JSON.

use lbm_core::multiphase::{Psi, ShanChen};
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
    Circle {
        cx: f64,
        cy: f64,
        r: f64,
    },
    Rect {
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
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
    /// Point time series of (ux, uy, rho).
    #[serde(rename_all = "camelCase")]
    Point { x: usize, y: usize, every: usize },
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
    warnings
}

fn edge_speeds(e: &EdgesSpec) -> [f64; 4] {
    [e.left, e.right, e.bottom, e.top].map(|s| match s {
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

/// Build the simulation (+ optional multiphase driver) from a scenario.
pub fn build(sc: &Scenario) -> Result<SimHandle, ConfigError> {
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

    let mp = sc
        .multiphase
        .as_ref()
        .map(|m| ShanChen::<T>::new(m.g).with_wall(m.g_wall));
    Ok((sim, mp))
}

// ---------------------------------------------------------------- presets

/// Built-in presets: (name, description, scenario JSON factory).
pub fn presets() -> Vec<(&'static str, &'static str, Scenario)> {
    let cavity = Scenario {
        version: 0,
        name: "cavity".into(),
        grid: Grid { nx: 128, ny: 128 },
        physics: Physics {
            nu: 0.02,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            precision: Precision::F64,
        },
        edges: EdgesSpec {
            left: EdgeSpec::BounceBack,
            right: EdgeSpec::BounceBack,
            bottom: EdgeSpec::BounceBack,
            top: EdgeSpec::MovingWall { u: [0.1, 0.0] },
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
        grid: Grid { nx: 440, ny: 164 },
        physics: Physics {
            nu: 0.04,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            precision: Precision::F64,
        },
        edges: EdgesSpec {
            left: EdgeSpec::VelocityInlet { u: [0.1, 0.0] },
            right: EdgeSpec::PressureOutlet { rho: 1.0 },
            bottom: EdgeSpec::BounceBack,
            top: EdgeSpec::BounceBack,
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
        grid: Grid { nx: 128, ny: 128 },
        physics: Physics {
            nu: 1.0 / 6.0,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            precision: Precision::F64,
        },
        edges: EdgesSpec {
            left: EdgeSpec::Periodic,
            right: EdgeSpec::Periodic,
            bottom: EdgeSpec::Periodic,
            top: EdgeSpec::Periodic,
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
    vec![
        ("cavity", "リッド駆動キャビティ（定常判定つき）", cavity),
        (
            "cylinder-karman",
            "円柱まわりのカルマン渦列 + 抗力プローブ",
            cylinder,
        ),
        ("two-phase-droplet", "Shan-Chen 二相液滴の平衡化", droplet),
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
}
