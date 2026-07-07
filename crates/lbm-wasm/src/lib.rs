//! wasm-bindgen bridge: exposes `lbm_core::compat::Simulation<f32>` behind
//! the TypeScript `Engine` interface (see `web/src/engine/types.ts`).
//!
//! Field access is zero-copy: the `*_ptr` methods return pointers into wasm
//! linear memory; the JS adapter wraps them in `Float32Array` views that are
//! valid until the next `step`/`init` call.

use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use serde::{Deserialize, Serialize};
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
type BridgeError = JsError;
#[cfg(not(target_arch = "wasm32"))]
type BridgeError = String;

#[cfg(target_arch = "wasm32")]
fn bridge_error(message: String) -> BridgeError {
    JsError::new(&message)
}

#[cfg(not(target_arch = "wasm32"))]
fn bridge_error(message: String) -> BridgeError {
    message
}

#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
enum JsEdge {
    Periodic,
    BounceBack,
    MovingWall { u: [f32; 2] },
    VelocityInlet { u: [f32; 2] },
    PressureOutlet { rho: f32 },
    Outflow,
    ConvectiveOutflow { u_conv: f32 },
}

impl JsEdge {
    fn to_core(&self) -> EdgeBC<f32> {
        match self {
            JsEdge::Periodic => EdgeBC::Periodic,
            JsEdge::BounceBack => EdgeBC::BounceBack,
            JsEdge::MovingWall { u } => EdgeBC::MovingWall { u: *u },
            JsEdge::VelocityInlet { u } => EdgeBC::VelocityInlet { u: *u },
            JsEdge::PressureOutlet { rho } => EdgeBC::PressureOutlet { rho: *rho },
            JsEdge::Outflow => EdgeBC::Outflow,
            JsEdge::ConvectiveOutflow { u_conv } => EdgeBC::ConvectiveOutflow { u_conv: *u_conv },
        }
    }
}

#[derive(Deserialize)]
struct JsEdges {
    left: JsEdge,
    right: JsEdge,
    bottom: JsEdge,
    top: JsEdge,
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct JsMultiphase {
    g: f64,
    #[serde(default)]
    g_wall: f64,
}

#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "camelCase")]
enum JsInit {
    Rest,
    #[serde(rename_all = "camelCase")]
    Droplet {
        cx: f64,
        cy: f64,
        r: f64,
        rho_liquid: f64,
        rho_vapor: f64,
    },
    #[cfg(test)]
    #[serde(rename_all = "camelCase")]
    TaylorGreen {
        amplitude: f64,
    },
}

#[derive(Deserialize)]
struct JsConfig {
    nx: usize,
    ny: usize,
    nu: f64,
    collision: String,
    edges: JsEdges,
    force: [f32; 2],
    #[serde(default)]
    multiphase: Option<JsMultiphase>,
    #[serde(default)]
    init: Option<JsInit>,
}

#[derive(Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct JsScenario {
    #[serde(default)]
    version: u32,
    name: String,
    grid: JsScenarioGrid,
    physics: JsScenarioPhysics,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    units: Option<serde_json::Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    compute: Option<serde_json::Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    wall: Option<serde_json::Value>,
    edges: JsScenarioEdges,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    inlet_profile: Option<JsScenarioInletProfile>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    obstacles: Vec<JsScenarioObstacle>,
    #[serde(default)]
    init: JsScenarioInit,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    multiphase: Option<JsMultiphase>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    rotor: Option<serde_json::Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    particles: Option<serde_json::Value>,
    run: JsScenarioRun,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    probes: Vec<serde_json::Value>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    outputs: Vec<serde_json::Value>,
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct JsScenarioGrid {
    nx: usize,
    ny: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    nz: Option<usize>,
}

impl JsScenarioGrid {
    fn nz(self) -> usize {
        self.nz.unwrap_or(1)
    }
}

#[derive(Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct JsScenarioPhysics {
    nu: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    collision: Option<JsScenarioCollision>,
    #[serde(default)]
    force: [f64; 2],
    #[serde(default, skip_serializing_if = "Option::is_none")]
    gravity: Option<[f64; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    precision: Option<JsPrecision>,
}

impl JsScenarioPhysics {
    fn precision(&self) -> JsPrecision {
        self.precision.unwrap_or_default()
    }

    fn collision(&self) -> Collision {
        match self.collision {
            Some(JsScenarioCollision::Bgk) => Collision::Bgk,
            Some(JsScenarioCollision::Trt { magic }) => Collision::Trt {
                magic: magic.unwrap_or(Collision::MAGIC_STD),
            },
            None => Collision::default(),
        }
    }
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
enum JsScenarioCollision {
    Bgk,
    #[serde(rename_all = "camelCase")]
    Trt {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        magic: Option<f64>,
    },
}

#[derive(Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
enum JsPrecision {
    F32,
    #[default]
    F64,
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
enum JsScenarioEdge {
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
    #[serde(rename_all = "camelCase")]
    ConvectiveOutflow {
        u_conv: f64,
    },
}

impl JsScenarioEdge {
    fn to_core_f32(self) -> EdgeBC<f32> {
        match self {
            JsScenarioEdge::Periodic => EdgeBC::Periodic,
            JsScenarioEdge::BounceBack => EdgeBC::BounceBack,
            JsScenarioEdge::MovingWall { u } => EdgeBC::MovingWall {
                u: [u[0] as f32, u[1] as f32],
            },
            JsScenarioEdge::VelocityInlet { u } => EdgeBC::VelocityInlet {
                u: [u[0] as f32, u[1] as f32],
            },
            JsScenarioEdge::PressureOutlet { rho } => EdgeBC::PressureOutlet { rho: rho as f32 },
            JsScenarioEdge::Outflow => EdgeBC::Outflow,
            JsScenarioEdge::ConvectiveOutflow { u_conv } => EdgeBC::ConvectiveOutflow {
                u_conv: u_conv as f32,
            },
        }
    }
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct JsScenarioEdges {
    left: JsScenarioEdge,
    right: JsScenarioEdge,
    bottom: JsScenarioEdge,
    top: JsScenarioEdge,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    front: Option<JsScenarioEdge>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    back: Option<JsScenarioEdge>,
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct JsScenarioInletProfile {
    edge: JsScenarioEdgeName,
    kind: JsScenarioProfileKind,
    umax: f64,
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
enum JsScenarioEdgeName {
    Left,
    Right,
    Bottom,
    Top,
}

impl JsScenarioEdgeName {
    fn to_core(self) -> Edge {
        match self {
            JsScenarioEdgeName::Left => Edge::Left,
            JsScenarioEdgeName::Right => Edge::Right,
            JsScenarioEdgeName::Bottom => Edge::Bottom,
            JsScenarioEdgeName::Top => Edge::Top,
        }
    }
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
enum JsScenarioProfileKind {
    Parabolic,
}

#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(tag = "shape", rename_all = "camelCase")]
enum JsScenarioObstacle {
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
    Sphere {
        cx: f64,
        cy: f64,
        cz: f64,
        r: f64,
    },
}

#[derive(Clone, Copy, Default, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "camelCase")]
enum JsScenarioInit {
    #[default]
    Rest,
    #[serde(rename_all = "camelCase")]
    Droplet {
        cx: f64,
        cy: f64,
        r: f64,
        rho_liquid: f64,
        rho_vapor: f64,
    },
    #[serde(rename_all = "camelCase")]
    Pool {
        height_frac: f64,
        rho_liquid: f64,
        rho_vapor: f64,
    },
}

#[derive(Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct JsScenarioRun {
    steps: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    stop_when_steady: Option<serde_json::Value>,
}

/// Browser-facing simulation handle.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub struct WasmSim {
    sim: Option<Simulation<f32>>,
    /// Shan-Chen driver when the config declares multiphase.
    multiphase: Option<ShanChen<f32>>,
    /// u8 mirror of the solid mask for cheap JS-side rendering.
    solid_u8: Vec<u8>,
    /// User-painted obstacles (kept separately so erasing can rebuild).
    painted: Vec<u8>,
    cfg_json: String,
}

fn build_sim(cfg_json: &str) -> Result<(Simulation<f32>, Option<ShanChen<f32>>), BridgeError> {
    let cfg: JsConfig = serde_json::from_str(cfg_json)
        .map_err(|e| bridge_error(format!("cannot parse config JSON: {e}")))?;
    let collision = match cfg.collision.as_str() {
        "bgk" => Collision::Bgk,
        _ => Collision::default(),
    };
    let mut sim = SimConfig {
        nx: cfg.nx,
        ny: cfg.ny,
        nu: cfg.nu,
        collision,
        edges: Edges {
            left: cfg.edges.left.to_core(),
            right: cfg.edges.right.to_core(),
            bottom: cfg.edges.bottom.to_core(),
            top: cfg.edges.top.to_core(),
        },
        force: cfg.force,
    }
    .build()
    .map_err(|e| bridge_error(format!("config error: {e}")))?;
    match cfg.init {
        Some(JsInit::Droplet {
            cx,
            cy,
            r,
            rho_liquid,
            rho_vapor,
        }) => {
            let r2 = r * r;
            sim.init_with(|x, y| {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let rho = if dx * dx + dy * dy <= r2 {
                    rho_liquid
                } else {
                    rho_vapor
                };
                (rho as f32, 0.0, 0.0)
            });
        }
        #[cfg(test)]
        Some(JsInit::TaylorGreen { amplitude }) => {
            init_taylor_green(&mut sim, amplitude as f32);
        }
        Some(JsInit::Rest) | None => {}
    }
    let mp = cfg
        .multiphase
        .map(|m| ShanChen::<f32>::new(m.g).with_wall(m.g_wall));
    Ok((sim, mp))
}

fn parse_scenario(scenario_json: &str) -> Result<JsScenario, BridgeError> {
    serde_json::from_str(scenario_json)
        .map_err(|e| bridge_error(format!("cannot parse scenario JSON: {e}")))
}

fn build_sim_from_scenario(
    scenario_json: &str,
) -> Result<(Simulation<f32>, Option<ShanChen<f32>>), BridgeError> {
    let sc = parse_scenario(scenario_json)?;
    if sc.grid.nz() != 1 {
        return Err(bridge_error(
            "wasm scenario bridge supports 2D f32 scenarios only; grid.nz must be omitted or 1"
                .to_string(),
        ));
    }
    if sc.physics.precision() != JsPrecision::F32 {
        return Err(bridge_error(
            "wasm scenario bridge supports physics.precision \"f32\" only; f64 is not available in the browser bridge".to_string(),
        ));
    }
    if sc.wall.is_some() {
        return Err(bridge_error(
            "wasm scenario bridge does not support scenario wall models".to_string(),
        ));
    }
    if sc.rotor.is_some() {
        return Err(bridge_error(
            "wasm scenario bridge does not support rotor scenarios".to_string(),
        ));
    }
    if sc.particles.is_some() {
        return Err(bridge_error(
            "wasm scenario bridge does not support particle scenarios".to_string(),
        ));
    }
    if let Some(gravity) = sc.physics.gravity {
        if gravity[2] != 0.0 {
            return Err(bridge_error(
                "wasm scenario bridge supports 2D gravity only; physics.gravity[2] must be 0"
                    .to_string(),
            ));
        }
    }

    let mut sim = SimConfig {
        nx: sc.grid.nx,
        ny: sc.grid.ny,
        nu: sc.physics.nu,
        collision: sc.physics.collision(),
        edges: Edges {
            left: sc.edges.left.to_core_f32(),
            right: sc.edges.right.to_core_f32(),
            bottom: sc.edges.bottom.to_core_f32(),
            top: sc.edges.top.to_core_f32(),
        },
        force: [sc.physics.force[0] as f32, sc.physics.force[1] as f32],
    }
    .build()
    .map_err(|e| bridge_error(format!("scenario config error: {e}")))?;

    if let Some(gravity) = sc.physics.gravity {
        sim.set_gravity([gravity[0] as f32, gravity[1] as f32]);
    }

    for obstacle in &sc.obstacles {
        match *obstacle {
            JsScenarioObstacle::Circle { cx, cy, r } => {
                let r2 = r * r;
                sim.set_solid_region(|x, y| {
                    let dx = x as f64 - cx;
                    let dy = y as f64 - cy;
                    dx * dx + dy * dy <= r2
                });
            }
            JsScenarioObstacle::Rect { x0, y0, x1, y1 } => {
                sim.set_solid_region(|x, y| x >= x0 && x <= x1 && y >= y0 && y <= y1);
            }
            JsScenarioObstacle::Sphere { .. } => {
                return Err(bridge_error(
                    "wasm scenario bridge supports 2D obstacles only; sphere requires a 3D engine"
                        .to_string(),
                ));
            }
        }
    }

    if let Some(profile) = sc.inlet_profile {
        match profile.kind {
            JsScenarioProfileKind::Parabolic => {
                let edge = profile.edge.to_core();
                let (nx, ny) = (sim.nx(), sim.ny());
                let len = match edge {
                    Edge::Left | Edge::Right => ny,
                    Edge::Bottom | Edge::Top => nx,
                };
                let h = (len - 2) as f64;
                let umax = profile.umax;
                let normal_sign: [f64; 2] = match edge {
                    Edge::Left => [1.0, 0.0],
                    Edge::Right => [-1.0, 0.0],
                    Edge::Bottom => [0.0, 1.0],
                    Edge::Top => [0.0, -1.0],
                };
                sim.set_inlet_profile(edge, move |c| {
                    if c == 0 || c as f64 >= h + 1.0 {
                        return [0.0, 0.0];
                    }
                    let yw = c as f64 - 0.5;
                    let mag = 4.0 * umax * yw * (h - yw) / (h * h);
                    [(mag * normal_sign[0]) as f32, (mag * normal_sign[1]) as f32]
                });
            }
        }
    }

    match sc.init {
        JsScenarioInit::Rest => {}
        JsScenarioInit::Droplet {
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
                (rho as f32, 0.0, 0.0)
            });
        }
        JsScenarioInit::Pool { .. } => {
            return Err(bridge_error(
                "wasm scenario bridge does not support pool initialization".to_string(),
            ));
        }
    }

    let mp = sc
        .multiphase
        .map(|m| ShanChen::<f32>::new(m.g).with_wall(m.g_wall));
    Ok((sim, mp))
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn normalize_scenario_json(scenario_json: &str) -> Result<String, BridgeError> {
    let sc = parse_scenario(scenario_json)?;
    serde_json::to_string(&sc)
        .map_err(|e| bridge_error(format!("cannot serialize scenario JSON: {e}")))
}

#[cfg(test)]
fn init_taylor_green(sim: &mut Simulation<f32>, amplitude: f32) {
    let nx = sim.nx() as f64;
    let ny = sim.ny() as f64;
    let kx = std::f64::consts::TAU / nx;
    let ky = std::f64::consts::TAU / ny;
    let u0 = amplitude as f64;
    sim.init_with(|x, y| {
        let x = x as f64;
        let y = y as f64;
        let rho = 1.0 - (3.0 * u0 * u0 / 4.0) * ((2.0 * kx * x).cos() + (2.0 * ky * y).cos());
        let ux = -u0 * (kx * x).cos() * (ky * y).sin();
        let uy = u0 * (kx * x).sin() * (ky * y).cos();
        (rho as f32, ux as f32, uy as f32)
    });
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
impl WasmSim {
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen(constructor))]
    pub fn new() -> WasmSim {
        WasmSim {
            sim: None,
            multiphase: None,
            solid_u8: Vec::new(),
            painted: Vec::new(),
            cfg_json: String::new(),
        }
    }

    /// (Re)initialise from an EngineConfig JSON string.
    pub fn init(&mut self, cfg_json: &str) -> Result<(), BridgeError> {
        let (sim, mp) = build_sim(cfg_json)?;
        let n = sim.nx() * sim.ny();
        self.cfg_json = cfg_json.to_string();
        self.painted = vec![0; n];
        self.sim = Some(sim);
        self.multiphase = mp;
        self.refresh_solid_mirror();
        Ok(())
    }

    /// Initialise from the shared scenario JSON schema. The browser bridge is
    /// intentionally restricted to the 2D f32 subset used by the GUI.
    pub fn init_scenario(&mut self, scenario_json: &str) -> Result<(), BridgeError> {
        let (sim, mp) = build_sim_from_scenario(scenario_json)?;
        let n = sim.nx() * sim.ny();
        self.cfg_json = scenario_json.to_string();
        self.painted = vec![0; n];
        self.sim = Some(sim);
        self.multiphase = mp;
        self.refresh_solid_mirror();
        Ok(())
    }

    fn refresh_solid_mirror(&mut self) {
        let solid: Vec<u8> = self
            .sim
            .as_ref()
            .map(|s| s.solid_field().iter().map(|&b| b as u8).collect())
            .unwrap_or_default();
        self.solid_u8 = solid;
    }

    pub fn step(&mut self, n: u32) {
        if let Some(sim) = self.sim.as_mut() {
            match &self.multiphase {
                Some(mp) => {
                    for _ in 0..n {
                        mp.update_force(sim);
                        sim.step();
                    }
                }
                None => sim.run(n as usize),
            }
        }
    }

    pub fn nx(&self) -> u32 {
        self.sim.as_ref().map_or(0, |s| s.nx() as u32)
    }
    pub fn ny(&self) -> u32 {
        self.sim.as_ref().map_or(0, |s| s.ny() as u32)
    }
    pub fn time(&self) -> f64 {
        self.sim.as_ref().map_or(0.0, |s| s.time() as f64)
    }

    pub fn rho_ptr(&self) -> *const f32 {
        self.sim
            .as_ref()
            .map_or(std::ptr::null(), |s| s.rho_field().as_ptr())
    }
    pub fn ux_ptr(&self) -> *const f32 {
        self.sim
            .as_ref()
            .map_or(std::ptr::null(), |s| s.ux_field().as_ptr())
    }
    pub fn uy_ptr(&self) -> *const f32 {
        self.sim
            .as_ref()
            .map_or(std::ptr::null(), |s| s.uy_field().as_ptr())
    }
    pub fn solid_ptr(&self) -> *const u8 {
        self.solid_u8.as_ptr()
    }

    /// Paint or erase an obstacle cell. Erasing rebuilds the simulation from
    /// the stored config (flow restarts) because removing walls from a live
    /// flow is not yet supported by the core.
    pub fn set_solid(&mut self, x: u32, y: u32, solid: bool) -> Result<(), BridgeError> {
        let (x, y) = (x as usize, y as usize);
        let Some(sim0) = self.sim.as_ref() else {
            return Ok(());
        };
        let nx = sim0.nx();
        let ny = sim0.ny();
        if x >= nx || y >= ny {
            return Ok(());
        }
        // Refuse to paint over open edges (the core would panic).
        let i = y * nx + x;
        if solid {
            if self.painted[i] == 1 || sim0.is_solid(x, y) {
                return Ok(());
            }
            // Only paint strictly interior cells (the outer ring is walls or
            // open faces), and never the cell directly inward from an open
            // face — set_solid panics there, because such a solid would
            // silently freeze the open BC's unknown populations (A-3).
            if x == 0 || y == 0 || x == nx - 1 || y == ny - 1 || !sim0.set_solid_allowed(x, y) {
                return Ok(());
            }
            self.painted[i] = 1;
            self.sim.as_mut().unwrap().set_solid(x, y);
            self.solid_u8[i] = 1;
        } else if self.painted[i] == 1 {
            self.painted[i] = 0;
            let painted = self.painted.clone();
            let (mut sim, mp) = build_sim(&self.cfg_json)
                .or_else(|_| build_sim_from_scenario(&self.cfg_json))
                .map_err(|_| bridge_error("rebuild failed".to_string()))?;
            let nx = sim.nx();
            sim.set_solid_region(|px, py| painted[py * nx + px] == 1);
            self.sim = Some(sim);
            self.multiphase = mp;
            self.refresh_solid_mirror();
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const NX: usize = 32;
    const NY: usize = 32;
    const STEPS: u32 = 100;
    const U0: f64 = 1.28 / NX as f64;

    fn tgv_json() -> String {
        serde_json::json!({
            "nx": NX,
            "ny": NY,
            "nu": 0.02,
            "collision": "bgk",
            "edges": {
                "left": { "type": "periodic" },
                "right": { "type": "periodic" },
                "bottom": { "type": "periodic" },
                "top": { "type": "periodic" }
            },
            "force": [0.0, 0.0],
            "init": { "kind": "taylorGreen", "amplitude": U0 }
        })
        .to_string()
    }

    #[cfg(target_arch = "wasm32")]
    fn native_after_100_steps() -> Simulation<f32> {
        let (mut sim, mp) = build_sim(&tgv_json()).unwrap();
        assert!(mp.is_none());
        sim.run(STEPS as usize);
        sim
    }

    fn sum_f32(values: &[f32]) -> f64 {
        values.iter().map(|&v| v as f64).sum()
    }

    fn tiny_scenario_json(precision: &str) -> String {
        serde_json::json!({
            "version": 0,
            "name": "wasm-native-smoke",
            "grid": { "nx": 8, "ny": 8 },
            "physics": {
                "nu": 0.08,
                "collision": { "type": "bgk" },
                "force": [1.0e-6, 0.0],
                "precision": precision
            },
            "edges": {
                "left": { "type": "periodic" },
                "right": { "type": "periodic" },
                "bottom": { "type": "bounceBack" },
                "top": { "type": "bounceBack" }
            },
            "obstacles": [
                { "shape": "rect", "x0": 3, "y0": 3, "x1": 3, "y1": 3 }
            ],
            "init": { "kind": "rest" },
            "run": { "steps": 4 },
            "probes": [],
            "outputs": [
                { "field": "rho", "format": "csv", "every": 0 }
            ]
        })
        .to_string()
    }

    #[cfg(not(target_arch = "wasm32"))]
    mod native_bridge {
        use super::*;

        fn err_text<T>(result: Result<T, BridgeError>) -> String {
            format!("{:?}", result.err().expect("expected bridge error"))
        }

        fn f32_slice(ptr: *const f32, len: usize) -> Vec<f32> {
            assert!(!ptr.is_null());
            unsafe { std::slice::from_raw_parts(ptr, len).to_vec() }
        }

        fn assert_round_trip_semantics(input: serde_json::Value) -> serde_json::Value {
            let input_json = input.to_string();
            let output_json = normalize_scenario_json(&input_json).unwrap();
            let output: serde_json::Value = serde_json::from_str(&output_json).unwrap();
            assert_eq!(output, input);
            output
        }

        #[test]
        fn native_scenario_entrypoint_accepts_steps_and_exposes_fields() {
            let mut wasm = WasmSim::new();
            wasm.init_scenario(&tiny_scenario_json("f32")).unwrap();

            assert_eq!(wasm.nx(), 8);
            assert_eq!(wasm.ny(), 8);
            assert_eq!(wasm.time(), 0.0);

            wasm.step(4);
            assert_eq!(wasm.time(), 4.0);

            let len = (wasm.nx() * wasm.ny()) as usize;
            let rho = f32_slice(wasm.rho_ptr(), len);
            let ux = f32_slice(wasm.ux_ptr(), len);
            let uy = f32_slice(wasm.uy_ptr(), len);
            let mass = sum_f32(&rho);

            assert!(rho.iter().all(|v| v.is_finite()));
            assert!(ux.iter().chain(uy.iter()).all(|v| v.is_finite()));
            assert!(mass > 0.0);
            assert!(
                ux.iter().any(|v| *v > 0.0),
                "body force should produce positive x velocity after stepping"
            );
        }

        #[test]
        fn native_scenario_entrypoint_rejects_invalid_config_with_error() {
            let invalid = serde_json::json!({
                "version": 0,
                "name": "bad-nu",
                "grid": { "nx": 8, "ny": 8 },
                "physics": {
                    "nu": 0.0,
                    "collision": { "type": "bgk" },
                    "force": [0.0, 0.0],
                    "precision": "f32"
                },
                "edges": {
                    "left": { "type": "periodic" },
                    "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" }
                },
                "init": { "kind": "rest" },
                "run": { "steps": 1 }
            })
            .to_string();

            let text = err_text(WasmSim::new().init_scenario(&invalid));
            assert!(
                text.contains("kinematic viscosity must be > 0"),
                "unexpected error: {text}"
            );
        }

        #[test]
        fn native_scenario_entrypoint_rejects_f64_precision() {
            let text = err_text(WasmSim::new().init_scenario(&tiny_scenario_json("f64")));
            assert!(
                text.contains("physics.precision \\\"f32\\\" only")
                    || text.contains("physics.precision \"f32\" only"),
                "unexpected error: {text}"
            );
        }

        #[test]
        fn schema_round_trip_preserves_2d_and_3d_semantics() {
            let two_d = serde_json::json!({
                "version": 0,
                "name": "round-trip-2d",
                "grid": { "nx": 24, "ny": 12 },
                "physics": {
                    "nu": 0.04,
                    "collision": { "type": "trt", "magic": 0.1875 },
                    "force": [0.000001, 0.0],
                    "gravity": [0.0, -0.000002, 0.0],
                    "precision": "f32"
                },
                "compute": { "backend": "cpu" },
                "edges": {
                    "left": { "type": "velocityInlet", "u": [0.03, 0.0] },
                    "right": { "type": "pressureOutlet", "rho": 1.0 },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" }
                },
                "inletProfile": { "edge": "left", "kind": "parabolic", "umax": 0.05 },
                "obstacles": [
                    { "shape": "circle", "cx": 7.5, "cy": 6.0, "r": 2.0 },
                    { "shape": "rect", "x0": 10, "y0": 4, "x1": 12, "y1": 5 }
                ],
                "init": {
                    "kind": "droplet",
                    "cx": 8.0,
                    "cy": 6.0,
                    "r": 2.5,
                    "rhoLiquid": 2.0,
                    "rhoVapor": 0.2
                },
                "multiphase": { "g": -5.0, "gWall": 0.5 },
                "run": {
                    "steps": 25,
                    "stopWhenSteady": { "epsilon": 1.0e-8, "checkEvery": 5 }
                },
                "probes": [
                    { "type": "point", "x": 5, "y": 6, "every": 2 }
                ],
                "outputs": [
                    { "field": "speed", "format": "png", "every": 10 },
                    { "field": "rho", "format": "csv", "every": 0 }
                ]
            });
            let two_d_out = assert_round_trip_semantics(two_d);
            assert_eq!(two_d_out["inletProfile"]["kind"], "parabolic");
            assert_eq!(two_d_out["multiphase"]["gWall"], 0.5);

            let three_d = serde_json::json!({
                "version": 0,
                "name": "round-trip-3d",
                "grid": { "nx": 8, "ny": 6, "nz": 4 },
                "physics": {
                    "nu": 0.03,
                    "collision": { "type": "bgk" },
                    "force": [0.0, 0.0],
                    "gravity": [0.0, 0.0, -0.000001],
                    "precision": "f64"
                },
                "compute": { "backend": "cpu" },
                "wall": "bouzidi",
                "edges": {
                    "left": { "type": "periodic" },
                    "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" },
                    "front": { "type": "periodic" },
                    "back": { "type": "periodic" }
                },
                "obstacles": [
                    { "shape": "sphere", "cx": 4.0, "cy": 3.0, "cz": 2.0, "r": 1.5 }
                ],
                "init": { "kind": "rest" },
                "run": { "steps": 5 },
                "probes": [
                    { "type": "point", "x": 3, "y": 2, "z": 1, "every": 1 }
                ],
                "outputs": [
                    { "field": "vorticityMag", "format": "vtk", "every": 5 }
                ]
            });
            let three_d_out = assert_round_trip_semantics(three_d);
            assert_eq!(three_d_out["grid"]["nz"], 4);
            assert_eq!(three_d_out["edges"]["front"]["type"], "periodic");
            assert_eq!(three_d_out["obstacles"][0]["shape"], "sphere");
        }
    }

    #[test]
    fn native_tgv_smoke_golden() {
        let mut initial = build_sim(&tgv_json()).unwrap().0;
        let mass0 = sum_f32(initial.rho_field());
        initial.run(STEPS as usize);
        let mass1 = sum_f32(initial.rho_field());
        let rel = ((mass1 - mass0) / mass0).abs();

        assert!(rel <= 1.0e-6, "relative mass drift {rel:e}");
        assert_eq!(initial.time(), STEPS as u64);
        assert_eq!(initial.rho(7, 11).to_bits(), 0x3f80_25b6);
        assert_eq!(initial.ux(7, 11).to_bits(), 0xbbb6_2bd2);
        assert_eq!(initial.uy(7, 11).to_bits(), 0xbc98_d05a);
    }

    #[cfg(target_arch = "wasm32")]
    mod wasm {
        use super::*;
        use wasm_bindgen::JsCast;
        use wasm_bindgen_test::*;

        fn f32_view(ptr: *const f32, len: usize) -> Vec<f32> {
            let memory = wasm_bindgen::memory()
                .dyn_into::<js_sys::WebAssembly::Memory>()
                .unwrap();
            let all = js_sys::Float32Array::new(&memory.buffer());
            all.subarray(ptr as u32 / 4, ptr as u32 / 4 + len as u32)
                .to_vec()
        }

        #[wasm_bindgen_test]
        fn wasm_tgv_smoke_matches_compat_f32() {
            let cfg = tgv_json();
            let (native_initial, _) = build_sim(&cfg).unwrap();
            let mass0 = sum_f32(native_initial.rho_field());
            let native = native_after_100_steps();

            let mut wasm = WasmSim::new();
            wasm.init(&cfg).unwrap();
            wasm.step(STEPS);

            assert_eq!(wasm.nx(), NX as u32);
            assert_eq!(wasm.ny(), NY as u32);
            assert_eq!(wasm.time(), STEPS as f64);

            let len = NX * NY;
            let rho = f32_view(wasm.rho_ptr(), len);
            let ux = f32_view(wasm.ux_ptr(), len);
            let uy = f32_view(wasm.uy_ptr(), len);
            let mass1 = sum_f32(&rho);
            let rel = ((mass1 - mass0) / mass0).abs();

            assert!(rel <= 1.0e-6, "relative mass drift {rel:e}");
            assert!(ux.iter().chain(uy.iter()).all(|v| v.is_finite()));

            for i in 0..len {
                assert_eq!(
                    rho[i].to_bits(),
                    native.rho_field()[i].to_bits(),
                    "rho[{i}]"
                );
                assert_eq!(ux[i].to_bits(), native.ux_field()[i].to_bits(), "ux[{i}]");
                assert_eq!(uy[i].to_bits(), native.uy_field()[i].to_bits(), "uy[{i}]");
            }
        }
    }
}

impl Default for WasmSim {
    fn default() -> Self {
        Self::new()
    }
}
