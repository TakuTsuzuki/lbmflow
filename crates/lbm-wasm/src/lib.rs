//! wasm-bindgen bridge: exposes `lbm_core::compat::Simulation<f32>` behind
//! the TypeScript `Engine` interface (see `web/src/engine/types.ts`).
//!
//! Field access is zero-copy: the `*_ptr` methods return pointers into wasm
//! linear memory; the JS adapter wraps them in `Float32Array` views that are
//! valid until the next `step`/`init` call.

use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use serde::Deserialize;
use wasm_bindgen::prelude::*;

#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
enum JsEdge {
    Periodic,
    BounceBack,
    MovingWall { u: [f32; 2] },
    VelocityInlet { u: [f32; 2] },
    PressureOutlet { rho: f32 },
    Outflow,
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

#[derive(Deserialize)]
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

/// Browser-facing simulation handle.
#[wasm_bindgen]
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

fn build_sim(cfg_json: &str) -> Result<(Simulation<f32>, Option<ShanChen<f32>>), JsError> {
    let cfg: JsConfig = serde_json::from_str(cfg_json)
        .map_err(|e| JsError::new(&format!("設定JSONを解釈できません: {e}")))?;
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
    .map_err(|e| JsError::new(&format!("設定エラー: {e}")))?;
    if let Some(JsInit::Droplet {
        cx,
        cy,
        r,
        rho_liquid,
        rho_vapor,
    }) = cfg.init
    {
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
    let mp = cfg
        .multiphase
        .map(|m| ShanChen::<f32>::new(m.g).with_wall(m.g_wall));
    Ok((sim, mp))
}

#[wasm_bindgen]
impl WasmSim {
    #[wasm_bindgen(constructor)]
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
    pub fn init(&mut self, cfg_json: &str) -> Result<(), JsError> {
        let (sim, mp) = build_sim(cfg_json)?;
        let n = sim.nx() * sim.ny();
        self.cfg_json = cfg_json.to_string();
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
    pub fn set_solid(&mut self, x: u32, y: u32, solid: bool) -> Result<(), JsError> {
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
            // set_solid panics on open edges; only paint strictly interior
            // cells for simplicity.
            if x == 0 || y == 0 || x == nx - 1 || y == ny - 1 {
                return Ok(());
            }
            self.painted[i] = 1;
            self.sim.as_mut().unwrap().set_solid(x, y);
            self.solid_u8[i] = 1;
        } else if self.painted[i] == 1 {
            self.painted[i] = 0;
            let painted = self.painted.clone();
            let (mut sim, mp) = build_sim(&self.cfg_json)
                .map_err(|_| JsError::new("再構築に失敗しました"))?;
            let nx = sim.nx();
            sim.set_solid_region(|px, py| painted[py * nx + px] == 1);
            self.sim = Some(sim);
            self.multiphase = mp;
            self.refresh_solid_mirror();
        }
        Ok(())
    }
}

impl Default for WasmSim {
    fn default() -> Self {
        Self::new()
    }
}
