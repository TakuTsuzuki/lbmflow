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
            let (mut sim, mp) =
                build_sim(&self.cfg_json).map_err(|_| JsError::new("再構築に失敗しました"))?;
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
