# WASM Bridge Design (Phase 5)

**Status (2026-07-07)**: **Landed**: `crates/lbm-wasm` exposes `WasmSim` over `compat::Simulation<f32>` and `web/src/engine/wasm.ts` adapts it to `Engine`.
**Landed with drift**: zero-copy field pointers and rebuild-on-erase are implemented; SCMP f32 config is also present.
**Superseded/current intent**: shared `lbm-scenario` config conversion and `set_inlet_profile_parabolic` are not the landed WASM API; performance targets were not remeasured in this sweep.

The layer that connects `lbm-core` to the GUI's (web/, TypeScript) `Engine` interface.

## Crate: crates/lbm-wasm

(landed 2026-07-07 — wasm-bindgen crate depends on `lbm-core` with default features disabled and emits the committed `web/src/engine/pkg` wrapper)

- `wasm-bindgen` + `lbm-core` (`default-features = false`, rayon disabled =
  single-threaded. Targets the scale that can run within a single browser frame)
- Precision is **fixed at f32** (halves memory, sufficient for WASM. Copy-free alignment
  with the JS-side Float32Array)
- Build: `wasm-pack build crates/lbm-wasm --target web --release`
  → outputs to `web/src/engine/pkg/`, wrapped by the `WasmEngine implements Engine` adapter

## Public API (1:1 with the TS Engine interface)

(landed with API drift 2026-07-07 — `WasmSim` exposes init/step/nx/ny/time/field pointers/solid editing; `set_inlet_profile_parabolic` is not implemented, and errors are English `JsError` strings)

```rust
#[wasm_bindgen]
pub struct WasmSim { inner: Option<Simulation<f32>>, cfg: ... }

#[wasm_bindgen]
impl WasmSim {
    #[wasm_bindgen(constructor)]
    pub fn new() -> WasmSim;
    /// cfg_json: the GUI's EngineConfig, JSON.stringify'd as-is.
    /// Errors are returned as JsError (Japanese message).
    pub fn init(&mut self, cfg_json: &str) -> Result<(), JsError>;
    pub fn step(&mut self, n: u32);
    pub fn nx(&self) -> u32;  pub fn ny(&self) -> u32;  pub fn time(&self) -> f64;
    /// Fields are returned as a view into wasm memory without copying (equivalent to
    /// Float32Array::view; document the caveat "valid only until the next step")
    pub fn rho_ptr(&self) -> *const f32;   // + len is nx*ny (view constructed on the JS side)
    pub fn ux_ptr(&self) -> *const f32;
    pub fn uy_ptr(&self) -> *const f32;
    pub fn solid_ptr(&self) -> *const u8;
    pub fn set_solid(&mut self, x: u32, y: u32, solid: bool);
    pub fn set_inlet_profile_parabolic(&mut self, edge: &str, umax: f32);
}
```

### EngineConfig(JSON) → SimConfig conversion

(partially superseded 2026-07-07 — WASM uses local `JsConfig`/`JsEdge`; GUI scenario export targets `lbm-scenario`, but conversion code is not shared)

- `collision: "bgk" | "trt"` → `Collision::Bgk | Trt{magic: 3/16}`
- edges' tagged union → `EdgeBC` (deserialized via serde, defined inside lbm-wasm)
- **This JSON representation is made identical in shape to the `edges`/`physics` sections of
  Agent Mode's scenario JSON (docs/AGENT_MODE_DESIGN.md)**, so the conversion code can be shared
  (serde types live in the shared crate `lbm-scenario` — both lbm-cli and lbm-wasm depend on it)

## Handling setSolid's "erase" operation

(landed via rebuild path 2026-07-07 — the planned live `clear_solid` core API is not implemented)

`Simulation` has no unset_solid (to protect the rim). The GUI's eraser works as follows:
- On the lbm-wasm side, separately maintain a "user drawing layer" (Vec<bool>)
- Erase operation = update the user layer → rebuild equivalent to `init(cfg)` + reapply the
  user layer (tens of ms; batched to a single rebuild during painting)
- Alternatively, in Phase 5 add `clear_solid(x,y)` to lbm-core (panics on open boundary/rim cells)
  and fill the surrounding cells' f with the local feq. **This is the planned approach** (an
  experience where editing doesn't stop the flow is more enjoyable for beginners). Add a
  robustness test to VALIDATION confirming "mass stays finite and no NaN appears after clear_solid."

## Performance target

(current intent 2026-07-07 — target retained as design guidance; no benchmark was run for this status sweep)

- 256×128 (32k cells) f32 single-threaded: target ≥ 30 MLUPS → ~15 step/frame at 60fps.
  GUI default is set to either 192×96 or 256×128.
- Prepare UI copy that directs users toward native (CLI/MCP) for larger grids, 3D, and heavy
  multiphase computation.
