/**
 * WASM bridge adapter for the Rust-based LBM core (lbm-core).
 *
 * Build: from the repository root, run
 *   `wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg`
 * (the generated pkg/ is committed, so a rebuild is not normally needed)
 *
 * Fields are zero-copy: returns views over wasm memory.
 * Views remain valid until the next step()/init() call (re-fetch every
 * frame when rendering).
 */
import init, { WasmSim } from "./pkg/lbm_wasm.js";
import type { Engine, EngineConfig } from "./types.ts";

export class WasmEngine implements Engine {
  private constructor(
    private sim: WasmSim,
    private mem: WebAssembly.Memory,
  ) {}

  static async create(): Promise<WasmEngine> {
    const out = await init();
    return new WasmEngine(new WasmSim(), out.memory);
  }

  init(cfg: EngineConfig): void {
    this.sim.init(JSON.stringify(cfg));
  }

  step(n: number): void {
    this.sim.step(n);
  }

  get nx(): number {
    return this.sim.nx();
  }
  get ny(): number {
    return this.sim.ny();
  }
  get time(): number {
    return this.sim.time();
  }

  rho(): Float32Array {
    return new Float32Array(this.mem.buffer, this.sim.rho_ptr(), this.nx * this.ny);
  }
  ux(): Float32Array {
    return new Float32Array(this.mem.buffer, this.sim.ux_ptr(), this.nx * this.ny);
  }
  uy(): Float32Array {
    return new Float32Array(this.mem.buffer, this.sim.uy_ptr(), this.nx * this.ny);
  }
  solidMask(): Uint8Array {
    return new Uint8Array(this.mem.buffer, this.sim.solid_ptr(), this.nx * this.ny);
  }

  setSolid(x: number, y: number, solid: boolean): void {
    this.sim.set_solid(x, y, solid);
  }
}
