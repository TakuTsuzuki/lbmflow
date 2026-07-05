/**
 * LBMFlow engine abstraction.
 *
 * The Rust-based WASM engine (wasm-bindgen) will eventually conform to this
 * interface and be swapped in. The UI side should depend only on the types
 * in this file.
 */

export type EdgeBC =
  | { type: "periodic" }
  | { type: "bounceBack" }
  | { type: "movingWall"; u: [number, number] }
  | { type: "velocityInlet"; u: [number, number] }
  | { type: "pressureOutlet"; rho: number }
  | { type: "outflow" };

export interface EngineConfig {
  nx: number;
  ny: number;
  nu: number; // kinematic viscosity (lattice units)
  collision: "bgk" | "trt";
  edges: { left: EdgeBC; right: EdgeBC; bottom: EdgeBC; top: EdgeBC };
  force: [number, number];
  /** Shan-Chen single-component multiphase (optional). Negative g causes cohesion (recommended -5.0) */
  multiphase?: { g: number; gWall?: number };
  /** Initial density field (uniform at rest if omitted). Used together with multiphase */
  init?: {
    kind: "droplet";
    cx: number;
    cy: number;
    r: number;
    rhoLiquid: number;
    rhoVapor: number;
  };
}

export interface Engine {
  init(cfg: EngineConfig): void;
  step(n: number): void; // advance by n steps
  readonly nx: number;
  readonly ny: number;
  readonly time: number; // elapsed step count
  rho(): Float32Array; // length nx*ny, index = y*nx+x (y=0 is the bottom edge; watch for vertical flip when rendering)
  ux(): Float32Array;
  uy(): Float32Array;
  solidMask(): Uint8Array; // 1 = obstacle
  setSolid(x: number, y: number, solid: boolean): void; // paint obstacles
}
