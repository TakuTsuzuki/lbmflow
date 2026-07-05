/**
 * LBMFlow エンジン抽象。
 *
 * 将来 Rust 製 WASM エンジン（wasm-bindgen）をこの interface に適合させて
 * 差し替える。UI 側はこのファイルの型のみに依存すること。
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
  nu: number; // 動粘性係数（格子単位）
  collision: "bgk" | "trt";
  edges: { left: EdgeBC; right: EdgeBC; bottom: EdgeBC; top: EdgeBC };
  force: [number, number];
  /** Shan-Chen 単成分多相（省略可）。g は負で凝集（推奨 -5.0） */
  multiphase?: { g: number; gWall?: number };
  /** 初期密度場（省略時は静止一様）。多相とセットで使う */
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
  step(n: number): void; // nステップ進める
  readonly nx: number;
  readonly ny: number;
  readonly time: number; // 経過ステップ数
  rho(): Float32Array; // 長さ nx*ny、index = y*nx+x（y=0が下端。描画時は上下反転に注意）
  ux(): Float32Array;
  uy(): Float32Array;
  solidMask(): Uint8Array; // 1 = 障害物
  setSolid(x: number, y: number, solid: boolean): void; // 障害物ペイント
}
