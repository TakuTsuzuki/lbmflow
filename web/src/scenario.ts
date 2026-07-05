/**
 * Agent モードのシナリオ JSON (v0) 書き出し。
 *
 * 書式は `crates/lbm-scenario/src/lib.rs` の `Scenario`（serde camelCase・
 * deny_unknown_fields）に一致させること。ここで出力した JSON はそのまま
 * `lbm run scenario.json` で実行できる。
 */

import type { EdgeBC, Engine, EngineConfig } from "./engine/types.ts";

/** Rust: `Obstacle::Rect`（両端を含む閉区間、格子座標。y=0 が下端） */
export interface RectObstacle {
  shape: "rect";
  x0: number;
  y0: number;
  x1: number;
  y1: number;
}

export interface ScenarioInitRest {
  kind: "rest";
}

export interface ScenarioInitDroplet {
  kind: "droplet";
  cx: number;
  cy: number;
  r: number;
  rhoLiquid: number;
  rhoVapor: number;
}

export type ScenarioField = "speed" | "ux" | "uy" | "rho" | "vorticity";

export interface ScenarioOutput {
  field: ScenarioField;
  format: "png" | "csv";
  /** N ステップごとにスナップショット（0 = 終了時のみ） */
  every: number;
}

/** Rust: `Scenario`（camelCase）。省略可のフィールドは付けない */
export interface Scenario {
  version: number;
  name: string;
  grid: { nx: number; ny: number };
  physics: {
    nu: number;
    collision: { type: "bgk" | "trt" };
    force: [number, number];
    precision: "f32" | "f64";
  };
  edges: { left: EdgeBC; right: EdgeBC; bottom: EdgeBC; top: EdgeBC };
  obstacles?: RectObstacle[];
  init: ScenarioInitRest | ScenarioInitDroplet;
  multiphase?: { g: number; gWall: number };
  run: { steps: number };
  outputs: ScenarioOutput[];
}

/**
 * solid マスクを重なりのない矩形の集合に分解する（貪欲法）。
 * 見つけたセルからまず右（+x）へ、次に上（+y）へ伸ばして矩形を取り出す。
 * 矩形は Rust 側と同じく閉区間 [x0..x1] × [y0..y1] で、合併はマスクと厳密に
 * 一致する（近似ではなく無損失の分解）。
 */
export function maskToRects(mask: Uint8Array, nx: number, ny: number): RectObstacle[] {
  const covered = new Uint8Array(nx * ny);
  const rects: RectObstacle[] = [];

  for (let y = 0; y < ny; y++) {
    const row = y * nx;
    for (let x = 0; x < nx; x++) {
      if (mask[row + x] !== 1 || covered[row + x] === 1) continue;

      // 右へ伸ばす
      let x1 = x;
      while (x1 + 1 < nx && mask[row + x1 + 1] === 1 && covered[row + x1 + 1] === 0) {
        x1++;
      }

      // 上へ伸ばす（行 [x..x1] 全体が solid かつ未カバーの間）
      let y1 = y;
      grow: while (y1 + 1 < ny) {
        const next = (y1 + 1) * nx;
        for (let xx = x; xx <= x1; xx++) {
          if (mask[next + xx] !== 1 || covered[next + xx] === 1) break grow;
        }
        y1++;
      }

      for (let yy = y; yy <= y1; yy++) {
        const r = yy * nx;
        for (let xx = x; xx <= x1; xx++) covered[r + xx] = 1;
      }
      rects.push({ shape: "rect", x0: x, y0: y, x1, y1 });
    }
  }
  return rects;
}

/**
 * 現在のエンジン状態 + UI 設定からシナリオ JSON を組み立てる。
 * 障害物は現在の solid マスクから矩形分解で取り出す。
 *
 * 注意: 壁型境界（bounceBack / movingWall）は lbm-core が外周 1 セルの
 * solid リムとして実現するため、エンジンのマスクには外周セルが含まれる。
 * このリムは `lbm run` 側が edges 指定から再生成する（さらに movingWall を
 * 障害物で上書きすると駆動壁が死ぬ）ので、書き出すのはユーザーが描ける
 * 内部セルのみとする（ペイントは外周 1 セルには置けない仕様）。
 */
export function buildScenario(engine: Engine, cfg: EngineConfig, name: string): Scenario {
  const { nx, ny } = engine;
  const mask = engine.solidMask();
  const interior = new Uint8Array(nx * ny);
  for (let y = 1; y < ny - 1; y++) {
    const row = y * nx;
    for (let x = 1; x < nx - 1; x++) {
      interior[row + x] = mask[row + x]!;
    }
  }
  const obstacles = maskToRects(interior, nx, ny);

  return {
    version: 0,
    name,
    grid: { nx: engine.nx, ny: engine.ny },
    physics: {
      nu: cfg.nu,
      collision: { type: cfg.collision },
      force: [cfg.force[0], cfg.force[1]],
      precision: "f64",
    },
    edges: {
      left: cfg.edges.left,
      right: cfg.edges.right,
      bottom: cfg.edges.bottom,
      top: cfg.edges.top,
    },
    ...(obstacles.length > 0 ? { obstacles } : {}),
    init: cfg.init
      ? {
          kind: "droplet",
          cx: cfg.init.cx,
          cy: cfg.init.cy,
          r: cfg.init.r,
          rhoLiquid: cfg.init.rhoLiquid,
          rhoVapor: cfg.init.rhoVapor,
        }
      : { kind: "rest" },
    ...(cfg.multiphase
      ? { multiphase: { g: cfg.multiphase.g, gWall: cfg.multiphase.gWall ?? 0 } }
      : {}),
    run: { steps: 20000 },
    outputs: [{ field: cfg.multiphase ? "rho" : "speed", format: "png", every: 0 }],
  };
}
