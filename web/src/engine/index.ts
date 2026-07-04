import type { Engine } from "./types.ts";
import { MockEngine } from "./mock.ts";

/**
 * エンジンの生成をここに集約する。
 *
 * WASM エンジンが完成したら、この関数の中身を
 * `return new WasmEngine(...)` に差し替えるだけで UI 全体が切り替わる。
 * （非同期ロードが必要になる場合は `Promise<Engine>` を返す形に変え、
 * main.ts の起動シーケンスで await する。）
 */
export function createEngine(): Engine {
  return new MockEngine();
}
