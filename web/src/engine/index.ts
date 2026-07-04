import type { Engine } from "./types.ts";
import { MockEngine } from "./mock.ts";

/**
 * エンジンの生成をここに集約する。
 *
 * 既定では Rust 製 WASM エンジン（本物の LBM）をロードし、失敗した場合や
 * `?engine=mock` クエリ付きのときはモックにフォールバックする。
 */
export async function createEngine(): Promise<Engine> {
  const params = new URLSearchParams(location.search);
  if (params.get("engine") === "mock") {
    console.info("LBMFlow: モックエンジンを使用します (?engine=mock)");
    return new MockEngine();
  }
  try {
    const { WasmEngine } = await import("./wasm.ts");
    const engine = await WasmEngine.create();
    console.info("LBMFlow: WASM エンジン (lbm-core) を使用します");
    return engine;
  } catch (err) {
    console.warn("LBMFlow: WASM エンジンのロードに失敗、モックへフォールバック:", err);
    return new MockEngine();
  }
}
