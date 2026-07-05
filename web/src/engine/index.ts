import type { Engine } from "./types.ts";
import { MockEngine } from "./mock.ts";

/**
 * Centralizes engine creation.
 *
 * By default, loads the Rust-based WASM engine (the real LBM); falls back
 * to the mock if that fails or if the `?engine=mock` query param is set.
 */
export async function createEngine(): Promise<Engine> {
  const params = new URLSearchParams(location.search);
  if (params.get("engine") === "mock") {
    console.info("LBMFlow: using the mock engine (?engine=mock)");
    return new MockEngine();
  }
  try {
    const { WasmEngine } = await import("./wasm.ts");
    const engine = await WasmEngine.create();
    console.info("LBMFlow: using the WASM engine (lbm-core)");
    return engine;
  } catch (err) {
    console.warn("LBMFlow: failed to load the WASM engine, falling back to mock:", err);
    return new MockEngine();
  }
}
