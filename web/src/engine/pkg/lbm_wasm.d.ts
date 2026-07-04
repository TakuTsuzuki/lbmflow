/* tslint:disable */
/* eslint-disable */

/**
 * Browser-facing simulation handle.
 */
export class WasmSim {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * (Re)initialise from an EngineConfig JSON string.
     */
    init(cfg_json: string): void;
    constructor();
    nx(): number;
    ny(): number;
    rho_ptr(): number;
    /**
     * Paint or erase an obstacle cell. Erasing rebuilds the simulation from
     * the stored config (flow restarts) because removing walls from a live
     * flow is not yet supported by the core.
     */
    set_solid(x: number, y: number, solid: boolean): void;
    solid_ptr(): number;
    step(n: number): void;
    time(): number;
    ux_ptr(): number;
    uy_ptr(): number;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_wasmsim_free: (a: number, b: number) => void;
    readonly wasmsim_init: (a: number, b: number, c: number) => [number, number];
    readonly wasmsim_new: () => number;
    readonly wasmsim_nx: (a: number) => number;
    readonly wasmsim_ny: (a: number) => number;
    readonly wasmsim_rho_ptr: (a: number) => number;
    readonly wasmsim_set_solid: (a: number, b: number, c: number, d: number) => [number, number];
    readonly wasmsim_solid_ptr: (a: number) => number;
    readonly wasmsim_step: (a: number, b: number) => void;
    readonly wasmsim_time: (a: number) => number;
    readonly wasmsim_ux_ptr: (a: number) => number;
    readonly wasmsim_uy_ptr: (a: number) => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
