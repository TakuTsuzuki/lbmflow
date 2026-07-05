# LBMFlow Web GUI

Browser GUI for LBMFlow, a Lattice Boltzmann Method (LBM) fluid simulator.
Implemented in Vite + TypeScript (vanilla, no framework, zero runtime dependencies).

It currently runs on a **mock engine** (a pure-TS analytic flow-field generator).
The architecture is designed so a Rust-based WASM engine can be plugged into the
same interface in the future.

## Getting started

```bash
cd web
npm install
npm run dev        # http://localhost:5173
```

Production build:

```bash
npm run build      # tsc(strict) → vite build; output goes to web/dist/
npm run preview    # check the dist/ build
```

## Usage

1. Pick a preset from the header (lid-driven cavity / flow around a cylinder /
   channel flow / free canvas)
2. Press ▶ Run (Space also works)
3. Drag on the canvas to draw obstacles (right-drag, or use "Erase" mode, to erase)
4. Adjust the visualized quantity (speed / vorticity / density) or parameters
   in the right panel

The simulation stops automatically when the tab is hidden.

## Directory layout

```
web/
├── index.html            # Static UI skeleton (English labels)
├── src/
│   ├── main.ts           # App wiring, RAF loop, obstacle painting
│   ├── style.css         # Dark theme (CSS variables, hand-written)
│   ├── presets.ts        # Preset definitions (EngineConfig + description + initial obstacles)
│   ├── colormap.ts       # viridis / RdBu LUTs (no external dependencies)
│   ├── render.ts         # Scalarize (|u| / vorticity / density) → LUT color → canvas transfer
│   └── engine/
│       ├── types.ts      # ★ Engine abstraction (wasm-bindgen contract)
│       ├── index.ts      # ★ Engine-creation swap point
│       └── mock.ts       # Mock engine (analytic flow-field generator)
└── vite.config.ts
```

## Engine swap design

The UI depends **only** on the `Engine` interface in `src/engine/types.ts`.

```ts
export interface Engine {
  init(cfg: EngineConfig): void;
  step(n: number): void;
  readonly nx: number;
  readonly ny: number;
  readonly time: number;
  rho(): Float32Array;   // length nx*ny, index = y*nx+x (y=0 is the bottom edge)
  ux(): Float32Array;
  uy(): Float32Array;
  solidMask(): Uint8Array;
  setSolid(x: number, y: number, solid: boolean): void;
}
```

Steps to migrate to the WASM engine:

1. On the wasm-bindgen side, expose a class matching the signature above
   (e.g. `WasmEngine`)
   - `rho()` and friends can return a `Float32Array` view into the WASM
     memory buffer, or a copy. The caller only ever holds onto it under the
     assumption that it's "valid until the next `step()`/`init()` call," so
     returning a view is fine
2. Rewrite `createEngine()` in `src/engine/index.ts` to return a `WasmEngine`
   - If asynchronous loading of `.wasm` is needed, change it to
     `createEngine(): Promise<Engine>` and `await` it in the startup
     sequence at the top of `main.ts` (the change stays confined to these
     2 files)
3. `mock.ts` can be kept around for demos and as a fallback

### Coordinate system convention

- `index = y * nx + x`, `y = 0` is the **bottom edge** (physics convention)
- At draw time, `render.ts` flips vertically when transferring to the canvas
  (where the top is y max)

## How the mock engine works (`src/engine/mock.ts`)

It doesn't solve real LBM; instead it analytically synthesizes the field as a
function of the elapsed step count `t`:

- The base flow is chosen from the boundary conditions (top-wall movingWall
  → cavity-like primary vortex; velocityInlet → uniform flow + alternating
  vortex shedding downstream of an obstacle, Kármán-vortex-street-like;
  body force → Poiseuille parabolic profile; all-periodic → decaying
  Taylor-Green vortex)
- Viscosity ν affects the vortex decay rate / vortex-street amplitude
  (the larger it is, the faster things settle down)
- `collision: "bgk"` adds a small amount of noise (a dramatization implying
  TRT is more stable — not real physics)
- Obstacle cells get u=0, ρ=1; neighboring cells are slowed down to look
  wall-like

## Known limitations

- The mock engine's flow field is a synthesized field optimized for
  appearance and is not physically correct
  (boundary conditions, ν, and the collision operator are only used for the
  "feel" of the dramatization)
- The `pressureOutlet` ρ setting is currently unused by the mock
- When the resolution changes, painted obstacles are carried over via
  nearest-neighbor sampling, so their outline becomes coarser
