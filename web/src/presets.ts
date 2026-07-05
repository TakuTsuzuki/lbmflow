import type { Engine, EngineConfig } from "./engine/types.ts";
import type { VisMode } from "./render.ts";

/** A preset = engine config + description + initial obstacle layout */
export interface Preset {
  id: string;
  /** Name shown in the select box */
  name: string;
  /** What happens / what to watch for (1-2 sentences, beginner-friendly) */
  description: string;
  /** Config at the reference resolution. The resolution select scales nx, ny */
  config: EngineConfig;
  /** Initial obstacle layout (optional). Draw using engine.nx / engine.ny */
  paintObstacles?: (engine: Engine) => void;
  /** Visualization mode to switch to when the preset is applied (optional) */
  defaultVis?: VisMode;
}

export const PRESETS: Preset[] = [
  {
    id: "cavity",
    name: "Lid-driven cavity",
    description:
      "Flow inside a box where only the top lid moves to the right. The fluid dragged along by the lid forms one large vortex circulating around the whole box. Watch how the vortex center drifts gradually.",
    config: {
      nx: 128,
      ny: 128,
      nu: 0.02,
      collision: "trt",
      edges: {
        left: { type: "bounceBack" },
        right: { type: "bounceBack" },
        bottom: { type: "bounceBack" },
        top: { type: "movingWall", u: [0.1, 0] },
      },
      force: [0, 0],
    },
  },
  {
    id: "cylinder",
    name: "Flow around a cylinder",
    description:
      "Fluid flowing in from the left hits a cylinder, shedding vortices alternately behind it (a Kármán vortex street). The smaller the viscosity ν, the more distinct the vortices become. Vorticity view is recommended.",
    config: {
      nx: 224,
      ny: 112,
      nu: 0.005,
      collision: "trt",
      edges: {
        left: { type: "velocityInlet", u: [0.08, 0] },
        right: { type: "outflow" },
        bottom: { type: "bounceBack" },
        top: { type: "bounceBack" },
      },
      force: [0, 0],
    },
    paintObstacles: (engine) => {
      paintCircle(engine, engine.nx * 0.24, engine.ny * 0.5, engine.ny * 0.11);
    },
  },
  {
    id: "poiseuille",
    name: "Channel flow (Poiseuille)",
    description:
      "A channel between two walls is driven by a constant force (standing in for a pump). The highlight is the parabolic velocity profile that forms: zero at the walls, maximum at the center.",
    config: {
      nx: 192,
      ny: 96,
      nu: 0.02,
      collision: "trt",
      edges: {
        left: { type: "periodic" },
        right: { type: "periodic" },
        bottom: { type: "bounceBack" },
        top: { type: "bounceBack" },
      },
      force: [2e-6, 0],
    },
  },
  {
    id: "droplet",
    name: "Two-phase fluid (droplet)",
    description:
      "A two-phase fluid simulation where liquid and gas (vapor) coexist (Shan-Chen model). Watch, in the density view, how the droplet is held round by surface tension as the interface stabilizes into a sharp boundary.",
    config: {
      nx: 128,
      ny: 128,
      nu: 1 / 6,
      collision: "trt",
      edges: {
        left: { type: "periodic" },
        right: { type: "periodic" },
        bottom: { type: "periodic" },
        top: { type: "periodic" },
      },
      force: [0, 0],
      multiphase: { g: -5.0 },
      init: {
        kind: "droplet",
        cx: 64,
        cy: 64,
        r: 26,
        rhoLiquid: 2.0,
        rhoVapor: 0.15,
      },
    },
    defaultVis: "density",
  },
  {
    id: "sandbox",
    name: "Free canvas",
    description:
      "A constant force drives flow through a channel between two walls. Use the brush to draw obstacles of any shape and see how the vortices and flow behind them change.",
    config: {
      nx: 192,
      ny: 112,
      nu: 0.006,
      collision: "trt",
      edges: {
        left: { type: "periodic" },
        right: { type: "periodic" },
        bottom: { type: "bounceBack" },
        top: { type: "bounceBack" },
      },
      force: [4e-6, 0],
    },
  },
];

/** Helper function that paints a circular obstacle */
export function paintCircle(
  engine: Engine,
  cx: number,
  cy: number,
  r: number,
): void {
  const x0 = Math.max(0, Math.floor(cx - r - 1));
  const x1 = Math.min(engine.nx - 1, Math.ceil(cx + r + 1));
  const y0 = Math.max(0, Math.floor(cy - r - 1));
  const y1 = Math.min(engine.ny - 1, Math.ceil(cy + r + 1));
  for (let y = y0; y <= y1; y++) {
    for (let x = x0; x <= x1; x++) {
      const dx = x - cx;
      const dy = y - cy;
      if (dx * dx + dy * dy <= r * r) engine.setSolid(x, y, true);
    }
  }
}
