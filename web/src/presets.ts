import type { Engine, EngineConfig } from "./engine/types.ts";

/** プリセット = エンジン設定 + 説明文 + 障害物の初期配置 */
export interface Preset {
  id: string;
  /** セレクトボックスに出す名前 */
  name: string;
  /** 何が起きるか・注目ポイント（1-2文、初学者向け） */
  description: string;
  /** 基準解像度での設定。解像度セレクトで nx, ny はスケールされる */
  config: EngineConfig;
  /** 初期障害物の配置（省略可）。engine.nx / engine.ny を使って描くこと */
  paintObstacles?: (engine: Engine) => void;
}

export const PRESETS: Preset[] = [
  {
    id: "cavity",
    name: "キャビティ流れ",
    description:
      "上のフタだけが右に動く箱の中の流れです。フタに引きずられた流体が箱全体を回る大きな渦を作ります。渦の中心が少しずつ動く様子に注目してください。",
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
    name: "円柱まわりの流れ",
    description:
      "左から流れてきた流体が円柱にぶつかり、後ろに渦が交互に放出されます（カルマン渦列）。粘性 ν を小さくするほど渦がはっきり現れます。渦度表示がおすすめです。",
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
    name: "チャネル流（ポアズイユ）",
    description:
      "上下の壁にはさまれた水路を、一定の力（ポンプの代わり）で押し流します。壁でゼロ・中央で最大の、放物線型の速度分布ができるのが見どころです。",
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
    id: "sandbox",
    name: "自由キャンバス",
    description:
      "全方向が周期境界（右端と左端がつながっている）の何もない空間に、渦模様が漂います。ブラシで好きな形の障害物を描いて、流れがどう変わるか試してみましょう。",
    config: {
      nx: 160,
      ny: 120,
      nu: 0.008,
      collision: "trt",
      edges: {
        left: { type: "periodic" },
        right: { type: "periodic" },
        bottom: { type: "periodic" },
        top: { type: "periodic" },
      },
      force: [0, 0],
    },
  },
];

/** 円形の障害物を塗る補助関数 */
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
