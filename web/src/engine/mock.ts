import type { Engine, EngineConfig } from "./types.ts";

/**
 * モックエンジン。
 *
 * 本物の LBM は解かず、時間 t の関数として「それらしい」2D 流れ場を
 * 解析的に合成する。UI 開発・デモ用。
 *
 * - 境界条件の組み合わせから基本流を選ぶ:
 *   - 上壁 movingWall            → キャビティ流れ風（せん断 + 主渦の旋回）
 *   - 左 velocityInlet           → 一様流 + 障害物下流のカルマン渦列風の渦放出
 *   - 外力あり + 上下 bounceBack → ポアズイユ放物線分布 + 進行波
 *   - すべて periodic            → 減衰しながら流れるテイラー・グリーン渦
 * - 粘性 ν は渦の減衰率・渦列の振幅に反映（大きいほど早く静まる）。
 * - collision="bgk" のときは微小なノイズを乗せ、"trt" はクリーン
 *   （「TRT の方が安定・高精度」というデモ上の演出。実物理ではない）。
 * - 障害物セルは u=0, ρ=1。周囲は距離に応じて減速させ、壁らしく見せる。
 */
export class MockEngine implements Engine {
  private cfg: EngineConfig = {
    nx: 128,
    ny: 96,
    nu: 0.01,
    collision: "trt",
    edges: {
      left: { type: "periodic" },
      right: { type: "periodic" },
      bottom: { type: "periodic" },
      top: { type: "periodic" },
    },
    force: [0, 0],
  };

  private _time = 0;
  private _rho = new Float32Array(0);
  private _ux = new Float32Array(0);
  private _uy = new Float32Array(0);
  private _solid = new Uint8Array(0);

  /** 障害物からの距離による減速係数（0=障害物内 … 1=十分遠い） */
  private damp = new Float32Array(0);
  private dampDirty = true;

  /** 障害物の重心と代表半径（渦放出の起点に使う） */
  private solidCx = 0;
  private solidCy = 0;
  private solidR = 0;
  private solidCount = 0;

  get nx(): number {
    return this.cfg.nx;
  }
  get ny(): number {
    return this.cfg.ny;
  }
  get time(): number {
    return this._time;
  }

  init(cfg: EngineConfig): void {
    this.cfg = {
      ...cfg,
      edges: { ...cfg.edges },
      force: [cfg.force[0], cfg.force[1]],
    };
    const n = cfg.nx * cfg.ny;
    this._rho = new Float32Array(n);
    this._ux = new Float32Array(n);
    this._uy = new Float32Array(n);
    this._solid = new Uint8Array(n);
    this.damp = new Float32Array(n);
    this.dampDirty = true;
    this._time = 0;
    this.recompute();
  }

  step(n: number): void {
    this._time += n;
    this.recompute();
  }

  rho(): Float32Array {
    return this._rho;
  }
  ux(): Float32Array {
    return this._ux;
  }
  uy(): Float32Array {
    return this._uy;
  }
  solidMask(): Uint8Array {
    return this._solid;
  }

  setSolid(x: number, y: number, solid: boolean): void {
    const { nx, ny } = this.cfg;
    if (x < 0 || y < 0 || x >= nx || y >= ny) return;
    const i = y * nx + x;
    const v = solid ? 1 : 0;
    if (this._solid[i] === v) return;
    this._solid[i] = v;
    this.dampDirty = true;
    if (solid) {
      // 一時停止中でも「壁になった」ことが見えるように即座に流速を消す
      this._ux[i] = 0;
      this._uy[i] = 0;
      this._rho[i] = 1;
    }
  }

  // ---------------------------------------------------------------- 内部実装

  /** 障害物マスクから減速マップと重心・半径を再計算する */
  private refreshSolidInfo(): void {
    const { nx, ny } = this.cfg;
    const solid = this._solid;
    const damp = this.damp;

    let count = 0;
    let sx = 0;
    let sy = 0;
    for (let i = 0; i < solid.length; i++) {
      if (solid[i] === 1) {
        count++;
        sx += i % nx;
        sy += (i / nx) | 0;
      }
    }
    this.solidCount = count;
    if (count > 0) {
      this.solidCx = sx / count;
      this.solidCy = sy / count;
      this.solidR = Math.max(2, Math.sqrt(count / Math.PI));
    }

    // 近傍に障害物があるセルを減速（チェビシェフ距離2セルの簡易判定）
    const R = 2;
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const i = y * nx + x;
        if (solid[i] === 1) {
          damp[i] = 0;
          continue;
        }
        let dmin = R + 1;
        for (let dy = -R; dy <= R && dmin > 0; dy++) {
          const yy = y + dy;
          if (yy < 0 || yy >= ny) continue;
          const row = yy * nx;
          for (let dx = -R; dx <= R; dx++) {
            const xx = x + dx;
            if (xx < 0 || xx >= nx) continue;
            if (solid[row + xx] === 1) {
              const d = Math.max(Math.abs(dx), Math.abs(dy));
              if (d < dmin) dmin = d;
            }
          }
        }
        damp[i] = dmin > R ? 1 : dmin / (R + 1);
      }
    }
    this.dampDirty = false;
  }

  /** 決定的な擬似乱数（セル座標と時間から）。BGK のノイズ演出用 */
  private static hashNoise(x: number, y: number, t: number): number {
    let h = (x * 374761393 + y * 668265263 + t * 2246822519) | 0;
    h = Math.imul(h ^ (h >>> 13), 1274126177);
    h ^= h >>> 16;
    return (h & 0xffff) / 0xffff - 0.5; // [-0.5, 0.5)
  }

  /** 現在時刻 _time における場をまるごと作り直す */
  private recompute(): void {
    if (this.dampDirty) this.refreshSolidInfo();

    const { nx, ny, nu, edges, force, collision } = this.cfg;
    const t = this._time;
    const rho = this._rho;
    const ux = this._ux;
    const uy = this._uy;
    const solid = this._solid;
    const damp = this.damp;

    // ---- 基本流の種類を判定
    const top = edges.top;
    const left = edges.left;
    const lidU = top.type === "movingWall" ? top.u[0] : 0;
    const inletU = left.type === "velocityInlet" ? left.u[0] : 0;
    const hasForce = force[0] !== 0 || force[1] !== 0;

    // ポアズイユ流の最大速度 Umax = f H^2 / (8 ν) を格子単位で（暴走防止に上限）
    const H = Math.max(2, ny - 2);
    const poisMax = hasForce
      ? Math.min(0.18, Math.abs((force[0] * H * H) / (8 * Math.max(nu, 1e-5))))
      : 0;
    const poisSign = force[0] >= 0 ? 1 : -1;

    // 粘性による全体減衰（テイラー・グリーン渦の厳密解 e^{-2νk²t} を流用）
    const k0 = (2 * Math.PI) / Math.max(nx, ny);
    const decay = Math.exp(-2 * nu * k0 * k0 * 4 * t);

    // ---- カルマン渦列風の放出渦（障害物 + 流入があるときだけ）
    // 各渦は誕生時刻から一様流に乗って下流へ流され、粘性で減衰する。
    type Vortex = { x: number; y: number; s: number; r: number };
    const vortices: Vortex[] = [];
    const advU = inletU !== 0 ? inletU : poisMax * poisSign;
    if (this.solidCount > 12 && Math.abs(advU) > 1e-4) {
      const period = Math.max(40, (this.solidR * 8) / Math.abs(advU) / 4);
      const amp = (Math.abs(advU) * 2.2) / (1 + nu * 260);
      const n0 = Math.floor(t / period);
      for (let k = n0; k > n0 - 9 && k >= 0; k--) {
        const age = t - k * period;
        const sgn = k % 2 === 0 ? 1 : -1;
        const vx = this.solidCx + this.solidR * 1.4 + advU * age;
        if (vx > nx + this.solidR || vx < -this.solidR) continue;
        vortices.push({
          x: vx,
          y: this.solidCy + sgn * this.solidR * 0.75,
          s: sgn * amp * Math.exp(-nu * age * 0.55),
          r: Math.max(3, this.solidR * 0.9),
        });
      }
    }

    // ---- キャビティ主渦（中心がゆっくり歳差運動する）
    let cav: Vortex | null = null;
    if (lidU !== 0) {
      const wob = 0.06 * Math.min(nx, ny);
      cav = {
        x: nx * 0.5 + wob * Math.cos(t * 0.004),
        y: ny * 0.62 + wob * 0.6 * Math.sin(t * 0.0031),
        s: -lidU * 1.6 * (1 - Math.exp(-t / 500)), // 徐々に発達
        r: Math.min(nx, ny) * 0.3,
      };
    }

    // ---- テイラー・グリーン渦（周期境界のとき）
    const isPeriodic =
      edges.left.type === "periodic" &&
      edges.right.type === "periodic" &&
      edges.top.type === "periodic" &&
      edges.bottom.type === "periodic";
    const tgAmp = isPeriodic && !hasForce ? 0.09 : 0;
    const kx = (2 * Math.PI * 2) / nx;
    const ky = (2 * Math.PI * 2) / ny;
    const driftX = 8e-3 * nx * 1e-2; // ゆっくり横に流れる
    const tgPhase = t * driftX * kx * 0.5;

    const noiseAmp = collision === "bgk" ? 0.02 : 0;

    for (let y = 0; y < ny; y++) {
      const fy = ny > 1 ? y / (ny - 1) : 0; // 0(下端) .. 1(上端)
      const row = y * nx;
      for (let x = 0; x < nx; x++) {
        const i = row + x;
        if (solid[i] === 1) {
          ux[i] = 0;
          uy[i] = 0;
          rho[i] = 1;
          continue;
        }

        let vx = 0;
        let vy = 0;

        // 一様流入
        if (inletU !== 0) {
          vx += inletU;
          if (left.type === "velocityInlet") vy += left.u[1];
        }

        // ポアズイユ放物線 + 見た目用の弱い進行波
        if (poisMax > 0) {
          const par = 4 * fy * (1 - fy);
          vx += poisSign * poisMax * par;
          vy +=
            0.06 *
            poisMax *
            par *
            Math.sin((2 * Math.PI * (x - poisSign * poisMax * t * 0.6)) / (nx * 0.5));
        }

        // 上壁駆動のせん断（壁近くほど強い）
        if (lidU !== 0) {
          vx += lidU * Math.pow(fy, 3) * decay0(t);
        }

        // テイラー・グリーン渦（減衰 + ドリフト）
        if (tgAmp > 0) {
          const px = kx * x - tgPhase;
          const py = ky * y;
          vx += tgAmp * decay * Math.cos(px) * Math.sin(py);
          vy += -tgAmp * decay * Math.sin(px) * Math.cos(py);
        }

        // キャビティ主渦
        if (cav) {
          const g = gaussVortex(x, y, cav);
          // 壁で消えるよう sin 包絡を掛ける
          const env =
            Math.sin((Math.PI * x) / Math.max(1, nx - 1)) *
            Math.sin((Math.PI * y) / Math.max(1, ny - 1));
          vx += g[0] * env;
          vy += g[1] * env;
        }

        // 放出渦（カルマン渦列風）
        for (let k = 0; k < vortices.length; k++) {
          const g = gaussVortex(x, y, vortices[k]!);
          vx += g[0];
          vy += g[1];
        }

        // bounceBack 壁の近くでは滑りなし風に減衰
        if (edges.bottom.type === "bounceBack" && fy < 0.06) {
          const w = fy / 0.06;
          vx *= w;
          vy *= w;
        }
        if (edges.top.type === "bounceBack" && fy > 0.94) {
          const w = (1 - fy) / 0.06;
          vx *= w;
          vy *= w;
        }
        const fx = nx > 1 ? x / (nx - 1) : 0;
        if (edges.left.type === "bounceBack" && fx < 0.06) {
          const w = fx / 0.06;
          vx *= w;
          vy *= w;
        }
        if (edges.right.type === "bounceBack" && fx > 0.94) {
          const w = (1 - fx) / 0.06;
          vx *= w;
          vy *= w;
        }

        // 障害物近傍の減速
        const d = damp[i]!;
        vx *= d;
        vy *= d;

        // BGK のノイズ演出
        if (noiseAmp > 0) {
          vx += noiseAmp * MockEngine.hashNoise(x, y, t) * 0.5;
          vy += noiseAmp * MockEngine.hashNoise(x + 7919, y, t) * 0.5;
        }

        // 速度上限（格子単位で |u| < 0.3 目安）
        const sp2 = vx * vx + vy * vy;
        if (sp2 > 0.09) {
          const f = 0.3 / Math.sqrt(sp2);
          vx *= f;
          vy *= f;
        }

        ux[i] = vx;
        uy[i] = vy;

        // 密度: ベルヌーイ的に速いところで低く、渦コアでさらに低く
        let r = 1 - 1.4 * sp2;
        if (noiseAmp > 0) {
          r += noiseAmp * 0.12 * MockEngine.hashNoise(x, y + 104729, t);
        }
        rho[i] = r;
      }
    }
  }
}

/** ガウス渦: 中心 (v.x, v.y)・強さ v.s・半径 v.r の回転速度場 */
function gaussVortex(
  x: number,
  y: number,
  v: { x: number; y: number; s: number; r: number },
): [number, number] {
  const dx = x - v.x;
  const dy = y - v.y;
  const r2 = dx * dx + dy * dy;
  const rr = v.r * v.r;
  if (r2 > rr * 16) return [0, 0];
  const f = (v.s * Math.exp(-r2 / (2 * rr))) / v.r;
  return [-dy * f, dx * f];
}

/** 起動直後の立ち上がり（キャビティのせん断が徐々に発達する演出） */
function decay0(t: number): number {
  return 1 - Math.exp(-t / 300);
}
