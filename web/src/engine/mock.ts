import type { Engine, EngineConfig } from "./types.ts";

let warnedMockActivation = false;

/**
 * Mock engine.
 *
 * Does not solve the real LBM; instead analytically synthesizes a
 * "plausible-looking" 2D flow field as a function of time t. For UI
 * development and demo purposes.
 *
 * - Picks a base flow from the combination of boundary conditions:
 *   - top movingWall            -> cavity-flow-like (shear + primary vortex swirl)
 *   - left velocityInlet        -> uniform flow + Karman-vortex-street-like
 *                                  vortex shedding downstream of the obstacle
 *   - external force + top/bottom bounceBack -> Poiseuille parabolic profile
 *                                  + traveling wave
 *   - all periodic              -> decaying Taylor-Green vortex
 * - Viscosity nu affects the vortex decay rate / vortex street amplitude
 *   (larger nu settles down faster).
 * - When collision="bgk", overlays a small amount of noise; "trt" stays
 *   clean (a demo flourish suggesting "TRT is more stable/accurate" —
 *   not real physics).
 * - Obstacle cells have u=0, rho=1. Surrounding cells are decelerated
 *   based on distance to make them look wall-like.
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

  /** Deceleration factor based on distance from the obstacle (0=inside obstacle ... 1=far enough away) */
  private damp = new Float32Array(0);
  private dampDirty = true;

  /** Obstacle centroid and characteristic radius (used as the vortex-shedding origin) */
  private solidCx = 0;
  private solidCy = 0;
  private solidR = 0;
  private solidCount = 0;

  constructor() {
    if (!warnedMockActivation) {
      console.warn(
        "LBMFlow: mock engine active. Output is synthetic UI fallback data, not a simulation result.",
      );
      warnedMockActivation = true;
    }
  }

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
      // Immediately zero out the flow velocity so that becoming a wall is
      // visible even while paused
      this._ux[i] = 0;
      this._uy[i] = 0;
      this._rho[i] = 1;
    }
  }

  // ---------------------------------------------------------------- internal implementation

  /** Recompute the deceleration map and centroid/radius from the obstacle mask */
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

    // Decelerate cells that have an obstacle nearby (simple check within
    // Chebyshev distance of 2 cells)
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

  /** Deterministic pseudo-random number (from cell coordinates and time). For the BGK noise flourish */
  private static hashNoise(x: number, y: number, t: number): number {
    let h = (x * 374761393 + y * 668265263 + t * 2246822519) | 0;
    h = Math.imul(h ^ (h >>> 13), 1274126177);
    h ^= h >>> 16;
    return (h & 0xffff) / 0xffff - 0.5; // [-0.5, 0.5)
  }

  /** Fully rebuild the field at the current time _time */
  private recompute(): void {
    if (this.dampDirty) this.refreshSolidInfo();

    const { nx, ny, nu, edges, force, collision } = this.cfg;
    const t = this._time;
    const rho = this._rho;
    const ux = this._ux;
    const uy = this._uy;
    const solid = this._solid;
    const damp = this.damp;

    // ---- Determine the base flow type
    const top = edges.top;
    const left = edges.left;
    const lidU = top.type === "movingWall" ? top.u[0] : 0;
    const inletU = left.type === "velocityInlet" ? left.u[0] : 0;
    const hasForce = force[0] !== 0 || force[1] !== 0;

    // Poiseuille flow max velocity Umax = f H^2 / (8 nu) in lattice units (capped to avoid blow-up)
    const H = Math.max(2, ny - 2);
    const poisMax = hasForce
      ? Math.min(0.18, Math.abs((force[0] * H * H) / (8 * Math.max(nu, 1e-5))))
      : 0;
    const poisSign = force[0] >= 0 ? 1 : -1;

    // Overall decay from viscosity (reuses the Taylor-Green vortex exact
    // solution e^{-2*nu*k^2*t})
    const k0 = (2 * Math.PI) / Math.max(nx, ny);
    const decay = Math.exp(-2 * nu * k0 * k0 * 4 * t);

    // ---- Karman-vortex-street-like shed vortices (only when there's an obstacle + inflow)
    // Each vortex is carried downstream by the uniform flow from its birth
    // time and decays via viscosity.
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

    // ---- Cavity primary vortex (center slowly precesses)
    let cav: Vortex | null = null;
    if (lidU !== 0) {
      const wob = 0.06 * Math.min(nx, ny);
      cav = {
        x: nx * 0.5 + wob * Math.cos(t * 0.004),
        y: ny * 0.62 + wob * 0.6 * Math.sin(t * 0.0031),
        s: -lidU * 1.6 * (1 - Math.exp(-t / 500)), // gradually develops
        r: Math.min(nx, ny) * 0.3,
      };
    }

    // ---- Taylor-Green vortex (for periodic boundaries)
    const isPeriodic =
      edges.left.type === "periodic" &&
      edges.right.type === "periodic" &&
      edges.top.type === "periodic" &&
      edges.bottom.type === "periodic";
    const tgAmp = isPeriodic && !hasForce ? 0.09 : 0;
    const kx = (2 * Math.PI * 2) / nx;
    const ky = (2 * Math.PI * 2) / ny;
    const driftX = 8e-3 * nx * 1e-2; // slowly drifts sideways
    const tgPhase = t * driftX * kx * 0.5;

    const noiseAmp = collision === "bgk" ? 0.02 : 0;

    for (let y = 0; y < ny; y++) {
      const fy = ny > 1 ? y / (ny - 1) : 0; // 0(bottom edge) .. 1(top edge)
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

        // Uniform inflow
        if (inletU !== 0) {
          vx += inletU;
          if (left.type === "velocityInlet") vy += left.u[1];
        }

        // Poiseuille parabola + weak traveling wave for visual effect
        if (poisMax > 0) {
          const par = 4 * fy * (1 - fy);
          vx += poisSign * poisMax * par;
          vy +=
            0.06 *
            poisMax *
            par *
            Math.sin((2 * Math.PI * (x - poisSign * poisMax * t * 0.6)) / (nx * 0.5));
        }

        // Top-wall-driven shear (stronger near the wall)
        if (lidU !== 0) {
          vx += lidU * Math.pow(fy, 3) * decay0(t);
        }

        // Taylor-Green vortex (decay + drift)
        if (tgAmp > 0) {
          const px = kx * x - tgPhase;
          const py = ky * y;
          vx += tgAmp * decay * Math.cos(px) * Math.sin(py);
          vy += -tgAmp * decay * Math.sin(px) * Math.cos(py);
        }

        // Cavity primary vortex
        if (cav) {
          const g = gaussVortex(x, y, cav);
          // Apply a sin envelope so it vanishes at the wall
          const env =
            Math.sin((Math.PI * x) / Math.max(1, nx - 1)) *
            Math.sin((Math.PI * y) / Math.max(1, ny - 1));
          vx += g[0] * env;
          vy += g[1] * env;
        }

        // Shed vortices (Karman-vortex-street-like)
        for (let k = 0; k < vortices.length; k++) {
          const g = gaussVortex(x, y, vortices[k]!);
          vx += g[0];
          vy += g[1];
        }

        // Decelerate near bounceBack walls to look no-slip-like
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

        // Deceleration near the obstacle
        const d = damp[i]!;
        vx *= d;
        vy *= d;

        // BGK noise flourish
        if (noiseAmp > 0) {
          vx += noiseAmp * MockEngine.hashNoise(x, y, t) * 0.5;
          vy += noiseAmp * MockEngine.hashNoise(x + 7919, y, t) * 0.5;
        }

        // Velocity cap (roughly |u| < 0.3 in lattice units)
        const sp2 = vx * vx + vy * vy;
        if (sp2 > 0.09) {
          const f = 0.3 / Math.sqrt(sp2);
          vx *= f;
          vy *= f;
        }

        ux[i] = vx;
        uy[i] = vy;

        // Density: Bernoulli-like, lower where the flow is fast, even
        // lower in the vortex core
        let r = 1 - 1.4 * sp2;
        if (noiseAmp > 0) {
          r += noiseAmp * 0.12 * MockEngine.hashNoise(x, y + 104729, t);
        }
        rho[i] = r;
      }
    }
  }
}

/** Gaussian vortex: rotational velocity field with center (v.x, v.y), strength v.s, radius v.r */
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

/** Startup ramp-up (a flourish making cavity shear develop gradually) */
function decay0(t: number): number {
  return 1 - Math.exp(-t / 300);
}
