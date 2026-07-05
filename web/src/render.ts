import type { Engine } from "./engine/types.ts";
import { RDBU, VIRIDIS, type Lut } from "./colormap.ts";

export type VisMode = "speed" | "vorticity" | "density";

export const VIS_MODE_LABEL: Record<VisMode, string> = {
  speed: "Speed |u|",
  vorticity: "Vorticity ω",
  density: "Density ρ",
};

export const VIS_MODE_HINT: Record<VisMode, string> = {
  speed: "Faster flow is shown in brighter colors.",
  vorticity: "Counterclockwise vortices are shown in red, clockwise vortices in blue.",
  density: "Areas denser than average are shown in red, sparser areas in blue.",
};

const SOLID_RGB: readonly [number, number, number] = [110, 118, 129];

export interface BrushPreview {
  /** Grid coordinates (x: 0..nx, y: 0..ny; y is 0 at the bottom edge) */
  gx: number;
  gy: number;
  radius: number;
  erase: boolean;
}

/** Display range (for the colorbar labels) */
export interface VisRange {
  lo: number;
  hi: number;
}

/**
 * Handles everything from scalarizing the field -> normalizing -> LUT
 * coloring -> offscreen canvas(nx*ny) -> upscaled blit to the main canvas.
 */
export class FieldRenderer {
  private off: HTMLCanvasElement;
  private offCtx: CanvasRenderingContext2D;
  private img: ImageData | null = null;
  private scalar = new Float32Array(0);

  /** Smooth the range with an exponential moving average to avoid flicker */
  private emaHi = 0;
  private emaMode: VisMode | null = null;

  constructor(private main: HTMLCanvasElement) {
    this.off = document.createElement("canvas");
    const ctx = this.off.getContext("2d");
    if (!ctx) throw new Error("Failed to get 2D context");
    this.offCtx = ctx;
  }

  /** Restart range learning, e.g. on reset */
  resetRange(): void {
    this.emaHi = 0;
    this.emaMode = null;
  }

  render(engine: Engine, mode: VisMode, brush: BrushPreview | null): VisRange {
    const nx = engine.nx;
    const ny = engine.ny;

    if (this.off.width !== nx || this.off.height !== ny || !this.img) {
      this.off.width = nx;
      this.off.height = ny;
      this.img = this.offCtx.createImageData(nx, ny);
      this.scalar = new Float32Array(nx * ny);
    }
    if (this.emaMode !== mode) {
      this.emaHi = 0;
      this.emaMode = mode;
    }

    const scalar = this.computeScalar(engine, mode);
    const range = this.updateRange(mode, scalar);
    const lut: Lut = mode === "speed" ? VIRIDIS : RDBU;

    this.paint(engine, scalar, lut, range);
    this.blit(engine, brush);
    return range;
  }

  // ------------------------------------------------------------ Scalarization

  private computeScalar(engine: Engine, mode: VisMode): Float32Array {
    const nx = engine.nx;
    const ny = engine.ny;
    const out = this.scalar;

    if (mode === "speed") {
      const ux = engine.ux();
      const uy = engine.uy();
      for (let i = 0; i < out.length; i++) {
        out[i] = Math.hypot(ux[i]!, uy[i]!);
      }
    } else if (mode === "density") {
      out.set(engine.rho());
    } else {
      // Vorticity ω = ∂uy/∂x − ∂ux/∂y (central difference, 0 at edges)
      const ux = engine.ux();
      const uy = engine.uy();
      out.fill(0);
      for (let y = 1; y < ny - 1; y++) {
        const row = y * nx;
        for (let x = 1; x < nx - 1; x++) {
          const i = row + x;
          const duy_dx = (uy[i + 1]! - uy[i - 1]!) * 0.5;
          const dux_dy = (ux[i + nx]! - ux[i - nx]!) * 0.5;
          out[i] = duy_dx - dux_dy;
        }
      }
    }
    return out;
  }

  // ---------------------------------------------------------- Range computation

  private updateRange(mode: VisMode, scalar: Float32Array): VisRange {
    let target: number;
    if (mode === "density") {
      let m = 0;
      for (let i = 0; i < scalar.length; i++) {
        const d = Math.abs(scalar[i]! - 1);
        if (d > m) m = d;
      }
      target = Math.max(1e-4, m);
    } else {
      let m = 0;
      for (let i = 0; i < scalar.length; i++) {
        const d = Math.abs(scalar[i]!);
        if (d > m) m = d;
      }
      target = Math.max(1e-4, m * (mode === "vorticity" ? 0.7 : 1));
    }

    // Fall back to the previous value so the range doesn't break even if
    // things diverge (Inf)
    if (!Number.isFinite(target)) target = this.emaHi > 0 ? this.emaHi : 1;

    this.emaHi = this.emaHi === 0 ? target : this.emaHi * 0.92 + target * 0.08;
    const hi = this.emaHi;

    if (mode === "speed") return { lo: 0, hi };
    if (mode === "vorticity") return { lo: -hi, hi };
    return { lo: 1 - hi, hi: 1 + hi }; // density
  }

  // ------------------------------------------------------------- Coloring

  private paint(
    engine: Engine,
    scalar: Float32Array,
    lut: Lut,
    range: VisRange,
  ): void {
    const nx = engine.nx;
    const ny = engine.ny;
    const solid = engine.solidMask();
    const img = this.img!;
    const px = img.data;
    const inv = 1 / (range.hi - range.lo);

    for (let y = 0; y < ny; y++) {
      const srcRow = y * nx;
      const dstRow = (ny - 1 - y) * nx; // Flip vertically (y=0 is the bottom edge)
      for (let x = 0; x < nx; x++) {
        const si = srcRow + x;
        const di = (dstRow + x) * 4;
        if (solid[si] === 1) {
          px[di] = SOLID_RGB[0];
          px[di + 1] = SOLID_RGB[1];
          px[di + 2] = SOLID_RGB[2];
          px[di + 3] = 255;
          continue;
        }
        let t = (scalar[si]! - range.lo) * inv;
        // !(t > 0) also catches NaN: guards against garbled rendering on divergence
        if (!(t > 0)) t = 0;
        else if (t > 1) t = 1;
        const li = ((t * 255) | 0) * 3;
        px[di] = lut[li]!;
        px[di + 1] = lut[li + 1]!;
        px[di + 2] = lut[li + 2]!;
        px[di + 3] = 255;
      }
    }
    this.offCtx.putImageData(img, 0, 0);
  }

  // ----------------------------------------------------------- Upscaled blit

  private blit(engine: Engine, brush: BrushPreview | null): void {
    const ctx = this.main.getContext("2d");
    if (!ctx) return;
    const w = this.main.width;
    const h = this.main.height;
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = "high";
    ctx.clearRect(0, 0, w, h);
    ctx.drawImage(this.off, 0, 0, w, h);

    if (brush) {
      const sx = w / engine.nx;
      const sy = h / engine.ny;
      const cx = brush.gx * sx;
      const cy = (engine.ny - brush.gy) * sy; // Flip vertically
      ctx.beginPath();
      ctx.ellipse(cx, cy, brush.radius * sx, brush.radius * sy, 0, 0, Math.PI * 2);
      ctx.strokeStyle = brush.erase ? "rgba(255,120,120,0.9)" : "rgba(255,255,255,0.9)";
      ctx.lineWidth = Math.max(1, w / 500);
      ctx.setLineDash([4, 4]);
      ctx.stroke();
      ctx.setLineDash([]);
    }
  }
}

/** Draw the LUT onto the panel's colorbar canvas */
export function drawColorbar(canvas: HTMLCanvasElement, mode: VisMode): void {
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const lut: Lut = mode === "speed" ? VIRIDIS : RDBU;
  const w = canvas.width;
  const h = canvas.height;
  const img = ctx.createImageData(w, h);
  for (let x = 0; x < w; x++) {
    const li = (((x / (w - 1)) * 255) | 0) * 3;
    for (let y = 0; y < h; y++) {
      const di = (y * w + x) * 4;
      img.data[di] = lut[li]!;
      img.data[di + 1] = lut[li + 1]!;
      img.data[di + 2] = lut[li + 2]!;
      img.data[di + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
}

/** Formatter for the colorbar's end labels */
export function formatRange(v: number): string {
  if (!Number.isFinite(v)) return "—";
  const a = Math.abs(v);
  if (a === 0) return "0";
  if (a >= 0.01 && a < 1000) return v.toPrecision(2);
  return v.toExponential(1);
}
