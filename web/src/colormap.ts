/**
 * Colormap LUT (256 steps x RGB).
 * Linearly interpolates anchor colors, generated once at startup. No external dependencies.
 */

export type Lut = Uint8ClampedArray; // length 256 * 3

/** Build a 256-step LUT from a sequence of anchor colors */
function buildLut(stops: ReadonlyArray<readonly [number, number, number]>): Lut {
  const lut = new Uint8ClampedArray(256 * 3);
  const nSeg = stops.length - 1;
  for (let i = 0; i < 256; i++) {
    const t = (i / 255) * nSeg;
    const s = Math.min(nSeg - 1, Math.floor(t));
    const f = t - s;
    const a = stops[s]!;
    const b = stops[s + 1]!;
    lut[i * 3] = a[0] + (b[0] - a[0]) * f;
    lut[i * 3 + 1] = a[1] + (b[1] - a[1]) * f;
    lut[i * 3 + 2] = a[2] + (b[2] - a[2]) * f;
  }
  return lut;
}

/** viridis (9 anchor points from matplotlib) */
export const VIRIDIS: Lut = buildLut([
  [68, 1, 84],
  [72, 40, 120],
  [62, 74, 137],
  [49, 104, 142],
  [38, 130, 142],
  [31, 158, 137],
  [53, 183, 121],
  [109, 205, 89],
  [180, 222, 44],
  [253, 231, 37],
]);

/** RdBu (ColorBrewer 11-class, diverging. Reversed so negative=blue / 0=white / positive=red) */
export const RDBU: Lut = buildLut([
  [5, 48, 97],
  [33, 102, 172],
  [67, 147, 195],
  [146, 197, 222],
  [209, 229, 240],
  [247, 247, 247],
  [253, 219, 199],
  [244, 165, 130],
  [214, 96, 77],
  [178, 24, 43],
  [103, 0, 31],
]);
