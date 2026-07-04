/**
 * カラーマップ LUT（256 段 × RGB）。
 * アンカー色を線形補間して起動時に一度だけ生成する。外部依存なし。
 */

export type Lut = Uint8ClampedArray; // 長さ 256 * 3

/** アンカー色列から 256 段の LUT を作る */
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

/** viridis（matplotlib 由来のアンカー 9 点） */
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

/** RdBu（ColorBrewer 11 クラス、発散型。負=青 / 0=白 / 正=赤 になるよう反転済み） */
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
