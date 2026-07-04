import "./style.css";
import { createEngine } from "./engine/index.ts";
import type { EngineConfig } from "./engine/types.ts";
import { PRESETS, type Preset } from "./presets.ts";
import {
  FieldRenderer,
  drawColorbar,
  formatRange,
  VIS_MODE_HINT,
  VIS_MODE_LABEL,
  type BrushPreview,
  type VisMode,
} from "./render.ts";

// ------------------------------------------------------------- DOM ヘルパ

function $<T extends HTMLElement>(id: string): T {
  const el = document.getElementById(id);
  if (!el) throw new Error(`要素が見つかりません: #${id}`);
  return el as T;
}

const presetSelect = $<HTMLSelectElement>("preset-select");
const btnRun = $<HTMLButtonElement>("btn-run");
const btnRunIcon = $<HTMLSpanElement>("btn-run-icon");
const btnRunLabel = $<HTMLSpanElement>("btn-run-label");
const btnReset = $<HTMLButtonElement>("btn-reset");
const canvas = $<HTMLCanvasElement>("sim-canvas");
const canvasWrap = $<HTMLDivElement>("canvas-wrap");
const presetDesc = $<HTMLParagraphElement>("preset-desc");
const visModeSelect = $<HTMLSelectElement>("vis-mode");
const visHint = $<HTMLParagraphElement>("vis-hint");
const colorbar = $<HTMLCanvasElement>("colorbar");
const colorbarMin = $<HTMLSpanElement>("colorbar-min");
const colorbarMax = $<HTMLSpanElement>("colorbar-max");
const brushDrawBtn = $<HTMLButtonElement>("brush-draw");
const brushEraseBtn = $<HTMLButtonElement>("brush-erase");
const brushSize = $<HTMLInputElement>("brush-size");
const brushSizeValue = $<HTMLSpanElement>("brush-size-value");
const btnClearSolids = $<HTMLButtonElement>("btn-clear-solids");
const nuSlider = $<HTMLInputElement>("nu-slider");
const nuValue = $<HTMLSpanElement>("nu-value");
const resSelect = $<HTMLSelectElement>("res-select");
const spfSlider = $<HTMLInputElement>("spf-slider");
const spfValue = $<HTMLSpanElement>("spf-value");
const collisionSelect = $<HTMLSelectElement>("collision-select");
const statusStep = $<HTMLElement>("status-step");
const statusSps = $<HTMLElement>("status-sps");
const statusGrid = $<HTMLElement>("status-grid");
const statusMode = $<HTMLElement>("status-mode");

// ------------------------------------------------------------- 状態

const engine = createEngine();
const renderer = new FieldRenderer(canvas);

let currentPreset: Preset = PRESETS[0]!;
let running = false;
let visMode: VisMode = "speed";
let stepsPerFrame = Number(spfSlider.value);
let brushErase = false;
let brushPreview: BrushPreview | null = null;
let painting = false;
let paintEraseOverride: boolean | null = null; // 右ドラッグは常に消しゴム
let lastPaint: { gx: number; gy: number } | null = null;

// steps/s 計測
let spsSteps = 0;
let spsT0 = performance.now();

// ---------------------------------------------------------- 粘性スライダー
// 対数スケール: slider 0..100 → ν ∈ [1e-4, 0.5]

const NU_LOG_MIN = Math.log10(1e-4);
const NU_LOG_MAX = Math.log10(0.5);

function sliderToNu(v: number): number {
  return 10 ** (NU_LOG_MIN + (v / 100) * (NU_LOG_MAX - NU_LOG_MIN));
}

function nuToSlider(nu: number): number {
  const t = (Math.log10(nu) - NU_LOG_MIN) / (NU_LOG_MAX - NU_LOG_MIN);
  return Math.round(Math.min(100, Math.max(0, t * 100)));
}

function currentNu(): number {
  return sliderToNu(Number(nuSlider.value));
}

function updateNuLabel(): void {
  nuValue.textContent = currentNu().toPrecision(2);
}

// ------------------------------------------------------------ リセット処理

function scaledDims(preset: Preset): { nx: number; ny: number } {
  const s = Number(resSelect.value);
  return {
    nx: Math.max(32, Math.round(preset.config.nx * s)),
    ny: Math.max(32, Math.round(preset.config.ny * s)),
  };
}

function buildConfig(preset: Preset): EngineConfig {
  const { nx, ny } = scaledDims(preset);
  return {
    ...preset.config,
    nx,
    ny,
    nu: currentNu(),
    collision: collisionSelect.value === "bgk" ? "bgk" : "trt",
    edges: { ...preset.config.edges },
    force: [preset.config.force[0], preset.config.force[1]],
  };
}

/**
 * シミュレーションを初期化する。
 * preserveSolids=true なら現在の障害物を（解像度が変わる場合は最近傍で
 * スケールして）引き継ぐ。false ならプリセットの初期障害物を配置する。
 */
function resetSim(preserveSolids: boolean): void {
  let oldMask: Uint8Array | null = null;
  let oldNx = 0;
  let oldNy = 0;
  if (preserveSolids && engine.nx > 0) {
    oldMask = new Uint8Array(engine.solidMask()); // コピー
    oldNx = engine.nx;
    oldNy = engine.ny;
  }

  engine.init(buildConfig(currentPreset));

  if (oldMask && oldNx > 0 && oldNy > 0) {
    for (let y = 0; y < engine.ny; y++) {
      const sy = Math.min(oldNy - 1, Math.floor((y / engine.ny) * oldNy));
      for (let x = 0; x < engine.nx; x++) {
        const sx = Math.min(oldNx - 1, Math.floor((x / engine.nx) * oldNx));
        if (oldMask[sy * oldNx + sx] === 1) engine.setSolid(x, y, true);
      }
    }
  } else {
    currentPreset.paintObstacles?.(engine);
  }

  renderer.resetRange();
  spsSteps = 0;
  spsT0 = performance.now();
  statusSps.textContent = "—";
  statusGrid.textContent = `${engine.nx}×${engine.ny}`;
  fitCanvas();
}

function applyPreset(preset: Preset): void {
  currentPreset = preset;
  presetDesc.textContent = preset.description;
  nuSlider.value = String(nuToSlider(preset.config.nu));
  collisionSelect.value = preset.config.collision;
  updateNuLabel();
  resetSim(false);
}

// ---------------------------------------------------------- 実行 / 停止

function setRunning(v: boolean): void {
  running = v;
  btnRunIcon.textContent = v ? "⏸" : "▶";
  btnRunLabel.textContent = v ? "停止" : "実行";
  btnRun.classList.toggle("btn-primary", !v);
  btnRun.classList.toggle("btn-running", v);
  if (!v) statusSps.textContent = "—";
  spsSteps = 0;
  spsT0 = performance.now();
}

// ------------------------------------------------------- キャンバスサイズ

/** ラッパー内にアスペクト比を保って収まるよう canvas の実寸を決める */
function fitCanvas(): void {
  const rect = canvasWrap.getBoundingClientRect();
  if (rect.width < 4 || rect.height < 4 || engine.nx === 0) return;
  const aspect = engine.nx / engine.ny;
  let w = rect.width;
  let h = w / aspect;
  if (h > rect.height) {
    h = rect.height;
    w = h * aspect;
  }
  canvas.style.width = `${w}px`;
  canvas.style.height = `${h}px`;
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const bw = Math.max(1, Math.round(w * dpr));
  const bh = Math.max(1, Math.round(h * dpr));
  if (canvas.width !== bw) canvas.width = bw;
  if (canvas.height !== bh) canvas.height = bh;
}

new ResizeObserver(() => fitCanvas()).observe(canvasWrap);

// ------------------------------------------------------------ ペイント

function eventToGrid(e: PointerEvent): { gx: number; gy: number } {
  const rect = canvas.getBoundingClientRect();
  const gx = ((e.clientX - rect.left) / rect.width) * engine.nx;
  const gy = engine.ny - ((e.clientY - rect.top) / rect.height) * engine.ny;
  return { gx, gy };
}

function paintDisk(gx: number, gy: number, radius: number, solid: boolean): void {
  const r2 = radius * radius;
  const x0 = Math.max(0, Math.floor(gx - radius));
  const x1 = Math.min(engine.nx - 1, Math.ceil(gx + radius));
  const y0 = Math.max(0, Math.floor(gy - radius));
  const y1 = Math.min(engine.ny - 1, Math.ceil(gy + radius));
  for (let y = y0; y <= y1; y++) {
    for (let x = x0; x <= x1; x++) {
      const dx = x + 0.5 - gx;
      const dy = y + 0.5 - gy;
      if (dx * dx + dy * dy <= r2) engine.setSolid(x, y, solid);
    }
  }
}

/** 前回位置から補間しながら塗る（速いドラッグでも途切れない） */
function paintStroke(gx: number, gy: number): void {
  const erase = paintEraseOverride ?? brushErase;
  const radius = Number(brushSize.value);
  if (lastPaint) {
    const dist = Math.hypot(gx - lastPaint.gx, gy - lastPaint.gy);
    const steps = Math.max(1, Math.ceil(dist / (radius * 0.5)));
    for (let s = 1; s <= steps; s++) {
      const f = s / steps;
      paintDisk(
        lastPaint.gx + (gx - lastPaint.gx) * f,
        lastPaint.gy + (gy - lastPaint.gy) * f,
        radius,
        !erase,
      );
    }
  } else {
    paintDisk(gx, gy, radius, !erase);
  }
  lastPaint = { gx, gy };
}

canvas.addEventListener("contextmenu", (e) => e.preventDefault());

canvas.addEventListener("pointerdown", (e) => {
  if (e.button !== 0 && e.button !== 2) return;
  e.preventDefault();
  canvas.setPointerCapture(e.pointerId);
  painting = true;
  paintEraseOverride = e.button === 2 ? true : null;
  lastPaint = null;
  const { gx, gy } = eventToGrid(e);
  paintStroke(gx, gy);
});

canvas.addEventListener("pointermove", (e) => {
  const { gx, gy } = eventToGrid(e);
  brushPreview = {
    gx,
    gy,
    radius: Number(brushSize.value),
    erase: paintEraseOverride ?? brushErase,
  };
  if (painting) paintStroke(gx, gy);
});

canvas.addEventListener("pointerup", (e) => {
  if (painting) canvas.releasePointerCapture(e.pointerId);
  painting = false;
  paintEraseOverride = null;
  lastPaint = null;
});

canvas.addEventListener("pointerleave", () => {
  brushPreview = null;
});

// ------------------------------------------------------------- UI 配線

for (const p of PRESETS) {
  const opt = document.createElement("option");
  opt.value = p.id;
  opt.textContent = p.name;
  presetSelect.appendChild(opt);
}

presetSelect.addEventListener("change", () => {
  const p = PRESETS.find((q) => q.id === presetSelect.value);
  if (p) applyPreset(p);
});

btnRun.addEventListener("click", () => setRunning(!running));
btnReset.addEventListener("click", () => resetSim(false));

visModeSelect.addEventListener("change", () => {
  visMode = visModeSelect.value as VisMode;
  statusMode.textContent = VIS_MODE_LABEL[visMode];
  visHint.textContent = VIS_MODE_HINT[visMode];
  drawColorbar(colorbar, visMode);
});

function setBrushErase(v: boolean): void {
  brushErase = v;
  brushDrawBtn.classList.toggle("active", !v);
  brushEraseBtn.classList.toggle("active", v);
}
brushDrawBtn.addEventListener("click", () => setBrushErase(false));
brushEraseBtn.addEventListener("click", () => setBrushErase(true));

brushSize.addEventListener("input", () => {
  brushSizeValue.textContent = brushSize.value;
});

btnClearSolids.addEventListener("click", () => {
  for (let y = 0; y < engine.ny; y++) {
    for (let x = 0; x < engine.nx; x++) {
      engine.setSolid(x, y, false);
    }
  }
});

nuSlider.addEventListener("input", updateNuLabel);
nuSlider.addEventListener("change", () => resetSim(true));
collisionSelect.addEventListener("change", () => resetSim(true));
resSelect.addEventListener("change", () => resetSim(true));

spfSlider.addEventListener("input", () => {
  stepsPerFrame = Number(spfSlider.value);
  spfValue.textContent = spfSlider.value;
});

// Space で実行 / 停止（フォーム要素にフォーカスがあるときは無効）
window.addEventListener("keydown", (e) => {
  if (e.code !== "Space") return;
  const t = e.target;
  if (
    t instanceof HTMLInputElement ||
    t instanceof HTMLSelectElement ||
    t instanceof HTMLTextAreaElement ||
    t instanceof HTMLButtonElement
  ) {
    return;
  }
  e.preventDefault();
  setRunning(!running);
});

// タブが隠れたら自動停止
document.addEventListener("visibilitychange", () => {
  if (document.hidden && running) setRunning(false);
});

// ------------------------------------------------------------ メインループ

function frame(): void {
  if (running) {
    engine.step(stepsPerFrame);
    spsSteps += stepsPerFrame;
    const now = performance.now();
    const dt = now - spsT0;
    if (dt >= 500) {
      statusSps.textContent = Math.round((spsSteps * 1000) / dt).toLocaleString("ja-JP");
      spsSteps = 0;
      spsT0 = now;
    }
  }

  const range = renderer.render(engine, visMode, brushPreview);
  colorbarMin.textContent = formatRange(range.lo);
  colorbarMax.textContent = formatRange(range.hi);
  statusStep.textContent = engine.time.toLocaleString("ja-JP");

  requestAnimationFrame(frame);
}

// ------------------------------------------------------------- 起動

drawColorbar(colorbar, visMode);
statusMode.textContent = VIS_MODE_LABEL[visMode];
visHint.textContent = VIS_MODE_HINT[visMode];
spfValue.textContent = spfSlider.value;
brushSizeValue.textContent = brushSize.value;
presetSelect.value = currentPreset.id;
applyPreset(currentPreset);
requestAnimationFrame(frame);
