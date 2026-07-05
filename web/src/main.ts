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
import { buildScenario } from "./scenario.ts";
import { showToast } from "./toast.ts";

// ------------------------------------------------------------- DOM helpers

function $<T extends HTMLElement>(id: string): T {
  const el = document.getElementById(id);
  if (!el) throw new Error(`Element not found: #${id}`);
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
const statusMlups = $<HTMLElement>("status-mlups");
const statusGrid = $<HTMLElement>("status-grid");
const statusMode = $<HTMLElement>("status-mode");
const colorbarMid = $<HTMLSpanElement>("colorbar-mid");
const btnExport = $<HTMLButtonElement>("btn-export");
const firstHint = $<HTMLDivElement>("first-hint");
const firstHintClose = $<HTMLButtonElement>("first-hint-close");

// ------------------------------------------------------------- State

const engine = await createEngine();
const renderer = new FieldRenderer(canvas);

let currentPreset: Preset = PRESETS[0]!;
let running = false;
let diverged = false;
let visMode: VisMode = "speed";
let stepsPerFrame = Number(spfSlider.value);
let brushErase = false;
let brushPreview: BrushPreview | null = null;
let painting = false;
let paintEraseOverride: boolean | null = null; // Right-drag is always erase
let lastPaint: { gx: number; gy: number } | null = null;

// steps/s measurement
let spsSteps = 0;
let spsT0 = performance.now();

// ---------------------------------------------------------- Viscosity slider
// Logarithmic scale: slider 0..100 -> ν ∈ [1e-4, 0.5]

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

// ------------------------------------------------------------ Reset handling

function scaledDims(preset: Preset): { nx: number; ny: number } {
  const s = Number(resSelect.value);
  return {
    nx: Math.max(32, Math.round(preset.config.nx * s)),
    ny: Math.max(32, Math.round(preset.config.ny * s)),
  };
}

function buildConfig(preset: Preset): EngineConfig {
  const { nx, ny } = scaledDims(preset);
  const scale = nx / preset.config.nx;
  return {
    ...preset.config,
    nx,
    ny,
    nu: currentNu(),
    collision: collisionSelect.value === "bgk" ? "bgk" : "trt",
    edges: { ...preset.config.edges },
    force: [preset.config.force[0], preset.config.force[1]],
    init: preset.config.init
      ? {
          ...preset.config.init,
          cx: preset.config.init.cx * scale,
          cy: preset.config.init.cy * scale,
          r: preset.config.init.r * scale,
        }
      : undefined,
  };
}

/**
 * Initialize the simulation.
 * If preserveSolids=true, carries over the current obstacles (scaled with
 * nearest-neighbor if the resolution changes). If false, places the preset's
 * initial obstacles.
 */
function resetSim(preserveSolids: boolean): void {
  let oldMask: Uint8Array | null = null;
  let oldNx = 0;
  let oldNy = 0;
  if (preserveSolids && engine.nx > 0) {
    oldMask = new Uint8Array(engine.solidMask()); // Copy
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
  diverged = false;
  spsSteps = 0;
  spsT0 = performance.now();
  statusSps.textContent = "—";
  statusMlups.textContent = "—";
  statusGrid.textContent = `${engine.nx}×${engine.ny}`;
  fitCanvas();
}

function applyPreset(preset: Preset): void {
  currentPreset = preset;
  presetDesc.textContent = preset.description;
  nuSlider.value = String(nuToSlider(preset.config.nu));
  collisionSelect.value = preset.config.collision;
  if (preset.defaultVis) {
    visMode = preset.defaultVis;
    visModeSelect.value = preset.defaultVis;
    visHint.textContent = VIS_MODE_HINT[visMode];
    statusMode.textContent = VIS_MODE_LABEL[visMode];
  }
  updateNuLabel();
  resetSim(false);
}

// ---------------------------------------------------------- Run / Stop

function setRunning(v: boolean): void {
  running = v;
  btnRunIcon.textContent = v ? "⏸" : "▶";
  btnRunLabel.textContent = v ? "Stop" : "Run";
  btnRun.classList.toggle("btn-primary", !v);
  btnRun.classList.toggle("btn-running", v);
  if (!v) {
    statusSps.textContent = "—";
    statusMlups.textContent = "—";
  } else {
    dismissFirstHint(true);
  }
  spsSteps = 0;
  spsT0 = performance.now();
}

// ------------------------------------------------------------ Divergence detection

/**
 * Treat the simulation as diverged if the fields show NaN / Inf / a
 * physically impossible velocity (|u| > 5). NaN is caught via isFinite
 * since comparison operators always evaluate to false for it.
 */
function fieldsDiverged(): boolean {
  const ux = engine.ux();
  const uy = engine.uy();
  const rho = engine.rho();
  for (let i = 0; i < ux.length; i++) {
    const a = ux[i]!;
    const b = uy[i]!;
    const sp2 = a * a + b * b;
    if (!Number.isFinite(sp2) || sp2 > 25 || !Number.isFinite(rho[i]!)) return true;
  }
  return false;
}

/** On divergence (or an engine exception): auto-stop and offer a way to recover */
function handleDivergence(fromError: boolean): void {
  if (diverged) return;
  diverged = true;
  setRunning(false);
  console.warn("LBMFlow: detected divergence", { fromError });
  showToast(
    "The simulation diverged. Raise the viscosity ν or lower the flow speed, then reset.",
    "danger",
    {
      label: "↺ Reset",
      onClick: () => {
        try {
          resetSim(true);
        } catch (err) {
          console.error("LBMFlow: reset failed", err);
          showToast("Could not recover. Please reload the page.", "danger");
        }
      },
    },
    15000,
  );
}

// ------------------------------------------------------- Canvas sizing

/** Determine the canvas's actual size to fit inside the wrapper while preserving aspect ratio */
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

// ------------------------------------------------------------ Painting

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

/** Paint while interpolating from the previous position (no gaps even on fast drags) */
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
  dismissFirstHint(true);
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

// ------------------------------------------------------- First-time hint

const HINT_STORAGE_KEY = "lbmflow.first-hint-dismissed";

function safeStorageGet(key: string): string | null {
  try {
    return localStorage.getItem(key);
  } catch {
    return null;
  }
}

function safeStorageSet(key: string, value: string): void {
  try {
    localStorage.setItem(key, value);
  } catch {
    // It's fine if saving fails, e.g. in private browsing mode
  }
}

function dismissFirstHint(persist: boolean): void {
  if (!firstHint.hidden) firstHint.hidden = true;
  if (persist) safeStorageSet(HINT_STORAGE_KEY, "1");
}

if (safeStorageGet(HINT_STORAGE_KEY) === null) firstHint.hidden = false;
firstHintClose.addEventListener("click", () => dismissFirstHint(true));

// --------------------------------------------------- Scenario JSON export

function downloadJson(filename: string, data: unknown): void {
  const blob = new Blob([JSON.stringify(data, null, 2) + "\n"], {
    type: "application/json",
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

btnExport.addEventListener("click", () => {
  try {
    const scenario = buildScenario(
      engine,
      buildConfig(currentPreset),
      `${currentPreset.id}-gui`,
    );
    downloadJson(`lbmflow-${currentPreset.id}.json`, scenario);
    showToast("Scenario JSON saved. You can run it as-is with the CLI's lbm run.", "success");
  } catch (err) {
    console.error("LBMFlow: failed to export scenario", err);
    showToast("Failed to export the scenario.", "danger");
  }
});

// ------------------------------------------------------------- UI wiring

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

// Space to run / stop (disabled while a form element has focus)
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

// Auto-stop when the tab is hidden
document.addEventListener("visibilitychange", () => {
  if (document.hidden && running) setRunning(false);
});

// ------------------------------------------------------------ Main loop

function frame(): void {
  if (running) {
    try {
      engine.step(stepsPerFrame);
      spsSteps += stepsPerFrame;
      const now = performance.now();
      const dt = now - spsT0;
      if (dt >= 500) {
        const sps = (spsSteps * 1000) / dt;
        statusSps.textContent = Math.round(sps).toLocaleString("en-US");
        // MLUPS = lattice point updates per second (in millions)
        statusMlups.textContent = ((sps * engine.nx * engine.ny) / 1e6).toFixed(1);
        spsSteps = 0;
        spsT0 = now;
      }
      if (fieldsDiverged()) handleDivergence(false);
    } catch (err) {
      console.error("LBMFlow: engine step failed", err);
      handleDivergence(true);
    }
  }

  // Rendering has a fallback that clamps NaN to the minimum color, so the
  // field keeps being displayed even after divergence. Avoid a blank white
  // screen even while the engine is throwing.
  try {
    const range = renderer.render(engine, visMode, brushPreview);
    colorbarMin.textContent = formatRange(range.lo);
    colorbarMid.textContent = formatRange((range.lo + range.hi) / 2);
    colorbarMax.textContent = formatRange(range.hi);
    statusStep.textContent = engine.time.toLocaleString("en-US");
  } catch (err) {
    if (!diverged) {
      console.error("LBMFlow: render failed", err);
      handleDivergence(true);
    }
  }

  requestAnimationFrame(frame);
}

// ------------------------------------------------------------- Startup

drawColorbar(colorbar, visMode);
statusMode.textContent = VIS_MODE_LABEL[visMode];
visHint.textContent = VIS_MODE_HINT[visMode];
spfValue.textContent = spfSlider.value;
brushSizeValue.textContent = brushSize.value;
presetSelect.value = currentPreset.id;
applyPreset(currentPreset);
requestAnimationFrame(frame);
