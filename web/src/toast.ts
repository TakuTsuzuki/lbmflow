/**
 * 画面下部に出す最小トースト通知。
 * 同一メッセージの多重表示は抑止し、アクションボタン 1 個まで対応する。
 */

export type ToastKind = "info" | "success" | "danger";

export interface ToastAction {
  label: string;
  onClick: () => void;
}

let region: HTMLDivElement | null = null;
const activeMessages = new Set<string>();

function ensureRegion(): HTMLDivElement {
  if (!region) {
    region = document.createElement("div");
    region.id = "toast-region";
    region.setAttribute("role", "status");
    region.setAttribute("aria-live", "polite");
    document.body.appendChild(region);
  }
  return region;
}

export function showToast(
  message: string,
  kind: ToastKind = "info",
  action?: ToastAction,
  durationMs = 6500,
): void {
  if (activeMessages.has(message)) return;
  activeMessages.add(message);

  const el = document.createElement("div");
  el.className = `toast toast-${kind}`;

  const text = document.createElement("span");
  text.className = "toast-text";
  text.textContent = message;
  el.appendChild(text);

  let closed = false;
  const close = (): void => {
    if (closed) return;
    closed = true;
    activeMessages.delete(message);
    el.classList.add("toast-out");
    el.addEventListener("transitionend", () => el.remove(), { once: true });
    // transition が発火しない環境（reduced-motion 等）向けの保険
    window.setTimeout(() => el.remove(), 400);
  };

  if (action) {
    const btn = document.createElement("button");
    btn.className = "toast-action";
    btn.type = "button";
    btn.textContent = action.label;
    btn.addEventListener("click", () => {
      try {
        action.onClick();
      } finally {
        close();
      }
    });
    el.appendChild(btn);
  }

  const closeBtn = document.createElement("button");
  closeBtn.className = "toast-close";
  closeBtn.type = "button";
  closeBtn.setAttribute("aria-label", "通知を閉じる");
  closeBtn.textContent = "✕";
  closeBtn.addEventListener("click", close);
  el.appendChild(closeBtn);

  ensureRegion().appendChild(el);
  window.setTimeout(close, durationMs);
}
