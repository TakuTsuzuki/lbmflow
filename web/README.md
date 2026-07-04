# LBMFlow Web GUI

格子ボルツマン法（LBM）流体シミュレータ LBMFlow のブラウザ GUI です。
Vite + TypeScript（vanilla、フレームワーク不使用・実行時依存ゼロ）で実装しています。

現在は **モックエンジン**（純 TS の解析的な流れ場生成器）で動作します。
将来、Rust 製の WASM エンジンを同じ interface に差し込む前提の構成です。

## 起動方法

```bash
cd web
npm install
npm run dev        # http://localhost:5173
```

本番ビルド:

```bash
npm run build      # tsc(strict) → vite build、成果物は web/dist/
npm run preview    # dist/ の動作確認
```

## 使い方

1. ヘッダーのプリセット（キャビティ流れ / 円柱まわりの流れ / チャネル流 / 自由キャンバス）を選ぶ
2. ▶実行 を押す（Space キーでも可）
3. キャンバスをドラッグすると障害物を描ける（右ドラッグ or「消す」モードで消去）
4. 右パネルで可視化する量（速さ / 渦度 / 密度）やパラメータを調整

タブが非表示になるとシミュレーションは自動停止します。

## ディレクトリ構成

```
web/
├── index.html            # UI の静的骨格（日本語ラベル）
├── src/
│   ├── main.ts           # アプリ配線・RAF ループ・障害物ペイント
│   ├── style.css         # ダークテーマ（CSS 変数、手書き）
│   ├── presets.ts        # プリセット定義（EngineConfig + 説明 + 初期障害物）
│   ├── colormap.ts       # viridis / RdBu の LUT（外部依存なし）
│   ├── render.ts         # スカラー化（|u|・渦度・密度）→ LUT 着色 → canvas 転写
│   └── engine/
│       ├── types.ts      # ★ エンジン抽象（wasm-bindgen 契約）
│       ├── index.ts      # ★ エンジン生成の差し替えポイント
│       └── mock.ts       # モックエンジン（解析的な流れ場生成器）
└── vite.config.ts
```

## エンジン差し替え設計

UI は `src/engine/types.ts` の `Engine` interface **のみ** に依存します。

```ts
export interface Engine {
  init(cfg: EngineConfig): void;
  step(n: number): void;
  readonly nx: number;
  readonly ny: number;
  readonly time: number;
  rho(): Float32Array;   // 長さ nx*ny、index = y*nx+x（y=0 が下端）
  ux(): Float32Array;
  uy(): Float32Array;
  solidMask(): Uint8Array;
  setSolid(x: number, y: number, solid: boolean): void;
}
```

WASM エンジンへの移行手順:

1. wasm-bindgen 側で上記シグネチャに対応するクラス（例 `WasmEngine`）を公開する
   - `rho()` などは WASM メモリ上のバッファを指す `Float32Array` ビュー、
     もしくはコピーを返す。呼び出し側は「次に `step()`/`init()` を呼ぶまで有効」
     という前提でしか保持しないので、ビュー返しで問題ない
2. `src/engine/index.ts` の `createEngine()` を `WasmEngine` を返すように書き換える
   - `.wasm` の非同期ロードが必要な場合は `createEngine(): Promise<Engine>` に変え、
     `main.ts` 冒頭の起動シーケンスで `await` する（変更点はこの 2 ファイルに閉じる)
3. `mock.ts` はデモ・フォールバック用として残してよい

### 座標系の約束

- `index = y * nx + x`、`y = 0` が **下端**（物理系の慣例）
- 描画時は `render.ts` が上下反転して canvas（上が y 最大）へ転写する

## モックエンジンの仕組み（`src/engine/mock.ts`）

本物の LBM は解かず、経過ステップ数 `t` の関数として場を解析的に合成します:

- 境界条件から基本流を選択（上壁 movingWall → キャビティ風の主渦、
  velocityInlet → 一様流 + 障害物下流の交互渦放出（カルマン渦列風）、
  外力 → ポアズイユ放物線分布、全周期 → 減衰テイラー・グリーン渦）
- 粘性 ν は渦の減衰率・渦列の振幅に反映（大きいほど早く静まる）
- `collision: "bgk"` では微小ノイズを付加（TRT が安定という演出。実物理ではない）
- 障害物セルは u=0・ρ=1、近傍セルは減速して壁らしく見せる

## 既知の制限

- モックエンジンの流れ場は見た目重視の合成場であり、物理的に正しくない
  （境界条件・ν・衝突演算子は「らしさ」の演出にのみ使われる）
- `pressureOutlet` の ρ 指定は現状モックでは未使用
- 解像度変更時、描いた障害物は最近傍サンプリングで引き継ぐため輪郭が粗くなる
