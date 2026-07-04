# WASM ブリッジ設計（Phase 5）

GUI（web/, TypeScript）の `Engine` インターフェースに `lbm-core` を接続する層。

## クレート: crates/lbm-wasm

- `wasm-bindgen` + `lbm-core`（`default-features = false`、rayon 無効 =
  シングルスレッド。ブラウザの 1 フレーム内で回せる規模が対象）
- 精度は **f32 固定**（メモリ半減・WASM では十分。JS 側の Float32Array と無コピー整合）
- ビルド: `wasm-pack build crates/lbm-wasm --target web --release`
  → `web/src/engine/pkg/` に出力し、`WasmEngine implements Engine` アダプタで包む

## 公開 API（TS の Engine インターフェースと 1:1）

```rust
#[wasm_bindgen]
pub struct WasmSim { inner: Option<Simulation<f32>>, cfg: ... }

#[wasm_bindgen]
impl WasmSim {
    #[wasm_bindgen(constructor)]
    pub fn new() -> WasmSim;
    /// cfg_json: GUI の EngineConfig をそのまま JSON.stringify したもの。
    /// エラーは JsError（日本語メッセージ）で返す。
    pub fn init(&mut self, cfg_json: &str) -> Result<(), JsError>;
    pub fn step(&mut self, n: u32);
    pub fn nx(&self) -> u32;  pub fn ny(&self) -> u32;  pub fn time(&self) -> f64;
    /// フィールドはコピーせず wasm メモリのビューを返す（Float32Array::view 相当、
    /// ドキュメントに「次の step まで有効」の注意書き）
    pub fn rho_ptr(&self) -> *const f32;   // + len は nx*ny（JS 側でビュー生成）
    pub fn ux_ptr(&self) -> *const f32;
    pub fn uy_ptr(&self) -> *const f32;
    pub fn solid_ptr(&self) -> *const u8;
    pub fn set_solid(&mut self, x: u32, y: u32, solid: bool);
    pub fn set_inlet_profile_parabolic(&mut self, edge: &str, umax: f32);
}
```

### EngineConfig(JSON) → SimConfig 変換

- `collision: "bgk" | "trt"` → `Collision::Bgk | Trt{magic: 3/16}`
- edges の tagged union → `EdgeBC`（serde でデシリアライズ、lbm-wasm 内に定義）
- **この JSON 表現は Agent モードのシナリオ JSON（docs/AGENT_MODE_DESIGN.md）の
  `edges`/`physics` 節と同一形**にし、変換コードを共有できるようにする
  （共有クレート `lbm-scenario` に serde 型を置く — lbm-cli と lbm-wasm 両方が依存）

## setSolid の「消す」対応

`Simulation` は unset_solid を持たない（リム保護のため）。GUI の消しゴムは:
- lbm-wasm 側で「ユーザー描画レイヤ」(Vec<bool>) を別途保持
- 消去操作 = ユーザーレイヤ更新 → `init(cfg)` 相当の再構築 + ユーザーレイヤ再適用
  （数十ms、ペイント中はまとめて 1 回）
- または Phase 5 で lbm-core に `clear_solid(x,y)`（開境界・リム上は panic）を追加し、
  周辺セルの f を局所 feq で埋める。**こちらを採用予定**（流れを止めずに編集できる
  体験のほうが初学者に楽しいため）。VALIDATION に「clear_solid 後も質量が有限で
  NaN が出ない」ロバスト性テストを追加する。

## パフォーマンス目安

- 256×128 (32k cells) f32 シングルスレッド: 目標 ≥ 30 MLUPS → 60fps で
  ~15 step/frame。GUI 既定は 192×96 か 256×128 に設定。
- 大きい格子・3D・混相の重い計算はネイティブ（CLI/MCP）へ誘導する UI 文言を用意。
