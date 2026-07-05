# LBMFlow — 格子ボルツマン法流体シミュレータ

商用グレードを目指す LBM シミュレータ。Rust コア + TypeScript GUI + Agent モード。
**必読**: [docs/PLAN.md](docs/PLAN.md)（フェーズ計画・体制）、
[docs/VALIDATION.md](docs/VALIDATION.md)（検証テスト仕様 = 受入基準）。

## ビルド・テスト

```bash
cargo build --workspace --release
cargo test --workspace --release          # 通常スイート（LBM は debug だと ~50x 遅い。必ず --release）
cargo test --release -- --include-ignored # 重いベンチ含むフル検証（~5分）
# WASM（web GUI 用。lbm-wasm はワークスペース外）:
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg
#   （生成後 pkg/.gitignore を削除して pkg をコミットする運用）
cd web && npm run build                   # GUI（tsc strict + vite）
./target/release/lbm presets run cavity   # CLI スモーク
```

## 体制・規約

- Fable が PM。実装は Opus/Sonnet サブエージェント / codex に委任。
  **検証テストは codex or Opus/Sonnet が仕様（VALIDATION.md）から敵対的に作成**し、実装と分離する。
- codex 実行例: `codex exec --sandbox workspace-write --skip-git-repo-check "<task>" < /dev/null`
  （モデル gpt-5.5。**`< /dev/null` 必須** — stdin が pipe だと EOF 待ちで永久にスタックする。
  進捗は `~/.codex/sessions/<date>/rollout-*.jsonl` の更新で確認できる）
- コード・識別子・コミットメッセージは英語、ドキュメント・ユーザー向け文言は日本語。
- 物理仕様を変更したら docs/PHYSICS.md に理由と実験結果を記録する。
- フェーズ完了ごとに git commit。テストが red のままコミットしない（WIP は例外、メッセージに明記）。

## コア設計の約束事（壊すと検証が全滅する）

- 単一コアは `crates/lbm-core`（旧 lbm-core2 = V2 アーキテクチャ。V1 は 2026-07-05 引退、
  等価性凍結値はブランチ履歴の `tests/v1_match.rs` ヘッダ）。旧 V1 API は
  `lbm_core::compat`（公開ファサード）が提供し、scenario / CLI / wasm の 2D 経路が使う。
- D2Q9 の方向順序は lattice.rs（`Lattice` trait 実装）の定義が唯一の正。0:(0,0), 1:(1,0),
  2:(0,1), 3:(-1,0), 4:(0,-1), 5:(1,1), 6:(-1,1), 7:(-1,-1), 8:(1,-1)。
- f の配置は q-major SoA（fields.rs、halo パッド付き）: `f[q*plane + cell]`、
  cell = z·(nx·ny) + y·nx + x。GPU コアレッシング前提と同一。公開 API には出ない。
- 1 ステップ = 衝突 → halo 交換 → streaming → 開放端 BC → 境界線 moments 修正
  （CpuSimd は collide+stream+moments を step_band で融合）。パス構造・格納順を変える
  改修は `tests/backend_simd_equiv.rs` と T13（分割不変）のビット/閾値ゲートを
  通してから入れる。
- 壁エッジは 1 セルのソリッドリム。壁面は half-way（リム中心と流体中心の中間）。
- 速度モーメントは Guo forcing の F/2 補正込み（`sim.ux()` などは物理速度）。
- tau = 3*nu + 0.5（cs² = 1/3）。
