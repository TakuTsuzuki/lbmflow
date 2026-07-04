# LBMFlow — 格子ボルツマン法流体シミュレータ

商用グレードを目指す LBM シミュレータ。Rust コア + TypeScript GUI + Agent モード。
**必読**: [docs/PLAN.md](docs/PLAN.md)（フェーズ計画・体制）、
[docs/VALIDATION.md](docs/VALIDATION.md)（検証テスト仕様 = 受入基準）。

## ビルド・テスト

```bash
cargo build --workspace
cargo test --release                      # 通常スイート（LBM は debug だと ~50x 遅い。必ず --release）
cargo test --release -- --include-ignored # 重いベンチ含むフル検証
```

## 体制・規約

- Fable が PM。実装は Claude サブエージェント / codex に委任。
  **検証テストは codex が仕様（VALIDATION.md）から敵対的に作成**し、実装と分離する。
- codex 実行例: `codex exec --sandbox workspace-write --skip-git-repo-check "<task>" < /dev/null`
  （モデル gpt-5.5。**`< /dev/null` 必須** — stdin が pipe だと EOF 待ちで永久にスタックする。
  進捗は `~/.codex/sessions/<date>/rollout-*.jsonl` の更新で確認できる）
- コード・識別子・コミットメッセージは英語、ドキュメント・ユーザー向け文言は日本語。
- 物理仕様を変更したら docs/PHYSICS.md に理由と実験結果を記録する。
- フェーズ完了ごとに git commit。テストが red のままコミットしない（WIP は例外、メッセージに明記）。

## コア設計の約束事（壊すと検証が全滅する）

- D2Q9 の方向順序は lattice.rs の定義が唯一の正。0:(0,0), 1:(1,0), 2:(0,1), 3:(-1,0),
  4:(0,-1), 5:(1,1), 6:(-1,1), 7:(-1,-1), 8:(1,-1)。
- f の配置は cell-major: `f[(y*nx+x)*9 + q]`。
- 壁エッジは 1 セルのソリッドリム。壁面は half-way（リム中心と流体中心の中間）。
- 速度モーメントは Guo forcing の F/2 補正込み（`sim.ux()` などは物理速度）。
- tau = 3*nu + 0.5（cs² = 1/3）。
