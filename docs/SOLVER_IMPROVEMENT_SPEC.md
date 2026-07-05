# ソルバー改善仕様書（全体レビュー 2026-07-05）

> **main 取り込み注記（2026-07-05 PM）**: 本書はレビューブランチ
> `claude/amazing-mirzakhani-4060d3`（V1 引退前のコミット 84abaa3 基準）で書かれた。
> main は V1 引退済みのため、以下のパス読み替えで読むこと:
> `crates/lbm-core2` → `crates/lbm-core`（改名）、旧 V1 の `crates/lbm-core/src/*` →
> `crates/lbm-core/src/compat/*`（ファサードに集約。V1 単体は削除済み）。
> **A-1（S0）は main では解消済み**（sync-tests.sh は perl 化ののち V1 引退で削除、
> 複製スイートは compat 直参照の正規テストに昇格済み）。B-4 の「compat 移行」も
> V1 引退で実施済み。実験は `scripts/spec-experiments/`（パス翻訳済み）で再実行可能 —
> E2/E7 は改名後 main 上で仕様書の数値と一致することを確認済み。

> **ステータス: v1（実験検証済み・実施可能版）**。v0 の全主張を §3 の実験 E1〜E10 で
> 検証済み。**反証された項目はゼロ**、記述修正 2 件（A-3 の症状の様相、A-6 の数値）。
> **A-1（S0）は検証と同時に本ブランチで実施済み**（sync-tests.sh の perl 化＋
> 事後条件ガード、再生成スイート全 green 確認済み）。
> 実験コードは `scripts/spec-experiments/`（`cargo run --release e2` 等）で再実行可能。

- 対象コミット: 84abaa3（main 相当、ブランチ claude/amazing-mirzakhani-4060d3）
- 対象範囲: `crates/lbm-core`（V1、凍結参照実装）、`crates/lbm-core2`（V2:
  kernels / lattice / fields / solver / backend / halo / subdomain / dist / gpu / compat）、
  検証・ベンチ体制（tests / VALIDATION.md / スクリプト / CI）
- レビュー方法: 6 系統の独立レビュー（V2 物理カーネル / V2 構造 / GPU / MPI / V1 コア /
  V&V 体制）を並列実施し、PM が S0/S1 級所見を実コードで裏取り。行番号は対象コミット時点。
- 重大度: **S0** = 正しさの誤り（テストが偽の安心を与える等）/ **S1** = 高リスク
  （潜在バグ・スケール障害・重要な設計欠陥）/ **S2** = 改善機会 / **S3** = 軽微。
- 工数: S = 数時間 / M = 1 日 / L = 数日。

---

## 0. エグゼクティブサマリ

**物理カーネルの正しさに S0 はない。** D3Q19 Zou–He は Hecht & Harting (2010) との
文献照合・手計算検証まで実施して一致。Guo forcing の TRT 分解、偏差格納の定数折込み、
ハロー交換の十分性（unknown 集合との厳密一致）、GPU push 型融合の単一ライター保証も
すべて確認済み。**課題は「入口・構造・運用・検証体制」に集中している**：

1. **検証体制に S0 が 1 件**: sync-tests.sh の sed が BSD sed 非対応（`\b`）で無音失敗し、
   「V2 検証スイート」57 テストが実際には V1 を再実行している（compat 層は物理検証を
   一度も受けていない = R5 未実証）。修正は未マージの `v1-retirement` ブランチに存在。
2. **入口ガードの欠如**: V2 ネイティブ経路（3D/MPI の本経路）に構成検証層がなく、
   「未被覆面」「periodic×open 同軸」「ν=0」等が無症状で非物理計算になる。V1 側も
   NaN 速度がバリデーションを素通しする（実測で確認済み）。
3. **M-E の型レベルブロッカー**: `Solver`/`MpiSolver` の `Fields = SoaFields<T>` 束縛により
   GPU バックエンドがオーケストレータに載らず、ステップ列が GpuSolver に二重化。
   マルチ GPU / MPI+GPU は現構造では着手不能（構造レビューと GPU レビューが独立に同定）。
4. **物理の三重実装**: Shan–Chen 力が V1 / compat / V2 native の 3 箇所に並走し、
   機能マトリクスが軸ごとに非対称（接触角・二成分は 2D/CPU 専用のまま固着リスク）。
5. **スケール障害**: MPI セットアップがグローバル配列の全ランク複製（10⁹ 格子で
   wall_u 約 24 GB/ランク）、通信オーバーラップなし、rank-0 直列 gather のみ、
   チェックポイントなし。加えて `MPI_THREAD_SINGLE` 宣言のまま rayon が起動し得る。
6. **検証の穴**: CI 不在（全品質主張が手動スナップショット）、MPI 経路は cargo test で
   0 行実行、GPU は CPU 相対等価のみ（絶対物理検証ゼロ）、f32×3D は製品経路なのに無検証。

---

## 1. レビュー済み範囲と確認済み事項（問題なし）

以下は改修不要と確認された土台。改修時の「壊してはならないもの」リストでもある。

- **D2Q9/D3Q19 テーブル**: C/W/OPP/TRT ペア/FACE_UNKNOWNS が全て C からの const fn 導出で
  ドリフト構造的に不可能。0–4 次モーメント等方性 1e-15、V1 テーブルとのビットロック検査済み。
- **D3Q19 面 Zou–He**: 閉包 ρ(1−u·n) = S0+2S⁻、NEBB 再構成、接線補正 N_k の 6 項の符号を
  文献と照合し一致。質量・法線運動量の厳密充足を導出確認。
- **Guo forcing**: TRT 整合の対称/反対称分解（cp=1−ω⁺/2, cm=1−ω⁻/2）、moments と
  Reduction の F/2 補正、V1 と定義同一。
- **偏差格納の定数折込み**: half-way BB 不変、probe の +2w 物理化、convective 質量ピン止めの
  定数相殺、outflow の weight 中立性 — 全経路で矛盾なし。
- **ハロー交換**: pull で halo から読まれる方向 = unknown 集合と厳密一致（全数検査テスト有）。
  x→y→z 二相フォワードのコーナー充足、相内ハザードなし。
- **MPI プロトコル**: 全 Irecv 先行 post → Isend → wait でデッドロックフリー、
  ワイルドカード受信なし、タグ一意、コミュニケータ分離、POD 生バイト転送。
  T13-MPI（2D コーナー・3D 2×2×2・ψ 交換込み）で場ビット一致の実証あり。
- **GPU push 型融合**: 全 (q,セル) スロット単一ライター、parity 管理の全経路整合、
  WGSL と kernels.rs の項単位一致（結合順含む）。f64 はコンパイル時拒否で silent degrade なし。
- **V1 融合パス**: バンド境界の排他・unsafe 行アクセス契約・copy_span の周期 wrap /
  開放端スロット温存を全 cx ケースで確認。compat 層のコピーに実質的 drift なし
  （差分は doc/import/unused_mut のみ）。
- **テストの決定論性**: 乱数・環境依存なし、フレーク源なし。T13 はフィールドの
  `assert_eq!(d, 0.0)` ビット一致水準。T14 圧力 BC 緩和には 1-ulp 対照テスト付き。

---

## 2. 改善項目

### WP-A: 正しさ・入口ガード（即時、全項目とも既存テスト無修正 green が前提）

#### A-1 [S0] sync-tests.sh の置換修正と再発ガード
- 対象: `scripts/sync-tests.sh:42`、生成先 `crates/lbm-core2/tests/`（17 ファイル）
- 現状: BSD sed で `\b` が無効 → 無置換コピー。生成 17 ファイル全てが `use lbm_core::`
  のままで、compat 参照 0 件（PM 実測で確認）。修正コミット 622bbb2（perl 化）は
  未マージの `v1-retirement` にのみ存在。
- 改善: 622bbb2 相当（perl 置換）を main 系へ適用し、スクリプトに事後条件を追加
  （生成物に `use lbm_core::` が残存 or 置換件数 0 で exit 1）。生成ヘッダ付きファイルが
  `lbm_core2::compat` を import していることを検査する静的ガードテストを追加。
- 受入: `grep -rl "AUTO-GENERATED" crates/lbm-core2/tests | xargs grep -L "lbm_core2::compat"`
  が空。再生成スイートが compat 経由で全 green（→ 実験 E1 が事前検証）。
- 工数: S ／ 実験: **E1 = 確定・実施済み**（perl 化＋事後条件ガードを適用し 16 ファイル
  再生成 → `cargo test -p lbm-core2 --release` **全 suite green（exit 0）**。
  接触角凍結値・RT・open BC 系を含む物理スイートが compat 経由で初めて実走・合格 =
  **R5 の初実証**。残作業は静的ガードテストの追加のみ）

#### A-2 [S1] V1+compat の設定検証を NaN セーフ化
- 対象: `crates/lbm-core/src/domain.rs:306-310`、`crates/lbm-core2/src/compat/domain.rs`（同一コード）
- 現状: `if s > MAX_SPEED` は s=NaN で false（素通し）。NaN inlet は 3 ステップで場が NaN、
  NaN MovingWall は rim 速度選択の比較失敗で**静止壁に化けて無症状**（レビューで実測）。
  `Trt { magic }` と `force` も無検証。直下の rho 検査は `!(x > 0.0)` の NaN セーフ形式であり、
  安全イディオム自体はコードベースに既在。
- 改善: 速度検査を `if !(s <= MAX_SPEED)` に反転。`magic > 0`・`force`/`u` の `is_finite()` を
  `validate()` に追加。V1 と compat に同一パッチ（テキスト同一性維持）。
- 受入: NaN/inf の u・force・magic≤0 の `build()` が全て Err。合法設定はビット不変。
- 工数: S ／ 実験: **E6 = 確定**（NaN inlet: build()=Ok のまま 3 步で非有限 rho 42 セル。
  NaN MovingWall: build()=Ok、場は静止壁とビット一致 — 無音の静止壁化を厳密に実証）

#### A-3 [S1] Outflow/ConvectiveOutflow×固体隣接の無音質量リークを構成拒否
- 対象: `crates/lbm-core/src/sim.rs:765-774,667-669`、`crates/lbm-wasm/src/lib.rs:240-242`
- 現状: 開放端セルの内側隣接が固体だと BC がスキップされ（`solid[i] || solid[j]` で
  コード確認）、未知スロットが初期値のまま永久凍結する。GUI の塗り絵（外周 1 セルのみ
  拒否）から到達可能な本番経路。
- 実験結果（E5/E5b で機構を実証、ただし症状の様相は形状依存で v0 の記述を修正）:
  質量ドリフトの符号・規模は形状に依存し単独では判定指標にならない。決定的なのは
  **静止系での定常非物理速度**: 静止した箱＋右 Outflow＋ポケット（初期 rho=2.0）で
  2000 步後、バグ経路のエッジセルは **ux = −0.115** の巨大な定常速度を持ち続ける
  （対照 = プラグ 1 セル内側では ux = +0.00000 で完全静止、物理どおり）。
  全セル有限のままなので NaN 監視では捕捉不能 — v0 の主張どおり「無音」である。
- 改善: 凍結方針に沿い最小修正 = `set_solid`（V1 `sim.rs:792` / compat / wasm）に
  「開放端エッジセルの内側直隣への固体設置」を拒否する assert を追加。
  BC フォールバック実装は V2 の課題（B-8 の設計ノートに記載）とする。
- 受入: E5b 形状で panic（GUI は塗り拒否）。合法形状は probe_state_hash ビット不変。
  E5b を回帰テスト化（拒否確認）。
- 工数: M ／ 実験: **E5/E5b = 機構確定・記述修正済み**

#### A-4 [S1] V2 ネイティブ構成検証層 `GlobalSpec::validate` の新設
- 対象: `crates/lbm-core2/src/solver.rs:286-357`、`crates/lbm-core2/src/params.rs:24-35`、
  利用側の重複検証 `crates/lbm-scenario/src/lib.rs:582-654`
- 現状: `Solver::build` は次元・配列長 assert のみ。(1) **未被覆面**（periodic でも open でも
  壁 rim でもない面）で stale 値が毎ステップ実データとして混入し無症状で非物理化
  （`GlobalSpec::default()` のまま D3Q19 を使うだけで到達）。(2) ν=0 が
  `omegas()` で omega_m=0 となり非物理な緩和のまま進行。(3) periodic×open 同軸の
  二重適用。(4) MAX_SPEED / rho>0 / u_conv 範囲 / 異軸 open 面（V1 は
  `AdjacentOpenEdges` で拒否、V2 は素通し = 3D エッジで Zou–He の前提崩壊）も未検査。
  同等検証が compat と lbm-scenario に**二重実装**されている。
- 改善: `GlobalSpec::validate() -> Result<(), SpecError>` を core2 に新設し
  `Solver::build` 冒頭で強制。検査項目: ν>0 / 全非周期面の被覆（open BC または全面 solid rim）/
  periodic×open 排他 / 異軸 open 拒否 / MAX_SPEED（NaN セーフ）/ rho_bc>0 / u_conv∈(0,1] /
  2D で force[2]==0 / open 面軸 extent≥3。`set_inlet_profile` にも MAX_SPEED 検査。
  scenario の手書き検査は validate 呼び出し＋エラー変換に置換。
- 受入: 上記不正構成が全て Err になる単体テスト。既存 T13/T14/T15/v1_match 無修正 green。
  scenario の 3D 検証テスト無修正 green で重複検査コードが消える。
- 工数: M ／ 実験: **E2 = 確定**（D3Q19・z 面未被覆・z 一様初期条件で 100 步:
  nonfinite=0 のまま質量ドリフト 2.7e-3、z 不変性破れ 1.9e-4、偽 uz 2.6e-3。
  被覆対照は全指標 0.0 — 「NaN を出さず静かに非物理」を定量実証）、
  **E3 = 確定**（`omegas(nu=0)` → TRT (ω₊,ω₋)=(2,0)・BGK (2,2)。Solver は
  構築・10 步ともエラーなし）

#### A-5 [S1] ハロー交換スコープ誤用の構築時拒否
- 対象: `crates/lbm-core2/src/halo.rs:47-62,246-268`、`crates/lbm-core2/src/solver.rs:262-284`
- 現状: `Subdomain::neighbors` はグローバル part id だが `exchange_f_generic` はローカル
  index として解決。`new_local_part` + `LocalPeriodic`/`InProcess` の誤用は doc 注意書きのみで、
  隣接 id が 0 の場合は自己ラップとして**無音で誤った物理**になる（id≥1 なら OOB panic）。
- 改善: `HaloExchange` に `const SCOPE: ExchangeScope { Local, Remote }` を追加し、
  `Solver::build(only=Some(part))` で `SCOPE == Remote` を要求（不一致は構築時 panic）。
- 受入: 誤用構成が構築時エラーになる回帰テスト。T13 / T13-MPI 無修正 green。
- 工数: S ／ 実験: **E4 = 確定**（part=1 of [2,1,1] 周期 x + LocalPeriodic が panic せず
  実行され、正しい 2 パート InProcess 実行と比べ所有ブロックの rho が最大 7.7e-2 乖離
  — 無音の誤った物理を実証）

#### A-6 [S2] MovingWall 法線成分の拒否
- 対象: `crates/lbm-core/src/domain.rs:284-327`、`crates/lbm-core/src/sim.rs:1421-1425`
- 現状: half-way BB の運動量注入は接線壁速度のみ整合。法線成分は発散せず質量を
  無音で注入/流出し続ける（符号は向きに依存）。
- 改善: `validate()` でエッジ法線成分を持つ MovingWall を `InvalidParameter` 拒否
  （V1+compat 同時）。doc に理由を追記。
- 受入: 法線成分入りが Err。cavity 系既存テスト green・ビット不変。
- 工数: S ／ 実験: **E7 = 確定**（32×32 閉箱・500 步: 接線 u=[0.05,0] は
  ドリフト +1.1e-13（厳密保存）、法線 u=[0,−0.05] は質量 900→395.5、**−56.1%**。
  エラーなし — 発散もしないため気づけない）

#### A-7 [S2] `init_with` の入力検証（V1+compat）
- 対象: `crates/lbm-core/src/sim.rs:830-910`
- 現状: rho=0 で即 NaN（0×inf）、速度の MAX_SPEED 検査なし（`set_inlet_profile` は
  検査ありで API 内非対称）。GUI の Droplet 初期化は JSON 値を無検証で流し込む。
- 改善: `assert!(r > 0 && r.is_finite())` + MAX_SPEED 検査（座標入りメッセージ）、
  doc に Panics 節。compat の `init_with` には「クロージャは純粋であること
  （近傍で最大 5 回再評価される）」を明記。
- 受入: 不正クロージャで座標付き panic。probe_state_hash ビット不変。
- 工数: S

#### A-8 [S2] `zou_he_face_3d` の D3Q19 専用ガードと ConvectiveOutflow 契約テスト
- 対象: `crates/lbm-core2/src/kernels.rs:494-508`（unknown 5 本ハードコード）、
  `kernels.rs:589-594` + `fields.rs:114-118`（stale-slot 暗黙契約）
- 現状: (1) unknown 集合を 5 本と決め打ち。D3Q27 追加時（Q_MAX=27 で計画内）に
  コンパイル・実行とも通ったまま 4 本を放置し検出不能に誤る。(2) ConvectiveOutflow の
  記憶項は「streaming が unknown slot を書かない」という 4 モジュール横断の暗黙契約に
  依存（GPU は edge stash で独自再実装）。in-place streaming（M-E 候補）で確実に壊れる。
- 改善: (1) 関数冒頭に `assert_eq!(L::unknowns(face).len(), 5)` ガード（中期は面方向テーブルの
  const 化で `dir_index` の panic 経路ごと解消）。(2) `Backend::stream` の契約に
  「open 面 unknown slot 不変」を明文化し、CPU/GPU 両方で pre-stream 値とのビット一致を
  検査する契約テストを追加。
- 受入: ガード assert のユニットテスト。契約テストが両バックエンド green。
- 工数: S＋S

#### A-9 [S2] V2 実行時 NaN watchdog
- 対象: `crates/lbm-core2/src/solver.rs:877-891`（`local_nonfinite_count` は手動呼び出しのみ）
- 現状: 発散検出が CLI 2D / CLI 3D / MPI / GPU で 4 通りにばらけ、GPU 経路は手段なし。
  カーネルは V1 等価維持のため無ガード（これは正しい）。
- 改善: 既存の f64 集約 `local_mass_partials` の finite 検査（NaN は総和に伝播）を利用した
  `Solver::run_guarded(steps, check_every) -> Result<(), Diverged{step}>` を標準化。
  CLI/MPI ドライバはこれを呼ぶだけにする。GPU は同 API の readback 経路で暫定対応。
- 受入: 1 セル NaN 注入 → N 步以内に step 番号付きで検出。オーバーヘッド <1%（512²）。
- 工数: S〜M

#### A-10 [S3] 小粒整合バンドル
- (a) `crates/lbm-core/src/multiphase.rs:310` の unused_mut 除去（compat との diff ノイズ解消、
  バイナリ不変）。(b) `sim.rs:801-806` の誤解を招くコメント修正（多相は固体 rho を読まない）。
- (c) `t15_3d.rs:455-470` の ±25%/±15% 表記不整合の解消（VALIDATION.md 側の旧数値削除）。
- (d) kernels.rs 冒頭 doc の「bit-for-bit」主張を実態（pre-fusion V1 とビット一致で出発、
  現行は ≤1e-11 拘束・実測 ~1.6e-14/50steps）に更新。(e) MCMP `update_forces` に
  solid/周期一致の debug_assert 追加。(f) `equilibrium()` と collide 内 feq の
  ビット一致 property test。
- 工数: S（一括）

### WP-B: 構造改修（M-E 着手前に完了させる前提整備）

#### B-1 [S1] Backend `Fields` 一般化と GpuSolver 統合（M-E 最重要）
- 対象: `crates/lbm-core2/src/solver.rs:218,242-243`、`dist.rs:280-284`、
  `gpu/solver.rs:40-50,142-169`、`gpu/backend.rs:810-814`、`halo.rs:33-40`
- 現状: `Solver`/`MpiSolver` が `Fields = SoaFields<T>` に固定され、`WgpuBackend`
  （`Fields = GpuFields`）が載らない。GpuSolver がステップ列を複製（既知乖離あり）。
  GPU 側の `stream` は `CellRange::full` を assert し two_pass 分割を拒否。
  マルチ GPU / MPI+GPU は型レベルで合成不可能。
- 改善（段階発注）:
  1. `Backend` に `stage_in/stage_out`（ホスト⇔デバイス転写。`WgpuBackend::upload` が既実体）
     を正式化し、`Solver` は `SoaFields` をホストステージングとして保持、編集境界でのみ転写
     （GpuSolver の `host_dirty`/`device_ahead` 機構の一般化）。
  2. gather/診断を `read_moments`/`reduce` 経由に統一し、`Solver<D2Q9, f32, WgpuBackend,
     LocalPeriodic>` を成立させ、GpuSolver の独自 step 列を削除。
  3. fused カーネルに band ディスパッチ（uniform で y 範囲指定）を追加し
     `stream(range)` assert を撤廃（オーバーラップと将来のマルチ GPU の前提）。
  4. `HaloExchange` を pack/unpack バッファ境界で `Backend::Fields` ジェネリックに
     （GPU の edge stash と dist.rs の面プロトコルが雛形）。
- 受入: T14 が同一オーケストレータ経由で green。GpuSolver 削除。bench_gpu の MLUPS
  退行 ≤3%（同一測定手順）。
- 工数: L ／ 実験: **E10**（現行 submit 粒度の実測 = 改修時の性能基準線）

#### B-2 [S1] Backend 同期点契約の整理（probe / moments / end_step）
- 対象: `crates/lbm-core2/src/backend.rs:94-110`、`gpu/backend.rs:803-831,855-868`
- 現状: `stream` が probe 力を同期返却する契約に対し GPU はゼロを返し（Solver に載ると
  `probed_force()` が黙って 0 を返す仕込み）、`update_moments` は submit フックに意味流用。
- 改善: probe 力を `stream` の返り値から外し `read_probed_force`（明示 readback）に正式化。
  `end_step` フックを trait に追加し submit を分離。`update_moments` は lazy 契約として明文化。
  2 パス非対応は capability メソッドで表明（B-1 の 3. で解消するまでの経過措置）。
  V2 の probed_force は band 部分和の固定順 fold にし、ラン間ビット決定化
  （rayon reduce 非決定の解消、状態ハッシュへの編入を可能に）。
- 受入: ゼロ返却・意味流用の残存なし。probe 付き T14 ケース追加で CPU/GPU 同一 API 一致。
  同一スレッド数 2 ラン間で probed_force ビット一致。
- 工数: M

#### B-3 [S1] Shan–Chen 実装の一本化と V2 native 多相
- 対象: `crates/lbm-core2/src/solver.rs:679-748`、`compat/multiphase.rs:134-374`、
  V1 `multiphase.rs`（scenario が本番 import）
- 現状: SC 力ステンシルが 3 系統 5 ループ（PM 確認: V1:195,351 / compat:196,352 /
  solver.rs:736）。壁付着・仮想壁密度・二成分は compat/V1 形のみ（= 2D/CPU 限定）、
  MPI/3D は neutral 単相のみ、GPU は多相なし。
- 改善: `Solver::update_shan_chen_force` に壁項（g_wall・wall_rho、compat 版 347-365 と
  同一累積順）を吸収 → compat `ShanChen` を薄い委譲に置換 → `MultiComponent` を
  「2 つの Solver + exchange_scalar」として V2 native 移植。GPU 多相は M-F 検討
  （本仕様では対象外、B-9 ノートに記録）。
- 受入: SC ステンシルループが lbm-core2 内 1 箇所。validation_contact_angle /
  multiphase / rt 無修正 green。MPI 経路の wall_rho 付き接触角 T13 拡張 1 本追加。
- 工数: M〜L

#### B-4 [S2] 2D 本番経路の compat 移行と V1 の真の凍結
- 対象: `crates/lbm-scenario/src/lib.rs:7-8`、`crates/lbm-cli/src/runner.rs:6-7`、
  `crates/lbm-wasm/Cargo.toml:12`
- 現状: compat を import する本番コードが 0 件。2D シナリオ・CLI・wasm/GUI は V1 で動き、
  V1 は「凍結」といいつつロードベアリング。V1 レビューで wasm の使用 API 全数が compat で
  カバー済み・ゼロコピー `*_ptr` 互換・`default-features = false` ビルド可を確認済み。
- 改善: scenario / CLI / wasm の `lbm_core::` を `lbm_core2::compat::` に置換。V1 は
  dev-dependency（等価テスト専用）に降格。`v1-retirement` ブランチの完全削除
  （33e130a）は本移行の実績を積んでから別途判断。
- 受入: workspace 全 green、wasm-pack ビルド成功、本番 crate の `lbm_core::` 参照 0 件、
  GUI プリセットのフィールドハッシュが V1 版と一致（wasm smoke 1 本追加 = D-11）。
- 工数: M ／ 依存: A-1（compat スイートが本物になっていること）

#### B-5 [S2] リスタート/状態注入 API（スナップショット）
- 対象: `crates/lbm-core2/src/solver.rs:948-965`（gather のみで load がない）
- 改善: `Solver::snapshot() -> StateV1 { f[Q], solid, wall_u, force_field, time }` /
  `restore()` を対で実装（内部でハロー充填 + update_moments）。MpiSolver は rank0 経由
  gather/scatter で同 API（C-8 の分散チェックポイントの土台）。
- 受入: 「N 步 → snapshot → restore → M 步」=「N+M 步連続」が f64 ビット一致。
- 工数: M

#### B-6 [S2] per-cell 緩和率の下地（LES 前提、M-F 直前の三方改修を回避）
- 対象: `crates/lbm-core2/src/params.rs:76-105`、`kernels.rs:91-173`
- 改善: `SoaFields` に `omega_field: Option<Vec<T>>` を追加し `collide_row` が Some 時のみ
  per-cell omega を使う。None 経路はビット同一（probe_state_hash で担保）。GPU は
  storage buffer 1 本のフラグ制御（GPU-8 の limit 引き上げが前提）。LES 本体は M-F。
- 受入: None で全既存テストビット同一。一様値 = スカラー指定一致テスト。
- 工数: M

#### B-7 [S2] 公開裏口の封鎖と診断の f64 統一
- (a) `fields_mut` を pub(crate) 化し、必要操作は dirty 自動管理の専用メソッドに
  （MPI では片ランク編集ミスが無症状ハングになる最悪故障モードの元）。
- (b) facade に `total_mass_f64()` を追加（f32 での診断量子化 ~0.06/10⁶ セルの解消）。
- (c) MpiSolver::step 冒頭に debug 時のみの dirty フラグ一致 Allreduce（1 byte）で
  fail-fast（リリースはゼロコスト）。
- 受入: マスク編集経路の dirty 自動化、片ランク編集 debug テストがハングでなく assert 失敗。
- 工数: S〜M

#### B-8 [S2] カーネル拡張点設計ノート（実装なし、docs 1 枚）
- per-cell omega 受け渡し規約（B-6）、MRT/cumulant カーネルの置き場（CollisionKind 分岐と
  変換行列の所在）、曲面境界（Bouzidi）の per-link 壁距離 sparse 構造、in-place streaming 時の
  ConvectiveOutflow 代替（GPU edge-stash 方式の一般化）、Outflow×固体隣接の
  BC フォールバック（A-3 の恒久解）。各拡張の DoD に「既存構成ビット不変」を固定。
- 工数: M（レビュー承認まで）

### WP-C: スケール・運用（MPI / GPU）

#### C-1 [S1] MPI セットアップのローカル化（グローバル配列複製の解消）
- 対象: `crates/lbm-core2/src/dist.rs:305-312`、`solver.rs:298-333,76-125,628-654`
- 現状: `MpiSolver::new` がグローバル compact 配列（solid/wall_u）を全ランクで要求、
  `build_wall_rims` も全域生成、`set_solid` はグローバル全セル×全ランク呼び出し。
  10⁹ 格子で wall_u ≈24 GB/ランク → 弱スケーリング構成で確実に OOM。
  **ランク数でなく格子サイズで顕在化する構造的ブロッカー**（現行テストは n≤8・小格子
  ゆえ未検出）。
- 改善: クロージャ受け `MpiSolver::new_with(solid: impl Fn(x,y,z)->bool, …)`（init_with と
  同型のローカル評価）+ `build_wall_rims` のローカル版 + バッチ `set_solids_where(pred)`。
  既存 API は小規模用に残置。
- 受入: ランクあたりピーク RSS が O(N/P)+定数（実測）。T13-MPI 新 API 経由で場ビット一致。
- 工数: M

#### C-2 [S2] 交換オーバーラップ（post/finish 分割と two_pass 接続）
- 対象: `dist.rs:171-189`、`solver.rs:372-429`
- 改善: `HaloExchange` に `post_f()/finish_f()` を追加（InProcess は即時完了）。
  step を collide → post → interior stream（既存 two_pass の interior 範囲）→ finish →
  boundary shell に再配線。第一段は x 位相のみオーバーラップで可。
  **前提修正**: two_pass の boundary_shells は幅 1 軸でシェルが重複し probe が二重計上される
  （`solver.rs:972-1007`、場は冪等で無害・T13 では不可視）。シェルを互いに素に修正して
  から接続する。
- 受入: 幅 1 軸で two_pass on/off の probed_force 一致（→ 実験 E8 が現状の二重計上を実証）。
  T13-MPI 全 PASS 維持。exchange 待ち占有率の計測可能な減少。
- 工数: L ／ 実験: **E8 = 確定**（[64,1,1]・障害物 1 セル・probe・20 步:
  on/off の probed_force 比 = **2.000**（厳密二重計上）、total_mass は両者一致
  — 場は無傷で probe だけ壊れる = T13 の場比較では検出不能、の主張どおり）

#### C-3 [S2] 並列 I/O（per-rank raw + manifest、rank0 全域バッファ排除）
- 対象: `dist.rs:504-574`
- 改善: 短期 = 各ランクが自ブロックを個別ファイル書き、rank0 は manifest のみ。
  中期 = MPI-IO subarray（rsmpi の File サポート要確認）。gather_* は検証用と明記して残す。
- 受入: 出力時に rank0 ピーク RSS が O(N/P)。出力時間がランク数に対し非増加。
- 工数: M〜L

#### C-4 [S2] probe Allreduce の遅延化
- 対象: `dist.rs:369-379`（probe 有効時に毎ステップ 3-double Allreduce）
- 改善: `probed_force()` を collective 化し `time` キーでキャッシュ（照会時のみ縮約）。
- 受入: probe 付きベンチで per-step collective が profile から消える。mpi_t13 PASS 維持。
- 工数: S

#### C-5 [S2] 交換バッファの永続化（GPU-aware MPI への布石）
- 対象: `dist.rs:125-126,175-181,199-229`、`solver.rs:691-711`（ψ plane 毎ステップ確保）
- 改善: `MpiExchange` に面×種別の送受バッファを保持し再利用。ψ plane / staging も
  フィールド化。バッファ所有を MpiExchange に寄せる（GPUDirect への差し替え点を一箇所に）。
- 受入: 定常ステップ中のヒープ確保ゼロ（計測）。bench_mpi 非退行。T13-MPI ビット一致維持。
- 工数: S〜M

#### C-6 [S2] ランク間 spec 整合検査
- 対象: `dist.rs:305-333`
- 現状: nu・faces・マスク内容の不一致は**メッセージ長が一致するため検出されず**、
  縫い目で不連続な「もっともらしい」場を出す（ジョブスクリプト事故クラス）。
- 改善: 構築時に spec 正規化バイト列 + マスク FNV ハッシュの Allreduce(min/max) 比較で
  不一致を項目名付き abort。
- 受入: nu だけ変えた 2 ランク注入テストが即座に明示エラー。正常系コスト測定不能。
- 工数: S

#### C-7 [S2] MPI スレッドレベルの Funneled 化
- 対象: `examples/mpi_t13.rs:391`、`examples/bench_mpi.rs:27`、`backend.rs:36,133-147`
- 現状: rsmpi 0.8.1 の `initialize()` = `Threading::Single`（PM がレジストリソースで確認）。
  一方 default feature `parallel` は 16,384 セル以上で rayon 起動 → 実サイズで
  MPI_THREAD_SINGLE 宣言下のマルチスレッド実行（MPI 規格違反。UCX/OFI 系で破損・
  ハングし得る。現テストは ≤6,144 セル/rank でたまたま直列）。
- 改善: `dist::init_mpi() -> Universe`（Funneled 要求、provided 不足＋parallel 有効なら
  明示エラー）を追加し、examples/ガイドを移行。
- 受入: 2 ランク × rayon 強制（parallel_min_cells 引き下げ）で T13-MPI 相当 PASS、
  provided ≥ Funneled をログ確認。
- 工数: S ／ 実験: **E9 = 確定（ソースレベル、§3 参照）**

#### C-8 [S2] 分散チェックポイント/リスタート
- B-5 の上に、collective な `MpiSolver::save(dir)/load(world, dir, backend)`（per-rank raw +
  rank0 manifest、spec ハッシュ・decomp 一致検証）。偏差格納 f の raw 保存で再開ビット一致。
- 受入: 「50 步 → save → load → 50 步」=「100 步連続」が場ビット一致。manifest 不一致は明示エラー。
- 工数: M ／ 依存: B-5, C-6

#### C-9 [S1] GPU submit チャンクの時間校正と device-lost の Result 化
- 対象: `gpu/backend.rs:187-188,226`（`submit_chunk: 200` 固定）、`:94-98,301-314`（expect panic）
- 現状: 200 步分（最大 ~1000 dispatch）を単一 submit。機構上は時間無制限 →
  低速 GPU では Windows TDR（既定 2 s）超過で device removed → **プロセス panic**。
  復旧経路なし。
- 改善: 初回チャンク実測から 1 submit ≈ 100–250 ms 目標に自動校正（上限 200 維持）。
  `wait_idle`/`map_staging` を `Result<_, GpuError>` 化し device lost を伝播。
  `set_device_lost_callback` で理由捕捉。
- 受入: bench_gpu の MLUPS 退行 ≤3%。校正ロジック単体テスト。poll 失敗で Err が返るテスト。
- 工数: M ／ 実験: **E10 = 確定**（本機 M5 Max/Metal 実測: 2048² TGV 5,719 MLUPS →
  200 步チャンク = **147 ms/submit**。最速級の consumer GPU でこの値なので、
  ~15 倍遅い GPU（数百 MLUPS 級 iGPU）で同格子が TDR 2 s を超える外挿が成立。
  参考: 512²=11,509 / 1024²=6,607 MLUPS、proto 比 −5.3〜−18.0%（±20% 受入線内））

#### C-10 [S2] GPU リソース上限の事前検証
- 対象: `gpu/backend.rs:650-692,66-91,236-241`
- 改善: `alloc` 冒頭で必要バイト数 vs `device.limits()`、`Q*n ≤ u32::MAX`（D3Q19 は
  2.26 億セルで溢れる）、dispatch 数 ≤65,535 を検査し、日本語の理由付き Err。
- 受入: 上限超過格子で `GpuSolver::new` が明示エラー。T14 無変更 green。
- 工数: S

#### C-11 [S2] GPU 診断経路の効率化
- 対象: `gpu/backend.rs:318-343,870-911`、`gpu/solver.rs:184-207`
- 改善: `f_cache` の Arc 化（クローン排除）、FluidCells の readback 不要化（host_solid で
  完結）、sync の 3 readback を 1 エンコーダ/1 wait に統合。GPU 側 2 段 reduction は
  M-E の高速モードとして追加（ホスト f64 経路は T14 用に維持）。
- 受入: 2048² で sync+診断 3 連が ≥3 倍高速。T14 診断値ビット同一。
- 工数: S（＋M）

#### C-12 [S2] FP16 配管（M-E 本体の前提）
- 対象: `gpu/backend.rs:80`（Features::empty 固定）、`gpu/wgsl.rs:216`、要素サイズ `*4` 散在
- 改善: `SHADER_F16` の条件要求、`generate::<L>(cfg: KernelCfg { storage: F32|F16 })` 化
  （変更点を「バッファ宣言 + load/store ラッパ」2 箇所に封じ込め、演算は f32 維持）、
  `GpuFields` に要素サイズ保持。非対応アダプタは明示エラー（silent fallback しない）。
- 受入: T16 新設（f16 格納の劣化を凍結帯で定量）。2048² TGV で f32 版比 MLUPS ≥1.5 倍。
- 工数: M ／ 依存: B-1（オーケストレータ統合後が手戻り最小）

#### C-13 [S2] GPU の契約層接続（scenario `backend: "gpu" | "auto"`）
- 対象: `crates/lbm-scenario/src/lib.rs:575-577`、`crates/lbm-cli/src/main.rs:190`
- 改善: feature gpu ビルドで解禁、能力チェック（f64/3D/未対応 BC は理由付き reject）を
  validate に実装、`auto` は実測閾値（例 n≥256²）で選択し結果をログ明示。
- 受入: cavity/cylinder プリセットが `backend:"gpu"` で完走し CPU との場の差が T14 許容内。
  `f64`+`gpu` が日本語理由付き拒否。
- 工数: M ／ 依存: B-1
- 公開ベンチ（M-E）は再現手順の公開が本体 — この項がその前提。

#### C-14 [S2] GPU メモリフットプリント削減
- 対象: `gpu/backend.rs:677-678,691,133`
- 改善: force_field/wall_u 不在時のダミーバッファ化、staging の遅延確保・right-sizing。
  現状は 2048² TGV で +284 MB（約 1.9 倍）— 3D 化で 8–16 GB 級 GPU の最大格子を直撃。
- 受入: 力場・固体なしで確保量 ≤ 2×f + moments + mask + O(perimeter)。T14 green。
- 工数: S

#### C-15 [S3] GPU 小粒バンドル
- (a) `max_storage_buffers_per_shader_stage` 等の limit 引き上げ（現状 step カーネルが
  既定上限 8 本ちょうど — 1 本追加で実行時 panic する地雷）。(b) naga による生成 WGSL の
  parse+validate 単体テスト（GPU 不要、golden file 付き）。(c) BcParams の Rust index ⇔
  WGSL フィールドの単一テーブル生成 or 照合テスト。(d) `GpuContext::new` の
  `Result<_, GpuInitError>` 化（adapter_info 付き。auto フォールバックの診断可能化）。
  (e) probe 力 f32 CAS 加算の非決定性を GPU_EVALUATION.md に 1 行明記。
- 工数: S×5

#### C-16 [S3] MPI 小粒バンドル
- (a) `choose_decomp`（表面積最小の自動分割）+ mpi_t13 の任意 n 対応（n=3,5,6・
  割り切れない dims をケース追加）。(b) bench_mpi の 3D/D3Q19・strong/weak モード化
  （現状 2D strip 固定では R3 本測定に使えない）。(c) 決定論的総和 `total_mass_deterministic()`
  （グローバル行単位部分和の固定順合成）のオプション追加。
- 工数: S＋S＋M

### WP-D: 検証・プロセス

#### D-1 [S1] CI の 3 段整備
- 現状 `.github/` 不在（PM 確認）。全品質主張が手動スナップショット。
- 改善: (1) push/PR で `cargo test --workspace --release`、(2) 夜間 `--include-ignored`、
  (3) GPU/MPI はセルフホストランナー（本機 M5 Max + ~/.local の Open MPI）でタグ実行。
  リモート未定の間は pre-merge ローカルフック + 実行ログの `docs/CI_LOG.md` 追記で開始。
- 受入: マージに default スイート green が機械的に強制される。週次 full+gpu+mpi 記録。
- 工数: M

#### D-2 [S1] MPI ロジックの cargo テスト化
- 対象: `dist.rs`（`#[test]` 0 件 — PM 確認）
- 改善: pack/unpack・phase plan を純関数として単体テスト（InProcess と同一バッファ内容、
  許容 ==0.0）。`test_mpi.sh` を夜間ジョブ化。1 ランク self-exchange smoke を cargo 圏内に。
- 受入: mpirun なしで dist.rs 主要ロジックがテストされる。T13-MPI 週次 PASS ログ。
- 工数: M

#### D-3 [S1] GPU の絶対物理検証と skip 機構
- 対象: `tests/t14_backend_equiv.rs`（CPU 相対等価 8 本のみ、adapter 無しで expect panic）
- 現状: CPU と GPU が**同じ向きに壊れる**バグ（スペック解釈共有ミス）は相対等価で検出不能。
- 改善: GPU 直の絶対テスト 2 本（TGV 収束次数 ≥1.7、キャビティ Ghia RMS ≤0.02U — f32 実測で
  校正して凍結）。adapter 不在は skip（`LBM_REQUIRE_GPU=1` で fail 昇格）。
  3D GPU・GPU 多相の非対応は VALIDATION.md の既知の制限節に明文化。
- 受入: `--features gpu` で絶対 2 本 green。GPU 無しホストで skip 終了。
- 工数: M

#### D-4 [S1] f32×3D の検証追加
- 対象: `crates/lbm-scenario/src/lib.rs:490-492`（`Sim3Handle::F32` = 製品経路）、
  `tests/t15_3d.rs`（f32 出現 0 — PM 確認）
- 改善: t15-1（z 不変 2D 退化、f32 相対 ≤1e-5）と t15-4（TGV3D 減衰率 ±2%）の f32 版 +
  質量ドリフト ≤1e-5/10³step。実測値を VALIDATION.md T15 に f32 行として凍結。
- 工数: S

#### D-5 [S1] 等価検証の地平延長と native 絶対検証
- 現状: V2 の絶対物理検証は「V1 と 500 步 ≤1e-11 → V1 は検証済み」の推移律のみ
  （Ghia 定常は ~99k 步で地平の 200 倍）。V2 native API（3D/scenario 経路）は連鎖の外。
- 改善: A-1 完了後、(1) 長時間等価 1 本（キャビティ Re=100 を 20k 步、≤1e-9、#[ignore] 可）、
  (2) native `Solver` 直の TGV 収束次数テスト 1 本（facade とコアのズレも検出可能に）。
- 工数: S ／ 依存: A-1

#### D-6 [S2] 受入基準原典の整合回復
- 対象: `docs/COMPETITIVE_SPEC.md:57-59`（R1: 球 ±5% 等 / R3: 弱スケ ≥85%）vs 実装
  （球 ±10%・D_h 正規化 / R3 は n=8 で 73.2% → コミットで「n≤4 局所」に縮小）
- 改善: R1/R3 を改定履歴付きで更新（±10% と D_h 定義の根拠リンク、n≤4 局所線と
  クラスタ条件未達の明記）。PLAN の「R1 達成」表記を「（3D キャビティ除く・±10% 改定版）」に
  訂正。T15.5（3D キャビティ Ghia 表）にバックログ位置を付与。
- 受入: 全「達成」宣言が原典現行版と 1:1 対応。
- 工数: S

#### D-7 [S2] VALIDATION.md への T13/T14 節追加
- T14 の「6 構成 1e-5 / 圧力 1e-4 + 1-ulp 対照」、T13 の「場 ==0.0 / 診断 1e-12」を
  仕様書へ一元化（codex 発注可能な形に）。テストヘッダは仕様参照に薄くする。
- 工数: S

#### D-8 [S2] codex 敵対テスト order #7（T14/T15 への追随）
- 現状: 敵対発注は T13（order #6）まで。GPU 等価・3D 物理は実装側自作のみ。
- 改善: D-7 の仕様改定から発注: T14 攻撃（境界面直上の初期不連続、probe が面に接する、
  u→MAX_SPEED 近傍）、T15 攻撃（z 退化を壊す摂動、極端アスペクト比、球オフセンター）。
- 受入: 1 巡 + triage 記録が TESTING_NOTES に残る。
- 工数: M（発注・triage 込み）

#### D-9 [S2] 性能回帰検出
- 改善: core2 版 bench_mlups（V2 CPU の性能主張が V1 Phase 9 数値の引用のままの解消）+
  `--check` モード（凍結値 ±25% 逸脱 fail）+ 夜間実行で JSON 履歴を `docs/bench_history/` に
  追記。probe_state_hash も夜間で強制。
- 受入: 25% 退行が 24h 以内に fail 顕在化。
- 工数: M ／ 実験: **E10** が GPU 側の基準値を提供

#### D-10 [S2] `lbm verify` と LIMITATIONS.md
- 改善: 検証サブセット（数分級）を実行し受入帯との比較表を出力する `lbm verify`。
  既知の制限（GPU=2D/f32、CP/リスタートなし、長時間検証上限、τ→0.5 ガイドライン等）を
  `docs/LIMITATIONS.md` に一元化しリリースと対で出す。
- 工数: L（verify）+ S（LIMITATIONS）

#### D-11 [S2] wasm smoke テスト
- `wasm-pack test --node` で TGV 100 步の質量保存 + native f64 一致 ≤1e-12。
  web 側は書き出し JSON のスキーマ round-trip のみ最小限。
- 工数: M ／ 依存: B-4 と同時実施が効率的

#### D-12 [S2] CpuScalar 性能ギャップの明示（M-E への引き継ぎ）
- 現状: V2 CPU はフェーズ分離 3 パス + per-cell 分岐で、V1 融合カーネル比の単核ギャップが
  未計測（弱スケーリング 97-99% の分母が甘い可能性）。
- 改善: M-E の CpuSimd（計画済み）で V1 step_band の q-major 移植。それまで公開ベンチには
  V1 単核比を併記。受入: CpuSimd が CpuScalar とビット一致かつ単核 ≥ V1 の 0.9x。
- 工数: L（M-E 本体）

---

## 3. 検証実験の結果（2026-07-05 実施、本機 M5 Max。再実行: `scripts/spec-experiments/`）

**10 実験すべて実施済み。反証ゼロ、記述修正 2 件**（E5→A-3 の症状の様相、
E7→A-6 の数値と符号）。実験ログの生出力は各実験の RESULT 行として下表に転記。

| ID | 検証対象 | 実測結果 | 判定 |
|---|---|---|---|
| E1 | A-1/B-4/D-5 | perl 化＋ガード適用→16 ファイル再生成→`cargo test -p lbm-core2 --release` **全 green（exit 0）** | **確定・実施済み**。R5 初実証。B-4 のブロッカー解消 |
| E2 | A-4 | 未被覆 z 面: nonfinite=0 のまま質量ドリフト 2.7e-3、z 不変性破れ 1.9e-4、偽 uz 2.6e-3（被覆対照は全指標 0.0） | **確定**（無音の非物理を定量実証） |
| E3 | A-4 | `omegas(0)` = TRT (2, 0)・BGK (2, 2)。Solver 構築・10 步エラーなし | **確定** |
| E4 | A-5 | part=1+[2,1,1]+LocalPeriodic: panic なし、正解比 rho 最大乖離 **7.7e-2** | **確定**（無音の誤物理） |
| E5/E5b | A-3 | 質量ドリフトは形状依存で判定不適（v0 記述を修正）。決定打: 静止箱でポケットのエッジセルが 2000 步後も **ux=−0.115** の定常速度を保持（対照は +0.00000） | **機構確定・記述修正** |
| E6 | A-2 | NaN inlet: build()=Ok→3 步で非有限 42 セル。NaN MovingWall: build()=Ok→静止壁と**ビット一致** | **確定** |
| E7 | A-6 | 接線: ドリフト +1.1e-13（厳密保存）。法線 u=[0,−0.05]: 質量 900→395.5（**−56.1%**、エラーなし）。v0 の「+115%」は向き依存の符号違いと判明 | **確定・数値修正** |
| E8 | C-2 | [64,1,1] 幅 1 軸: two_pass on/off の probed_force 比 = **2.000**、total_mass は一致（場は無傷） | **確定** |
| E9 | C-7 | rsmpi 0.8.1 `environment.rs:268-270` で `initialize()` = `initialize_with_threading(Threading::Single)` を直接確認。PARALLEL_MIN_CELLS=16,384 に対し実運用格子（512²/rank=262,144 セル）で rayon 起動は算術的に確定。実行時デモは省略（共有メモリ BTL では症状が出にくく、判定に不要） | **確定（ソースレベル）** |
| E10 | C-9/B-1/D-9 | 2048² TGV **5,719 MLUPS** → 200 步チャンク = **147 ms/submit**（512²=11,509 / 1024²=6,607。proto 比 −5.3〜−18.0%、±20% 線内） | **確定**（TDR 外挿成立、校正初期値の根拠） |

補助（コード検査で確定済み・実験不要）: `Fields=SoaFields` 束縛（B-1）、GPU probe ゼロ返却
（B-2）、SC 5 ループ（B-3）、`MpiSolver::new` のグローバル配列シグネチャと 24 GB/ランク算術
（C-1）、`submit_chunk: 200`（C-9）、dist.rs `#[test]` 0 件（D-2）、`.github` 不在（D-1）、
t15 の f32 出現 0（D-4）、rsmpi `initialize()`→`Threading::Single`（C-7、レジストリソース確認）。

---

## 4. 実施順序（確定）

実験依存はすべて解消済み。各項目は独立に発注可能な粒度で記述してある。

- **実施済み（本ブランチ）**: A-1 の本体（スクリプト修正＋ガード＋再生成＋全 green 確認）。
  残: 静的ガードテスト 1 本（R-Phase 1 に含める）。
- **R-Phase 1（即時・~2 日）**: A-2〜A-10 + A-1 残作業 + D-6/D-7（文書整合）。
  すべて S〜M。既存テスト無修正 green・合法設定ビット不変が共通 DoD。
- **R-Phase 2（M-E 前提整備・~1.5 週）**: B-1 → B-2 →（並行）B-3, B-5〜B-8 /
  C-9〜C-11 / D-1〜D-5。B-4（compat 移行）は E1 で実証済みのためいつでも着手可
  （D-11 の wasm smoke と同時実施を推奨）。
- **R-Phase 3（スケール・M-E と並走可）**: C-1, C-4〜C-7, C-16 →（C-2, C-3, C-8）/
  C-12〜C-15 は B-1 後 / D-8〜D-10。
- M-E 本体（FP16 実装・マルチ GPU・公開ベンチ・CpuSimd = D-12）は本仕様の
  B-1/B-2/C-9/C-12/C-13/D-9 を前提とする。

体制: 実装は Opus/Sonnet サブエージェント / codex へ WP 単位で発注。検証テスト
（D-3/D-4/D-8 と各受入テスト）は codex が VALIDATION.md 改定版（D-7）から敵対的に作成し、
実装と分離する（従来プロトコル維持）。
