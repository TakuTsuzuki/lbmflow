# MPI 分散ガイド（M-D, 2026-07-05）

lbm-core の feature `mpi` は HaloExchange の MPI 実装（`dist::MpiExchange`）と
1ランク=1サブドメインのドライバ（`dist::MpiSolver`）を提供する。
設計は docs/ARCHITECTURE_V2.md §2.3 / docs/HPC_SCALING.md 段階計画 3 に対応。

## 現状スコープ（正直な現在地）

- **検証済み**: 単一ノード内マルチランク（Open MPI 5.0.9 / arm64 macOS、
  共有メモリ BTL 経由）。T13-MPI で分散実行 ≡ 単一ランク実行を確認済み
  （場はビット一致、診断は f64 再結合差のみ。下記）。
- **未対応**:
  - GPU バックエンドとの併用（`--features mpi,gpu` はビルドは通るが
    `MpiSolver` は `CpuScalar` 系 `SoaFields` バックエンド前提。
    device-resident ハロー転送・GPUDirect は M-E 以降）。
  - 通信と計算のオーバーラップ（`exchange_f` は各軸位相で同期完結。
    two-pass ストリーミングの縫い目は Solver に既設なので、
    Isend 発行→内部計算→wait→境界計算への差し替えが M-E 候補）。
  - 並列 I/O（rank-0 gather 経由の全場復元まで。VTK 並列/HDF5 は未着手）。
  - マルチノード実測（クラスタアクセス待ち。§クラスタ測定リスト参照）。

## ビルド

rsmpi（crate `mpi` 0.8）はビルド時に `mpicc` を探す。**arm64 ネイティブの
MPI が PATH の先頭に必要**（/usr/local に x86_64 版 MPI が居る環境では
PATH 順で事故る。`file $(which mpirun)` で arm64 を確認すること）。

```bash
# ソースビルドした Open MPI（例: $HOME/.local/openmpi）を使う
export PATH=$HOME/.local/openmpi/bin:$PATH
file $(which mpirun)   # → Mach-O 64-bit executable arm64 を確認

cargo build -p lbm-core --release --features mpi
```

Open MPI のソースビルド手順（参考。5.0.9 / arm64、Fortran 無効で ~15分）:

```bash
mkdir -p ~/.local/src && cd ~/.local/src
curl -sL https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.9.tar.bz2 | tar xj
cd openmpi-5.0.9
./configure --prefix=$HOME/.local/openmpi CC=clang CXX=clang++ --disable-mpi-fortran
make -j$(sysctl -n hw.ncpu) && make install
```

既定ビルド（feature 無し）は rsmpi に一切依存しない。`cargo test --workspace`
は MPI 環境なしで従来どおり通る。

## 実行

```bash
# T13-MPI 検証（-n 1,2,4: 2D 4ケース / -n 8: 3D TGV 2x2x2。非ゼロ exit で失敗）
./scripts/test_mpi.sh

# 弱スケーリング（ランクあたり 512^2、ranks {1,2,4,8}、表出力）
./scripts/bench_mpi.sh          # LOCAL=512 STEPS=200 RANKS="1 2 4 8" で調整可
```

API 最小例（1ランク=1パート。デカルト分割はランク数と一致させる）:

```rust
use lbm_core::dist::MpiSolver;
use lbm_core::prelude::*;

let universe = mpi::initialize().unwrap();
let world = universe.world();
let spec = GlobalSpec::<f64> { dims: [1024, 512, 1], ..Default::default() };
let mut s: MpiSolver<D2Q9, f64, CpuScalar> =
    MpiSolver::new(&world, &spec, &[], &[], [world.size() as usize, 1, 1],
                   CpuScalar::default());
s.init_with(|x, y, _| (1.0, [0.0, 0.0, 0.0]));
s.run(1000);                       // step/診断/gather は全て collective
let mass = s.total_mass();         // Allreduce（全ランクで同値）
let rho = s.gather_rho();          // rank 0 のみ Some(全体場)
drop(s);                           // 複製コミュニケータを finalize 前に解放
```

**collective 契約**: `step` / `init_with` / `update_shan_chen_force` /
診断（`total_mass` / `total_momentum` / `nonfinite_count`）/ `gather_*` /
マスク編集は全ランクが同じ順序で呼ぶこと。`set_solid` は**全ランクが同じ
座標列で呼ぶ**（所有ランクが格納、他ランクはハロー再交換の予約だけ行う）。
`MpiSolver` は複製コミュニケータを保持するため、`mpi::initialize()` の
Universe を drop（= MPI_Finalize）する**前に** solver を drop すること。

## 交換プロトコル（実装メモ）

- InProcess と**同一の x → y → z 位相・同一の pack/unpack**（`halo.rs` の
  共有ヘルパを両実装が呼ぶ）。コーナー/エッジは面隣接経由の前送で、位相ごとに
  先行軸のハロー込み層を転送（MPI でも面リンク 6 本のみ）。
- 各軸位相で両側 2 面の層を Irecv → Isend 発行 → 全完了待ち → unpack。
  受信面 `F` 宛メッセージのタグは `base + F.index()`（f: 100, ψ: 200,
  マスク: 300/400, gather: 500）。周期軸で decomp=1 の自己ラップは MPI を
  介さずローカルコピー。
- 転送はスカラ型の生バイト（f64/f32 とも可逆）なので、分割実行の場は
  単一ランク実行と**ビット一致**する。診断のみ rank 部分和 → Allreduce の
  f64 再結合差を許容（T13 流儀: atol + rtol、1e-11）。
- ランク配置は `solver::partition` のデカルト分割そのまま
  （part id = `(pz·dy+py)·dx+px` = rank）。MPI_Cart は不使用。

## T13-MPI 実測（2026-07-05, M5 Max / Open MPI 5.0.9）

| ケース | -n | decomp | 場 max\|Δ\| | 診断 max\|Δ\| |
|---|---|---|---|---|
| 2D TGV 96×64 | 1/2/4 | 1×1 / 2×1 / 2×2 | **0.0**（ビット一致） | ≤ 3.3e-14 |
| キャビティ 64×64（蓋が縫い目跨ぎ） | 1/2/4 | 同上 | **0.0** | ≤ 2.3e-14 |
| 円柱+力プローブ（縫い目上）+ 放物線流入 | 1/2/4 | 同上 | **0.0** | ≤ 9.1e-13 |
| Shan-Chen 液滴（ψ を exchange_scalar、2×2 コーナー） | 1/2/4 | 同上 | **0.0** | ≤ 4.5e-11* |
| 3D TGV 24³ (D3Q19) | 8 | 2×2×2 | **0.0** | ≤ 4.6e-15 |

\* 液滴の診断差は total_mass ≈ 1.5e3 に対する再結合差（相対 ~3e-14）。
合否は `atol + rtol·|ref|`（両 1e-11）で判定し全 PASS。

## 弱スケーリング（単一ノード実測、2026-07-05）

ランクあたり 512²（D2Q9 f64 TGV、ランク内は**直列**バックエンド、
decomp [n,1,1]、200 step 計測・20 step ウォームアップ）:

| ranks | time | MLUPS 合計 | MLUPS/rank | 効率 |
|---|---|---|---|---|
| 1 | 1.304 s | 40.2 | 40.2 | 100% |
| 2 | 1.313 s | 79.9 | 40.0 | 99.4% |
| 4 | 1.346 s | 155.9 | 39.0 | 97.0% |
| 8 | 1.781 s | 235.5 | 29.4 | 73.2% |

**読み方（重要）**: この測定は Open MPI の**共有メモリ経由**であり、
インターコネクトの実測ではない。さらに測定機（M5 Max）は 6 Super + 12
Performance の異種コア構成で、対照実験（通信ゼロの独立 1 ランクジョブ ×8
並走）でも 33.7 MLUPS/rank（= 84% 相当）までしか出ない。つまり n=8 の
73.2% の内訳は「ハード（異種コア+帯域競合）の天井 84%」×「MPI 化による
残り ~87.5%（ロックステップで最遅ランクに同期するジッタ結合が主、
通信量自体は ~50 KB/step/rank で無視できる）」。均質コア内に収まる
n≤4 は R3 のローカル合格線 ≥85% を満たす（97-99%）。
**真の弱スケーリングはクラスタ実測が必要**（下記）。

## クラスタでやるべき測定リスト（R3 完了条件）

1. **弱スケーリング本測定**: ランクあたり 3D 128³（D3Q19）を 1→64 ランク
   （ノード内→複数ノード）で。R3 合格線: 64 ランクで効率 ≥80%。
   2D 512² 版も比較用に（本ガイドの単一ノード表と接続する）。
2. **強スケーリング**: 固定 1024³ を 8→512 ランクで（表面積/体積比の劣化点）。
3. **通信/計算比の実測**: MPI_T プロファイル or mpiP で exchange_f の
   占有率。>10% なら two-pass オーバーラップ（M-E）を前倒し。
4. **ノード間 BTL/MTL 確認**: UCX/OFI の選択、eager/rendezvous 閾値
   （層メッセージは ~10-200 KB 帯）とタグマッチング競合の有無。
5. **ランク×スレッドのハイブリッド最適点**: 本ガイドの測定はランク内直列。
   ノードあたり「ランク数 × rayon スレッド数」の格子を振る
   （`CpuScalar::parallel_min_cells` 閾値はそのまま使える）。
6. **プロセス配置/バインド**: `--map-by`/`--bind-to`（macOS では不可、
   Linux クラスタで必須）。NUMA ノード跨ぎの層 pack/unpack 帯域も確認。
7. **正しさの再確認**: scripts/test_mpi.sh をクラスタの MPI 実装
  （Open MPI 以外に MPICH/Cray MPICH）でも全 PASS させる
  （rsmpi は両対応。ビット一致要件は実装非依存のはず）。
8. **大規模での診断コスト**: Allreduce（診断）と rank-0 gather（出力）の
   スケール限界。出力は並列 I/O（HDF5/ADIOS2 系）への移行判断材料を取る。

## 既知の罠（今回踏んだもの）

- **x86_64 MPI との同居**: /usr/local の Homebrew (x86_64) 版 mpicc が
  PATH 先頭にあると rsmpi のプローブが x86_64 フラグを拾いリンクに失敗するか、
  Rosetta 経由の mpirun で起動が壊れる。常に arm64 版を PATH 先頭に。
- **MPI_Finalize 後の MPI_Comm_free**: `MpiSolver`（と `MpiExchange`）は
  複製コミュニケータを Drop で解放する。`mpi::initialize()` の Universe より
  **先に** drop されるスコープ設計にすること（examples 参照）。
- **マスク編集の collective 性**: `set_solid` を所有ランクだけで呼ぶと
  ハロー再交換（collective）の呼び出し回数がランク間でずれてデッドロックする。
  `MpiSolver::set_solid` は非所有ランクでも dirty マークだけ立てる設計。
  `Solver` を直接使う場合は `mark_masks_dirty()` を全ランクで呼ぶこと。
- **プローブの Allreduce タイミング**: probed_force は step 毎の Allreduce。
  プローブ未設定時は省略される（ベンチに余計な collective を入れない）。

## run_guarded (R-Phase 1, A-9)

`MpiSolver::run_guarded(steps, check_every)` is a **collective** call: every rank
must invoke it with the same arguments. Each check is a 2-double Allreduce over the
mass partials (NaN propagates through the sum), and the divergence branch —
`Err(Diverged{step})` — is taken uniformly on all ranks, so there is no
divergence-induced deadlock. Overhead at 512²/rank with check_every=100 is <0.5%.
The trajectory is bit-identical to plain `run`.
