# 公開ベンチ比較ドラフト（M-E: 性能タイトル向け）

**ステータス**: ドラフト（Web公開情報の収集結果、2026-07-05 取得）
**ルール**: 出典URLのない数字は載せない。見つからなかった定量は「公開定量なし/未確認」と明記する。
自陣の数字は共有負荷下の暫定値（アイドル再測定前）。**対外公開は §8 チェックリスト完了後**。

対象は COMPETITIVE_SPEC.md §1 の競合: FluidX3D / M-Star CFD / waLBerla / Palabos / OpenLB。

---

## 1. サマリ（読み方の注意つき）

- **単GPU 3D (D3Q19) の公開王者は FluidX3D**。A100 PCIe 40GB で FP32格納 8,526 MLUPS、
  FP16S格納 16,035 MLUPS。ハイエンド（MI300X 41,327 / H100 NVL 32,922、いずれもFP16S）は桁が違う。
- **うちは 3D GPU が未実装**なのでこの土俵にまだ乗れない。2D D2Q9 GPU の 5,857〜11,365 MLUPS は
  データ移動量が D3Q19 の約半分の別競技であり、**3D 表に混ぜてはいけない**（§6）。
- **CPU 3D は健闘**: M5 Max 18C の 260 MLUPS は、公開値のある 64〜128 コア級サーバCPUの
  実測（204〜330 MLUPS）と同じ桁（条件差あり、§4）。
- **M-Star は公開テキストとしての MLUPS 定量なし**（チャートのみ、§3）。
- **マルチノードは waLBerla / OpenLB の独壇場**（兆セル・TLUPS級、§5)。うちは M-D 未達で実績ゼロ。

---

## 2. FluidX3D 公開ベンチ（単GPU、D3Q19）

出典: [FluidX3D GitHub README](https://github.com/ProjectPhysX/FluidX3D)
（生ファイル: [raw README.md](https://raw.githubusercontent.com/ProjectPhysX/FluidX3D/master/README.md)、2026-07-05 取得）

測定条件（README 記載の原文要旨）: 「D3Q19 SRT、拡張機能なし（暗黙 mid-grid bounce-back 境界のみの
純LBM）、空の立方体ボックス、十分なサイズ（典型 256³）」。演算は常に FP32、格納精度を
FP32 / FP16S / FP16C で切替。README は算術強度 2.37（FP32/FP32）、5.27（FP32/FP16S）、
16.56（FP32/FP16C）FLOPs/Byte を明記し「性能はメモリ帯域のみで制限される」と述べる。
メモリは Esoteric-Pull + FP16 圧縮で 55 Bytes/cell（従来 FP64 LBM は約 344 Bytes/cell）。

| GPU | FP32/FP32 | FP32/FP16S | FP32/FP16C | 出典 |
|---|---:|---:|---:|---|
| AMD MI300X | 22,867 | **41,327** | 31,670 | [README](https://github.com/ProjectPhysX/FluidX3D) |
| NVIDIA H100 NVL | 20,303 | **32,922** | 18,424 | 同上 |
| NVIDIA H100 SXM5 | 17,602 | **29,561** | 20,227 | 同上 |
| NVIDIA RTX 5090 | 9,522 | 18,459 | **19,141** | 同上 |
| NVIDIA A100 PCIe 80GB | 9,657 | **17,896** | 10,817 | 同上 |
| NVIDIA A100 PCIe 40GB | 8,526 | **16,035** | 11,088 | 同上 |
| NVIDIA RTX 4090 | 5,624 | 11,091 | **11,496** | 同上 |
| NVIDIA RTX 3090 | 5,418 | **10,732** | 10,215 | 同上 |
| Apple M2 Ultra (76-CU) | 4,629 | **8,769** | 7,972 | 同上 |
| Apple M3 Ultra (60-CU) | 4,438 | **8,174** | 8,086 | 同上 |
| Apple M1 Ultra (64-CU) | （最速モード値 8,418） | | | 同上 |
| Apple M2 Max (38-CU) | 2,405 | **4,641** | 2,444 | 同上 |
| Apple M1 Max (24-CU) | 2,369 | **4,496** | 2,777 | 同上 |
| Apple M5 (10-CU) | （最速モード値 1,613） | | | 同上 |

補足（README 由来、2026-07-05 時点）:
- A100 の最速モード値は変種で異なる: SXM4 80GB 18,448 / PCIe 80GB 17,896 / PCIe 40GB 16,035 / SXM4 40GB 16,013。
- **Apple M4 世代・M5 Max/Pro/Ultra のエントリは表に存在しない**（= うちの M5 Max と同一機種の
  公開比較値は現状ない。最も近い公開値は M2 Max / M2・M3 Ultra）。
- README に専用のマルチGPUベンチ表は見当たらず（今回取得分）。マルチGPU対応自体はあり、マルチノード
  （MPI）は非対応。ライセンスは非商用向け無償。

## 3. M-Star CFD（本命比較対象）: 公開定量の状況

**結論: テキストとして公表された MLUPS 数値は見つからなかった（＝公開定量なし扱い）。
ただし公式ドキュメントにスケーリングチャート（SVG）と条件記載はある。**

確認できた事実（すべて公式ソース）:

| 項目 | 内容 | 出典 |
|---|---|---|
| スケーリングベンチ | v3.3.123 の結果。ケース2種: 撹拌槽（Rushton 1枚+粒子、格子 1M〜512M 点）/ バッフル付きパイプ（静的形状、2.6M〜970M 点）。1/2/4/8 GPU の「シミュレーション平均 MLUPS」をチャートで公表 | [Scaling Performance](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html) |
| 測定プラットフォーム | AWS p3.8xlarge（8× V100 SXM2 16GB, CUDA 11.5, peer access 制限）/ GCE a2-highgpu-8g（8× A100 SXM4 40GB, CUDA 11.2, full peer access） | 同上 |
| チャートの軸レンジ | 8×A100 構成のチャート縦軸目盛は撹拌槽ケースで最大 16,000 MLUPS、パイプケースで最大 25,000 MLUPS（データ点の数値ラベルはなく、正確な値は図の目視読み取りが必要） | [agitated チャートSVG](https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_agitated.svg) / [pipe チャートSVG](https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_pipe.svg) |
| 再実行手段 | ベンチマークパッケージは v3.3.140+ で利用可能（顧客が自環境で再実行できる） | [Scaling Performance](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html) |
| サイジング経験則 | 「GPU あたり 30〜60M 格子点以上が目安」「GPU RAM 1GB ≈ 2〜4M 格子点 + 1M 粒子」 | 同上 / [Hardware Guide](https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html) |
| ハード要件 | NVIDIA GPU 前提（推奨: GeForce 40/50 系〜RTX 6000 Ada / RTX PRO 6000〜H100/B100 SXM + NVLINK/NVSWITCH）。「帯域はメモリ容量・演算性能の次に比較すべき仕様」 | [Hardware Guide](https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html) |
| マーケ上の主張 | 「millions and billions of lattice grid points で動くよう構築」「詳細で正確なプロセスシミュレーションを数分で」 | [mstarcfd.com/software](https://mstarcfd.com/software/) |

注意: 今回 SVG からのデータ点自動抽出も試みたが、座標読み取りが信頼できない
（軸最大値を超える読み値が出る等）ため**本ドラフトにはデータ点数値を載せない**。
正確な値が必要なら §8 のとおり図の目視読み取りで別途確定させる。
比較の含意: M-Star のベンチは「撹拌槽+粒子」等の**フル物理ケース**であり、FluidX3D の
「空箱カーネル」数字と直接比較してはならない（§6.5）。

## 4. クロスコード比較表: 同一デバイス（NVIDIA A100）での D3Q19/Q27 単GPU

各コードの公開値を A100 に揃えた表。**それでもケース・格子・精度・ストリーミング実装が
異なるため「同条件」ではない**。条件列を必ず読むこと。

| コード | A100 変種 | 格子/ケース | 精度(演算/格納) | MLUPS | 種別 | 出典 |
|---|---|---|---|---:|---|---|
| FluidX3D | PCIe 40GB | D3Q19 SRT・空箱 典型256³ | FP32/FP32 | 8,526 | 実測(公表表) | [README](https://github.com/ProjectPhysX/FluidX3D) |
| FluidX3D | PCIe 40GB | 同上 | FP32/FP16S | 16,035 | 実測(公表表) | 同上 |
| Palabos (GPU port, C++ stdpar) | SXM4 40GB | D3Q19 BGK・Taylor-Green、L=590 | FP32 | 理論ピーク 9,481 の 75〜85% ≈ **7,100〜8,060**（換算） | 論文記載%からの換算 | [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1) |
| Palabos (GPU port) | SXM4 40GB | D3Q19 BGK、L=480 | FP64 | 理論ピーク 4,921 の同効率帯 | 論文記載 | 同上 |
| waLBerla (waLBerla-wind) | JUWELS Booster | **D3Q27** cumulant・風車つきフルソルバ | FP32 | 1,677（ルーフライン上限 7,513 の 22.3%） | 実測(論文) | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.5 | 4× A100 ノード | D3Q19 BGK・1000³ キャビティ・Periodic Shift | FP32 | ノードあたり 24,800 → **GPUあたり ≈6,200**（換算） | 実測(公式)+換算 | [OpenLB 1.5 release](https://www.openlb.net/news/openlb-release-1-5-available-for-download/) |
| **LBMFlow（うち）** | —（A100 実測環境なし） | — | — | **未測定**（3D GPU 未実装 + NVIDIA 実機なし） | — | COMPETITIVE_SPEC.md §5 |

参考（他GPUのクロスコード点）: STLBM（Palabos 系研究コード）は D3Q19 **FP64** キャビティ N=128 で
GTX 1080 Ti ≈820 / RTX 2080 Ti ≈1,100 / V100 PCIe ≈2,300 MLUPS（AA-pattern、論文図由来の近似値）
— [PLOS ONE 10.1371/journal.pone.0250306](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306)。
また Palabos GPU 論文は関連研究として「waLBerla CUDA バックエンドは A100-SXM4 40GB で理論ピーク約85%」
と記述 — [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1)。

## 5. CPU（単ノード）3D 比較

| コード | CPU | 格子/ケース | 精度 | MLUPS | 出典 |
|---|---|---|---|---:|---|
| **LBMFlow（うち）** | Apple M5 Max 18C | D3Q19 | f32 | **260**（共有負荷下・暫定） | 本リポジトリ実測（PERFORMANCE.md 系） |
| STLBM | AMD EPYC 64コア | D3Q19 キャビティ N=128 | FP64 | ≈300（AA-pattern, SoA） | [PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306) |
| STLBM | Intel Xeon 48コア | 同上 | FP64 | ≈330（swap-AoS） | 同上 |
| waLBerla-wind | AMD EPYC 7763 128コア(1ノード) | **D3Q27** cumulant（風車なし） | FP32 | 204（ルーフライン上限 461） | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.3 (2018) | Magnus 32,784コア | D3Q19 | 原文参照 | 総計 142,479（≈4.3 MLUPS/コア） | [openlb.net/performance](https://www.openlb.net/performance/) |

正直な読み: 18 コアのラップトップ SoC で 64〜128 コア級サーバ CPU の公開値と同じ桁に居るのは
強い材料。ただし STLBM は FP64（帯域約2倍消費 → MLUPS 半減相当）、waLBerla-wind は
D3Q27（データ移動 ≈27/19 倍）なので、**精度・ステンシルを補正すると「同等」ではなく
「フェアに見て互角圏、条件次第」**。断言は正式測定後にする。

## 6. 2D（D2Q9）: うちの GPU 数字の置き場所

| コード | デバイス | 格子 | 精度 | MLUPS | 出典 |
|---|---|---|---|---:|---|
| **LBMFlow（うち）** | M5 Max GPU (Metal, wgpu) | D2Q9 | f32 | **5,857〜11,365**（共有負荷下・暫定、条件により幅） | 本リポジトリ実測 |
| **LBMFlow（うち）** | M5 Max 18C CPU | D2Q9 | f32 | **1,183**（同上） | 本リポジトリ実測 |

**競合の公開 2D 数字は今回見つからなかった**: FluidX3D は 3D 専用、M-Star は 3D 製品、
waLBerla/OpenLB/Palabos の代表公開ベンチも 3D。つまり 2D 表は現状「単独走」であり、
対外的な速度主張の主戦場にはならない。D2Q9 は 1 更新あたりのデータ移動が D3Q19 の約半分
（9 vs 19 分布関数）なので、**D2Q9 の MLUPS を D3Q19 の表に並べると約2倍下駄を履く**。
対外資料では必ず分離する。

## 7. マルチノード / スケーリングの公開実績

| コード | マシン / 規模 | 実績 | 出典 |
|---|---|---|---|
| waLBerla | JUQUEEN（BG/Q）458,752 コア / 1.8M スレッド | 1兆セル超、最大 1.93 兆セル更新/s（アブストラクト記載値。ACM ページ直接取得不可のため要目視確認）。SuperMUC 32,768 コアで強スケーリング実証 | [SC13 DOI:10.1145/2503210.2503273](https://dl.acm.org/doi/10.1145/2503210.2503273) |
| waLBerla | JUQUEEN | 「最大シミュレーションは1兆セル超」「40万コア超への良好なスケーラビリティ」（本文から直接引用確認済み） | [arXiv:1511.07261 (ar5iv)](https://ar5iv.labs.arxiv.org/html/1511.07261) |
| waLBerla-wind | JUWELS Booster 30ノード/120 A100 | 弱スケーリングで「GPUあたり性能はほぼ一定」（17.5M セル/GPU、平均 74.46 steps/s ≈ GPUあたり約1,300 MLUPS 換算） | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.5 | HoreKa 128ノード/512 A100 | 総計 ≈1.33 TLUPS（D3Q19 FP32 キャビティ）。64→128 GPU 強スケーリング効率 0.64〜0.81（格子 575³〜2300³）。LES付き乱流ノズル実ケースでベンチ性能の 92%（224 GPU） | [openlb.net/performance](https://www.openlb.net/performance/) |
| OpenLB 1.5 | HoreKa 2ノード/8 A100 | 1000³ FP32 キャビティ 42.2 GLUPS（1ノード4×A100 は 24.8 GLUPS、CPU 2ノード AVX-512 は 2.7 GLUPS → GPU スピードアップ 15.6×） | [OpenLB 1.5 release](https://www.openlb.net/news/openlb-release-1-5-available-for-download/) |
| OpenLB 1.9 | Aurora 1,000ノード（系の約10%） | ピーク 21,120 GLUPS、4兆セル（D3Q19 FP32） | [openlb.net/performance](https://www.openlb.net/performance/) |
| OpenLB (2026) | HoreKa 異種混成（CPU/GPU 3パーティション） | 最大18Gセル、強スケーリング効率 0.66〜0.91（区間別）、単GPUノードで〜1e9 セル可 | [arXiv:2506.21804](https://arxiv.org/html/2506.21804v1) |
| Palabos GPU | DGX 4× A100 40GB | 弱スケーリング理想の 80〜90%、強スケーリング 65〜80% | [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1) |
| M-Star | 8× V100 / 8× A100（クラウド単ノード） | 1/2/4/8 GPU スケーリングチャート公表（数値は図、§3） | [Scaling Performance](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html) |
| **LBMFlow（うち）** | — | **実績なし**（M-D で MPI 予定。単ノード内マルチランク弱スケーリング ≥85% が受入基準） | COMPETITIVE_SPEC.md §3 R3 |

## 8. 条件差の注意書き（比較表を読む前に必須）

1. **MLUPS の定義は共通だが測り方が違う**。定義は「百万格子点更新/秒」
   （waLBerla の論文は「MLUP/s per core = 1コアが1秒に更新するセル数」と定義
   — [ar5iv:1511.07261](https://ar5iv.labs.arxiv.org/html/1511.07261)）。ただし
   何を1回の「更新」に含めるか（境界処理・出力・通信）はコードごとに異なる。
2. **ステンシルで割り引く**: D3Q27 は D3Q19 より1更新のデータ移動が大きく（27 vs 19 DDF）、
   帯域律速では MLUPS が構造的に低く出る。waLBerla-wind（D3Q27 cumulant）の 1,677 MLUPS を
   FluidX3D（D3Q19 SRT）の 8,526 と並べて「waLBerla は遅い」と読むのは誤り。
   2D D2Q9 はさらに別枠（§6）。
3. **精度で割り引く**: 帯域律速 LBM では格納精度が半分になると MLUPS はほぼ倍
   （FluidX3D: RTX 4090 で FP32 5,624 → FP16S 11,091 — [README](https://github.com/ProjectPhysX/FluidX3D)。
   Palabos: FP32 理論ピーク 9,481 vs FP64 4,921 GLUPS — [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1)）。
   精度面は Lehmann らが「FP64 と FP32 の精度差はほぼ全ケースで無視できる」「多くのケースで
   16-bit でも十分」と報告 — [arXiv:2112.08926 / Phys. Rev. E 106, 015308](https://arxiv.org/abs/2112.08926)。
   比較表には必ず「演算精度/格納精度」を併記する。
4. **格子サイズ依存**: 小さい格子は性能が出ない。FluidX3D は「十分なサイズ（典型 256³）」の
   空箱で測る。Palabos GPU 論文は「性能はメッシュ解像度の増加関数」と明記
   （[arXiv:2506.09242](https://arxiv.org/html/2506.09242v1)）。M-Star も「GPU あたり 30〜60M
   格子点以上」をスケーリング効率の目安にする（[docs](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html)）。
   → うちの正式測定も格子サイズスイープで飽和値と到達サイズを両方出す。
5. **カーネル単体 vs フルソルバ**: 同じコード・同じ GPU でも、空箱カーネルとフル物理では
   4〜5 倍違う（waLBerla-wind: ルーフライン 7,513 → 風車つき実測 1,677 MLUPS = 22.3%
   — [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1)。OpenLB: 実ケース乱流ノズルで
   ベンチ性能の 92% という好例もある — [openlb.net/performance](https://www.openlb.net/performance/)）。
   FluidX3D の表は前者、M-Star のチャートは後者（撹拌槽+粒子）。**この2つを直接比較しない**。
6. **平均区間・ウォームアップ**: M-Star は「シミュレーション全体の平均 MLUPS」を公表
   （[docs](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html)）。
   初期化・JIT・キャッシュ効果を含むかで数字が動く。うちの正式測定はウォームアップ除外を明記する。
7. **実装最適化の幅そのものが大きい**: 同一ハードでもデータレイアウト・ストリーミング実装で
   性能は大きく変わる（LBM ベンチカーネル研究の趣旨 — [arXiv:1711.11468](https://arxiv.org/abs/1711.11468)。
   実例: 素朴 211 → 最適化 550 MLUPS — [multiphase code-gen 論文](https://journals.sagepub.com/doi/full/10.1177/10943420211016525)）。
   「コード対コード」の差はハード差・条件差と切り分けて主張する。
8. **ベンダー公表帯域とルーフラインの扱い**: FluidX3D は算術強度を公開しルーフラインで
   説明する（[README](https://github.com/ProjectPhysX/FluidX3D)）。うちも「実測帯域→理論上限→実測 MLUPS」
   の3点セットで公表すると検証文化と整合する。

## 9. うちが正直に負けている点（現時点）

1. **3D GPU が存在しない**。競合の主戦場（単GPU D3Q19）で比較可能な数字を出せていない。
   R2 の受入基準（単GPU D3Q19 f32 ≥1,500 MLUPS）を満たして初めて表に乗れる。
2. **FP16 格納モード未実装**。FluidX3D は FP16S で FP32 比ほぼ2倍（例: 4090 11,091 vs 5,624）を
   実証済み（[README](https://github.com/ProjectPhysX/FluidX3D)）。柱4の実装待ち。
3. **ハイエンド NVIDIA/AMD の数字に届く手段がない**。H100 NVL 32,922 / MI300X 41,327 MLUPS
   （FP16S、[README](https://github.com/ProjectPhysX/FluidX3D)）級は、仮に 3D GPU 実装が完璧でも
   Apple Silicon の帯域では物理的に不可能。CUDA/HIP バックエンドと実機アクセスが必要（SPEC §5）。
4. **マルチノード実績ゼロ**。waLBerla は兆セル・40万コア超（[ar5iv:1511.07261](https://ar5iv.labs.arxiv.org/html/1511.07261)）、
   OpenLB は 512 GPU で 1.33 TLUPS / Aurora 4兆セル（[openlb.net](https://www.openlb.net/performance/)）。
   うちの R3 目標（64ランク弱スケーリング ≥80%）はこれらの何桁も下の初期目標にすぎない。
5. **フル物理のベンチがない**。M-Star のチャートは撹拌槽+粒子という「売り物のワークロード」で
   測っている。うちは空箱系カーネルの数字しかなく、LES・移動境界・スカラー輸送を積んだときの
   性能低下率（waLBerla-wind の例では -78%）を まだ知らない。
6. **自陣の数字自体が暫定**。共有負荷下の測定であり、公開したら検証文化（全主張実測紐づけ）に
   自分で違反する。§10 完了までは対外に出さない。

## 10. うちだけの点（比較表の脚注ではなく本文で主張する差別化）

- **検証スイート同梱・ワンコマンド再実行**（敵対的検証 56+ テスト、Ghia/Schäfer-Turek/RT/
  等変性 4e-16）。競合の検証は「論文・公開ベンチ集」形式が主で、顧客環境での再実行可能性を
  製品仕様にしているのはうちだけ（M-Star はスケーリング用ベンチパッケージを v3.3.140+ で
  提供しており「性能の再実行」は部分的に可能 — [docs](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html)。
  **物理精度の検証再実行**を同梱する点がうちの差別化）。
- **エージェントネイティブ**（JSON Schema 自己記述 + MCP。M-Star の Python API は人間向け、
  FluidX3D は C++ セットアップ編集）。
- **ポータビリティ**: wgpu (Metal/Vulkan/DX12) + WASM。M-Star は NVIDIA 前提
  （[Hardware Guide](https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html)）。
  FluidX3D は OpenCL で全ベンダー対応だが非商用ライセンス・MPI なし。
  → 「商用可・全ベンダー・ブラウザまで」の組合せはうちだけ。
- **精度の透明性**: FP16 導入時も「何を失うか」を検証スイートで数値化して出す
  （Lehmann らの知見 [arXiv:2112.08926](https://arxiv.org/abs/2112.08926) を製品仕様に落とす）。

## 11. アイドル機での正式測定チェックリスト（公開前必須）

環境:
- [ ] 共有負荷なしのアイドル実機（本 M5 Max）。電源接続・熱定常化・他アプリ終了を記録
- [ ] OS / wgpu / ドライバ / コンパイラのバージョン、メモリ構成・公称帯域を記録
- [ ] 実測メモリ帯域（STREAM 相当）を取り、ルーフライン上限を先に計算して併記

測定プロトコル:
- [ ] ウォームアップ N ステップを除外した計測窓、5 回以上の中央値、分散を記録
- [ ] 格子サイズスイープ（2D: 512²〜4096²、3D: 128³〜メモリ上限）で飽和カーブを公開
- [ ] MLUPS 定義を明記（総セル数×ステップ/実時間、出力 I/O 除外の旨）
- [ ] 「カーネル単体（空箱・境界最小）」と「代表シナリオ（障害物+出力あり）」を別掲

比較条件合わせ（FluidX3D と直接対決するとき）:
- [ ] FluidX3D を**同一マシン**でビルド・実行し（Apple Silicon の OpenCL で可）、README 値では
      なく同一環境実測で比較する（先方の測定条件: D3Q19 SRT・空箱・典型 256³・FP32 演算）
- [ ] ステンシル（D3Q19 vs D2Q9）と格納精度（f32 vs FP16S/C）をラベルに必ず併記
- [ ] うちの偏差格納 f32 と先方 FP32 の違い（精度検証済みである旨）を脚注に

公開物:
- [ ] 本ドラフトの自陣数字を正式値に差し替え、「暫定」注記を削除
- [ ] M-Star チャートの目視読み取り値を確定させる場合は、読み取り方法と誤差幅を明記
- [ ] waLBerla SC13 アブストラクト数値（1.93 兆セル更新/s・1.8M スレッド）の一次確認
      （ACM ページ or 論文 PDF 目視）
- [ ] 3D GPU 実装後: 単GPU D3Q19 f32 の M5 Max 実測を表2の Apple 行（M2 Max/M2 Ultra）と比較
- [ ] FP16 格納実装後: FP16S 相当モードで表2 と再比較 + 検証スイートの劣化定量を同時公開

## 12. 出典一覧

| # | ソース | 用途 |
|---|---|---|
| 1 | https://github.com/ProjectPhysX/FluidX3D | FluidX3D ベンチ表・測定条件・算術強度・メモリ/セル |
| 2 | https://raw.githubusercontent.com/ProjectPhysX/FluidX3D/master/README.md | 同上（生データ、2026-07-05 取得） |
| 3 | https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html | M-Star スケーリングベンチ条件・プラットフォーム・経験則 |
| 4 | https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_agitated.svg / …_pipe.svg | M-Star チャート実体（軸レンジ確認） |
| 5 | https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html | M-Star NVIDIA 要件・VRAM 経験則 |
| 6 | https://mstarcfd.com/software/ | M-Star マーケ主張（定性） |
| 7 | https://arxiv.org/html/2402.13171v1 | waLBerla-wind: A100 D3Q27 実測 1,677 / ルーフライン 7,513、EPYC 7763 実測、弱スケーリング |
| 8 | https://ar5iv.labs.arxiv.org/html/1511.07261 | waLBerla 兆セル・40万コア超の直接引用、MLUP/s 定義 |
| 9 | https://dl.acm.org/doi/10.1145/2503210.2503273 | waLBerla SC13（1.93兆更新/s 等はアブストラクト由来、要目視確認） |
| 10 | https://www.openlb.net/performance/ | OpenLB 512×A100 1.33 TLUPS、Aurora 21,120 GLUPS、Magnus 142,479 MLUPS、強スケーリング効率 |
| 11 | https://www.openlb.net/news/openlb-release-1-5-available-for-download/ | OpenLB 1000³ キャビティ 42.2/24.8/2.7 GLUPS、GPU 15.6× |
| 12 | https://arxiv.org/html/2506.21804v1 | OpenLB 異種混成 HPC、18G セル、効率 0.66〜0.91 |
| 13 | https://arxiv.org/html/2506.09242v1 | Palabos GPU: A100 理論ピーク 9.481/4.921 GLUPS・実測 75〜85%・スケーリング、waLBerla 85% 言及 |
| 14 | https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306 | STLBM: GPU/CPU の D3Q19 FP64 MLUPS |
| 15 | https://arxiv.org/abs/2112.08926 | FP64/FP32/16bit 精度影響（Phys. Rev. E 106, 015308） |
| 16 | https://arxiv.org/abs/1711.11468 | LBM ベンチカーネル（実装差の影響、方法論） |
| 17 | https://journals.sagepub.com/doi/full/10.1177/10943420211016525 | 最適化幅の実例（211→550 MLUPS） |
| 18 | 本リポジトリ実測（PERFORMANCE.md 系譜） | 自陣数字（共有負荷下・暫定） |

**未確認として今回掲載を見送った数字**: OpenLB 単一 A100 の「8.3 GLUPS（Periodic Shift, D3Q19 FP32）」
— 検索スニペット上は [ScienceDirect の第三者論文アブストラクト](https://www.sciencedirect.com/science/article/abs/pii/S0010465522003228)
に出現するが、当該ページ・一次論文（[Wiley cpe.7509](https://onlinelibrary.wiley.com/doi/full/10.1002/cpe.7509)）とも
直接取得不可（403/402）のため不採用。waLBerla の「889,602 MLUPS」も出典機種の帰属が
確認できず不採用（1.93 兆更新/s 系の数字と整合しないため）。
