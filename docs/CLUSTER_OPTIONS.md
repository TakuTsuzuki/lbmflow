# クラスタ/クラウド HPC 選択肢調査 — R3 マルチノード実測に向けて（2026-07-05）

docs/MPI_GUIDE.md「クラスタでやるべき測定リスト」8 項目（R3 完了条件）を実測するための
計算資源の選択肢調査と意思決定メモ。対象読者: 日本在住の個人〜小規模組織。

**免責**: 価格・仕様は各節に出典 URL と取得日を明記（全て 2026-07-05 取得、断りない限り税抜/税込は出典表記に従う）。
クラウド価格は変動するため発注前に再確認のこと。円換算は **1 USD = 160 円**（2026-07-03 実勢 161.27 円
[Trading Economics](https://tradingeconomics.com/japan/currency)）、**1 EUR = 185 円（概算仮定）** で計算した概算。
「要確認」と書いた箇所は一次情報を確認できていない。

---

## 0. 必要リソースの見積り

MPI_GUIDE の 8 項目から逆算した所要量（推定はリポジトリ実測値 M5 Max 40 MLUPS/rank・D2Q9 から外挿）:

| 測定 | 規模 | メモリ所要（推定） | コア/ノード所要 |
|---|---|---|---|
| 1. 弱スケーリング 1→64 ランク（3D 128³/rank D3Q19 f64） | 64 ランク | f ダブルバッファ 0.59 GiB + モーメント等 ≈ **0.7〜0.8 GiB/rank**（計 ~50 GiB） | 64 物理コア以上、**2 ノード以上に分散必須** |
| 2. 強スケーリング 固定 1024³、8→512 ランク | 512 ランク | 全体 **~340〜370 GiB**（8 ランク時 ~45 GiB/rank） | **512 物理コア** + 集約メモリ ~400 GiB |
| 3. 通信/計算比（mpiP / MPI_T） | 上記に相乗り | — | プロファイラのビルド権限 |
| 4. BTL/MTL（UCX/OFI、eager/rendezvous） | 2 ノード以上 | — | **RDMA 系ファブリック（IB/EFA 等）があるほど価値大** |
| 5. ランク×スレッドのハイブリッド格子 | 1〜2 ノード | — | ノード内コア数が多いほど格子を広く振れる |
| 6. --map-by/--bind-to、NUMA | 1 ノード〜 | — | Linux 必須（macOS 不可）、2 NUMA 以上が望ましい |
| 7. 他 MPI 実装（MPICH 系）での test_mpi.sh 全 PASS | 小規模 | — | MPICH をビルドできれば可 |
| 8. 診断 Allreduce / rank-0 gather のスケール限界 | 64〜512 ランク | gather は rank0 に全場（1024³ で ~170 GiB/場は非現実→縮小問題で傾向を取る） | ランク数が多いほど良い |

- 3D 128³/rank の 1 ステップは ~0.1〜0.3 秒/rank（8〜20 MLUPS/rank 想定）。220 ステップ計測 ≈ 25〜60 秒/構成。
  弱スケーリング曲線 + ハイブリッド格子 + 強スケーリング曲線 + プロファイル取得で**正味計算 1〜2 時間**、
  環境構築・試運転込みで **1 セッション = 8〜16 ノード × 4〜8 時間（= 約 50〜150 ノード時間）** を campaign の基準単位とする。
- 注意: 1024³ の 8 ランク起点は 45 GiB/rank 必要 → ノードメモリ 128 GiB なら 2 ランク/ノード × 4 ノードに分散すれば可。
  富岳系（32 GiB/ノード）では 8 ランク起点は物理的に不可（16 ランク起点 1 ランク/ノード ≈ 23 GiB/rank から、または 768³ に縮小）。

---

## 1. 総覧比較表

費用は上記「1 キャンペーン（≈50〜150 ノード時間）」の概算。◎○△× は 8 項目をどこまで測れるか（§5 に詳細）。

| 選択肢 | キャンペーン費用目安 | リードタイム | 測れる範囲 | 手間 | 一言 |
|---|---|---|---|---|---|
| **AWS hpc7g ×8〜16（ParallelCluster + EFA）** | **1.3〜3.5 万円**（us-east-1） | 2〜5 日（クォータ申請） | **8/8 ◎** | 中 | 本命。arm64 で手元と同 ISA |
| Azure HBv4 ×3〜4（InfiniBand NDR 400G） | OD 2.8 万円 / Spot 約 0.5 万円 | 3〜10 日（HB クォータ審査） | **8/8 ◎**（項目 4 は UCX/IB で最良） | 中〜高 | IB 実測の対照実験に最適 |
| GCP H3 ×8 / H4D ×4 | H3 3.8 万円 / H4D 5.3 万円 | 2〜5 日 | H3 7/8（RDMA なし）/ H4D 8/8 | 中 | 日本近傍は H4D シンガポールのみ |
| 富岳 試行課題（無償枠） | **0 円** | **2〜4 週**（審査 1〜2 週 + HPCI 手続） | 8/8（項目 4 は Tofu 固有） | 中〜高 | 無償・随時受付。A64FX 移植が必要 |
| ABCI 3.0（rt_HF マルチノード） | 実消費 ~17 万円 / **最低購入 22 万円** | 3〜4 週 | 8/8 | 低〜中 | CPU 測定に H200×8 ノード課金は過剰 |
| TSUBAME4.0（産業利用・成果公開） | **最低 1 口 11 万円**（実消費 ~1.5 万円分） | 数週間（要確認） | 8/8 | 低〜中 | 192 コア/ノード + IB。GPU 展開の布石に |
| 東大 Wisteria/Miyabi（トライアル） | 数千円/セット〜（企業枠は要確認） | 数週間（要確認） | 8/8（Odyssey は富岳同系） | 中〜高 | 学術系最安。企業はトライアル枠経由 |
| FOCUS スパコン（産業界専用） | **1〜3 万円**（年会費 1 万 + 従量、初年度無料枠 1 万円分） | 1〜3 週（要問合せ） | 6.5/8（512 ランク不可） | 低〜中 | 国内・安価・InfiniBand あり。規模が小さい |
| Hetzner 専用サーバ ×4 + 手動 MPI | **約 15〜18 万円/月 + 初期費 ~6 万円**（月契約） | 3〜10 日 | 4/8（≥80% 合格線はネットワーク律速で困難） | **高** | 常設ミニクラスタ用途なら意味あり |
| 単一大インスタンス（c8g/c7a.48xlarge, 192 vCPU） | **0.3〜1.3 万円** | **即日** | 4.5/8（ノード間なし） | **低** | 今日できる前哨戦。R3 単独達成は不可 |

---

## 2. クラウド HPC

### 2.1 AWS — hpc7g + EFA + ParallelCluster（本命）

- **hpc7g.16xlarge**: Graviton3E **64 物理 vCPU / 128 GiB / EFA 200 Gbps**。
  オンデマンド **$1.6832/h（us-east-1）、$2.1117/h（東京 ap-northeast-1、+25%）**。**Spot 非対応**（HPC 系は全て）。
  出典: [Vantage](https://instances.vantage.sh/aws/ec2/hpc7g.16xlarge)・[aws-pricing.com](https://aws-pricing.com/hpc7g.16xlarge.html)（2026-07-05）。
  提供リージョンは バージニア北部 / 東京 / アイルランド / GovCloud（[AWS ニュース 2023-09](https://aws.amazon.com/about-aws/whats-new/2023/09/amazon-ec2-hpc7g-instances-additional-regions/)）。
  2026-07-05 時点で Graviton4 世代の hpc8g は未確認。
- **hpc7a.96xlarge**（AMD EPYC 9R14 **192 コア / 768 GiB / EFA 300 Gbps**）: $9.0793/h、Spot 非対応
  （[Vantage](https://instances.vantage.sh/aws/ec2/hpc7a.96xlarge)、2026-07-05）。ノード数を減らしたい場合の代替。
- **構成例と費用**（us-east-1、1 USD=160 円）:
  - hpc7g ×8 × 6 h = 48 NH × $1.6832 ≈ **$81 ≈ 1.3 万円**（512 物理コア: 弱 64 ランクも強 512 ランクも成立）
  - hpc7g ×16 × 8 h = 128 NH ≈ **$215 ≈ 3.5 万円**（余裕をみた上限）
  - ヘッドノード（c7g.large 等）+ EBS 共有は数百円オーダー。東京リージョンなら +25%。
- **MPI 対応**: EFA は libfabric(OFI) 経由。ParallelCluster が EFA ドライバ + Open MPI + Slurm を自動構築
  （[ParallelCluster docs](https://docs.aws.amazon.com/parallelcluster/latest/ug/slurm-workload-manager-v3.html)、
  hpc7g 構築例: [Sean Smith のガイド](https://swsmith.cc/posts/hpc7g-parallelcluster.html)）。
  rsmpi は mpicc プローブで素直に通る見込み（Open MPI 系）。項目 4 は **OFI/EFA 側面**を実測でき、
  `FI_PROVIDER=tcp` との比較も可能。UCX/InfiniBand 側面は測れない（→ Azure か国内 IB 機で補完）。
- **セットアップの手間（中）**: ① アカウント + **サービスクォータ「Running On-Demand HPC instances」の引き上げ申請**
  （新規アカウントは既定値が小さい/0 のことがある。申請から 1〜3 営業日、要ユースケース記入）
  ② `pip install aws-parallelcluster` → YAML 1 枚で Slurm クラスタ生成 ③ placement group（cluster）指定。
  マネージド Slurm の AWS PCS もあるが、単発キャンペーンには ParallelCluster で十分。
- **注意**: HPC インスタンスはリージョン内 AZ が限られる。実行後は `pcluster delete-cluster` で確実に破棄（課金停止）。

### 2.2 Azure — HBv4（InfiniBand NDR 400 Gbps）

- **HB176rs_v4**: AMD EPYC Genoa-X **176 物理コア（SMT 無効）/ 768 GiB / NDR InfiniBand 400 Gb/s**
  （[Microsoft Learn](https://learn.microsoft.com/en-us/azure/virtual-machines/hbv4-series-overview)）。
  オンデマンド **$7.20/h、Spot $1.331/h**（基準リージョン、[Vantage](https://instances.vantage.sh/azure/vm/hb176)、2026-07-05）。
  **HPC VM なのに Spot 可**なのが最大の魅力（中断リスクは短時間ベンチなら許容しやすい）。
- リージョン: East US 系・West Europe・**Southeast Asia・Korea Central** など。**東日本には無い**
  （[Spare Cores](https://sparecores.com/server/azure/Standard_HB176rs_v4)、2026-07-05）。
- 構成例: ×4 ノード（704 コア）× 6 h = OD $173 ≈ **2.8 万円** / Spot $32 ≈ **5 千円**。
  1024³ 強スケーリングは 3 ノードで足りる（528 コア・2.3 TiB）。
- MPI: Azure HPC イメージ（AlmaLinux/Ubuntu-HPC）に Mellanox OFED + HPC-X（Open MPI/UCX）+ Intel MPI が同梱。
  **項目 4 の UCX/IB 側面（eager/rendezvous、タグマッチング）を最も HPC らしい形で測れる**。
- 手間（中〜高）: HB ファミリの vCPU クォータ申請が必要で、**新規/従量課金サブスクリプションでは却下されることがある**
  （サポートリクエスト経由、数日〜）。クラスタ組みは CycleCloud か手動（同一 PPG + IB 確認）。ここが AWS より一段面倒。
- 次世代 HBv5（HBM 搭載 EPYC・800 Gbps IB）は発表済みだが本調査では価格未確認（要確認）。

### 2.3 GCP — H3 / H4D + compact placement

- **h3-standard-88**: Sapphire Rapids **88 vCPU（SMT 無効）/ 352 GB / 200 Gbps**、compact placement 対応。
  **$4.9236/h（us-central1）、Spot 非対応**。提供は us-central1 / europe-west4 / northamerica-northeast1 の 3 リージョンのみ
  （[gcloud-compute.com](https://gcloud-compute.com/h3-standard-88.html)、2026-07-05）。
  **ノード間は RDMA なし（gVNIC/TCP）** → 項目 4 の実測価値が下がるのが難点。
  ×8 × 6 h = $236 ≈ **3.8 万円**。
- **h4d-standard-192**: EPYC Turin **192 vCPU（SMT 無効）/ 720 GB / Cloud RDMA (Falcon)**。
  **2026-03 GA**、us-central1-a / europe-west4-b / **asia-southeast1-a（シンガポール）**
  （[Google Cloud ブログ](https://cloud.google.com/blog/products/compute/h4d-vms-now-ga)）。
  参考価格 **~$13.74/h 相当**（$10,033/月、[CloudPrice](https://cloudprice.net/gcp/compute/instances/h4d-standard-192)、2026-07-05。
  ブログには DWS Flex Start で「3 セント/コア時間から」の記載 → 条件次第で半額以下）。×4 × 6 h ≈ $330 ≈ **5.3 万円**。
- 手間（中）: Cluster Toolkit で Slurm 構築可。クォータ + ゾーン在庫の確認が必要。
  RDMA 対応 MPI の設定（Intel MPI / Open MPI の対応バージョン）は GCP ドキュメント準拠で一手間。
- 判定: 価格・情報量とも AWS/Azure に対する優位が薄い。**GCP に既存資産がある場合のみ**推奨。

### 2.4 参考 — OCI（未調査）

Oracle Cloud の BM.HPC 系はベアメタル + RDMA クラスタネットワーキングで安価という評判があるが、
今回は価格・在庫を未調査（要調査。候補に残す価値はある）。

---

## 3. 国内学術・公的系

### 3.1 富岳 — 試行課題（無償）が第一候補

- **試行課題（一般/産業）: 無償・随時受付・審査結果 1〜2 週間・最長 6 ヶ月・上限 10 万 NH**。
  さらに簡易な「ファーストタッチオプション」は 1,000 NH 固定・最長 3 ヶ月
  （[HPCI 富岳試行課題](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_trial)、2026-07-05）。
  1,000 NH でも本キャンペーン（50〜150 NH）には十分すぎる。
- 有償に進む場合: 従量制 **98.64 円/NH（成果非公開）、49.32 円/NH（公開）**
  （[HPCI 料金ページ](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_price)、2026-07-05）。
  例: 128 ノード × 10 h = 1,280 NH ≈ 12.6 万円（非公開）と、有償でもクラウドと同水準。
- ノード: A64FX 48 コア / **32 GiB HBM2** / Tofu-D（公称仕様: [R-CCS](https://www.r-ccs.riken.jp/fugaku/)）。
  - 弱スケーリング 64 ランク: 8 ノード × 8 ランク等で成立。512 ランクも 11 ノード〜で余裕。
  - **1024³ の 8 ランク起点は 32 GiB/ノード制約で不可**（§0 注意参照。16 ランク起点 or 768³ に変更）。
- 手間（中〜高）: ① HPCI アカウント・電子証明書の手続きが加わる ② **aarch64 は良いが、MPI は富士通 MPI
  （Open MPI ベース、ラッパー名 mpifcc 等）**。rsmpi の mpicc プローブとの相性・bindgen まわりで一手間の可能性（要検証）。
  ③ ジョブは pjsub（富士通 TCS）。項目 7 の「実装非依存のビット一致」を Open MPI 系以外に広げる意味では価値が高い。
- 資格: 募集要領上、企業・大学等を排除していないが、**無所属の個人での応募可否は要確認**（申請書に所属・研究計画が必要）。

### 3.2 ABCI 3.0（産総研）

- 料金: **1 ポイント = 220 円（税込）、最低購入 1,000 pt = 22 万円**。
  マルチノード MPI に必要な **rt_HF（フルノード）は 16 pt/h = 3,520 円/ノード時**（Spot/On-demand 区分）
  （[ABCI 料金 2026 年度](https://abci.ai/ja/how_to_use/tariffs.html)、2026-07-05）。
  例: 8 ノード × 6 h = 768 pt ≈ **16.9 万円分**（最低購入 22 万円の枠内）。
- ノード（H）: Xeon Platinum 8558 48c ×2 = **96 物理コア / 2 TB / InfiniBand NDR200 ×8** + **H200 ×8**。
  マルチノードは rt_HF のみ、**最大 128 ノード**（[ABCI 3.0 docs: システム概要](https://docs.abci.ai/v3/ja/system-overview/)・
  [ジョブ実行](https://docs.abci.ai/v3/ja/job-execution/)、2026-07-05）。
- リードタイム: 利用申請の審査 ~10 営業日 + ポイント付与まで最大 10 日 → **実質 3〜4 週間**
  （[ご利用の流れ](https://abci.ai/ja/how_to_use/)・[ポイント申請案内](https://abci.ai/news/2025/01/23/ja_news_Point_Application.html)、2026-07-05）。
- 判定: CPU の R3 測定だけなら **H200×8 が遊ぶノードに 3,520 円/h は割高**で、最低購入 22 万円も重い。
  ただし M-E（GPU/CUDA 展開、GPUDirect）まで見据えるなら国内最有力なので、**「R3 は別で済ませ、M-E の時に申請」**が合理的。
  資格: 約款上は法人所属が前提の運用（無所属個人は実質困難。個人事業主は要問合せ）。

### 3.3 TSUBAME4.0（東京科学大）

- 料金（学外・従量）: **成果公開 275 円/ノード時、非公開 1,100 円/ノード時**。
  **学外は 1 口 = 400 ノード時 = 11 万円（公開）/ 44 万円（非公開）が最低購入単位**。ポイントは年度末（3/31）失効
  （[TSUBAME4 利用料の概略](https://www.t4.cii.isct.ac.jp/fare_overview)、2026-07-05）。
- ノード: EPYC 9654 ×2 = **192 物理コア / 768 GiB / InfiniBand NDR200 ×4 / H100 ×4**（公称構成は
  [TSUBAME4 サイト](https://www.t4.cii.isct.ac.jp/)参照）。キャンペーン実消費は ~40 NH ≈ 1.1 万円分だが最低 1 口 11 万円。
- 判定: H100 込みノードが 275 円/h（公開条件）は破格で、**GPU 展開の布石も兼ねるなら ABCI より安い**。
  ただし成果公開義務（利用報告書）と年度失効、申請リードタイム（数週間、要確認）を織り込むこと。
  トライアルユース（無償枠）の現行有無は料金ページには記載がなく要確認。

### 3.4 東大情報基盤センター — Wisteria/BDEC-01・Miyabi

- 通常利用（トライアル）FY2026: **Wisteria/BDEC-01 1 セット 2,250 円（720 トークン）、Miyabi 1 セット 7,500 円（720 トークン）**、
  最大 12 セット（[利用負担金ページ](https://www.cc.u-tokyo.ac.jp/guide/application/charge_trial.php)、2026-07-05。
  トークン→ノード時間の換算係数は機種別なので同ページで要確認）。
  **企業向けには「企業利用トライアル」枠（Wisteria）**が別にある（[案内](https://www.cc.u-tokyo.ac.jp/guide/trial/company.php)）。
- Wisteria-Odyssey は A64FX + Tofu-D（富岳同系）→ 富岳試行の前哨・代替になる。Miyabi は GH200 主体（GPU 布石側）。
- 判定: 金額は最安級だが、**学術利用が主対象で企業・個人は枠と審査の確認が必須**。リードタイムは数週間（要確認）。

### 3.5 FOCUS スパコン（計算科学振興財団、神戸）— 産業界専用

- 料金（従量）: **F システム（CPU）300 円/ノード時**（利用ノード数で 150 円まで逓減）、
  **X システム（A64FX）80 円/ノード時**、S システム（EPYC 9654 192c、8 コア VM 単位）60 円/VM 時など。
  アカウント料 **1 万円/従事者・年度**、**初年度は従量 1 万円分の無料枠**。無償の**試行利用**制度あり（ポーティング・ベンチ用途）
  （[料金](https://www.j-focus.or.jp/focus/fee.html)・[利用形態](https://www.j-focus.or.jp/focus/form.html)・
  [試行利用](https://www.j-focus.or.jp/focus/free-trial.html)、2026-07-05）。
- システム（[利用案内](https://www.j-focus.jp/user_guide/ug0001000000/)、2026-07-05）:
  F = Xeon E5-2698v4 40c/128GB/**InfiniBand FDR 56G**（12 ノード規模・2016 年整備）、
  X = **A64FX 48c/32GB/InfiniBand EDR 100G**、Z = Xeon 40c/**EDR 100G**、S = EPYC 9654 192c/**100GbE**。
  各システムのノード台数は要問合せ（unyo@j-focus.or.jp）。
- 判定: **国内・安価・実 InfiniBand・産業利用専用（中小企業歓迎）**で、項目 1/3/4/5/6/7 は測れる。
  ただし規模が小さく **512 ランクの強スケーリング（項目 2）は不可**（F は最大 480 コア）。
  X（A64FX+EDR）は 80 円/h で富岳移植のリハーサルにもなる。個人事業主の可否・リードタイムは要問合せ（目安 1〜3 週）。

---

## 4. お手軽段（すぐ・安く）

### 4.1 単一の大コア数インスタンスでランク数だけ稼ぐ

- 候補（192 vCPU 級、2026-07-05 取得）:
  - **c8g.48xlarge**（Graviton4 192 vCPU/384 GiB）: OD $7.657/h、**Spot $2.727/h**（[Vantage](https://instances.vantage.sh/aws/ec2/c8g.48xlarge)）。
    arm64 なので手元 M5 Max・hpc7g と ISA が揃うのが利点。
  - **c7a.48xlarge**（EPYC 192 vCPU/384 GiB）: OD $9.853/h、**Spot $3.312/h**（[Vantage](https://instances.vantage.sh/aws/ec2/c7a.48xlarge)）。
  - GCP c3d/c4d の 360 vCPU 級、Azure HBv4 単騎（176c、Spot $1.33/h）も同用途に使える。
- 費用: **Spot なら 8 時間で 3〜4 千円、オンデマンドでも 1 万円前後。即日開始可**（大型サイズの vCPU クォータだけ注意）。
- できること: 項目 1 の「単一ノード版 64 ランク」（MPI_GUIDE の既存表を均質コア・Linux で取り直し、n≤64 まで延長）、
  **項目 5（ハイブリッド格子）・6（map-by/bind-to、2 NUMA）・7（MPICH を入れて共有メモリ経路の再検証）・8（〜192 ランクの傾向）**。
- **限界（明確に）**: 全ランクが共有メモリ（vader/xpmem）で繋がるため、**ノード間ネットワーク（項目 4）と
  「マルチノード弱スケーリング ≥80%」という R3 の本丸は原理的に測れない**。
  また 384 GiB では 1024³（~370 GiB）が載らないため、項目 2 は 768³ への縮小か r8g/m8g 等（768 GiB〜1.5 TiB、価格要確認）が必要。
  位置づけは「クラスタ本番前の較正・デバッグ・単一ノード基準線の確定」。

### 4.2 Hetzner 等のベアメタル複数台 + 手動 MPI

- 価格例: **AX162-R**（EPYC 9454P 48c/96t、256 GB DDR5）**€199/月 + 初期費 €79**
  （[Hetzner プレス](https://www.hetzner.com/pressroom/new-ax162/)。2026-06 の値上げで €238/月 という情報もあり:
  [Northflank まとめ](https://northflank.com/blog/hetzner-cloud-server-price-increases)、いずれも 2026-07-05 取得）。
  4 台で **月額 €796〜952 ≈ 15〜18 万円 + 初期費 €316 ≈ 6 万円**（月契約。時間貸しではない）。
  Hetzner Cloud の CCX63（48 vCPU = 24 物理コア、€1.37/h、[Spare Cores](https://sparecores.com/server/hcloud/ccx63)）は
  時間貸しだが物理コアが半分でネットワーク保証もなく、MPI 用途には勧めない。
- ネットワーク: 標準 1 GbE、vSwitch/10G アップリンクはオプション（構成・追加費は要確認）。**RDMA はない（TCP のみ）**。
- 定量的な見立て: 3D 128³/rank の面ハロー ≈ 0.6〜2.5 MB/面/step。16 ランク/ノードでノード出入り数十 MB/step、
  ステップ 0.1〜0.3 秒 → **1 GbE（~120 MB/s）は確実に飽和、10 GbE でも TCP レイテンシ込みで効率 ≥80% は厳しい**。
  R3 の合格線を「コードのせいでなくネットワークのせいで」落とす結果になりがち。
- 判定: 項目 4 の「TCP BTL での挙動」と項目 1 の下界は取れるが、**R3 達成の主戦場にはならない**。
  月額固定でクラスタを常設したい（CI 的に回す）段階になったら再検討。

---

## 5. 8 項目カバレッジ対応表

◎=そのまま測れる ○=測れる（注記付き） △=部分的 ×=不可

| # | 測定（MPI_GUIDE §測定リスト） | AWS hpc7g | Azure HBv4 | GCP H3/H4D | 富岳(試行) | ABCI 3.0 | TSUBAME4 | FOCUS | Hetzner | 単一大ノード |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 弱スケーリング 1→64（マルチノード） | ◎ | ◎ | ○/◎ | ◎ | ◎ | ◎ | ○(F/X/Z) | △(効率低) | △(ノード内のみ) |
| 2 | 強スケーリング 1024³ 8→512 | ◎(8 ノード) | ◎(3 ノード) | ○(H3 は TCP 律速) | ○(16 ランク起点) | ◎ | ◎ | ×(≤480 コア) | × | ×(192 まで・メモリ不足) |
| 3 | 通信/計算比（mpiP/MPI_T） | ◎ | ◎ | ◎ | ○(富士通 MPI 流儀) | ◎ | ◎ | ◎ | ◎ | ○(共有メモリ比) |
| 4 | BTL/MTL・閾値（UCX/OFI） | ○(OFI/EFA 面) | ◎(UCX/IB 面) | △(H3 TCP)/○(H4D RDMA) | ○(Tofu 固有) | ◎(IB NDR) | ◎(IB NDR) | ◎(IB FDR/EDR) | △(TCP のみ) | ×
| 5 | ランク×スレッド最適点 | ◎(64c) | ◎(176c) | ◎ | ○(48c/32GB) | ◎(96c) | ◎(192c) | ○ | ○ | ◎(192c) |
| 6 | map-by/bind-to・NUMA | ◎ | ◎ | ◎ | ○(pjsub 流儀) | ◎ | ◎ | ◎ | ◎ | ◎(2 NUMA) |
| 7 | 他 MPI 実装で test_mpi.sh | ◎(MPICH 追加可) | ◎(HPC-X/Intel MPI 同梱) | ◎ | ◎(富士通 MPI=別系統) | ◎ | ◎ | ○ | ◎ | ○(共有メモリ経路のみ) |
| 8 | 診断 Allreduce/gather の限界 | ◎(512 ランク) | ◎ | ◎ | ◎(それ以上も) | ◎ | ◎ | △(〜480) | △ | △(〜192) |

---

## 6. 推奨シナリオ

### 案 A — 最安で 60 点（総額 数千円〜1.5 万円、今週から）

1. **今日〜明日**: AWS **c8g.48xlarge（Spot）1 台**を数時間借り、
   項目 5・6・7（共有メモリ範囲）・8（〜192 ランク）と「単一ノード 64 ランクの基準線」を消化する。
   Linux/Graviton なので手元 arm64 ビルドがほぼそのまま通る。費用 3〜10 千円。
   （1024³ を触るなら 768³ に縮小するか、メモリの大きい r8g 系を選ぶ）
2. **同時に**: **富岳 試行課題（ファーストタッチ 1,000 NH、無償）を申請**しておく（審査 1〜2 週 + HPCI 手続）。
   通れば残りの本丸（マルチノード弱スケーリング・項目 2・4）を**無償**で測れる。
- 到達点: 富岳が通るまでは R3 の「64 ランクで ≥80%」のマルチノード実測だけが残る（=60 点）。
  富岳到着後に 100 点まで伸ばせるが、**A64FX への移植・富士通 MPI との接続確認という技術リスクと数週間の待ち**を抱える。

### 案 B — しっかり 90 点（総額 2〜5 万円、リードタイム 2〜5 日）

1. AWS アカウントで **HPC インスタンスのクォータ引き上げを申請**（512〜1,024 vCPU、1〜3 営業日）。
2. **ParallelCluster で hpc7g.16xlarge × 8〜16（us-east-1、EFA、cluster placement group）** の Slurm クラスタを半日だけ立てる。
   8 項目を 1 セッションで全て実測（項目 4 は OFI/EFA 側面）。費用 1.3〜3.5 万円 + ヘッドノード少額。
3. （+5 点のオプション）**Azure HBv4 Spot ×3〜4 を数時間**（約 5 千円）追加し、
   UCX/InfiniBand 側面（項目 4）と「EFA と IB で効率がどう変わるか」の対照を取る → 95 点。
   ※ HB クォータ審査に日数がかかる・却下リスクがあるため、こちらは「取れたらやる」扱い。
- 到達点: R3 完了条件を自前のスケジュールで満たせる。測定は再現可能（YAML + スクリプトを repo に残せる）。

### 併走の推奨（どちらの案でも）

- **富岳試行課題は無償なので、案 B でも出しておいて損がない**（項目 7 の「Open MPI 以外」の最有力サンプル）。
- **ABCI 3.0 / TSUBAME4 は R3 用ではなく M-E（GPU/CUDA、GPUDirect）用の布石**として位置づける。
  CPU だけの R3 に最低購入 11〜22 万円を払うのは割に合わない。
- 国内で完結させたい/請求書払いが必須なら **FOCUS（試行利用→従量）**が現実解（ただし項目 2 の 512 ランクは断念）。

---

## 7. 次のアクション

### ユーザー（tsuzuki さん）がやること

1. **方針の決定**: 案 A / 案 B / 併走（推奨は「案 B + 富岳申請の併走」）と予算上限の確認。
2. AWS アカウント準備と**クォータ申請**（案 B: us-east-1 の「Running On-Demand HPC instances」を 512〜1,024 vCPU、
   案 A: Standard/Spot の 192 vCPU）。申請文面はこちらで下書き可能。
3. **富岳試行課題の申請**（HPCI アカウント取得含む）: 所属・課題概要（LBM ソルバの弱/強スケーリング実測、数百字）が必要。
   個人事業主としての応募可否は HPCI ヘルプデスクに事前確認。
4. （オプション）Azure サブスクリプション作成と HB ファミリのクォータ申請（却下されたら諦めて AWS のみ）。
5. FOCUS を使う場合: unyo@j-focus.or.jp へ利用資格（小規模事業者）・X/Z システムのノード数・試行利用の可否を問い合わせ。

### 私たち（リポジトリ側）で準備しておけること

1. **bench_mpi の 3D 対応**: 現行 `crates/lbm-core2/examples/bench_mpi.rs` は 2D 512²/rank・帯分割固定。
   128³/rank・D3Q19・任意デカルト分割・RESULT 行の共通フォーマット（ranks/nodes/threads/MLUPS/効率）に拡張する。
2. **Slurm 版ドライバ**: `scripts/bench_mpi.sh` の sbatch テンプレート化
   （ranks×nodes×threads の格子実行、`--map-by`/`--bind-to` の系統振り、結果 CSV 集計、MPI_GUIDE の表形式への整形）。
3. **ParallelCluster 設定 YAML の草案**（head: c7g.large、compute: hpc7g×16・EFA on・placement group、共有 /home、
   ジョブ完了後の自動スケールダウン）と実行手順書。
4. **項目 3 用のプロファイル手順**: mpiP のビルド手順メモ or Open MPI の `MPI_T`/OSU マイクロベンチによる代替計測スクリプト。
5. **項目 7 用の MPI 実装マトリクス**: test_mpi.sh を Open MPI / MPICH / Intel MPI で回す手順（rsmpi の再ビルド込み）。
6. **富岳向け移植チェックリスト**: aarch64 ビルド（済: M5 Max で常用）、富士通 MPI の mpicc ラッパー対応、pjsub ジョブ雛形。

---

## 8. 出典一覧（すべて 2026-07-05 取得）

- AWS: [hpc7g 製品ページ](https://aws.amazon.com/ec2/instance-types/hpc7g/) /
  [Vantage hpc7g.16xlarge](https://instances.vantage.sh/aws/ec2/hpc7g.16xlarge) /
  [aws-pricing.com hpc7g.16xlarge（リージョン別）](https://aws-pricing.com/hpc7g.16xlarge.html) /
  [hpc7g 追加リージョン](https://aws.amazon.com/about-aws/whats-new/2023/09/amazon-ec2-hpc7g-instances-additional-regions/) /
  [Vantage hpc7a.96xlarge](https://instances.vantage.sh/aws/ec2/hpc7a.96xlarge) /
  [Vantage c8g.48xlarge](https://instances.vantage.sh/aws/ec2/c8g.48xlarge) /
  [Vantage c7a.48xlarge](https://instances.vantage.sh/aws/ec2/c7a.48xlarge) /
  [ParallelCluster docs](https://docs.aws.amazon.com/parallelcluster/latest/ug/slurm-workload-manager-v3.html) /
  [hpc7g + ParallelCluster 構築例](https://swsmith.cc/posts/hpc7g-parallelcluster.html) /
  [EFA docs](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/efa.html)
- GCP: [gcloud-compute.com h3-standard-88](https://gcloud-compute.com/h3-standard-88.html) /
  [H3 発表ブログ](https://cloud.google.com/blog/products/compute/new-h3-vm-instances-are-optimized-for-hpc) /
  [H4D GA ブログ](https://cloud.google.com/blog/products/compute/h4d-vms-now-ga) /
  [CloudPrice h4d-standard-192](https://cloudprice.net/gcp/compute/instances/h4d-standard-192)
- Azure: [HBv4 シリーズ概要（Microsoft Learn）](https://learn.microsoft.com/en-us/azure/virtual-machines/hbv4-series-overview) /
  [Vantage HB176rs_v4](https://instances.vantage.sh/azure/vm/hb176) /
  [Spare Cores HB176rs_v4（リージョン）](https://sparecores.com/server/azure/Standard_HB176rs_v4)
- ABCI: [料金（2026 年度）](https://abci.ai/ja/how_to_use/tariffs.html) / [ご利用の流れ](https://abci.ai/ja/how_to_use/) /
  [ポイント申請受付](https://abci.ai/news/2025/01/23/ja_news_Point_Application.html) /
  [ABCI 3.0 システム概要](https://docs.abci.ai/v3/ja/system-overview/) / [ジョブ実行](https://docs.abci.ai/v3/ja/job-execution/)
- 富岳: [試行課題（随時募集）](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_trial) /
  [有償課題の利用料金](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_price) /
  [富岳（R-CCS）](https://www.r-ccs.riken.jp/fugaku/)
- TSUBAME4.0: [利用料の概略](https://www.t4.cii.isct.ac.jp/fare_overview)
- 東大: [通常利用（トライアル）負担金](https://www.cc.u-tokyo.ac.jp/guide/application/charge_trial.php) /
  [企業利用（トライアル）](https://www.cc.u-tokyo.ac.jp/guide/trial/company.php)
- FOCUS: [利用料金](https://www.j-focus.or.jp/focus/fee.html) / [利用形態](https://www.j-focus.or.jp/focus/form.html) /
  [試行利用（無償）](https://www.j-focus.or.jp/focus/free-trial.html) /
  [システム概要（利用案内）](https://www.j-focus.jp/user_guide/ug0001000000/)
- Hetzner: [AX162 プレスリリース](https://www.hetzner.com/pressroom/new-ax162/) /
  [AX162-R 製品ページ](https://www.hetzner.com/dedicated-rootserver/ax162-r/) /
  [2026-06 値上げまとめ（Northflank）](https://northflank.com/blog/hetzner-cloud-server-price-increases) /
  [Spare Cores CCX63](https://sparecores.com/server/hcloud/ccx63)
- 為替: [Trading Economics USD/JPY](https://tradingeconomics.com/japan/currency)（2026-07-03: 161.27）
