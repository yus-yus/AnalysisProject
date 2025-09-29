# AnalysisProject

## プロジェクトの概要

AnalysisProjectは、細胞シミュレーションデータのクラスタリング・時系列解析・可視化を自動化するPythonプロジェクトです。  
XML形式のシミュレーション出力から細胞のクラスタリングや各種統計量を算出し、可視化まで一連の流れをサポートします。

---

## ディレクトリ構成

```
AnalysisProject/
├── analysis/
│   ├── modules/         # 解析用モジュール（snapshot, timeseriesなど）
│   ├── utils/           # ユーティリティ関数
│   ├── config/          # 設定ファイル（config_analysis.iniなど）
│   ├── plot_cell_count.py
│   ├── plot_cluster_lifespan.py
│   ├── plot_cluster_tree.py
│   ├── plot_snapshot.py
├── examples/            # サンプルデータや出力例（任意）
├── requirements.txt     # 必要なPythonパッケージ
├── README.md            # このファイル
```

---

## 使い方

### 1. 環境構築

```sh
python -m venv .venv
.venv\Scripts\activate  # Windowsの場合
pip install -r requirements.txt
```

### 2. 解析の実行例

#### クラスタリング前後のスナップショット可視化

```sh
python analysis/plot_snapshot.py
```
- プロンプトに従い、解析結果フォルダを入力してください（例: `C:\Users\user\sim_results\run1`）。
- 描画したいXMLファイル名（未入力なら `final.xml`）を入力してください。

#### 細胞数の時系列プロット

```sh
python analysis/plot_cell_count.py
```
- 解析対象ディレクトリとステップ間隔を入力してください。

#### クラスタ寿命の可視化

```sh
python analysis/plot_cluster_lifespan.py
```

---

## 出力例

- クラスタリング前後の細胞分布を並べて可視化した画像
- 細胞数やクラスタ寿命の時系列グラフ
- クラスタツリー構造の可視化

（例：`examples/` フォルダや `results/` フォルダに画像ファイルが保存されます）

---

必要に応じて、設定ファイル（`analysis/config/config_analysis.ini`）を編集してパラメータを調整してください。
