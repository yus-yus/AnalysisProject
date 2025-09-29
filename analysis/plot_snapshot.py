import matplotlib.pyplot as plt
import pandas as pd
import os

from modules.snapshot import analyze_snapshot
from utils.common_utils import load_cell_dataframe, load_config

def plot_snapshot(df_cell, ax, config_ini_path=None, before_clustering=False):
    """スナップショットの細胞分布をプロット

    Args:
        df_cell (pd.DataFrame): 解析結果を含むDataFrame
        ax: Axes object
        config_ini_path (str): 設定ファイルのパス
        before_clustering (bool): クラスタリング前ならTrue

    Returns:
        ax: Axes object
    """
    # 設定ファイルの読み込み
    if config_ini_path is None:
        config_ini_path = os.path.join(os.path.dirname(__file__), "config", "config_analysis.ini")
        config_ini_path = os.path.abspath(config_ini_path)
    config = load_config(config_ini_path)

    print(config)

    # 色の設定
    cell_colors = {k: v for k, v in config["CELL_COLOR"].items()}
    cmap = plt.get_cmap("tab20")

    # 生きている細胞のみplot
    df_cell = df_cell[df_cell["dead"] == False]

    if before_clustering:
        # クラスタリング前はCELL_COLORで色分け
        for cell_type, color in cell_colors.items():
            sub_df = df_cell[df_cell["cell_type"] == cell_type]
            ax.scatter(
                sub_df["position_x"], sub_df["position_y"],
                label=cell_type, color=color, edgecolor="black", linewidth=0.5, s=15
            )
    else:
        # クラスタリング後はクラスタごとに色分け（適当な色でOK）
        if "new_cluster" in df_cell.columns:
            # capillaryのみ先に表示
            capillary_df = df_cell[df_cell["cell_type"] == "capillary"]
            if not capillary_df.empty:
                ax.scatter(
                    capillary_df["position_x"], capillary_df["position_y"],
                    label="capillary", color=cell_colors.get("capillary", "yellow"),
                    edgecolor="black", linewidth=0.5, s=15
                )
            # noise（capillary以外のみ）
            noise_df = df_cell[(df_cell["cell_type"] != "capillary") & (df_cell["new_cluster"] == -1)]
            if not noise_df.empty:
                ax.scatter(
                    noise_df["position_x"], noise_df["position_y"],
                    label="noise", color="gray", edgecolor="black", linewidth=0.5, s=15
                )
            # クラスタごとに色分け（capillary以外のみ）
            unique_clusters = [c for c in df_cell["new_cluster"].unique() if c != -1]
            for idx, cluster_id in enumerate(unique_clusters):
                sub_df = df_cell[(df_cell["new_cluster"] == cluster_id) & (df_cell["cell_type"] != "capillary")]
                if not sub_df.empty:
                    ax.scatter(
                        sub_df["position_x"], sub_df["position_y"],
                        label=f"cluster{cluster_id}", color=cmap(idx % cmap.N),
                        edgecolor="black", linewidth=0.5, s=15
                    )
        else:
            # クラスタ情報がない場合はグレー
            ax.scatter(
                df_cell["position_x"], df_cell["position_y"],
                label="cell", color="gray", edgecolor="black", linewidth=0.5, s=15
            )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xlim(-1000, 1000)
    ax.set_ylim(-1000, 1000)
    ax.autoscale(False)
    ax.legend()

    return ax


def plot_snapshot_comparison(df_before, df_after, config_ini_path=None):
    """クラスタリング前後のスナップショットを比較プロット

    Args:
        df_before (pd.DataFrame): クラスタリング前のDataFrame
        df_after (pd.DataFrame): クラスタリング後のDataFrame
        config_ini_path (str): 設定ファイルのパス
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    ax1 = plot_snapshot(df_before, axes[0], config_ini_path, before_clustering=True)
    ax1.set_title("Before Clustering")

    ax2 = plot_snapshot(df_after, axes[1], config_ini_path, before_clustering=False)
    ax2.set_title("After Clustering")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # 例: コマンドライン引数やinputでパスを取得
    dir_path = input("解析対象のoutput.xmlがあるディレクトリパスを入力してください: ").strip()
    xml_name = input("描画したいXMLファイル名（未入力ならfinal.xml）: ").strip()
    if not xml_name:
        xml_name = "final.xml"

    xml_path = os.path.join(dir_path, xml_name)
    print(f"使用するXML: {xml_path}")

    # クラスタリング前のデータを読み込み
    df_before = load_cell_dataframe(xml_path)   
    # クラスタリング後のデータを解析
    df_after = analyze_snapshot(xml_path, config_ini_path=None)

    plot_snapshot_comparison(df_before, df_after, config_ini_path=None)