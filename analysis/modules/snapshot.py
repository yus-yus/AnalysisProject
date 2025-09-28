import os

from utils.common_utils import load_config, load_cell_dataframe
from utils.clustering_utils import calculate_neighbors, cluster_cells, label_clusters, merge_clusters, count_capillaries


def analyze_snapshot(xml_path, config_ini_path):
    """スナップショットを解析する関数

    Args:
        xml_path (str): XMLファイルのパス
        config_ini_path (str): 設定ファイルのパス

    Returns:
        DataFrame: 解析結果を含むDataFrame
    """
    if config_ini_path is None:
        # このファイル（snapshot.py）から見たconfigの絶対パスを自動設定
        config_ini_path = os.path.join(os.path.dirname(__file__), "..", "config", "config_analysis.ini")
        config_ini_path = os.path.abspath(config_ini_path)
    config = load_config(config_ini_path)

    df_cell = load_cell_dataframe(xml_path)
    df_cell = calculate_neighbors(df_cell, config)
    df_cell = cluster_cells(df_cell, config)
    df_cell = label_clusters(df_cell, config)
    df_cell = merge_clusters(df_cell, config)
    df_cell = count_capillaries(df_cell)

    # 必要なカラムだけ返す
    return df_cell[[
        "ID", "position_x", "position_y", "status", "dead", "new_cluster",
        "fibroblast_count", "macrophage_count", "capillary_count_new"
    ]]
