# クラスターの時系列変化をツリー構造で可視化

import networkx as nx
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

from modules.timeseries import analyze_timeseries
from utils.common_utils import get_sorted_output_xml_files, build_cluster_graph_from_global_ids, get_ancestor_nodes

def visualize_cluster_graph(G):
    """
    クラスタ時系列グラフを可視化
    G: networkx.DiGraph
    ノード属性: size, label, fibro, macro, cap
    """


    label_color_map = {
        "autocrine": "skyblue",
        "paracrine": "lightcoral",
        "noise": "gray"
    }

    steps = sorted(set([n[0] for n in G.nodes()]))
    nodes_by_step = {step: [] for step in steps}
    for n in G.nodes():
        nodes_by_step[n[0]].append(n)

    pos = {}
    for step in steps:
        nodes = nodes_by_step[step]
        n_nodes = len(nodes)
        x_positions = [i / (n_nodes - 1) if n_nodes > 1 else 0.5 for i in range(n_nodes)]
        for i, node in enumerate(nodes):
            pos[node] = (x_positions[i], -step)

    node_sizes = [G.nodes[n]["size"] * 10 for n in G.nodes()]
    node_colors = [label_color_map.get(G.nodes[n]["label"], "gray") for n in G.nodes()]
    labels = {
        node: (
            f"{node[1]}\n{G.nodes[node]['size']}({G.nodes[node]['fibro']}, {G.nodes[node]['macro']}, {G.nodes[node]['cap']})"
        )
        for node in G.nodes()
    }

    plt.figure(figsize=(12, 8))
    nx.draw(
        G, pos,
        with_labels=True,
        labels=labels,
        node_size=node_sizes,
        node_color=node_colors,
        edge_color="gray",
        arrows=True,
        arrowstyle='-|>',
        arrowsize=12,
        font_size=10
    )

    legend_elements = [Patch(facecolor=color, label=label) for label, color in label_color_map.items()]
    plt.legend(handles=legend_elements, loc="upper left")
    plt.title("Cluster Time-Series Graph (Step-aligned layout)")
    plt.axis('off')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # 出力XMLファイルのリストを取得
    dir_path = input("解析対象のoutput.xmlがあるディレクトリパスを入力してください: ").strip()
    output_xml_list = get_sorted_output_xml_files(dir_path)

    # ステップ間隔を入力
    try:
        step_interval = int(input("解析に使用するステップ間隔を入力してください（例: 10）: "))
        if step_interval < 1:
            raise ValueError("ステップ間隔は1以上の整数にしてください。")
    except ValueError as e:
        print(f"入力エラー: {e}")
        exit(1)


    step_clusters_with_global_ids = analyze_timeseries(output_xml_list, step_interval)
    G = build_cluster_graph_from_global_ids(step_clusters_with_global_ids, step_interval=1)

    print("=== Step List ===")
    print(sorted(step_clusters_with_global_ids.keys()))

    print("=== Global Cluster ID Presence ===")
    for step, df in step_clusters_with_global_ids.items():
        if df["global_cluster_ID"].isnull().any():
            print(f"[Warning] Null global_cluster_ID at step {step}")

    # === 全クラスタ描画 ===
    visualize_cluster_graph(G)

    # === 再描画オプション: サイズ上位クラスタの祖先のみ ===
    try:
        top_n = int(input("再描画する上位クラスタ数（最終ステップのサイズ順）を入力してください（例: 5）: "))
        if top_n < 1:
            raise ValueError("1以上を指定してください。")
    except ValueError as e:
        print(f"入力エラー: {e}")
        exit(1)

    final_step = max(step_clusters_with_global_ids.keys())
    df_final = step_clusters_with_global_ids[final_step]
    top_clusters = df_final.sort_values("cell_counts", ascending=False).head(top_n)

    # 対象ノード（step, global_id）
    target_nodes = [(final_step, gid) for gid in top_clusters["global_cluster_ID"]]
    ancestor_nodes = get_ancestor_nodes(G, target_nodes)

    # サブグラフを作成
    G_sub = G.subgraph(ancestor_nodes).copy()

    print(f"[Subgraph] Nodes: {G_sub.number_of_nodes()}, Edges: {G_sub.number_of_edges()}")
    visualize_cluster_graph(G_sub)

 
