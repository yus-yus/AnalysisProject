import configparser
import os
from pathlib import Path
import pcdl
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path
import re
from collections import defaultdict
import pandas as pd

# 設定ファイルを読み込むユーティリティ関数
def load_config(config_ini_path):
    config = configparser.ConfigParser()
    config.read(config_ini_path, encoding='utf-8')
    return config

# セルデータの読み込み
def load_cell_dataframe(xml_path):
    final_xml_path = Path(xml_path)
    df_cell = pcdl.TimeStep(final_xml_path.as_posix(), physiboss=False, verbose=False, microenv=False).get_cell_df()
    df_cell["status"] = "noise"
    df_cell = df_cell.reset_index()
    return df_cell

# svgファイルから時間情報を抽出する関数
def extract_time_from_svg(svg_path):
    try:
        with open(svg_path, "r", encoding="utf-8") as f:
            text = f.read()

        # 「Current time: ...」形式を優先して抽出
        match = re.search(r"Current time:\s+(\d+)\s+days,\s+(\d+)\s+hours,\s+and\s+([\d\.]+)\s+minutes", text)
        if match:
            days = int(match.group(1))
            hours = int(match.group(2))
            minutes = float(match.group(3))
            total_hours = days * 24 + hours + minutes / 60
            return total_hours

        # fallback: 「0 days, 8 hours, 44 minutes, and 34.1991 seconds」形式にも対応
        match = re.search(r"(\d+)\s+days,\s+(\d+)\s+hours,\s+(\d+)\s+minutes,\s+and\s+([\d\.]+)\s+seconds", text)
        if match:
            days = int(match.group(1))
            hours = int(match.group(2))
            minutes = int(match.group(3))
            seconds = float(match.group(4))
            total_hours = days * 24 + hours + minutes / 60 + seconds / 3600
            return total_hours

    except Exception as e:
        print(f"[SVG read error] {svg_path}: {e}")

    return None

# SVGファイルから時間情報を抽出し、step_to_time辞書を作成する関数
def build_step_to_time(xml_files, output_dir):
    """
    xml_files: List[Path] - output*.xmlファイルのリスト
    output_dir: str or Path - snapshot*.svgがあるディレクトリ
    extract_time_from_svg: 関数(svg_path) -> float or None
    return: {step_index: time_val}
    """
    step_to_time = {}
    for step_index, xml_path in enumerate(xml_files):
        svg_name = xml_path.name.replace("output", "snapshot").replace(".xml", ".svg")
        svg_path = os.path.join(output_dir, svg_name)
        time_val = extract_time_from_svg(svg_path)
        if time_val is not None:
            step_to_time[step_index] = time_val
        else:
            print(f"[Warning] SVGから時間抽出できず: {svg_path}")
    return step_to_time

# outputXXXX.xml ファイルをステップ順に取得
def get_sorted_output_xml_files(target_dir_path):
    target_dir = Path(target_dir_path)
    xml_files = list(target_dir.glob("*.xml"))
    numbered_files = []

    for f in xml_files:
        match = re.match(r"output(\d+)\.xml$", f.name)
        if match:
            num = int(match.group(1))
            numbered_files.append((num, f))

    numbered_files.sort(key=lambda x: x[0])
    return [f for _, f in numbered_files]

# 各ステップのDataFrameを辞書形式で取得
def extract_clusters_per_step(df_by_step):
    step_clusters = {}
    for step, df in df_by_step.items():
        clustered = df[df["new_cluster"] != -1]
        if clustered.empty:
            step_clusters[step] = pd.DataFrame(columns=[
                "cluster_ID", "cell_ids", "position_x", "position_y", "label", "cell_counts",
                "fibroblast_count", "macrophage_count", "capillary_count"
            ])
            continue

        cluster_summary = (
            clustered.groupby("new_cluster")
            .agg({
                "ID": lambda x: list(x),
                "position_x": "mean",
                "position_y": "mean",
                "status": "first",
                "fibroblast_count": "max",
                "macrophage_count": "max",
                "capillary_count_new": "max"
            })
            .rename(columns={
                "ID": "cell_ids",
                "status": "label",
                "fibroblast_count": "fibroblast_count",
                "macrophage_count": "macrophage_count",
                "capillary_count_new": "capillary_count"
            })
            .reset_index()
            .rename(columns={"new_cluster": "cluster_ID"})
        )

        cluster_summary["cell_counts"] = cluster_summary["cell_ids"].apply(len)
        step_clusters[step] = cluster_summary

    return step_clusters


# エッジベースでグローバルクラスタIDを割当て
def assign_global_cluster_ids_edge_based(step_clusters, iou_threshold=0.1):
    global_id_counter = 0
    global_clusters = defaultdict(dict)
    prev_clusters = None
    prev_global_map = {}

    for step in sorted(step_clusters.keys()):
        current = step_clusters[step]
        global_map = {}

        print(f"\n=== Step {step} ===")
        if prev_clusters is None:
            # 初回ステップは単純に新規ID割当
            print("Initial step: assigning new global IDs")
            for _, row in current.iterrows():
                global_map[row["cluster_ID"]] = global_id_counter
                print(f"Assign global ID {global_id_counter} to cluster_ID {row['cluster_ID']}")
                global_id_counter += 1
        else:
            # 1. 前後クラスタ間のIoUベースのエッジ（つながり）を全部収集
            edges = []
            for _, row_curr in current.iterrows():
                curr_cid = row_curr["cluster_ID"]
                curr_cells = set(row_curr["cell_ids"])
                for _, row_prev in prev_clusters.iterrows():
                    prev_cid = row_prev["cluster_ID"]
                    prev_cells = set(row_prev["cell_ids"])
                    intersection = len(curr_cells & prev_cells)
                    union = len(curr_cells | prev_cells)
                    iou = intersection / union if union > 0 else 0
                    if iou >= iou_threshold:
                        edges.append((prev_cid, curr_cid, iou))
                        print(f"Edge detected: Prev cluster {prev_cid} -> Curr cluster {curr_cid} | IoU = {iou:.3f}")

            # 2. 辞書構造に変換（親→子、子→親）
            prev_to_curr = defaultdict(list)
            curr_to_prev = defaultdict(list)
            for p, c, iou in edges:
                prev_to_curr[p].append((c, iou))
                curr_to_prev[c].append((p, iou))

            print(f"Total edges found: {len(edges)}")
            print(f"Parent to current mapping: {dict(prev_to_curr)}")
            print(f"Current to parent mapping: {dict(curr_to_prev)}")

            # 3. currentの各クラスタにID割当て
            for _, row_curr in current.iterrows():
                curr_cid = row_curr["cluster_ID"]
                parents = curr_to_prev.get(curr_cid, [])

                if len(parents) == 0:
                    # 親なし（新規出現）
                    global_map[curr_cid] = global_id_counter
                    print(f"New cluster {curr_cid} has no parent, assigned new global ID {global_id_counter}")
                    global_id_counter += 1
                elif len(parents) == 1:
                    # 親1つ → その親のglobal IDを継続
                    parent_cid = parents[0][0]
                    parent_gid = prev_global_map.get(parent_cid, None)
                    global_map[curr_cid] = parent_gid
                    print(f"Cluster {curr_cid} has one parent {parent_cid}, continuing global ID {parent_gid}")
                else:
                    # 複数親（統合） → とりあえず最もIoUが高い親のIDを継続
                    parents.sort(key=lambda x: x[1], reverse=True)
                    parent_cid = parents[0][0]
                    parent_gid = prev_global_map.get(parent_cid, None)
                    global_map[curr_cid] = parent_gid
                    print(f"Cluster {curr_cid} has multiple parents { [p[0] for p in parents] }, continuing global ID {parent_gid} from parent {parent_cid}")

            # 4. 親クラスタで繋がらないクラスタは新規ID割当て（未割当チェック）
            assigned_cids = set(global_map.keys())
            for _, row_curr in current.iterrows():
                cid = row_curr["cluster_ID"]
                if cid not in assigned_cids:
                    global_map[cid] = global_id_counter
                    print(f"Unassigned cluster {cid} assigned new global ID {global_id_counter}")
                    global_id_counter += 1

        print(f"Global map for step {step}: {global_map}")

        global_clusters[step] = global_map
        prev_clusters = current
        prev_global_map = global_map

    # 最後に各stepのDataFrameにglobal_cluster_ID列を追加
    for step, df in step_clusters.items():
        df["global_cluster_ID"] = df["cluster_ID"].map(global_clusters[step])
    return step_clusters

# グローバルクラスタIDベースでクラスタ時系列グラフを構築
def build_cluster_graph_from_global_ids(step_clusters, step_interval=1):
    G = nx.DiGraph()
    steps = sorted(step_clusters.keys())

    for step in steps:
        df = step_clusters[step]
        for _, row in df.iterrows():
            node = (step, row["global_cluster_ID"])
            G.add_node(node,
                        size=row["cell_counts"],
                        label=row["label"],
                        fibro=row["fibroblast_count"],
                        macro=row["macrophage_count"],
                        cap=row["capillary_count"])



    for i in range(len(steps) - step_interval):
        step = steps[i]
        next_step = steps[i + step_interval]
        df_curr = step_clusters[step]
        df_next = step_clusters[next_step]

        cluster_map_curr = {row["global_cluster_ID"]: set(row["cell_ids"]) for _, row in df_curr.iterrows()}
        cluster_map_next = {row["global_cluster_ID"]: set(row["cell_ids"]) for _, row in df_next.iterrows()}

        for gid1, cells1 in cluster_map_curr.items():
            for gid2, cells2 in cluster_map_next.items():
                shared = len(cells1 & cells2)
                if shared > 0:
                    node1 = (step, gid1)
                    node2 = (next_step, gid2)
                    G.add_edge(node1, node2, weight=shared)

    print(f"[build_cluster_graph_from_global_ids] Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    return G

def get_ancestor_nodes(G, target_nodes):
    visited = set()
    queue = list(target_nodes)

    while queue:
        node = queue.pop()
        if node in visited:
            continue
        visited.add(node)
        parents = list(G.predecessors(node))
        queue.extend(parents)

    return visited