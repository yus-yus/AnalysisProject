# 時系列解析メインスクリプト
from utils.common_utils import extract_clusters_per_step, assign_global_cluster_ids_edge_based
from modules.snapshot import analyze_snapshot


def analyze_timeseries(output_xml_list, step_interval=1):
    df_by_step = {}
    for i, f in enumerate(range(len(output_xml_list))):
        if i % step_interval != 0:
            continue
        xml_path = output_xml_list[i]
        df_by_step[i] = analyze_snapshot(xml_path)

    step_clusters = extract_clusters_per_step(df_by_step)
    step_clusters_with_global_ids = assign_global_cluster_ids_edge_based(step_clusters, iou_threshold=0.1)

    return step_clusters_with_global_ids

    