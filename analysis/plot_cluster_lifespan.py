import matplotlib.pyplot as plt
from collections import defaultdict

from modules.timeseries import analyze_timeseries
from utils.common_utils import get_sorted_output_xml_files, build_step_to_time


def plot_cluster_lifespan(step_clusters_with_global_ids, step_to_time=None):
    """クラスタの寿命（存続ステップ数）をラベルごとに集計・可視化

    Args:
        step_clusters_with_global_ids (dict): {step: DataFrame}, 各DataFrameに 'global_cluster_ID' 列がある前提
        step_to_time (dict): {step: 時間値} ステップ番号と時間の対応辞書
    """
    # クラスタごとに出現ステップを記録
    cluster_lifespan = defaultdict(list)  # {label: [lifespan, ...]}
    cluster_steps = defaultdict(lambda: defaultdict(list))  # {label: {gid: [step, ...]}}

    for step, df in step_clusters_with_global_ids.items():
        for _, row in df.iterrows():
            gid = row['global_cluster_ID']
            label = row['label']
            if label == "unknown":
                continue  # unknown labelはスキップ
            cluster_steps[label][gid].append(step)

    for label, gid_dict in cluster_steps.items():
        for gid, steps in gid_dict.items():
            if step_to_time:
                times = [step_to_time[s] for s in steps if s in step_to_time]
                if times:
                    lifespan = max(times) - min(times)  # 存続時間
                    cluster_lifespan[label].append(lifespan)
            else:
                lifespan = max(steps) - min(steps)  # 存続ステップ数
                cluster_lifespan[label].append(lifespan)

    # 可視化
    labels = sorted(cluster_lifespan.keys())
    data = [cluster_lifespan[label] for label in labels if cluster_lifespan[label]]

    plt.figure(figsize=(10, 6))
    plt.boxplot(data, labels=labels, patch_artist=True)
    plt.title('Cluster Lifespan by Label')
    plt.xlabel('Label')
    plt.ylabel('Lifespan (steps)')
    plt.grid(True)
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

    # ステップから時間へのマッピングを構築
    step_to_time = build_step_to_time(output_xml_list, dir_path)

    # 結果のplot
    plot_cluster_lifespan(step_clusters_with_global_ids, step_to_time)
