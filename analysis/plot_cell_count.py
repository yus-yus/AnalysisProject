from collections import defaultdict
import matplotlib.pyplot as plt

from modules.timeseries import analyze_timeseries
from utils.common_utils import build_step_to_time, get_sorted_output_xml_files

def plot_cell_count_timecourse(step_clusters_with_global_ids, step_to_time):
	"""
	クラスタごとの細胞数の時系列を可視化
	step_clusters_with_global_ids: {step: DataFrame}, 各DataFrameに 'label' と 'cell_counts' 列がある前提
	step_to_time: {step: 時間値}
	"""
	label_to_timecourse = defaultdict(list)
	steps = sorted(step_clusters_with_global_ids.keys())

	for step in steps:
		time = step_to_time.get(step)
		if time is None:
			continue
		df = step_clusters_with_global_ids[step]
		grouped = df.groupby('label')['cell_counts'].sum()
		for label in grouped.index:
			if label != 'unknown':
				label_to_timecourse[label].append((time, grouped[label]))

	plt.figure(figsize=(10, 6))
	for label, values in label_to_timecourse.items():
		values.sort()
		x_vals = [t for t, _ in values]
		y_vals = [v for _, v in values]
		plt.plot(x_vals, y_vals, label=str(label), marker='o')

	plt.title('Timecourse of Cell Count by Label')
	plt.xlabel('Time (hours)')
	plt.ylabel('Cell Count')
	plt.grid(True)
	plt.legend(title='Label')
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
    step_to_time = build_step_to_time(output_xml_list, dir_path)

    # print(step_to_time)

    plot_cell_count_timecourse(step_clusters_with_global_ids, step_to_time)  