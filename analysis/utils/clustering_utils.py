# cluster_utils.py
import pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import cdist

# DBSCANクラスタリング
def cluster_cells(df_cell, config):
    eps = int(config.get('DBSCAN', 'eps'))
    min_samples = int(config.get('DBSCAN', 'min_samples'))
    df_cell_temp = df_cell[(df_cell['cell_type'].isin(['macrophage', 'fibroblast'])) & (df_cell['dead'] == False)]
    X = df_cell_temp[['position_x', 'position_y']].to_numpy()
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    clusters = dbscan.fit_predict(X)
    df_cell['cluster'] = -1
    df_cell.loc[df_cell_temp.index, 'cluster'] = clusters
    df_cell['cluster'] = df_cell['cluster'].apply(lambda x: int(x) if not pd.isna(x) and x >= 0 else -1)
    return df_cell

# クラスタラベル付け処理
def label_clusters(df_cell, config):
    cluster_labels = []
    for cluster_id in df_cell['cluster'].unique():
        if cluster_id == -1:
            continue
        macrophage_count = df_cell[(df_cell['cluster'] == cluster_id) & (df_cell['cell_type'] == 'macrophage')].shape[0]
        fibroblast_count = df_cell[(df_cell['cluster'] == cluster_id) & (df_cell['cell_type'] == 'fibroblast')].shape[0]
        if macrophage_count >= 1 and fibroblast_count >= 1:
            label = 'paracrine'
        elif fibroblast_count >= 1:
            label = 'autocrine'
        else:
            label = 'unknown'
        cluster_labels.append((label, macrophage_count, fibroblast_count, cluster_id))
    for label, macrophage_count, fibroblast_count, cluster_id in cluster_labels:
        df_cell.loc[df_cell['cluster'] == cluster_id, 'status'] = label
        df_cell.loc[df_cell['cluster'] == cluster_id, 'macrophage_count'] = macrophage_count
        df_cell.loc[df_cell['cluster'] == cluster_id, 'fibroblast_count'] = fibroblast_count
    # noiseにも必ずカラムを追加
    df_cell.loc[df_cell['cluster'] == -1, 'status'] = 'noise'
    df_cell['macrophage_count'] = df_cell.get('macrophage_count', 0)
    df_cell['fibroblast_count'] = df_cell.get('fibroblast_count', 0)
    return df_cell

# 近傍セルの計算
def calculate_neighbors(df_cell, config):
    distance_threshold = int(config.get('NEIGHBOR', 'distance_threshold'))
    positions = df_cell[['position_x', 'position_y']].to_numpy()
    distances = cdist(positions, positions)
    neighbor_mask = (distances <= distance_threshold) & (distances > 0)
    neighbor_ids = [
        df_cell.loc[neighbor_mask[i], 'ID'].tolist() for i in range(len(df_cell))
    ]
    df_cell['neighbor'] = [','.join(map(str, ids)) for ids in neighbor_ids]
    return df_cell

# クラスタ統合処理（Union-Findによるクラスタ結合）
def merge_clusters(df_cell, config):
    n = int(config.get('CLUSTER', 'min_common_cells'))
    class UnionFind:
        def __init__(self, n):
            self.parent = list(range(n))
            self.rank = [0] * n
        def find(self, x):
            if self.parent[x] != x:
                self.parent[x] = self.find(self.parent[x])
            return self.parent[x]
        def union(self, x, y):
            root_x = self.find(x)
            root_y = self.find(y)
            if root_x != root_y:
                if self.rank[root_x] > self.rank[root_y]:
                    self.parent[root_y] = root_x
                elif self.rank[root_x] < self.rank[root_y]:
                    self.parent[root_x] = root_y
                else:
                    self.parent[root_y] = root_x
                    self.rank[root_x] += 1
    unique_clusters = df_cell['cluster'].dropna().unique()
    cluster_map = {cluster: i for i, cluster in enumerate(unique_clusters)}
    reverse_cluster_map = {i: cluster for cluster, i in cluster_map.items()}
    uf = UnionFind(len(unique_clusters))
    for i, row in df_cell.iterrows():
        neighbors = row.get('neighbor', None)
        if pd.isna(neighbors):
            continue
        neighbors = [int(n) for n in str(neighbors).split(',') if n.strip() != '']
        cluster1 = row['cluster']
        if cluster1 == -1:
            continue
        status1 = row.get('status', None)
        for neighbor in neighbors:
            neighbor_row = df_cell.loc[df_cell['ID'] == neighbor]
            if neighbor_row.empty:
                continue
            cluster2 = neighbor_row['cluster'].values[0]
            if cluster2 == -1:
                continue
            status2 = neighbor_row['status'].values[0] if 'status' in neighbor_row else None
            common_cells = len(df_cell[(df_cell['cluster'] == cluster1) & (df_cell['ID'].isin(neighbors))])
            if status1 == status2 and common_cells >= n:
                uf.union(cluster_map[cluster1], cluster_map[cluster2])
    new_cluster_labels = {cluster: uf.find(cluster_map[cluster]) for cluster in unique_clusters}
    new_cluster_labels = {cluster: reverse_cluster_map[root] for cluster, root in new_cluster_labels.items()}
    df_cell['new_cluster'] = df_cell['cluster'].map(new_cluster_labels)
    return df_cell

# クラスタごとにcapillary数をカウント
def count_capillaries(df_cell):
    capillary_counts = {}
    cluster_col = 'new_cluster' if 'new_cluster' in df_cell.columns else 'cluster'
    for cluster_id in df_cell[cluster_col].unique():
        if cluster_id == -1:
            continue
        cluster_cells = df_cell[df_cell[cluster_col] == cluster_id]['ID'].tolist()
        capillary_ids = set()
        for cell_id in cluster_cells:
            neighbors = df_cell.loc[df_cell['ID'] == cell_id, 'neighbor'].iloc[0]
            if pd.isna(neighbors):
                continue
            neighbors = [int(n) for n in str(neighbors).split(',') if n.strip() != '']
            for neighbor_id in neighbors:
                neighbor_cell_type = df_cell.loc[df_cell['ID'] == neighbor_id, 'cell_type'].values[0]
                if neighbor_cell_type == 'capillary':
                    capillary_ids.add(neighbor_id)
        capillary_counts[cluster_id] = len(capillary_ids)
    # capillary_count_newカラムで追加
    df_cell['capillary_count_new'] = df_cell[cluster_col].map(capillary_counts)
    return df_cell
