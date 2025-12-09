import pandas as pd
import os

def load_cluster_file(file_path):
    df = pd.read_csv(file_path, sep="\t", header=None, names=["node", "cluster_id"])
    return df


def concatenate_clusters(thresholds, output_file):
    result_df = pd.DataFrame()

    all_nodes = set()

    for threshold in thresholds:
        cluster_file = f"Louvain_{threshold}_cluster.out"
        if os.path.exists(cluster_file):
            cluster_df = load_cluster_file(cluster_file)
            s = cluster_df.groupby('cluster_id', as_index=False).agg({'node':'count'})
            d = dict(zip(s['cluster_id'], s['node']))
            sorted_d = dict(sorted(d.items(), key=lambda item: item[1], reverse=True))
            r = {j:i for i, j  in enumerate(sorted_d.keys())}
            
            cluster_df = cluster_df.set_index("node")
            cluster_df = cluster_df.rename(columns={"cluster_id": f"cluster_{threshold}"})
            
            all_nodes.update(cluster_df.index)

            if result_df.empty:
                result_df = cluster_df
                result_df[f"index_{threshold}"] =  result_df[f"cluster_{threshold}"].map(r)
            else:
                result_df = result_df.join(cluster_df, how="outer")
                result_df[f"index_{threshold}"] =  result_df[f"cluster_{threshold}"].map(r)
    
    result_df = result_df.loc[sorted(all_nodes)]

    sorted_cols = [f"index_{t}" for t in thresholds]
    print(sorted_cols)
    sorted_result_df =  result_df.sort_values(by=sorted_cols, ascending=[True]*len(sorted_cols))
    sorted_result_df.to_csv(output_file, sep="\t", header=True, index=True)
    return sorted_result_df


def main():
    thresholds = list(range(87, 100))
    output_file = "Louvain_87-99_concatenated_clusters_sorted.txt"

    
    result_df = concatenate_clusters(thresholds, output_file)
    
if __name__ == "__main__":
    main()
