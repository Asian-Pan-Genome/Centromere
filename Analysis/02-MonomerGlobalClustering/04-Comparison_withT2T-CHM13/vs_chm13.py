import os
from sklearn.metrics import jaccard_score
from sklearn.metrics import f1_score
import sys
import pandas as pd

def stat_chm13_anno():
    mn_anno = {}
    with open("CHM13.anno.txt", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            mn = tokens[0]
            anno = tokens[0].split('::')[0]
            if 'd' not in anno:
            #if 'L' in anno:
            #if ('L' not in anno) and ('d' not in anno):
                mn_anno[mn] = anno
    return mn_anno

def stat_chm13_cluster():
    mn_cluster = {}
    with open("CHM13.cluster.txt", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            mn_cluster[tokens[0]] = tokens[1]
    return mn_cluster

def stat_linked(infile):
    i = 0
    oid_mid = {}
    with open(infile, 'r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            for j in tokens[2].split(','):
                oid_mid[j] = i
            i += 1
    return oid_mid


def rename_cluster(mn_cluster, oid_mid):
    mn_mcluster = {}
    for mn, cluster in mn_cluster.items():
        mid = oid_mid[cluster]
        mn_mcluster[mn] = mid
    return mn_mcluster


def compare_input(mn_anno, mn_mcluster, outfile):
    with open(outfile, 'w') as outf:
        for mn, anno in mn_anno.items():
            outf.write(f"{mn}\t{anno}\t{mn_mcluster[mn]}\n")
        
def get_jaccard_index(outfile):
    df = pd.read_csv(outfile, sep="\t", names=['mn', 'anno', 'cluster'])
    sorted_df = df.sort_values(by = ['anno', 'cluster'], ascending=[True, True])
    unique_anno = sorted_df['anno'].unique()
    unique_cluster = sorted_df['cluster'].unique()

    anno_mapping = {value: idx for idx, value in enumerate(unique_anno)}
    cluster_mapping = {value: idx for idx, value in enumerate(unique_cluster)}

    sorted_df['anno_mapped'] = sorted_df['anno'].map(anno_mapping)
    sorted_df['cluster_mapped'] = sorted_df['cluster'].map(cluster_mapping)

    anno_labels = sorted_df['anno_mapped'].tolist()
    cluster_labels = sorted_df['cluster_mapped'].tolist()
    
    #jaccard_macro = jaccard_score(anno_labels, cluster_labels, average="macro")
    #jaccard_micro = jaccard_score(anno_labels, cluster_labels, average="micro")

    #print("Jaccard Index (macro):", jaccard_macro)
    #print("Jaccard Index (micro):", jaccard_micro)

    f1 = f1_score(anno_labels, cluster_labels, average='micro')  # micro 平均用于多类数据
    print("F1 Score:", f1)

def calculate_jaccard(outfile):
    df = pd.read_csv(outfile, sep="\t", names=['mn', 'anno', 'cluster'])
    #sorted_df = df.sort_values(by = ['anno', 'cluster'], ascending=[True, True])
    jaccard_scores = 0
    max_possible_score = 0
    
    anno_group = df.groupby('anno')['mn'].apply(set).to_dict()  # A dictionary: anno_class -> set of mn
    cluster_group = df.groupby('cluster')['mn'].apply(set).to_dict()  # A dictionary: cluster_class -> set of mn

    jaccard_scores = 0
    max_possible_score = 0
    
    # Iterate over each combination of 'anno_class' and 'cluster_class'
    for anno_class, anno_set in anno_group.items():
        for cluster_class, cluster_set in cluster_group.items():
            overlap = len(anno_set.intersection(cluster_set))

            jaccard = overlap / (len(anno_set) + len(cluster_set) - overlap)
            
            anno_size = len(anno_set)
            cluster_size = len(cluster_set)
            weight_factor = (overlap / anno_size) * (overlap / cluster_size)  if anno_size and cluster_size else 0
            
            weighted_jaccard = jaccard * weight_factor
            jaccard_scores += weighted_jaccard
            
            max_possible_score += weight_factor
    
    if max_possible_score > 0:
        final_score = jaccard_scores / max_possible_score
    else:
        final_score = 0
    
    return final_score


def main():
    mn_anno = stat_chm13_anno()
    mn_cluster = stat_chm13_cluster()

    infile = sys.argv[1]
    outfile = sys.argv[2]

    oid_mid = stat_linked(infile)
    mn_mcluster = rename_cluster(mn_cluster, oid_mid)
    #get_jaccard_index(mn_anno, mn_cluster)
    compare_input(mn_anno, mn_mcluster, outfile)
    score = calculate_jaccard(outfile)
    print(score)


if __name__ == "__main__":
    main()
