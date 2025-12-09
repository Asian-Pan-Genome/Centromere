import numpy as np
from sklearn.metrics import silhouette_score
import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import networkx as nx
from community import community_louvain 
import igraph as ig
import leidenalg  
import random

def build_graph(infile,n):
    G = nx.Graph()
    used = set()
    with open(infile, "r") as inf1:
        for line in inf1:
            tokens = line.strip().split('\t')
            i = tokens[0]
            j = tokens[1]
            identity = round(float(tokens[2]),2)
            G.add_edge(i, j, weight=identity * 100)
            used.add(i)
            used.add(j)
    for t in range(0,n):
        if str(t) not in used:
            G.add_node(str(t))
    return G

def louvain_clustering(G):
    partition = community_louvain.best_partition(G)
    modularity = community_louvain.modularity(partition, G)
    return partition, modularity

def get_silhouette_score(G, partition):
    labels = list(partition.values())
    #X = nx.to_numpy_matrix(G, weight="weight")
    X = nx.to_numpy_array(G, weight="weight")
    silhouette = silhouette_score(X, labels, metric="euclidean")
    return silhouette


def leiden_clustering(G):
    ig_G = ig.Graph.TupleList(G.edges(data=False))
    partition = leidenalg.find_partition(ig_G, leidenalg.ModularityVertexPartition)
    return {v['name']: part for v, part in zip(ig_G.vs, partition.membership)}

def plot_network_with_clusters(G, partition, outdir,title=""):
    pos = nx.spring_layout(G, seed=42)  # 使用 spring_layout 布局
    plt.figure(figsize=(10, 10))
    colors = [partition.get(node, random.randint(0, 10)) for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_color=colors, cmap=plt.cm.tab20, node_size=50)
    nx.draw_networkx_edges(G, pos, alpha=0.3)
    plt.title(title)
    plt.savefig(os.path.join(outdir,f"{title}.png"), dpi=300)
    #plt.show()

def main(threshold, n):
    os.makedirs("cluster_out", exist_ok=True)

    #threshold = 87
    G = build_graph(f"0.{threshold}.pairwise.txt", n)

    # Louvain clustering
    louvain_partition, louvain_modularity = louvain_clustering(G)
    #plot_network_with_clusters(G, louvain_partition, "cluster_out", title=f"Louvain_{threshold}_threthold")
    with open(os.path.join("cluster_out", f"Louvain_{threshold}_cluster.out"), "w") as outf:
        for node, cluster in louvain_partition.items():
            #print(f"Node {node}: Cluster {cluster}")
            outf.write(f"{node}\t{cluster}\n")
    
    #sscore = get_silhouette_score(G, louvain_partition)
    sscore = 0
    with open(os.path.join("cluster_out", f"Louvain_{threshold}_cluster.stat"), "w") as statout:
        statout.write(f"{threshold}\t{louvain_modularity}\t{sscore}\n")
    # Leiden clustering
    #leiden_partition = leiden_clustering(G)
    #plot_network_with_clusters(G, leiden_partition, "cluster_out", title=f"Leiden_{threshold}_threthold")
    #with open(os.path.join("cluster_out", f"Leiden_{threshold}_cluster.out"), "w") as outf:
    #    for node, cluster in leiden_partition.items():
    #        outf.write(f"{node}\t{cluster}\n")


if __name__ == "__main__":
    threshold=sys.argv[1]
    n = int(sys.argv[2]) #clusternum
    main(threshold, n)
