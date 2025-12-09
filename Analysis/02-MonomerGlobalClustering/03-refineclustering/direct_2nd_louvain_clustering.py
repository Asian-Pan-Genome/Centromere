import sys
import os
from collections import defaultdict
import networkx as nx
from community import community_louvain


def load_inner_similarity_file(infile):
    pairwise_centroid_identity = []
    with open(infile, 'r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            pairwise_centroid_identity.append([tokens[0], int(tokens[1]), tokens[2], int(tokens[3]), float(tokens[4])])
    return pairwise_centroid_identity


def build_graph(pairwise_centroid_identity):

    G = nx.Graph()
    for i, isize, j, jsize, sim in pairwise_centroid_identity:
        if sim >= 0.95:
            G.add_edge(i, j, weight=sim * 100)
    return G

def louvain_clustering(G):
    partition = community_louvain.best_partition(G)
    return partition


def main(infile, nid, outfile):
    #infile = os.path.join('tmp', f"{nid}.inner.txt")
    pairwise_centroid_identity = load_inner_similarity_file(infile)
    G = build_graph(pairwise_centroid_identity)
    partition=louvain_clustering(G)
    p_dict = defaultdict(list)
    for node, part in partition.items():
        p_dict[part].append(node)

    with open(outfile, 'w') as outf:
        for p, t in p_dict.items():
            outf.write(f"{nid}.{p}\t{len(t)}\t{','.join(t)}\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3]) 
