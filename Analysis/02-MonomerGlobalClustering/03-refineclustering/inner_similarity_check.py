import sys
import statistics
import numpy as np
import os
from collections import defaultdict
import networkx as nx
from community import community_louvain


def stat_merge_file(infile):
    nid_mids = {}
    with open(infile, "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            mids = tokens[2].split(',')
            nid_mids[tokens[0]] = mids
    return nid_mids

def load_inner_similarity_file(infile):
    pairwise_centroid_identity = []
    with open(infile, 'r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            pairwise_centroid_identity.append([tokens[0], int(tokens[1]), tokens[2], int(tokens[3]), float(tokens[4])])
    return pairwise_centroid_identity

def load_centroid_size():
    c_size = {}
    with open("centroid.size", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            c_size[tokens[0]] = int(tokens[1])
    #print(list(c_size.keys())[:10])
    return c_size

def estimate_global_identity(pairwise_centroid_identity):
    estimate_global = []
    allcentroids = set()

    intra_class_calculated = defaultdict(bool)

    for i, isize, j, jsize, sim in pairwise_centroid_identity:
        if not intra_class_calculated[i]:
            estimate_global.extend([0.96] * (isize * (isize - 1) // 2))
            intra_class_calculated[i] = True

        if not intra_class_calculated[j]:
            estimate_global.extend([0.96] * (jsize * (jsize - 1) // 2))
            intra_class_calculated[j] = True

        estimate_global.extend([sim] * (isize * jsize))

    estimate_global = np.array(estimate_global)

    min_sim = np.min(estimate_global)
    max_sim = np.max(estimate_global)
    mean_sim = np.mean(estimate_global)
    median_sim = np.median(estimate_global)
    std_dev = np.std(estimate_global)
    q1 = np.percentile(estimate_global, 25)
    q3 = np.percentile(estimate_global, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr

    min_percent = np.sum(estimate_global == min_sim) / len(estimate_global)
    return min_sim, mean_sim, median_sim, max_sim, std_dev, q1, q3, iqr, lower_bound, upper_bound, min_percent 

def inner_similarity(nid_mids, c_size):

    for nid, mids in nid_mids.items():        
        if len(mids) > 2:
            #statout = open(f"tmp/{nid}.stat", 'w')
            #statout.write(f"index\tmin\tmean\tmedian\tmax\tstd_dev\tq1\tq3\tiqr\tlower_bound\tupper_bound\tmin_percent\n")
            tmp = sorted(list(map(int, mids)))
            sorted_mids = list(map(str,tmp))
            #print(mids, sorted_mids)

            pairwise_centroid_identity = []
            for i in range(len(sorted_mids)):
                ref = sorted_mids[i]
                queries = sorted_mids[i+1:]

                #print(ref, queries)
                with open(f"../step-2/tmp/centroid_{ref}.pairwise.txt", "r") as inf:
                    for line in inf:
                        tokens = line.strip().split('\t')
                        if tokens[1] in queries:
                            pairwise_centroid_identity.append([ref, c_size[ref], tokens[1], c_size[tokens[1]], float(tokens[2])])

            #print(pairwise_centroid_identity)
            with open(f"tmp/{nid}.inner.txt", 'w') as outf:
                for info in pairwise_centroid_identity:
                    outf.write(f"{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}\t{info[4]}\n")

            #min_sim, mean_sim, median_sim, max_sim, \
            #std_dev, q1, q3, iqr, lower_bound, upper_bound, min_percent = estimate_global_identity(pairwise_centroid_identity)
            #statout.write(f"{nid}\t{min_sim}\t{mean_sim}\t{median_sim}\t{max_sim}\t{std_dev}\t{q1}\t{q3}\t{iqr}\t{lower_bound}\t{upper_bound}\t{min_percent}\n")
            #statout.close()

def size90(pairwise_centroid_identity):
    all_size = {}
    for i, isize, j, jsize, identity in pairwise_centroid_identity:
        all_size[i] = isize
        all_size[j] = jsize

    sorted_all_size = dict(sorted(all_size.items(), key=lambda item: item[1], reverse=True))
    #print(sorted_all_size)
    large_clusters = []
    total_size = sum(sorted_all_size.values())
    cum = 0
    for k, ksize in sorted_all_size.items():
        cum += ksize
        large_clusters.append(k)
        if cum / total_size >= 0.9:
            break
    #print(large_clusters)
    lsize = {x: all_size[x] for x in large_clusters}
    ssize = {x: all_size[x] for x in all_size if x not in large_clusters}
    return all_size, lsize, ssize
    
def build_graph(lsize, ssize, pairwise_centroid_identity):

    G = nx.Graph()
    updated_large_clusters = set()
    for i, isize, j, jsize, sim in pairwise_centroid_identity:
        if (i in lsize) and (j in lsize) and (sim >= 0.94):
            G.add_edge(i, j, weight=sim * 100)
            updated_large_clusters.add(i)
            updated_large_clusters.add(j)

    updated_lsize = {}
    updated_ssize = {}

    for x in lsize.keys():
        if x in updated_large_clusters:
            updated_lsize[x] = lsize[x]
        else:
            updated_ssize[x] = lsize[x]
    for x in ssize.keys():
        updated_ssize[x] = ssize[x] 
    return G, updated_lsize, updated_ssize

def louvain_clustering(G):
    partition = community_louvain.best_partition(G)
    return partition


def recluster(lsize, ssize, updated_lsize, updated_ssize, partition, pairwise_centroid_identity):
    print('lsize', lsize, 'ssize', ssize, 'updated_lsize', updated_lsize, 'updated_ssize', updated_ssize)
    p_dict = defaultdict(list)
    node_part = {}
    for node, part in partition.items():
        p_dict[part].append(node)
        node_part[node] = part

    if len(updated_lsize) == 0:
        for i, lid in enumerate(lsize.keys()):
            p_dict[i].append(lid)
            node_part[lid] = i
    else:
        lsize = updated_lsize
        ssize = updated_ssize

    s_sims = defaultdict(list)
    for i, isize, j, jsize, identity in pairwise_centroid_identity:
        if (i in lsize) and (j in ssize) and (identity >= 0.95):
            s_sims[j].append([i,identity])
        elif (i in ssize) and (j in lsize) and (identity >= 0.95):
            s_sims[i].append([j, identity])

        else:
            pass

    s_p = {}
    for s, sims in s_sims.items():
        identity_lst = [i[1] for i in sims]
        lids = [i[0] for i in sims]
        max_idex = identity_lst.index(max(identity_lst))
        max_lid = lids[max_idex]
        s_part = node_part[max_lid]
        s_p[s] = s_part
        p_dict[s_part].append(s)

    num_p = len(p_dict)

    for sid in ssize.keys():
        if sid not in s_p:
            p_dict[num_p+1].append(sid)
            num_p += 1

    return p_dict


def whether_recluster(lsize, pairwise_centroid_identity):
    min_identity = min([identity for i, isize, j, jsize, identity in pairwise_centroid_identity])
    if (min_identity <= 0.93) and (len(lsize) > 1) :
        flag = True
    else:
        flag = False
    return flag    

def main(nid, outfile):
    outf = open(outfile, 'w')
    infile = os.path.join('tmp', f"{nid}.inner.txt")
    pairwise_centroid_identity = load_inner_similarity_file(infile)
    all_size, lsize, ssize = size90(pairwise_centroid_identity)
    flag = whether_recluster(lsize, pairwise_centroid_identity)
    print(str(flag))
    if flag:
        G, updated_lsize, updated_ssize = build_graph(lsize, ssize, pairwise_centroid_identity)
        partition = louvain_clustering(G)
        p_dict = recluster(lsize, ssize, updated_lsize, updated_ssize, partition, pairwise_centroid_identity)
        print(p_dict)
        y = sum([len(t) for p, t in p_dict.items()])
        if y == len(all_size):
            for p, t in p_dict.items():
                outf.write(f"{nid}.{p}\t{len(t)}\t{','.join(t)}\n")
    else:
        ids = list(all_size.keys())
        outf.write(f"{nid}\t{len(ids)}\t{','.join(ids)}\n")
    

if __name__ == "__main__":
    #nid_mids = stat_merge_file(sys.argv[1])
    #c_size = load_centroid_size()
    #inner_similarity(nid_mids, c_size) 
    main(sys.argv[1], sys.argv[2])
