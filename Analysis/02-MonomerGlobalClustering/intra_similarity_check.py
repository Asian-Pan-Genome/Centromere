import sys
from collections import defaultdict
def stat_merge_file(infile):
    nid_mids = {}
    with open(infile, "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            mids = tokens[2].split(',')
            nid_mids[tokens[0]] = mids
    return nid_mids

def load_centroid_size():
    c_size = {}
    with open("centroid.size", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            c_size[tokens[0]] = int(tokens[1])
    #print(list(c_size.keys())[:10])
    return c_size

def intra_similarity(target, nid_mids, c_size):
    targets = target.split('_')
    targets_mids = []

    for t in targets:
        targets_mids.extend(nid_mids[t])
    tmp = sorted(list(map(int, targets_mids)))
    sorted_mids = list(map(str,tmp))
    print(len(sorted_mids))

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

    with open(f"tmp/{target}.intra.txt", 'w') as outf:
        for info in pairwise_centroid_identity:
            outf.write(f"{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}\t{info[4]}\n")


if __name__ == "__main__":
    linkedfile=sys.argv[1]
    target = sys.argv[2]
    nid_mids = stat_merge_file(linkedfile)
    c_size = load_centroid_size()
    intra_similarity(target, nid_mids, c_size)

