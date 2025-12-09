import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def stat_linked(infile):
    i = 0
    oid_mid = {}
    with open(infile, 'r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            print(i,tokens[0])
            for j in tokens[2].split(','):
                oid_mid[j] = i
            i += 1
    return oid_mid


def load_inner_similarity_file(infile):
    pairwise_centroid_identity = []
    with open(infile, 'r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            pairwise_centroid_identity.append([tokens[0], int(tokens[1]), tokens[2], int(tokens[3]), float(tokens[4])])
    return pairwise_centroid_identity


def load_compares(targets):
    compares = [int(i) for i in targets.split(',')]
    pairwise = []
    for k in range(len(compares)):
        for t in range(k+1, len(compares)):
            pairwise.append((compares[k], compares[t]))
    return compares, pairwise


def plot_boxplot(compares, pairwise, oid_mid, pairwise_centroid_identity):
    out = defaultdict(list)
    for i, isize, j, jsize, identity in pairwise_centroid_identity:
        imid = oid_mid[i]
        jmid = oid_mid[j]
        if (imid in compares) and (jmid in compares):
            if imid == jmid:
                out[imid].extend([identity])
            else:
                if (imid, jmid) in pairwise:
                    out[(imid, jmid)].append(identity)
                elif (jmid, imid) in pairwise:
                     out[(jmid, imid)].append(identity)

                else:
                    print(jmid, imid, "please check")
    data = [{"key": str(k), "value": v} for k, values in out.items() for v in values]
    df = pd.DataFrame(data)

    plt.figure(figsize=(10, 6))
    sns.boxplot(x="key", y="value", data=df)
    plt.xlabel("Key", fontsize=12)
    plt.ylabel("Value Distribution", fontsize=12)
    plt.title("Boxplot of Values Grouped by Key", fontsize=14)
    plt.xticks(rotation=45) 
    plt.show()

def main():
    linkedfile=sys.argv[1]
    innerfile=sys.argv[2]
    targets=sys.argv[3]
    oid_mid = stat_linked(linkedfile)
    pairwise_centroid_identity = load_inner_similarity_file(innerfile)
    compares, pairwise = load_compares(targets)
    print(compares, pairwise)
    plot_boxplot(compares, pairwise, oid_mid, pairwise_centroid_identity)
if __name__ == "__main__":
    main() 
        

