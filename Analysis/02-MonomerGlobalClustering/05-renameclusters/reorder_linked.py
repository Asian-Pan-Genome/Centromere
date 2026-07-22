from collections import defaultdict
nid_line = {}

i=0
with open("final.linked.txt", "r") as inf:
    for line in inf:
        line = line.strip()
        nid_line[i] = line
        i += 1


order_nid = defaultdict(list)
with open("fasttree.order", "r") as inf:
    for line in inf:
        tokens = line.strip().split('\t')
        nid = int(tokens[0])
        order = int(tokens[1])
        order_nid[order].append(nid)


with open("final.linked.ordered.v2.txt", "w") as outf:
    for order, nid in order_nid.items():
        if len(nid) == 1:
            line = nid_line[nid[0]]
            outf.write(f"{order}\t{line}\n")
        else:
            tmp = []
            num = []
            oids = []
            for i in nid:
                line = nid_line[i]
                tokens = line.split('\t')
                num.append(int(tokens[1]))
                tmp.append(tokens[2])
                oids.append(tokens[0])
            print(tmp)
            print(num)
            print(oids)
            outf.write(f"{order}\t{'#'.join(oids)}\t{sum(num)}\t{','.join(tmp)}\n")



