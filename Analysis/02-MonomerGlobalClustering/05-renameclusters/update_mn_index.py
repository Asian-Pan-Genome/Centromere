

def load_centroid_size():
    c_size = {}
    with open("centroid.size", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            c_size[tokens[0]] = int(tokens[1])
    #print(list(c_size.keys())[:10])
    return c_size


def load_final_linked(c_size):
    oid_nid = {}
    nid_oids = {}
    outf = open("final.linked.updated.v2.txt", "w") 

    with open("final.linked.ordered.v2.txt", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            nid = tokens[0]
            mark = tokens[1]
            clusternum = tokens[2]
            oids = tokens[3].split(',')
            nid_oids[nid] = oids
            clustersize = sum([c_size[i] for i in oids])
            for i in oids:
                oid_nid[i] = nid
            outf.write(f"{nid}\t{mark}\t{clusternum}\t{clustersize}\t{tokens[3]}\n")
    outf.close()
    return oid_nid, nid_oids


def update_mn_index(nid_oids):
    outf = open("all.mn.updated.v2.index", "w")

    for nid, oids in nid_oids.items():
        for i in oids:
            idfile = f"../step-0/original_seq/{i}.fai"
            with open(idfile, 'r') as inf:
                for line in inf:
                    tokens = line.strip().split('\t')
                    outf.write(f"{tokens[0]}\t{nid}\n")
    outf.close()

def main():
    c_size = load_centroid_size()
    oid_nid, nid_oids = load_final_linked(c_size)
    update_mn_index(nid_oids)


if __name__ == "__main__":
    main()

                    
