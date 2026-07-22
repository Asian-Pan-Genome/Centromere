import pysam

def stat_linked(infile):
    i = 0
    mid_oids = {}
    with open(infile, 'r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            mid_oids[i] = tokens[2].split(',')
            i += 1
    return mid_oids


def load_centroid_size():
    c_size = {}
    with open("../centroid.size", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            c_size[tokens[0]] = int(tokens[1])
    print(list(c_size.keys())[:10])
    return c_size


def load_centroid_seqname():
    oid_seqname={}
    with open("../all.centroid.txt", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            oid_seqname[tokens[1]] = tokens[8]
    return oid_seqname


def select_largest_original_cluster(mid_oids, c_size, oid_seqname):

    outf = open("marker_seq.fa", "w")

    for mid, oids in mid_oids.items():
        oid_osize = { i: c_size[i] for i in oids}
        sorted_oid_osize = dict(sorted(oid_osize.items(), key=lambda item: item[1], reverse=True))

        max_oid = list(sorted_oid_osize.keys())[0]
        max_centroid_seqname = oid_seqname[max_oid]
        fa = pysam.FastaFile("../../step-0/original_seq/" + max_oid)
        max_centroid_seq = fa.fetch(max_centroid_seqname)
    
        outf.write(f">{mid}@{max_centroid_seqname}\n{max_centroid_seq}\n")

    outf.close()


def main():
    mid_oids = stat_linked("final.linked.txt")
    c_size = load_centroid_size()
    oid_seqname = load_centroid_seqname()
    select_largest_original_cluster(mid_oids, c_size, oid_seqname)

if __name__ == "__main__":
    main()
            
        
        
        
        
    
