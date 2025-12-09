import sys
from get_continuous_block import *
import os

def search_start(t):
    min_value = min(t)
    min_indices = [i for i, x in enumerate(t) if x == min_value]
    min_index = min_indices[0]
    for i in min_indices:
        for j in range(1, len(t) - i):
            if i + j >= len(t) or min_index + j >= len(t):
                break
            if t[i + j] < t[min_index + j]:
                min_index = i
                break
            elif t[i + j] > t[min_index + j]:
                break
    return min_index

def unified_order(ihor):
    ihor_lst = ihor.split(',')
    t = []
    for k in ihor_lst:
        if k != '-':
            t.append(int(k))
        else:
            t.append(float('inf'))
    min_index = search_start(t)
    update_ihor_lst = ihor_lst[min_index:] + ihor_lst[:min_index]
    update_ihor = ",".join(list(map(str, update_ihor_lst)))         
    return update_ihor
    

def graph_hordecomposition_stat(cblocks, seqid_clustid,  outdir):
    #cblocks = CountinuousBlock(asatbed, distance)
    #b_hors = {}
    statfile = os.path.join(outdir, "graph", "graph.HORdecomposition.stat.xls")
    outf = open(statfile, 'w')

    for i, iblock in enumerate(cblocks):
        chrom = iblock[1]
        start = iblock[2]
        end = iblock[3]
        blen = iblock[4]
        num = iblock[6]
        
        horfile = os.path.join(outdir, "graph",  str(i) +  ".HORdecomposition.tsv")
        if os.path.exists(horfile):
            with open(horfile, 'r') as inf:
                for line in inf:
                    line = line.strip()
                    tokens = line.split('\t')
                    if ',' in tokens[6]:
                        update_hor = unified_order(tokens[6])
                        n_mer = len(update_hor.split(','))
                    else:
                        update_hor = tokens[6]
                        n_mer = 1
                    outf.write(f"{line}\t{update_hor}\t{n_mer}\n")
        else:
            mns = iblock[5].split(',')
            index=0
            for k, mn in enumerate(mns):
                ichr, mn_start, mn_end, strand = regex_mn(mn)
                cid = seqid_clustid.get(mn,'-')
                outf.write(f"{i}\t{index+k}\t{index+k+1}\t{chrom}\t{mn_start}\t{mn_end}\t{cid}\t{cid}\t1\n")
#if __name__ == "__main__":
#    main(sys.argv[1], sys.argv[2], sys.argv[3]) 
