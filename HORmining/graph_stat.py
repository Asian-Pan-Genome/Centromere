import sys
#from get_continuous_block import *
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

def unified_order(horfile):
    hors = {}
    with open(horfile, 'r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            hors[tokens[0]] = tokens[1]
    update_hors = {}
    for m, ihor in hors.items():
        ihor_lst = ihor.split(',')[:-1]
        t = []
        for k in ihor_lst:
            if k != '-':
                t.append(int(k))
            else:
                t.append(float('inf'))
        min_index = search_start(t)
        update_ihor_lst = ihor_lst[min_index:] + ihor_lst[:min_index]
        update_ihor = ",".join(list(map(str, update_ihor_lst)))
        update_hors[m] = update_ihor
            
    return update_hors
    

def graph_hor_stat(cblocks, outdir):
    #cblocks = CountinuousBlock(asatbed, distance)
    #b_hors = {}
    statfile = os.path.join(outdir, "graph", "graph.hor.stat.xls")
    outf = open(statfile, 'w')

    for i, iblock in enumerate(cblocks):
        chrom = iblock[1]
        start = iblock[2]
        end = iblock[3]
        blen = iblock[4]
        num = iblock[6]
        
        horfile = os.path.join(outdir, "graph",  str(i) +  ".HORs.tsv")
        if os.path.exists(horfile):
            update_hors = unified_order(horfile)
            #b_hors[i] = update_hors
            for m, update_ihor in update_hors.items():
                horlen = len(update_ihor.split(','))
                outf.write(f"{i}\t{chrom}\t{start}\t{end}\t{blen}\t{num}\t{m}\t{update_ihor}\t{horlen}\n")
        else:
            #b_hors[i] = {}
            outf.write(f"{i}\t{chrom}\t{start}\t{end}\t{blen}\t{num}\tNa\tNa\tNa\n")

#if __name__ == "__main__":
#    main(sys.argv[1], sys.argv[2], sys.argv[3]) 
