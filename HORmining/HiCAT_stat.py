from collections import OrderedDict, defaultdict
import sys
from get_continuous_block import *
import os

def is_subsequence(str1, str2):
    '''check whether str2 is a hor-stv in str1'''
    str1_parts = str1.split('_')
    str2_parts = str2.split('_')

    extended_str1 = str1_parts + str1_parts
    
    it = iter(extended_str1)
    return all(x in it for x in str2_parts)

def is_combination(str1, str2):
    '''check whether str2 is a combination of str1 hor-stvs'''
    str1_parts = str1.split('_')
    str2_parts = str2.split('_')
    
    for i in range(1, len(str2_parts)):
        first_part = str2_parts[:i]
        second_part = str2_parts[i:]
        
        if is_subsequence(str1, '_'.join(first_part)) and is_subsequence(str1, '_'.join(second_part)):
            return True
        
    return False

def belong_to_same(str1, str2):
    if is_subsequence(str1, str2):
        return True
    elif is_combination(str1, str2):
        return True
    else:
        return False

def simplified_hor(horfile):
    hors_dict = OrderedDict()
    with open(horfile, 'r') as inf:
        for line in inf:
            line = line.strip()
            if line.startswith('#'):
                continue
            else:
                tokens = line.split('\t')
                ihor = tokens[0]
                if ihor not in hors_dict:
                    hors_dict[ihor] = 1
                else:
                    hors_dict[ihor] += 1
    hors = list(hors_dict.keys())
    out = [hors[0]]
    print("out:", out)
    for i in range(len(hors)):
        for j in range(len(hors)):
            if i < j:
                flag = belong_to_same(hors[i], hors[j])
                #print(hors[i], hors[j], flag)
                if not flag:
                    out.append(hors[j])
       
    return out

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


def stat_alllayer(chrom, mns, f_all):

    starts = []
    ends = []

    mns = mns.split(',')
    for mn in mns:
        ichr, start, end, strand = regex_mn(mn)
        starts.append(start)
        ends.append(end)

    outf = open(f_all+".reorder.xls", "w") 
    hors = defaultdict(int)
    with open(f_all, 'r') as inf:
        for line in inf:
            line = line.strip()
            tokens = line.split('\t')

            start_index = int(tokens[0])
            end_index = int(tokens[1])

            start_pos = starts[start_index]
            try:
                end_pos = ends[end_index]
            except:
                print(f_all, chrom, len(ends), end_index)

            hor = tokens[3]
            count = int(tokens[2])
            if hor != "-":
                ihor_lst = hor.split('_')
                t = []
                for k in ihor_lst:
                    if k != '-':
                        t.append(int(k))
                    else:
                        t.append(float('inf'))
                
                #initial order
                initial_min_index = t.index(min(t))
                initial_ihor_lst = ihor_lst[initial_min_index:] + ihor_lst[:initial_min_index]
                initial_t = t[initial_min_index:] + t[:initial_min_index]
                
                min_index = search_start(initial_t)
                update_ihor_lst = initial_ihor_lst[min_index:] + initial_ihor_lst[:min_index]
                update_ihor = "_".join(list(map(str, update_ihor_lst)))
                outf.write(chrom + "\t" + start_pos + "\t" + end_pos + "\t" + line + "\t" + update_ihor + "\n")
                hors[update_ihor] += count
    sorted_hors = dict(sorted(hors.items(), key=lambda item: item[1], reverse=True))
    outf.close()
    return sorted_hors
    
def hicat_hor_stat(outdir, asatbed,  distance=10000):
    cblocks = CountinuousBlock(asatbed, distance)
    #b_hors = {}
    statfile = os.path.join(outdir, "HiCAT", "HiCAT.hor.stat.xls")
    summaryfile = os.path.join(outdir, "HiCAT", "HiCAT.hor.summary.xls")
    outf = open(statfile, 'w')
    outs = open(summaryfile, 'w')

    for i, iblock in enumerate(cblocks):
        chrom = iblock[1]
        start = iblock[2]
        end = iblock[3]
        blen = iblock[4]
        mns = iblock[5]
        num = iblock[6]

        f_all = os.path.join(outdir, "HiCAT", str(i) +  ".all_layer.xls")
        #horfile = os.path.join(outdir, str(i) +  ".hor.xls")
        if os.path.exists(f_all):
            hors = stat_alllayer(chrom, mns, f_all)
            
            for hi, ihor in enumerate(list(hors.keys())):
                horlen = len(ihor.split('_'))
                icount = hors[ihor]
                outf.write(f"{i}\t{chrom}\t{start}\t{end}\t{blen}\t{num}\t{hi}\t{ihor}\t{horlen}\t{icount}\n")
        else:
            outf.write(f"{i}\t{chrom}\t{start}\t{end}\t{blen}\t{num}\tNa\tNa\tNa\tNa\n")

        f_all_reorder = os.path.join(outdir, "HiCAT", str(i) +  ".all_layer.xls.reorder.xls")
        if os.path.exists(f_all_reorder) and os.path.getsize(f_all_reorder) > 0:
            with open(f_all_reorder, "r") as inf:
                for line in inf:
                    line = line.strip()
                    outs.write(f"{i}\t{start}\t{end}\t{blen}\t{num}\t{line}\n")
        else:
            tmp = "\t".join(['Na']*9)
            outs.write(f"{i}\t{start}\t{end}\t{blen}\t{num}\t{tmp}\n")

    outf.close()
    outs.close()


if __name__ == "__main__":
    hicat_hor_stat(sys.argv[1], sys.argv[2])

