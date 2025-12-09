import pybedtools
import re

def CountinuousBlock(hmmer_file, distance):
    cblocks = []
    bed = pybedtools.BedTool(hmmer_file)
    merged = bed.merge(d=distance, c=[4], o='collapse')
    for index, interval in enumerate(merged):
        filter_mns = FilterMn(interval[3])
        print(index, len(interval[3].split(',')), len(filter_mns.split(',')))
        cblocks.append([index, interval[0], interval[1], interval[2], int(interval[2])-int(interval[1]), filter_mns, len(filter_mns.split(','))])
    #print(cblocks[0][0:5], cblocks[0][-1])
    #print(cblocks[1][0:5], cblocks[1][-1])
    return cblocks

def OutFilterMn(cblocks, seqid_clustid, outfile):
    with open(outfile, 'w') as outf:
        for iblock in cblocks:
            mns = iblock[5].split(',')
            for mn in mns:
                mn = mn.strip()
                ichr, start, end, strand = regex_mn(mn)
                cid = seqid_clustid.get(mn,'-')
                outf.write(f"{ichr}\t{start}\t{end}\t{cid}\t0\t{strand}\n")
                
def FilterMn(mnchain):
    ###filter monomers overlaped with previous one larger than 80%###
    chain_lst = mnchain.split(',')
    pos = {}
    for index, mn in enumerate(chain_lst):
        ichr, start, end, strand = regex_mn(mn)
        start = int(start)
        end = int(end)
        pos[index] = [min([start, end]), max([start, end])]

    out = [chain_lst[0]]
    for i in range(1,len(pos)):
        pre_s, pre_e = pos[i-1][0], pos[i-1][1]
        sub_s, sub_e = pos[i][0], pos[i][1]
        overlap_s = max(pre_s, sub_s)
        overlap_e = min(pre_e, sub_e)
        if overlap_s <= overlap_e:
            overlap_len = overlap_e - overlap_s
            ratio = overlap_len / (pre_e - pre_s)
            if ratio < 0.8:
                out.append(chain_lst[i])
            else:
                print("mn0: ", chain_lst[i-1], "mn1: ", chain_lst[i], pre_s, pre_e, sub_s, sub_e, overlap_s, overlap_e, "overlap ratio: ", ratio)
        else:
            out.append(chain_lst[i])
    return ",".join(out)

def regex_mn(mnid):
    if "::" in mnid:
        regex = re.search(r"(.+)::(.+):(.+)-(.+)\((.+)\)$", mnid)
        ichr = regex.group(2)
        start = regex.group(3)
        end = regex.group(4)
        strand = regex.group(5)
    else:
        regex = re.search(r"(.+):(.+)-(.+)\((.+)\)$", mnid)
        ichr = regex.group(1)
        start = regex.group(2)
        end = regex.group(3)
        strand = regex.group(4)
    return ichr, start, end, strand


def MonomerCluster(clustid_file):
    seqid_clustid = {}
    with open(clustid_file, 'r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            seqid_clustid[tokens[0]] = tokens[1]
    return seqid_clustid
