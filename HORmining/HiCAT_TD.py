from collections import Counter
from get_continuous_block import regex_mn
from HiCAT_HOR import miningMonomerTDPattern, buildingHor
import os

def SingleBlockClusterStat(iblock, seqid_clustid):
    mns = iblock[5].split(',')
    cids = []
    for mn in mns:
        mn = mn.strip()
        cid = seqid_clustid.get(mn,'-')
        cids.append(cid)
    count = Counter(cids)
    return cids, count

def detect_HOR_TD(cids, mns, max_hor_len):
    new_monomer_sequence, top_layer, all_layer = miningMonomerTDPattern(cids, max_hor_len)
    mod_mns = []
    mns = mns.split(',')
    for mnid in mns:
        ichr, start, end, strand = regex_mn(mnid)
        mod_mns.append(ichr + '_' + start + '_' + end + '_' + strand)
    HORs = buildingHor(mod_mns, all_layer)
    return top_layer, all_layer, HORs


def filter_HOR_TD(HORs, cids):
    filter_HORs = {}
    for i in HORs.keys():
        pattern = i
        database = HORs[i][1]
        filter_HOR = []
        init_sequences = [0] * len(cids)
        sort_database = sorted(database, key=lambda x: x[1] - x[0])
        for j in sort_database:
            start = j[0]
            end = j[1]
            jump = 0
            for k in range(start, end + 1):
                if init_sequences[k] == 1:
                    jump = 1
                    break
            if jump == 1:
                continue
            else:
                for k in range(start, end + 1):
                    init_sequences[k] = 1
                filter_HOR.append(j)
        if len(filter_HOR) > 0:
            filter_HORs[pattern] = sorted(filter_HOR, key=lambda x: x[0])

    filter_HORs_list = {}
    for i in filter_HORs.keys():
        if '-' in i:
            continue
        pattern = i.split('_')
        pattern_database = filter_HORs[i]
        in_flag = 0
        for j in range(len(pattern)):
            prefix_pattern = pattern[j:]
            suffix_pattern = pattern[:j]
            loop_pattern = prefix_pattern + suffix_pattern
            s_loop_pattern = '_'.join(loop_pattern)
            if s_loop_pattern in filter_HORs_list.keys():
                in_flag = 1
                for k in pattern_database:
                    filter_HORs_list[s_loop_pattern].append([k[0], k[1], k[2], k[3], k[4]])
                break
        if in_flag == 0:
            s_pattern = '_'.join(pattern)
            filter_HORs_list[s_pattern] = []
            for j in pattern_database:
                filter_HORs_list[s_pattern].append([j[0], j[1], j[2], j[3], j[4]])
    return filter_HORs_list


def score_HOR_TD(filter_HORs_list, cids):
    pattern_score = []
    for i in filter_HORs_list.keys():
        cov = 0
        database = filter_HORs_list[i]
        for j in database:
            cov += (j[1] + 1 - j[0])
        cov_rate = cov / len(cids)
        pattern = i.split('_')
        repeat_number = 0
        for j in database:
            repeat_number += j[4]
        pattern_rate = repeat_number / (cov / len(pattern))

        head_tail_num = 0
        if len(database) > 1:
            for k in range(0, len(database) - 1):
                pre_end = database[k][1]
                sub_start = database[k + 1][0]
                if sub_start - pre_end == 1:
                    head_tail_num += 1
            head_tail_ratio = head_tail_num / (len(database) - 1)
        else:
            head_tail_ratio = 0.1
        score0 = cov_rate * pattern_rate
        score1 = cov_rate * pattern_rate * head_tail_ratio
        pattern_score.append([score0, score1, cov_rate, pattern_rate, head_tail_ratio, i])
    pattern_score = sorted(pattern_score, key=lambda x: x[0], reverse=True)

    marker_sequence = [0] * len(cids)
    filter_final_HORs = {}
    for i in pattern_score:
        database = filter_HORs_list[i[5]]
        for j in database:
            start = j[0]
            end = j[1]
            for k in range(start, end + 1):
                marker_sequence[k] = 1
        filter_final_HORs[i[5]] = database

    cov = sum(marker_sequence)
    cov_rate = cov / len(cids)
    max_cov_pattern = pattern_score[0][5].split('_')
    max_cov_pattern_monomer = set(max_cov_pattern)
    max_cov_rate = pattern_score[0][2]
    max_pattern_score = pattern_score[0][0]
    max_score1 = pattern_score[0][1]
    max_head_tail_ratio = pattern_score[0][4]
    hornum = len(filter_final_HORs)

    if len(filter_final_HORs.keys()) == 0:
        CE_rate = -1
    else:
        CE_rate = len(max_cov_pattern_monomer) / len(max_cov_pattern)

    statistics = [cov_rate, len(max_cov_pattern), max_cov_rate, CE_rate, max_pattern_score, max_score1,
                  max_head_tail_ratio, hornum]

    return pattern_score, filter_final_HORs, statistics

def output_block(outdir, iblock, filter_final_HORs, blockindex, statistics, top_layer, all_layer):

    with open(os.path.join(outdir, str(blockindex)+".hor.xls"), "w") as out_hor:
        out_hor.write("#" + ",".join(map(str, iblock[:5])) + "," + str(iblock[-1]) + "\n")
        for ihor, infos in filter_final_HORs.items():
            for info in infos:
                s = info[0]
                e = info[1]
                strand = info[2]
                true_pattern = "_".join(info[3])
                num = info[4]
                out_hor.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (ihor, s, e, strand, true_pattern, num))

    with open(os.path.join(outdir, str(blockindex)+".hor.stat.xls"), 'w') as out_stat:
        out_stat.write("#" + ",".join(map(str, iblock[:5])) + "," + str(iblock[-1]) + "\n")
        out_stat.write("total_coverage" + "\t" + "max_hor_length" + "\t" + "max_hor_coverage" + "\t" + "CE_rate" + "\t" + "max_hor_pattern_score" + "\t" + "max_score1" + "\t" + "max_hor_htratio" + "\t" + "hor_count" + "\n")
        out_stat.write("\t".join([str(i) for i in statistics]) + "\n" )

    top_layer_set = set()
    sorted_top_layer = sorted(top_layer,key=lambda x:x[0])
    layers = {}
    layers_info = {}

    with open(os.path.join(outdir, str(blockindex)+".top_layer.xls"), 'w') as out_top:
        for j in sorted_top_layer:
            start = j[0]
            end = j[1]
            pattern = '_'.join(j[2])
            repeat_number = j[3]
            top_layer_set.add(pattern+'_'+str(start)+'_'+str(end))
            out_top.write(str(start)+'\t'+str(end)+'\t'+str(repeat_number)+'\t'+str(pattern)+'\n')
            layers[pattern +'_'+str(start)+'_'+str(end)] = []
            layers_info[pattern +'_'+str(start)+'_'+str(end)] = [start,end,repeat_number,pattern]
    # fix same start region of top and cov
    sorted_all_layer = sorted(all_layer, key=lambda x: x[0])
    with open(os.path.join(outdir, str(blockindex)+".all_layer.xls"), 'w') as out_all:
        for j in sorted_all_layer:
            start = j[0]
            end = j[1]
            pattern = '_'.join(j[2])
            repeat_number = j[3]
            if pattern + '_' + str(start) + '_' + str(end) in layers.keys():
                pass
            else:
                for l in layers.keys():
                    top_start = int(l.split('_')[-2])
                    top_end = int(l.split('_')[-1])
                    if start >= top_start and end <= top_end:
                        layers[l].append([start,end,repeat_number,pattern])

        for j in layers_info.keys():
            start = layers_info[j][0]
            end = layers_info[j][1]
            repeat_number = layers_info[j][2]
            pattern = layers_info[j][3]
            out_all.write(str(start) + '\t' + str(end) + '\t' + str(repeat_number) + '\t' + str(pattern) + '\t' + 'top' + '\n')
            for k in layers[j]:
                sub_start = k[0]
                sub_end = k[1]
                sub_repeat_number = k[2]
                sub_pattern = k[3]
                out_all.write(str(sub_start) + '\t' + str(sub_end) + '\t' + str(sub_repeat_number) + '\t' + str(sub_pattern) + '\t' + 'cover' + '\n')

def main_HiCAT(cblocks, seqid_clustid, outdir):
    hicat_outdir = os.path.join(outdir, "HiCAT")
    os.makedirs(hicat_outdir, exist_ok=True)

    for i, iblock in enumerate(cblocks):
        cids, count = SingleBlockClusterStat(iblock, seqid_clustid)

        with open(os.path.join(hicat_outdir, str(i) + ".chain.xls"), 'w') as outf1:
            outf1.write("#" + ",".join(map(str, iblock[:5])) + "," + str(iblock[-1]) + "\n")
            outf1.write(",".join(cids) + "\n")

        with open(os.path.join(hicat_outdir, str(i) + ".chain.stat.xls"), 'w') as outf2:
            outf2.write("#" + ",".join(map(str, iblock[:5])) + "," + str(iblock[-1]) + "\n")
            sort_count = dict(sorted(count.items(), key=lambda item: item[1], reverse=True))
            for mnid, c in sort_count.items():
                outf2.write(mnid + "\t" + str(c) + "\n")

        mns = iblock[5]
        if iblock[1] != "chrY":
            top_layer, all_layer, HORs = detect_HOR_TD(cids, mns, 40)
        else:
            top_layer, all_layer, HORs = detect_HOR_TD(cids, mns, 50)
        filter_HORs_list = filter_HOR_TD(HORs, cids)
        if len(filter_HORs_list) > 0:
            pattern_score, filter_final_HORs, statistics =  score_HOR_TD(filter_HORs_list, cids)
            output_block(hicat_outdir, iblock, filter_final_HORs, i, statistics, top_layer, all_layer)
            print(i, "th block hor identification done")
        else:
            print(i, "th block has no hor")
        
