from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import alignment_v2 as alignment

def load_uniqhor():
    hor_nmer = {}
    hor_maxrepeat = {}
    with open("all.HiCAT.hor.summary.final.HOR_maxrepeat.xls", "r") as inf:
        for line in inf:
            if line.startswith('reorder_hor'):
                continue
            tokens = line.strip().split('\t')
            hor = tokens[0]
            maxrepeat = int(tokens[1])
            nmer = len(hor.split('_'))
            hor_nmer[hor] = nmer
            hor_maxrepeat[hor] = maxrepeat
    sorted_hor_nmer = dict(sorted(hor_nmer.items(), key=lambda item: item[1], reverse=True))
    return sorted_hor_nmer, hor_maxrepeat

def pairwise_overlap(ihor, jhor):
    imnset = set(ihor.split('_'))
    jmnset = set(jhor.split('_'))
    overlap = imnset.intersection(jmnset)
    return overlap

def pairwise_align(ihor, jhor):
    ihor_lst = ihor.split('_')
    jhor_lst = jhor.split('_')
    align1_1, align2_1, ed, identity = alignment.needle(ihor_lst, jhor_lst)
    align1_2, align2_2, ed2, identity2 = alignment.needle(jhor_lst, ihor_lst)
    if len(ihor_lst) > len(jhor_lst):
        rev_ihor_lst = ihor_lst[::-1]
        align1_3, align2_3, ed3, identity3 = alignment.needle(rev_ihor_lst, jhor_lst)
        align1_4, align2_4, ed4, identity4 = alignment.needle(jhor_lst, rev_ihor_lst)
    else:
        rev_jhor_lst = jhor_lst[::-1]
        align1_3, align2_3, ed3, identity3 = alignment.needle(rev_jhor_lst, ihor_lst)
        align1_4, align2_4, ed4, identity4 = alignment.needle(ihor_lst, rev_jhor_lst)
    max_match = max(identity, identity2, identity3, identity4)
    min_len = min(len(ihor_lst), len(jhor_lst))
    #print(ihor)
    #print(jhor)
    #print(max_match, min_len)
    #print('-----------------')
    return max_match / min_len


def same_layer_hors_group(hors, k=0):
    hor_group = {}
    if len(hors) == 1:
        hor_group[hors[0]] = k
        return hor_group
    else:
        for i in range(len(hors)):
            for j in range(i+1, len(hors)):
                ihor = hors[i]
                jhor = hors[j]
                overlap = pairwise_overlap(ihor, jhor)
                #print(i,j, len(overlap),overlap)
                if len(overlap) > 0:
                    match_ratio = pairwise_align(ihor, jhor)
                else:
                    match_ratio = 0
                #print(ihor,jhor, len(overlap),overlap, match_ratio )
                if match_ratio >= 0.5:
                    if ( ihor not in hor_group ) and (jhor not in hor_group):
                        hor_group[ihor] = k
                        hor_group[jhor] = k
                        k += 1
                    elif (ihor in hor_group) and (jhor not in hor_group):
                        hor_group[jhor] = hor_group[ihor]
                    elif (ihor in hor_group) and (jhor in hor_group):
                        if hor_group[ihor] != hor_group[jhor]:
                            min_idex = min(hor_group[ihor], hor_group[jhor])
                            max_idex = max(hor_group[ihor], hor_group[jhor])
                            for hor, group in hor_group.items():
                                if group == max_idex:
                                    hor_group[hor] = min_idex
                else:
                    if (ihor not in hor_group) and (jhor not in hor_group):
                        hor_group[ihor] = k 
                        k += 1
                        hor_group[jhor] = k
                        k += 1
                    elif (ihor in hor_group) and (jhor not in hor_group):
                        hor_group[jhor] = k
                        k += 1
        return hor_group
                

def layer_hors_group(hor_group, hors):
    initial_group_count = len(hor_group)
    print(initial_group_count)

    group_hors = {}
    for hor, group in hor_group.items():
        if group not in group_hors:
            group_hors[group] = []
        group_hors[group].append(hor)
    current_max_k = max(group_hors.keys())

    new_hor_flag = {}
    for ihor in hors:
        group_bh = {}
        for group, tophors in group_hors.items():
            tmp = []
            for thor in tophors:
                overlap = pairwise_overlap(ihor, thor)
                match_ratio = pairwise_align(ihor, thor) if len(overlap) > 0 else 0
                tmp.append(match_ratio)
            bh = max(tmp)
            group_bh[group] = bh
        new_hor_flag[ihor] = group_bh

    no_match_hors = []
    for ihor, group_dict in new_hor_flag.items():
        sorted_group_dict = dict(sorted(group_dict.items(), key=lambda item: item[1], reverse=True))
        print(ihor, sorted_group_dict)
        bh_group =list(sorted_group_dict.keys())[0]
        bh_flag = list(sorted_group_dict.values())[0]
        if bh_flag < 0.5:
            no_match_hors.append(ihor)
        else:
            hor_group[ihor] = bh_group
    no_match_hor_group = same_layer_hors_group(no_match_hors, current_max_k+1)
    updated_hor_group = hor_group | no_match_hor_group
    for ihor, group in updated_hor_group.items():
        print(ihor, group)
    assert initial_group_count + len(hors) == len(updated_hor_group), \
        f"Mismatch in expected and actual HOR group assignments: " \
        f"initial_group_count={initial_group_count}, hors={len(hors)}, updated={len(updated_hor_group)} " \
        f"exclude_hors: {[i for i in hors if i not in updated_hor_group]} "\
        f"no_match_hors: {no_match_hors}"
    return updated_hor_group
              

def layer_comparison(sorted_hor_nmer):
    nmer_hors = {}
    for hor, nmer in sorted_hor_nmer.items():
        if nmer not in nmer_hors:
            nmer_hors[nmer] = []
        nmer_hors[nmer].append(hor)

    for i, (nmer, hors) in enumerate(nmer_hors.items()):
        #print(nmer, len(hors))
        if i == 0:
            hor_group = same_layer_hors_group(hors, 0)
            print(hor_group)
            print('--------------intial done-------------------')
            #break
        else:
            print('--------------',i, nmer, '------------------')
            new_hor_group = layer_hors_group(hor_group, hors)
            hor_group = new_hor_group
            #break
    outfile = "all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.xls"
    #outfile = "test.xls"
    outf = open(outfile, 'w') 
    sorted_hor_group = dict(sorted(hor_group.items(), key=lambda item: item[1], reverse=False))
    for ihor, group in sorted_hor_group.items():
        outf.write(f"{ihor}\t{group}\n")
    outf.close()
    return sorted_hor_group

def load_v0_class():
    group_hors = {}
    with open("all.HiCAT.hors.final.updated.uniqhor.v0.class.xls","r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            ihor= tokens[0]
            group= tokens[1]
            if group not in group_hors:
                group_hors[group] = []
            group_hors[group].append(ihor)
    return group_hors


def split_class(targets, group_hors):
    for target in targets:
        #target = targets[0]
        hors = group_hors[target]
        total_mns = list(set([mn for ihor in hors for mn in ihor.split('_')]))
        data = pd.DataFrame(index=hors, columns=total_mns, dtype=int).fillna(0)
        for ihor in hors:
            for mn in ihor.split('_'):
                data.loc[ihor, mn] += 1

        print(data.iloc[:5,:5])

        import seaborn as sns
        import matplotlib.pyplot as plt

        g = sns.clustermap(data, metric="euclidean", method="ward", cmap="viridis")
        reordered_data = g.data2d
        reordered_data.to_csv(f"{target}.matrix.xls", sep="\t", index=True, header=True)

        g.savefig(f"{target}.split.png", dpi=300)
        #plt.show()
    
if __name__ == "__main__":
    sorted_hor_nmer, hor_maxrepeat = load_uniqhor()
    layer_comparison(sorted_hor_nmer)  
    #group_hors = load_v0_class()
    #split_class(['113'], group_hors)
