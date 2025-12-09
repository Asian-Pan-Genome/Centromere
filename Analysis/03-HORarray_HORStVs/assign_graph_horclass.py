import pandas as pd
import alignment_v2 as alignment


def graph_horclass():
    df = pd.read_csv("all.graph.hordecomposition.summary.final.pan_hor_dr_type.xls", sep="\t", header=0)
    df['nmer'] = df['reorder_hor'].apply(lambda x: len(x.split('_')))
    graph_hors_nmer = dict(zip(df['reorder_hor'], df['nmer']))
    sorted_graph_hors_nmer = dict(sorted(graph_hors_nmer.items(), key=lambda item: item[1], reverse=True))

    print("graph_uniq_hors", len(sorted_graph_hors_nmer))

    nmer_hors = {}
    for hor, nmer in sorted_graph_hors_nmer.items():
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

    outfile = "all.graph.hordecomposition.summary.final.horclass.05.xls"
    #outfile = "test.xls"
    outf = open(outfile, 'w') 
    sorted_hor_group = dict(sorted(hor_group.items(), key=lambda item: item[1], reverse=False))
    for ihor, group in sorted_hor_group.items():
        outf.write(f"{ihor}\t{group}\n")
    outf.close()
    return sorted_hor_group


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


def unified_hicat_graph_horclass():
    graph_horclass_df = pd.read_csv("all.graph.hordecomposition.summary.final.horclass.05.xls", sep="\t", names=['hor', 'group'])

    hicat_df = pd.read_csv("all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.horindex.update.xls", sep="\t",header=0)
    hicat_hor_horclass = dict(zip(hicat_df['reorder_hor'], hicat_df['horclass']))

    graph_horclass_df['hicat_horclass'] = graph_horclass_df['hor'].map(hicat_hor_horclass)

    def process_group(group):
        # Exclude NaN values
        non_nan_classes = group['hicat_horclass'].dropna().unique()
        
        if len(non_nan_classes) == 1:  # If all non-NaN values are the same
            group['hicat_horclass'] = group['hicat_horclass'].fillna(non_nan_classes[0])
        # If there are >1 unique non-NaN values, leave NaN unchanged
        return group

    graph_horclass_df = graph_horclass_df.groupby('group').apply(process_group)
    graph_horclass_df.to_csv("all.graph.hordecomposition.summary.final.horclass.05.hicat_horclass.xls", sep="\t", index=False, header=True)


def vs_CHM13():
    clusterdf = pd.read_csv("../../cluster.info.xls", sep="\t", header=0)
    clusterdf['Nid'] = clusterdf['Nid'].astype(str)
    nid_anno = dict(zip(clusterdf['Nid'], clusterdf['VS_CHM13_anno_bh']))
    
    graph_horclass_df = pd.read_csv("all.graph.hordecomposition.summary.final.horclass.05.hicat_horclass_manual_check.xls", sep="\t", header=0)
    def assign_anno(row):
        hor = row['hor'].split('_')
        anno = '_'.join([nid_anno[i] for i in hor])
        return anno

    graph_horclass_df['anno'] = graph_horclass_df.apply(assign_anno, axis=1)
    graph_horclass_df.to_csv("all.graph.hordecomposition.summary.final.horclass.05.hicat_horclass_manual_check.anno.xls", sep="\t", index=False, header=True)
    
if __name__ == "__main__":
    # graph_horclass() 
    # unified_hicat_graph_horclass()
    vs_CHM13()
