import pandas as pd

def load_order():
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-11/cumsum/cumsum.order.xls", sep="\t", header=0)

    ##sorted##
    df_sorted = df.sort_values(by='cumsum_Order')
    df_apg_first = df.sort_values(by='ref_apg_hprc_hgsvc')
    df_hprc_first = df.sort_values(by='ref_hprc_hgsvc_apg')

    order_sorted = list(df_sorted['haploid'])
    order_apg_first = list(df_apg_first['haploid'])
    order_hprc_first = list(df_hprc_first['haploid'])

    out = [order_sorted, order_apg_first, order_hprc_first]
    return out


def get_hicat_cumsum_data(order, outfile):
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-11/cumsum/all.HiCAT.hor.summary.final.xls", sep="\t", header=0)
    stvdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/cumsum/20250519_hicat_hor_final_summary.xls", sep="\t", header=0)
    filtereddf = stvdf[~stvdf['HORarray_reported_by_CHM13'].isin(['Novel(single_mn)', 'Novel(dimeric_mns)'])].copy()
    print(filtereddf.shape)
    hors = list(filtereddf['reorder_hor'])

    finaldf = df[df['reorder_hor'].isin(hors)].copy()

    sample_hap_horset = finaldf.groupby('sample_hap')['reorder_hor'].apply(set).to_dict()

    outf = open(outfile, "w")

    for i, sh in enumerate(order):
        ihorset = sample_hap_horset[sh]
        if i == 0:
            union_out = ihorset
            interset_out = ihorset
        else:
            union_out = union_out | ihorset
            interset_out = interset_out & ihorset
        outf.write(f"{i}\t{sh}\t{len(union_out)}\t{len(interset_out)}\n")

    outf.close()
    return hors
    
def get_hicat_cumsum_pantype(order, outfile):
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-11/cumsum/all.HiCAT.hor.summary.final.xls", sep="\t", header=0)
    stvdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/cumsum/20250519_hicat_hor_final_summary.xls", sep="\t", header=0)
    filtereddf = stvdf[~stvdf['HORarray_reported_by_CHM13'].isin(['Novel(single_mn)', 'Novel(dimeric_mns)'])].copy()
    print(filtereddf.shape)
    hors = list(filtereddf['reorder_hor'])

    hor_pantype = dict(zip(filtereddf['reorder_hor'], filtereddf['pan_hor_type']))

    finaldf = df[df['reorder_hor'].isin(hors)].copy()
    print("finaldf reorder_hor num :", finaldf['reorder_hor'].nunique())

    finaldf['pan_hor_type_update'] = finaldf['reorder_hor'].map(hor_pantype)
    pantype = list(finaldf['pan_hor_type_update'].unique())
    print(pantype)

    checkdf = finaldf[finaldf['pan_hor_type_update'].isna()].copy()
    print(checkdf.head())

    sample_hap_horset = finaldf.groupby(['sample_hap', 'pan_hor_type_update'])['reorder_hor'].apply(set).to_dict()


    pan_type_union_dict = {pt : set() for pt in pantype}
    pan_type_intersect_dict = {pt : set() for pt in pantype}

    outf = open(outfile, "w")

    for i, sh in enumerate(order):
        for pt in pantype:
            if (sh, pt) in sample_hap_horset:
                ihorset = sample_hap_horset[(sh, pt)]
            else:
                outf.write(f"{i}\t{sh}\t{pt}\t{len(pan_type_union_dict[pt])}\t{len(pan_type_intersect_dict[pt])}\n")
                continue

            if i == 0:
                pan_type_union_dict[pt] = ihorset
                pan_type_intersect_dict[pt] = ihorset
            else:
                pan_type_union_dict[pt] = pan_type_union_dict[pt] | ihorset
                pan_type_intersect_dict[pt] = pan_type_intersect_dict[pt] & ihorset

            outf.write(f"{i}\t{sh}\t{pt}\t{len(pan_type_union_dict[pt])}\t{len(pan_type_intersect_dict[pt])}\n")
    outf.close()
    return hors


def get_graph_cumsum_data(order, outfile):
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-11/cumsum/all.graph.hordecomposition.summary.final.pan_hor_dr_type.xls", sep="\t", header=0)
    stvdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/cumsum/20250519_graph_hor_final_summary.xls", sep="\t", header=0)
    filtereddf = stvdf[~stvdf['HORarray_reported_by_CHM13'].isin(['Novel(single_mn)', 'Novel(dimeric_mns)'])].copy()
    print(filtereddf.shape)
    hors = list(filtereddf['hor'])
    
    finaldf = df[df['reorder_hor'].isin(hors)].copy()

    sample_hap_horset = finaldf.groupby('sample_hap')['reorder_hor'].apply(set).to_dict()
    outf = open(outfile, "w")

    for i, sh in enumerate(order):
        ihorset = sample_hap_horset[sh]
        if i == 0:
            union_out = ihorset
            interset_out = ihorset
        else:
            union_out = union_out | ihorset
            interset_out = interset_out & ihorset
        outf.write(f"{i}\t{sh}\t{len(union_out)}\t{len(interset_out)}\n")

    outf.close()
    return hors

def get_graph_cumsum_pantype(order, outfile):
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-11/cumsum/all.graph.hordecomposition.summary.final.pan_hor_dr_type.xls", sep="\t", header=0)
    stvdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/cumsum/20250519_graph_hor_final_summary.xls", sep="\t", header=0)
    filtereddf = stvdf[~stvdf['HORarray_reported_by_CHM13'].isin(['Novel(single_mn)', 'Novel(dimeric_mns)'])].copy()
    print(filtereddf.shape)
    hors = list(filtereddf['reorder_hor'])
    pantype = list(df['pan_hor_type'].unique())
    print(pantype)
    # hor_pantype = dict(zip(filtereddf['reorder_hor'], filtereddf['pan_hor_type']))

    finaldf = df[df['reorder_hor'].isin(hors)].copy()

    sample_hap_horset = finaldf.groupby(['sample_hap', 'pan_hor_type'])['reorder_hor'].apply(set).to_dict()


    pan_type_union_dict = {pt : set() for pt in pantype}
    pan_type_intersect_dict = {pt : set() for pt in pantype}

    outf = open(outfile, "w")

    for i, sh in enumerate(order):
        for pt in pantype:
            if (sh, pt) in sample_hap_horset:
                ihorset = sample_hap_horset[(sh, pt)]
            else:
                outf.write(f"{i}\t{sh}\t{pt}\t{len(pan_type_union_dict[pt])}\t{len(pan_type_intersect_dict[pt])}\n")
                continue

            if i == 0 or (pan_type_union_dict[pt] == set() and pan_type_intersect_dict[pt] == set()):
                pan_type_union_dict[pt] = ihorset
                pan_type_intersect_dict[pt] = ihorset
            else:
                pan_type_union_dict[pt] = pan_type_union_dict[pt] | ihorset
                pan_type_intersect_dict[pt] = pan_type_intersect_dict[pt] & ihorset

            outf.write(f"{i}\t{sh}\t{pt}\t{len(pan_type_union_dict[pt])}\t{len(pan_type_intersect_dict[pt])}\n")
    outf.close()
    return hors


def get_hicat_venn_data(hors):
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-11/cumsum/cumsum.order.xls", sep="\t", header=0)
    project_haploid = df.groupby('project')['haploid'].apply(list).to_dict()
    print(project_haploid['REF'])

    hicat_df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-10/statistics/all.HiCAT.hor.summary.final.xls", sep="\t", header=0)
    hicat_df = hicat_df[hicat_df['reorder_hor'].isin(hors)].copy()

    sample_hap_horset = hicat_df.groupby('sample_hap')['reorder_hor'].apply(set).to_dict()
    for project, haploids in project_haploid.items():
        tmp_set = set()
        for hap in haploids:
            tmp_set = tmp_set | sample_hap_horset[hap]
        
        with open(f"hicat_{project}_uniq_hors.xls", "w") as outf:
            for ihor in list(tmp_set):
                outf.write(f"{ihor}\n")

def get_graph_venn_data(hors):
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-11/cumsum/cumsum.order.xls", sep="\t", header=0)
    project_haploid = df.groupby('project')['haploid'].apply(list).to_dict()
    print(project_haploid['REF'])

    graph_df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-10/statistics/all.graph.hordecomposition.summary.final.pan_hor_dr_type.xls", sep="\t", header=0)
    graph_df = graph_df[graph_df['reorder_hor'].isin(hors)].copy()
    sample_hap_horset = graph_df.groupby('sample_hap')['reorder_hor'].apply(set).to_dict()

    for project, haploids in project_haploid.items():
        tmp_set = set()
        for hap in haploids:
            tmp_set = tmp_set | sample_hap_horset[hap]
        
        with open(f"graph_{project}_uniq_hors.xls", "w") as outf:
            for ihor in list(tmp_set):
                outf.write(f"{ihor}\n")



def main_curve():
    orders = load_order()

    for i in range(len(orders)):
        order = orders[i]
        if i == 0:
            # hicat_outfile = "hicat_cumsum_ordered_count.xls"
            # graph_outfile = "graph_cumsum_ordered_count.xls"
            continue

        elif i == 1:
            # hicat_outfile = "hicat_apg_first_cumsum_ordered_count.xls"
            # graph_outfile = "graph_apg_first_cumsum_ordered_count.xls"
            continue

        else:
            hicat_outfile = "hicat_hprc_first_cumsum_ordered_count_pantype.xls"
            graph_outfile = "graph_hprc_first_cumsum_ordered_count_pantype.xls"


        # hicat_hors = get_hicat_cumsum_data(order, hicat_outfile)
        # hormon_hors = get_graph_cumsum_data(order, graph_outfile)
        hicat_pantype_hors = get_hicat_cumsum_pantype(order, hicat_outfile)
        hormon_pantype_hors = get_graph_cumsum_pantype(order, graph_outfile)

    
    get_hicat_venn_data(hicat_pantype_hors)
    get_graph_venn_data(hormon_pantype_hors)


if __name__ == "__main__":
    main_curve()
