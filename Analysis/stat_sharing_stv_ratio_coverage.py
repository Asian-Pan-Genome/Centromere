import pandas as pd

def stat_sharing_array(infile, method="hicat"):
    df = pd.read_csv(infile, sep='\t', header=0)
    filter_types = ["Novel(single_mn)", "Novel(dimeric_mns)"]
    filterdf = df[~df['HORarray_reported_by_CHM13'].isin(filter_types)].copy()

    filter_stat= filterdf.groupby("horclass").agg(
        unique_chroms=('chrom_manual_check', lambda x: "_".join(sorted(set(x))))
    ).reset_index()

    print(filter_stat.head())

    shared_tmp = filter_stat[filter_stat['unique_chroms'].str.contains('_')].copy()
    shared_arrays = shared_tmp['horclass'].unique().tolist()
    print("Total shared HORarray count:", len(shared_arrays))
    print(shared_arrays)

    shareddf = filterdf[filterdf['horclass'].isin(shared_arrays)].copy()
    shareddf.loc[:, 'shared_flag'] = filterdf['chrom_manual_check'].apply(
        lambda x: 1 if '_' in x  else 0
    )
    print(shareddf.head())
    sharedstvs = shareddf[shareddf['shared_flag'] == 1]['reorder_hor'].unique().tolist()



    centdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", sep="\t", header=0)
    complete_df = centdf[centdf['filterflag'] == 0].copy()

    ###stat number of haploid 
    if method == "hormon":
        horclass_length_file = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/horclass_length/all.HORmon.horclass.length.xls"
    else:
        horclass_length_file = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/horclass_length/all.HiCAT.horclass.length.xls"
    
    horclass_length_df = pd.read_csv(horclass_length_file, sep="\t", names = ['sample_hap_chrom', 'horclass', 'total_length'])
    shared_array_info = horclass_length_df.merge(centdf[['sample_hap_chrom','sample_hap', 'filterflag']], on='sample_hap_chrom', how='left')
    shared_array_info = shared_array_info[(shared_array_info['filterflag'] == 0) &
                                          (shared_array_info['horclass'].isin(shared_arrays))].copy()
    shared_array_haploids = shared_array_info.groupby('horclass').agg(
        shared_haploids_num =('sample_hap', 'nunique')).reset_index()
    shared_array_haploids.to_csv(f"all_{method}_shared_array_haploids_num.xls", sep="\t", header=True, index=False)




    # all_shared_stv_stat_dfs = []
    # for index, row in complete_df.iterrows():
    #     sample = row['sample']
    #     hap = row['hap']
    #     chrom = row['chrom'] 
    #     project = row['project']

    #     if project == "HGSVC" or project == "HPRC":
    #         sample_hap_chrom = f"{sample}_{hap}_{chrom}"
    #     elif sample == "CHM13":
    #         sample_hap_chrom = chrom
    #     else:
    #         sample_hap_chrom = f"{sample}#{hap}_{chrom}"

    #     if method == "hormon":
    #         stv_file = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HORmon_split/{sample}_{hap}.graph.hordecomposition.final.xls"
    #         stv_df = pd.read_csv(stv_file, sep="\t", names = ['sample_hap_chrom', 'start', 'end', 'hor', 'horclass', 'horclass_color', 'hor_color'])
    #     else:
    #         stv_file = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HiCAT_bhlayer_update_split/{sample}_{hap}.HiCAT.horstv.bed"
    #         stv_df = pd.read_csv(stv_file, sep="\t", names = ['sample_hap_chrom', 'start', 'end', 'hor', 'horclass', 'horclass_color', 'stv_index', 'hor_color'])

    #     shared_array_df = stv_df[stv_df['horclass'].isin(shared_arrays) & (stv_df['sample_hap_chrom'] == sample_hap_chrom)].copy()
    #     shared_stv_df = stv_df[stv_df['hor'].isin(shared_horstvs) & (stv_df['sample_hap_chrom'] == sample_hap_chrom) & (stv_df['horclass'] != "-")].copy()
    #     shared_stv_df.loc[:, "length"] = shared_stv_df['end'] - shared_stv_df['start']

    #     shared_stv_stat = shared_stv_df.groupby(["sample_hap_chrom", "horclass"]).agg(
    #         stv_count=('hor', 'nunique'),
    #         shared_length=('length', 'sum')
    #     ).reset_index()  




def stat_sharing_stv(infile, method="hicat"):
    df = pd.read_csv(infile, sep='\t', header=0)

    filter_types = ["Novel(single_mn)", "Novel(dimeric_mns)"]
    # target_shared_chroms = ["chr13", "chr14", "chr21", "chr22"]
    # target_shared_chroms = ["chr1", "chr5", "chr19", "chr16"]
    # target_shared_chroms = ['chr18', 'chr20']

    # chrom_regex = r'\b(?:' + '|'.join(target_shared_chroms) + r')(?:_(?:' + '|'.join(target_shared_chroms) + r'))*\b'
    # print(chrom_regex)

    filterdf = df[ (~df['HORarray_reported_by_CHM13'].isin(filter_types)) & 
                   (df['chrom_manual_check'].str.contains(chrom_regex, na=False)) 
                 ].copy()

    filterdf.loc[:, 'shared_flag'] = filterdf['chrom_manual_check'].apply(
        lambda x: 0 if x in target_shared_chroms else 1
    )
    print(filterdf['chrom_manual_check'].unique().tolist())

    shared_df = filterdf[filterdf['shared_flag'] == 1].copy()
    print("Total shared HORstv count:", len(shared_df))
    shared_horstvs = shared_df['reorder_hor'].unique().tolist()

    filter_stat = filterdf.groupby("horclass").apply(
        lambda group: pd.Series({
            "horstv_count": group["reorder_hor"].nunique(),
            "shared_count": group["shared_flag"].sum()
        }), include_groups=False
    ).reset_index()

    shared_stat = filter_stat[filter_stat['shared_count'] != 0].copy()
    shared_stat.loc[:, 'shared_ratio']  = shared_stat['shared_count'] / shared_stat['horstv_count']

    # shared_stat.to_csv(f"chr13-14-21-22_{method}_shared_horstv_ratio.xls", sep="\t", header=True, index=False)
    # shared_stat.to_csv(f"chr1-5-19-16_{method}_shared_horstv_ratio.xls", sep="\t", header=True, index=False)
    shared_stat.to_csv(f"chr18-20_{method}_shared_horstv_ratio.xls", sep="\t", header=True, index=False)
    print(shared_stat)

    all_shared_stv_stat_dfs = []

    centdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", sep="\t", header=0)
    complete_df = centdf[centdf['chrom'].isin(target_shared_chroms) & (centdf['filterflag'] == 0)].copy()

    for index, row in complete_df.iterrows():
        sample = row['sample']
        hap = row['hap']
        chrom = row['chrom'] 
        project = row['project']

        if project == "HGSVC" or project == "HPRC":
            sample_hap_chrom = f"{sample}_{hap}_{chrom}"
        elif sample == "CHM13":
            sample_hap_chrom = chrom
        else:
            sample_hap_chrom = f"{sample}#{hap}_{chrom}"

        if method == "hormon":
            stv_file = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HORmon_split/{sample}_{hap}.graph.hordecomposition.final.xls"
            stv_df = pd.read_csv(stv_file, sep="\t", names = ['sample_hap_chrom', 'start', 'end', 'hor', 'horclass', 'horclass_color', 'hor_color'])
        else:
            stv_file = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HiCAT_bhlayer_update_split/{sample}_{hap}.HiCAT.horstv.bed"
            stv_df = pd.read_csv(stv_file, sep="\t", names = ['sample_hap_chrom', 'start', 'end', 'hor', 'horclass', 'horclass_color', 'stv_index', 'hor_color'])


        shared_stv_df = stv_df[stv_df['hor'].isin(shared_horstvs) & (stv_df['sample_hap_chrom'] == sample_hap_chrom) & (stv_df['horclass'] != "-")].copy()
        shared_stv_df.loc[:, "length"] = shared_stv_df['end'] - shared_stv_df['start']

    #     shared_stv_df.to_csv(f"chr1-5-19-16/{sample}_{hap}_{chrom}_{method}_shared_stv.bed", sep="\t", header=True, index=False)
    #     shared_stv_df.to_csv(f"chr13-14-21-22/{sample}_{hap}_{chrom}_{method}_shared_stv.bed", sep="\t", header=True, index=False)
        shared_stv_df.to_csv(f"chr18-20/{sample}_{hap}_{chrom}_{method}_shared_stv.bed", sep="\t", header=True, index=False)

        shared_stv_stat = shared_stv_df.groupby(["sample_hap_chrom", "horclass"]).agg(
            stv_count=('hor', 'nunique'),
            shared_length=('length', 'sum')
        ).reset_index()  

        all_shared_stv_stat_dfs.append(shared_stv_stat)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

    all_shared_stv_stat_df = pd.concat(all_shared_stv_stat_dfs, ignore_index=True)
    all_shared_stv_stat_df['horclass'] = all_shared_stv_stat_df['horclass'].astype(int)
    
    if method == "hormon":
        horclass_length_file = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/horclass_length/all.HORmon.horclass.length.xls"
    else:
        horclass_length_file = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/horclass_length/all.HiCAT.horclass.length.xls"
    horclass_length_df = pd.read_csv(horclass_length_file, sep="\t", names = ['sample_hap_chrom', 'horclass', 'total_length'])

    all_shared_stv_stat_df = pd.merge(all_shared_stv_stat_df, horclass_length_df, on=['sample_hap_chrom', 'horclass'], how='left')
    all_shared_stv_stat_df.loc[:, 'shared_ratio'] = all_shared_stv_stat_df['shared_length'] / all_shared_stv_stat_df['total_length']
    # all_shared_stv_stat_df.to_csv(f"chr1-5-19-16_{method}_shared_stv_length.xls", sep="\t", header=True, index=False)
    # all_shared_stv_stat_df.to_csv(f"chr13-14-21-22_{method}_shared_stv_length.xls", sep="\t", header=True, index=False)
    all_shared_stv_stat_df.to_csv(f"chr18-20_{method}_shared_stv_length.xls", sep="\t", header=True, index=False)

def stat_number_of_shared_haploids(infile, summary_out, final_out):

    
    centdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", sep="\t", header=0)
    complete_df = centdf[centdf['filterflag'] == 0].copy()
    print(complete_df.head())

    shared_stat_df = pd.read_csv(infile, sep="\t", header=0)
    shared_stat_df['sample_hap_chrom'] = shared_stat_df['sample_hap_chrom'].apply(
        lambda x: f"CHM13#{x}" if str(x).startswith('chr') else x
    )
    print(shared_stat_df.head())

    merged_df = pd.merge(shared_stat_df, complete_df, on='sample_hap_chrom', how='left')
    print(merged_df.head())

    summary_df = merged_df.groupby(['sample_hap', 'horclass']).agg(
        unique_chroms=('chrom', lambda x: "_".join(sorted(set(x))))
    ).reset_index()
    print(summary_df.head())
    summary_df.to_csv(summary_out, sep="\t", header=True, index=False)

    final_df = summary_df.groupby(['horclass', 'unique_chroms']).agg(
        shared_haploids_num =('sample_hap', 'nunique')).reset_index()
    print(final_df.head())
    final_df.to_csv(final_out, sep="\t", header=True, index=False)



if __name__ == "__main__":
    stat_sharing_array('all_hicat_horstv_summary_20250807.xls', method="hicat")
    stat_sharing_array('all_graph_horstv_summary_20250807.xls', method="hormon")
    # stat_sharing_stv('all_hicat_horstv_summary_20250807.xls', method="hicat")
    # stat_sharing_stv('all_graph_horstv_summary_20250807.xls', method="hormon")

    # methods = ["hormon", "hicat"]
    # for i in methods:
        # stat_number_of_shared_haploids(f"chr13-14-21-22_{i}_shared_stv_length.xls", f"chr13-14-21-22_{i}_shared_haploids_info.xls", f"chr13-14-21-22_{i}_shared_haploids_stat.xls")
        # stat_number_of_shared_haploids(f"chr1-5-19-16_{i}_shared_stv_length.xls", f"chr1-5-19-16_{i}_shared_haploids_info.xls", f"chr1-5-19-16_{i}_shared_haploids_stat.xls")
        # stat_number_of_shared_haploids(f"chr18-20_{i}_shared_stv_length.xls", f"chr18-20_{i}_shared_haploids_info.xls", f"chr18-20_{i}_shared_haploids_stat.xls")


