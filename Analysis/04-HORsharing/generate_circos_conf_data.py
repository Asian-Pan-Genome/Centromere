import pandas as pd
import numpy as np

def get_hormon_karyotype():
    outf = open("hormon_karyotype.txt", "w")
    outf2 = open("hormon_horclass_highlight.txt", 'w')
    outf3 = open("hormon_horclass_sf_highlight.txt", 'w')
    
    infile = "/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/circos/hormon_hor_num_each_chrom.xls"
    chroms = []

    sf_df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/circos/hormon_horclass_sf.xls", sep="\t", header=0)
    sf_df['hicat_horclass'] = sf_df['hicat_horclass'].astype(str)
    horclass_sfcolor = dict(zip(sf_df['hicat_horclass'], sf_df['SF_color']))

    with open(infile, 'r') as inf:
        for line in inf:
            if line.startswith('chrom_manual_check'):
                continue
            tokens = line.strip().split("\t")
            ichr = tokens[0]
            if ichr not in chroms:
                if chroms != []:
                    outf.write(f"chr\t-\ths{name}\t{name}\t0\t{prev_end}\tchr{name}\n")

                start = 0
                end = int(tokens[2])
                name = ichr.replace("chr", "")
                anno = tokens[3]
                if "Novel" in anno:
                    outf2.write(f"hs{name}\t{start}\t{end}\tfill_color={tokens[1]}_color,r0=0.55r,r1=0.6r\n")
                else:
                    outf2.write(f"hs{name}\t{start}\t{end}\tfill_color={tokens[1]}_color,r0=0.5r,r1=0.55r\n")
                outf3.write(f"hs{name}\t{start}\t{end}\tfill_color={horclass_sfcolor[tokens[1]]},r0=0.61r,r1=0.63r\n")
                prev_end = end
                chroms.append(ichr)
            else:
                start = prev_end
                end = start + int(tokens[2])
                name = ichr.replace("chr", "")
                anno = tokens[3]
                if "Novel" in anno:
                    outf2.write(f"hs{name}\t{start}\t{end}\tfill_color={tokens[1]}_color,r0=0.55r,r1=0.6r\n")
                else:
                    outf2.write(f"hs{name}\t{start}\t{end}\tfill_color={tokens[1]}_color,r0=0.5r,r1=0.55r\n")
                outf3.write(f"hs{name}\t{start}\t{end}\tfill_color={horclass_sfcolor[tokens[1]]},r0=0.61r,r1=0.63r\n")
                prev_end = end
            
        outf.write(f"chr\t-\ths{name}\t{name}\t0\t{prev_end}\tchr{name}\n")

    outf.close()
    outf2.close()
    outf3.close()


def get_hormon_link():

    hormon_hor_df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/20250519_graph_hor_final_summary.xls", sep="\t", header=0)
    hors = list(hormon_hor_df[(hormon_hor_df['HORarray_reported_by_CHM13'] != "Novel(dimeric_mns)")]['hor'])
    hors_manual_check_chroms = hormon_hor_df.set_index('hor')['chrom_manual_check'].str.split('_').to_dict()
    print(hors_manual_check_chroms['400_6649_6649_22781_23272_6435_23278_23057_5288_23060_5368_22856_4168_23076_6421_22505_17496_6421_22508_17499_23053_550'])

    
    decomposition_cols = ['sample', 'hap', 'chrom', 'blockid', 'horindex_start', 'horindex_end', 'sample_hap_chrom', 'start', 'end', 'hor', 'reorder_hor', 'nmer']
    decomposition_df  = pd.read_csv("~/human/alpha/vsearch/20241022/step-10/statistics/all.graph.hordecomposition.summary.v2.xls", names = decomposition_cols, sep='\t')

    target_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
    print(target_chroms)
    decomposition_hor_df = decomposition_df.loc[decomposition_df['reorder_hor'].isin(hors) & (decomposition_df['chrom'].isin(target_chroms))].copy()

    hor_repeat_in_each_haploid_df = decomposition_hor_df.groupby(['sample', 'hap', 'sample_hap_chrom', 'chrom', 'reorder_hor']).size().reset_index(name='count')
    print("hor_repeat_in_each_haploid_df round1: ", hor_repeat_in_each_haploid_df.shape)
    print(hor_repeat_in_each_haploid_df.head())

    ## select chroms for manual check##
    hor_repeat_in_each_haploid_df = hor_repeat_in_each_haploid_df[
    hor_repeat_in_each_haploid_df.apply(
        lambda row: row['chrom'] in hors_manual_check_chroms.get(row['reorder_hor'], []),
        axis=1
    )]
    print("hor_repeat_in_each_haploid_df round2: ", hor_repeat_in_each_haploid_df.shape)
    print(hor_repeat_in_each_haploid_df.head())

    ## add hicat_horclass ##
    hor_repeat_in_each_haploid_df = hor_repeat_in_each_haploid_df.merge(
        hormon_hor_df[['hor', 'hicat_horclass']],
        how='left',
        left_on='reorder_hor',
        right_on='hor'
    )
    print("hor_repeat_in_each_haploid_df round3: ", hor_repeat_in_each_haploid_df.shape)
    print(hor_repeat_in_each_haploid_df.head())

    # add population ##
    popdf = pd.read_csv("~/human/anno/statistics/centhap/populaion.xls", sep="\t", header=0)
    popdf['superpopulation'] = popdf['superpopulation'].str.replace("EAS-APG", "EAS")

    hor_repeat_in_each_haploid_df = hor_repeat_in_each_haploid_df.merge(
        popdf[['sample', 'superpopulation']],
        how='left',
        on = "sample"
    )
    print(hor_repeat_in_each_haploid_df.head())



    # hor_maxrepeat_in_chrom_df = hor_repeat_in_each_haploid_df.groupby(['chrom', 'reorder_hor', 'hicat_horclass'])['count'].max().reset_index(name = 'max_count')
    # hormon_order_df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/circos/hormon_hor_num_each_chrom.xls", sep="\t", header=0)
    # hormon_order_df['sort_order'] = range(len(hormon_order_df))
    # hormon_order_df.to_csv("hormon_hor_num_each_chrom_order.xls", sep="\t", header=True, index=False)

    # hor_maxrepeat_in_chrom_df = hor_maxrepeat_in_chrom_df.merge(
    # hormon_order_df[['chrom_manual_check', 'hicat_horclass', 'sort_order']],
    # left_on=['chrom', 'hicat_horclass'],
    # right_on=['chrom_manual_check', 'hicat_horclass'],
    # how='left')

    # hor_maxrepeat_in_chrom_df = hor_maxrepeat_in_chrom_df.sort_values(
    # by=['sort_order', 'max_count'],  # 排序列
    # ascending=[True, False]          # sort_order 升序, max_count 降序
    # )

    # hor_maxrepeat_in_chrom_df['index'] = hor_maxrepeat_in_chrom_df.groupby('chrom').cumcount() + 1
    # hor_maxrepeat_in_chrom_df = hor_maxrepeat_in_chrom_df[
    # hor_maxrepeat_in_chrom_df.apply(
    #     lambda row: row['chrom'] in hors_manual_check_chroms.get(row['reorder_hor'], []),
    #     axis=1
    # )]
    # hor_maxrepeat_in_chrom_df.to_csv("hormon_hor_index_max_count_in_each_chrom.xls", sep="\t", header=True, index=False)

    # hor_repeat_in_each_haploid_df = hor_repeat_in_each_haploid_df.merge(
    # hor_maxrepeat_in_chrom_df[['reorder_hor', 'chrom', 'index']],
    # on=['reorder_hor', 'chrom'],
    # how='left')
    # print("hor_repeat_in_each_haploid_df round4: ", hor_repeat_in_each_haploid_df.shape)
    
    # hor_repeat_in_each_haploid_df = hor_repeat_in_each_haploid_df.sort_values(by=['index', 'count'], ascending=[True, False])
    # hor_repeat_in_each_haploid_df.to_csv("hormon_hor_index_count_in_haploid.xls", sep="\t", header=True, index=False )

    # for i, (horclass, subdf) in enumerate(hor_maxrepeat_in_chrom_df.groupby('hicat_horclass')):

    #     if subdf['chrom'].nunique() == 1:
    #         continue
        

    #     outpair = open(f"hormon.{horclass}.link.txt", 'w')

    #     for j, (reorder_hor, ssubdf)  in enumerate(subdf.groupby('reorder_hor')):
    #         if len(ssubdf) != 1:
    #             # print(ssubdf.head())
    #             for k in range(len(ssubdf)):
    #                 for t in range(k + 1, len(ssubdf)):
    #                     row_i = ssubdf.iloc[k]
    #                     row_j = ssubdf.iloc[t]
    #                     chrom_i = row_i['chrom'].replace("chr", "")
    #                     chrom_j = row_j['chrom'].replace("chr", "")
    #                     index_i = row_i['index']
    #                     index_j = row_j['index']

    #                     outpair.write(f"hs{chrom_i}\t{index_i}\t{index_i}\ths{chrom_j}\t{index_j}\t{index_j}\n")
    #     outpair.close()

    ##generate haploid shared ##
    # hor_haploid_shared_each_chrom = hor_repeat_in_each_haploid_df.groupby(['chrom', 'reorder_hor', 'index'])['sample_hap_chrom'].nunique().reset_index(name='shared_count')
    # hor_haploid_shared_each_chrom['chrom'] = hor_haploid_shared_each_chrom['chrom'].str.replace("chr", "hs")
    # hor_haploid_shared_each_chrom['log_shared_count'] = np.log10(hor_haploid_shared_each_chrom['shared_count']) + 0.1
    # hor_haploid_shared_each_chrom['end'] = hor_haploid_shared_each_chrom['index']
    # out_haploid_shared_each_chrom = hor_haploid_shared_each_chrom[['chrom', 'index', 'end', 'log_shared_count']]
    # out_haploid_shared_each_chrom.to_csv("hormon_hor_haploid_shared_each_chrom.xls", sep="\t", header=False, index=False)


    ##generate horstv repeat number##
    # hor_repeat_scatter_df = hor_repeat_in_each_haploid_df[['chrom', 'index', 'count']].copy()
    # hor_repeat_scatter_df['chrom'] = hor_repeat_scatter_df['chrom'].str.replace("chr", "hs")
    # hor_repeat_scatter_df['end'] = hor_repeat_scatter_df['index']

    # max_count_df = hor_repeat_scatter_df.groupby(['chrom', 'index', 'end'])['count'].max().reset_index()
    # max_count_df.rename(columns={'count': 'max_count'}, inplace=True)
    # max_count_df['log_max_count'] = np.log10(max_count_df['max_count']) + 0.1
    # out_max_count_df = max_count_df[['chrom', 'index', 'end', 'log_max_count']]
    # out_max_count_df.to_csv("hormon_hor_max_repeat_scatter.xls", sep="\t", header=False, index=False)

    # median_count_df = hor_repeat_scatter_df.groupby(['chrom', 'index', 'end'])['count'].median().reset_index()
    # median_count_df.rename(columns={'count': 'median_count'}, inplace=True)
    # median_count_df['log_median_count'] = np.log10(median_count_df['median_count']) + 0.1
    # out_median_count_df = median_count_df[['chrom', 'index', 'end', 'log_median_count']]
    # out_median_count_df.to_csv("hormon_hor_median_repeat_scatter.xls", sep="\t", header=False, index=False)

    # min_count_df = hor_repeat_scatter_df.groupby(['chrom', 'index', 'end'])['count'].min().reset_index()
    # min_count_df.rename(columns={'count': 'min_count'}, inplace=True)
    # min_count_df['log_min_count'] = np.log10(min_count_df['min_count']) + 0.1
    # out_min_count_df = min_count_df[['chrom', 'index', 'end', 'log_min_count']]
    # out_min_count_df.to_csv("hormon_hor_min_repeat_scatter.xls", sep="\t", header=False, index=False)

    ##generate horstv superpopulation frequncy##
    # Group by chrom, reorder_hor, and superpopulation to calculate the total count for each group
    # superpop_counts = hor_repeat_in_each_haploid_df.groupby(['chrom', 'reorder_hor', 'superpopulation', 'index'])['sample_hap_chrom'].nunique().reset_index(name='superpop_shared_count')
    # print(superpop_counts.head())
    # hor_repeat_in_each_haploid_df['sample_hap'] = hor_repeat_in_each_haploid_df['sample'] + "_" + hor_repeat_in_each_haploid_df['hap']
    # total_superpop_counts = hor_repeat_in_each_haploid_df.groupby(['superpopulation'])['sample_hap'].nunique().reset_index(name ='total_superpop_count')

    # # Merge the two DataFrames to calculate the percentage
    # percentage_df = superpop_counts.merge(total_superpop_counts, on='superpopulation', how='left')
    # percentage_df['percentage'] = percentage_df['superpop_shared_count'] / percentage_df['total_superpop_count']
    # percentage_df.to_csv("hormon_hor_superpopulation_percentage.xls", sep="\t", header=True, index=False)

    ##compare EAS vs AFR
    # wide_percentage_df = percentage_df.pivot(index=['chrom', 'reorder_hor', 'index'], columns='superpopulation', values='percentage').fillna(0)
    # wide_percentage_df.reset_index(inplace=True)
    # wide_percentage_df['EAS_vs_AFR'] = wide_percentage_df['EAS'] - wide_percentage_df['AFR']
    # wide_percentage_df.to_csv("hormon_hor_AFR_vs_EAS.xls", sep="\t", header=True, index=False)
    # wide_percentage_df['chrom'] = wide_percentage_df['chrom'].str.replace("chr", "hs")
    # wide_percentage_df['end'] = wide_percentage_df['index']
    # out_wide_percentage_df = wide_percentage_df[['chrom', 'index', 'end', 'EAS_vs_AFR']]
    # out_wide_percentage_df.to_csv("hormon_hor_AFR_vs_EAS_scatter.xls", sep="\t", header=False, index=False)

    

    ##compare EAS vs non-EAS
    # percentage_df['superflag'] = percentage_df['superpopulation'].apply(lambda x: "EAS" if x == "EAS" else "non-EAS")
    # vs_percentage_df = (
    # percentage_df
    # .groupby(['chrom', 'reorder_hor', 'index', 'superflag'])[['superpop_shared_count', 'total_superpop_count']]
    # .sum()
    # .reset_index()
    # )
    # print(vs_percentage_df.head())

    # vs_percentage_df['EAS_percent_vs'] = vs_percentage_df['superpop_shared_count'] / vs_percentage_df['total_superpop_count']
    # wide_vs_percentage_df = vs_percentage_df.pivot(index=['chrom', 'reorder_hor', 'index'], columns='superflag', values='EAS_percent_vs').fillna(0)
    # wide_vs_percentage_df.reset_index(inplace=True)
    # print(wide_vs_percentage_df.head())
    # print(wide_vs_percentage_df.columns)
    # wide_vs_percentage_df['EAS_vs_non-EAS'] =  wide_vs_percentage_df['EAS'] - wide_vs_percentage_df['non-EAS']
    # wide_vs_percentage_df.to_csv("hormon_hor_EAS_vs_nonEAS.xls", sep="\t", header=True, index=False)
    # wide_vs_percentage_df['chrom'] = wide_vs_percentage_df['chrom'].str.replace("chr", "hs")
    # wide_vs_percentage_df['end'] = wide_vs_percentage_df['index']
    # out_wide_vs_percentage_df = wide_vs_percentage_df[['chrom', 'index', 'end', 'EAS_vs_non-EAS']]
    # out_wide_vs_percentage_df.to_csv("hormon_hor_EAS_vs_nonEAS_scatter.xls", sep="\t", header=False, index=False)


    ##out##
    # for superpop, subdf in percentage_df.groupby('superpopulation'):
    #     subdf['chrom'] = subdf['chrom'].str.replace("chr", "hs")
    #     subdf['end'] = subdf['index']
    #     out_superpop_df = subdf[['chrom', 'index', 'end', 'percentage']]
    #     out_superpop_df.to_csv(f"hormon_hor_superpopulation_{superpop}_percentage_heatmap.xls", sep="\t", header=False, index=False)
    
    ##generate horstv superpopulation repeat##
    # superpop_repeat = hor_repeat_in_each_haploid_df.groupby(['chrom', 'reorder_hor', 'superpopulation', 'index'])['count'].sum().reset_index(name='superpop_sum_repeat')
    # superpop_repeat = superpop_repeat.merge(total_superpop_counts, on='superpopulation', how='left')
    # superpop_repeat['superpop_mean_repeat'] = superpop_repeat['superpop_sum_repeat'] / superpop_repeat['total_superpop_count']
    ##compare EAS vs AFR
    # wide_superpop_repeat_df = superpop_repeat.pivot(index=['chrom', 'reorder_hor', 'index'], columns='superpopulation', values='superpop_mean_repeat').fillna(0)
    # wide_superpop_repeat_df.reset_index(inplace=True)
    # wide_superpop_repeat_df['EAS_vs_AFR'] = wide_superpop_repeat_df['EAS'] - wide_superpop_repeat_df['AFR']
    # wide_superpop_repeat_df.to_csv("hormon_hor_AFR_vs_EAS_repeat.xls", sep="\t", header=True, index=False)
    # wide_superpop_repeat_df['chrom'] = wide_superpop_repeat_df['chrom'].str.replace("chr", "hs")
    # wide_superpop_repeat_df['end'] = wide_superpop_repeat_df['index']
    

    # def process_eas_vs_afr(value):
    #     if value == 1:
    #         return 0.1
    #     elif value == -1:
    #         return -0.1
    #     elif value == 0:
    #         return 0
    #     elif value > 1:
    #         return np.log10(value)
    #     elif 0 < value < 1 or -1 < value < 0:
    #         return value
    #     elif value < -1:
    #         return -np.log10(abs(value))
    #     else:
    #         return value  

    # wide_superpop_repeat_df['log_EAS_vs_AFR'] = wide_superpop_repeat_df['EAS_vs_AFR'].apply(process_eas_vs_afr)
    # out_wide_superpop_repeat_df = wide_superpop_repeat_df[['chrom', 'index', 'end', 'log_EAS_vs_AFR']]
    # out_wide_superpop_repeat_df.to_csv("hormon_hor_AFR_vs_EAS_repeat_scatter.xls", sep="\t", header=False, index=False)


    ##compare EAS vs non-EAS
    # hor_repeat_in_each_haploid_df['superflag'] = hor_repeat_in_each_haploid_df['superpopulation'].apply(lambda x: "EAS" if x == "EAS" else "non-EAS")
    # vs_superpop_repeat = hor_repeat_in_each_haploid_df.groupby(['chrom', 'reorder_hor', 'index', 'superflag'])['count'].sum().reset_index(name='superpop_sum_vs')
    # total_superpop_counts['superflag'] = total_superpop_counts['superpopulation'].apply(lambda x: "EAS" if x == "EAS" else "non-EAS")
    # total_eas_noneas_counts = total_superpop_counts.groupby('superflag')['total_superpop_count'].sum().reset_index(name='total_superpop_count')
    # vs_superpop_repeat = vs_superpop_repeat.merge(total_eas_noneas_counts, on='superflag', how='left')
    # vs_superpop_repeat['superpop_mean_vs'] = vs_superpop_repeat['superpop_sum_vs'] / vs_superpop_repeat['total_superpop_count']

    # wide_vs_superpop_repeat_df = vs_superpop_repeat.pivot(index=['chrom', 'reorder_hor', 'index'], columns='superflag', values='superpop_mean_vs').fillna(0)
    # wide_vs_superpop_repeat_df.reset_index(inplace=True)
    # wide_vs_superpop_repeat_df['EAS_vs_nonEAS'] =  wide_vs_superpop_repeat_df['EAS'] - wide_vs_superpop_repeat_df['non-EAS']
    # wide_vs_superpop_repeat_df.to_csv("hormon_hor_EAS_vs_nonEAS_repeat.xls", sep="\t", header=True, index=False)
    # wide_vs_superpop_repeat_df['chrom'] = wide_vs_superpop_repeat_df['chrom'].str.replace("chr", "hs")
    # wide_vs_superpop_repeat_df['end'] = wide_vs_superpop_repeat_df['index']
    # out_wide_vs_superpop_repeat_df = wide_vs_superpop_repeat_df[['chrom', 'index', 'end', 'EAS_vs_nonEAS']]
    # out_wide_vs_superpop_repeat_df.to_csv("hormon_hor_EAS_vs_nonEAS_repeat_scatter.xls", sep="\t", header=False, index=False)


def generate_highlight_superpopulation_specific():

    specific_df = pd.read_csv("hormon_hor_superpopulation_specific.xls", sep="\t", header = 0)
    index_df = pd.read_csv("hormon_hor_index_max_count_in_each_chrom.xls", sep="\t", header= 0)
    print(specific_df.shape)

    merged_df = index_df.merge(specific_df, how="left", on="reorder_hor")
    merged_df = merged_df.dropna(subset=['EAS_vs_AFR_repeat_num_p-value'])
    print(merged_df.shape)
    print(merged_df.head())


    merged_df['hs'] = merged_df['chrom'].str.replace("chr", "hs")
    merged_df['hs_start'] = merged_df['index']
    merged_df['hs_end']  = merged_df['index']
    merged_df['hs_anno']  = merged_df['EAS_vs_AFR_repeat_num_p-value'].apply(
        lambda x: "fill_color=210,46,119,r0=0.62r,r1=0.64r" if x == "EAS-specific" else "fill_color=49,155,98,r0=0.66r,r1=0.68r"
    )
    print(merged_df.head())

    # print("h3-10mer:", index_df[index_df['reorder_hor'] == "2975_24578_3036_25494_3379_3072_24569_3042_25511_3379"])
    # print("h1-11mer:", index_df[index_df['reorder_hor'] == "2975_24578_3036_25494_3059_24578_3072_24569_3042_25511_3379"])
    # print("h2-11mer:", index_df[index_df['reorder_hor'] == "2975_24578_3036_25494_3059_24579_3072_24569_3042_25511_3379"])
    # print("h4-7mer:",  index_df[index_df['reorder_hor'] == "2975_24578_3072_24569_3042_25511_3379"])

    print("dimer:", index_df[index_df['reorder_hor'] == "21663_29601"])
    print("2-dimer:", index_df[index_df['reorder_hor'] == "21663_29601_21663_29759"])
    print("2-dimer-hormon:", index_df[index_df['reorder_hor'] == "21663_29759"])
    print("D-E-F-E:", index_df[index_df['reorder_hor'] == "21663_29601_21663_29716"])
    print("a2-F:", index_df[index_df['reorder_hor'] == "21391_29375_21533_29716_21663_29601"])
    print("12-mer:", index_df[index_df['reorder_hor'] == "21377_29254_21663_29179"])
    print("8-mer:",  index_df[index_df['reorder_hor'] == "21391_29375_21533_29716_21663_29601_21663_29453"])




    outdf = merged_df[['hs', 'hs_start', 'hs_end', 'hs_anno']].copy()
    outdf.to_csv("hormon_hor_superpopulation_specific_circos.txt", sep="\t", header=False, index=False)
    print(outdf.head())



if __name__ == "__main__":
    # get_hormon_link()
    # get_hormon_karyotype()
    generate_highlight_superpopulation_specific()
