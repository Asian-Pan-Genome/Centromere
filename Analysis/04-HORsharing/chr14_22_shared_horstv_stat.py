import pandas as pd


##hicat##
def shared_hicat_stat():
    stv_cols = ['chrom', 'start', 'end', 'hor', 'horclass', 'horclass_color', 'horstv', 'horstv_color']
    chr14_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr14.HiCAT.horstv.bed"
    chr22_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr22.HiCAT.horstv.bed"


    chr14_df = pd.read_csv(chr14_stv_file, sep="\t", names=stv_cols)
    chr14_active_df = chr14_df[chr14_df['horclass'] == '94'].copy()
    print(chr14_active_df.head())


    chr22_df = pd.read_csv(chr22_stv_file, sep="\t", names=stv_cols)
    chr22_active_df = chr22_df[chr22_df['horclass'] == '94'].copy()
    print(chr22_active_df.head())


    chr14_active_horstvs = set(chr14_active_df['horstv'].tolist())
    chr22_active_horstvs = set(chr22_active_df['horstv'].tolist())
    shared_horstv = list(chr14_active_horstvs.intersection(chr22_active_horstvs))
    print(f"total horstv in chr14: {len(chr14_active_horstvs)}")
    print(f"total horstv in chr22: {len(chr22_active_horstvs)}")
    print(f"shared horstv: {len(shared_horstv)}")

    chr14_shared_df  = chr14_active_df[chr14_active_df['horstv'].isin(shared_horstv)].copy()
    chr22_shared_df  = chr22_active_df[chr22_active_df['horstv'].isin(shared_horstv)].copy()

    chr14_shared_stat = chr14_shared_df.groupby(['chrom', 'hor', 'horstv']).size().reset_index(name='count')
    chr22_shared_stat = chr22_shared_df.groupby(['chrom', 'hor', 'horstv']).size().reset_index(name='count')
    print(chr14_shared_stat.head())
    print(chr22_shared_stat.head())

    merged_stat = pd.concat([chr14_shared_stat, chr22_shared_stat], ignore_index=True)
    merged_stat.loc[merged_stat['chrom'] == "chr1", 'chrom'] = "CHM13#chr1"
    merged_stat.columns = ['sample_hap_chrom', 'hor', 'horstv', 'count']


    ##load cent info##
    df_cen = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", sep='\t', header=0)
    merged_stat = merged_stat.merge(df_cen, on="sample_hap_chrom", how='left')
    print(merged_stat.head())


    ##convert to wide##
    wide_stat = merged_stat.pivot_table(index=['sample_hap', 'hor', 'horstv', 'project', 'filterflag'], columns='chrom', values='count', fill_value=0).reset_index()
    wide_stat.to_csv("chr14_chr22.HiCAT.shared.horstv.stat.xls", sep='\t', index=False, header=True)
    print(wide_stat.head())



def shared_hormon_stat():
    stv_cols = ['chrom', 'start', 'end', 'hor', 'horclass', 'horclass_color',  'horstv_color']
    chr13_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr14.graph.hordecomposition.final.xls"
    chr21_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr22.graph.hordecomposition.final.xls"


    chr13_df = pd.read_csv(chr13_stv_file, sep="\t", names=stv_cols)
    chr13_active_df = chr13_df[chr13_df['horclass'] == '94'].copy()
    print(chr13_active_df.head())


    chr21_df = pd.read_csv(chr21_stv_file, sep="\t", names=stv_cols)
    chr21_active_df = chr21_df[chr21_df['horclass'] == '94'].copy()
    print(chr21_active_df.head())

    chr13_active_horstvs = set(chr13_active_df['hor'].tolist())
    chr21_active_horstvs = set(chr21_active_df['hor'].tolist())
    shared_horstv = list(chr13_active_horstvs.intersection(chr21_active_horstvs))
    print(f"total horstv in chr13: {len(chr13_active_horstvs)}")
    print(f"total horstv in chr21: {len(chr21_active_horstvs)}")
    print(f"shared horstv: {len(shared_horstv)}")

    chr13_shared_df  = chr13_active_df[chr13_active_df['hor'].isin(shared_horstv)].copy()
    chr21_shared_df  = chr21_active_df[chr21_active_df['hor'].isin(shared_horstv)].copy()

    chr13_shared_stat = chr13_shared_df.groupby(['chrom', 'hor']).size().reset_index(name='count')
    chr21_shared_stat = chr21_shared_df.groupby(['chrom', 'hor']).size().reset_index(name='count')
    print(chr13_shared_stat.head())
    print(chr21_shared_stat.head())

    merged_stat = pd.concat([chr13_shared_stat, chr21_shared_stat], ignore_index=True)
    merged_stat.loc[merged_stat['chrom'] == "chr1", 'chrom'] = "CHM13#chr1"
    merged_stat.columns = ['sample_hap_chrom', 'hor', 'count']


    ##load cent info##
    df_cen = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", sep='\t', header=0)
    merged_stat = merged_stat.merge(df_cen, on="sample_hap_chrom", how='left')
    print(merged_stat.head())


    ##convert to wide##
    wide_stat = merged_stat.pivot_table(index=['sample_hap', 'hor', 'project', 'filterflag'], columns='chrom', values='count', fill_value=0).reset_index()
    wide_stat.to_csv("chr14_chr22.hormon.shared.horstv.stat.xls", sep='\t', index=False, header=True)
    print(wide_stat.head())




def main():
    # shared_hicat_stat()
    shared_hormon_stat()


if __name__ == "__main__":
    main()

