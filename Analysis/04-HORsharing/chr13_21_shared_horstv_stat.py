import pandas as pd
from collections import Counter

##hicat##
def shared_hicat_stat():
    stv_cols = ['chrom', 'start', 'end', 'hor', 'horclass', 'horclass_color', 'horstv', 'horstv_color']
    chr13_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr13.HiCAT.horstv.bed"
    chr21_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr21.HiCAT.horstv.bed"
    chr14_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr14.HiCAT.horstv.bed"
    chr22_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr22.HiCAT.horstv.bed"


    shared_horclass = ['8', '42', '49', '62', '65', '94', '133', '134', '149', '165', '200', '245', '254', '265', '266']

    chr13_df = pd.read_csv(chr13_stv_file, sep="\t", names=stv_cols)
    # chr13_active_df = chr13_df[chr13_df['horclass'] == '42'].copy()
    chr13_active_df = chr13_df[chr13_df['horclass'].isin(shared_horclass)].copy()
    print(chr13_active_df.head())


    chr21_df = pd.read_csv(chr21_stv_file, sep="\t", names=stv_cols)
    # chr21_active_df = chr21_df[chr21_df['horclass'] == '42'].copy()
    chr21_active_df = chr21_df[chr21_df['horclass'].isin(shared_horclass)].copy()
    print(chr21_active_df.head())

    chr14_df = pd.read_csv(chr14_stv_file, sep="\t", names=stv_cols)
    # chr21_active_df = chr21_df[chr21_df['horclass'] == '42'].copy()
    chr14_active_df = chr14_df[chr14_df['horclass'].isin(shared_horclass)].copy()
    print(chr14_active_df.head())

    chr22_df = pd.read_csv(chr22_stv_file, sep="\t", names=stv_cols)
    # chr21_active_df = chr21_df[chr21_df['horclass'] == '42'].copy()
    chr22_active_df = chr22_df[chr22_df['horclass'].isin(shared_horclass)].copy()
    print(chr22_active_df.head())


    chr13_active_horstvs = set(chr13_active_df['horstv'].tolist())
    chr21_active_horstvs = set(chr21_active_df['horstv'].tolist())
    chr14_active_horstvs = set(chr14_active_df['horstv'].tolist())
    chr22_active_horstvs = set(chr22_active_df['horstv'].tolist())

    sets = [chr13_active_horstvs, chr21_active_horstvs, chr14_active_horstvs, chr22_active_horstvs]
    element_counts = Counter()
    for s in sets:
        element_counts.update(s)
    shared_horstv = [element for element, count in element_counts.items() if count >= 2]
    # shared_horstv = list(chr13_active_horstvs.intersection(chr21_active_horstvs))
    print(f"total horstv in chr13: {len(chr13_active_horstvs)}")
    print(f"total horstv in chr21: {len(chr21_active_horstvs)}")
    print(f"total horstv in chr14: {len(chr14_active_horstvs)}")
    print(f"total horstv in chr22: {len(chr22_active_horstvs)}")
    print(f"shared horstv: {len(shared_horstv)}")

    chr13_shared_df  = chr13_active_df[chr13_active_df['horstv'].isin(shared_horstv)].copy()
    chr21_shared_df  = chr21_active_df[chr21_active_df['horstv'].isin(shared_horstv)].copy()
    chr14_shared_df  = chr14_active_df[chr14_active_df['horstv'].isin(shared_horstv)].copy()
    chr22_shared_df  = chr22_active_df[chr22_active_df['horstv'].isin(shared_horstv)].copy()

    chr13_shared_df['activeflag']  = chr13_shared_df['horclass'].apply(lambda x: 'active' if x in ['42', '94'] else 'inactive')
    chr21_shared_df['activeflag']  = chr21_shared_df['horclass'].apply(lambda x: 'active' if x in ['42', '94'] else 'inactive')
    chr14_shared_df['activeflag']  = chr14_shared_df['horclass'].apply(lambda x: 'active' if x in ['42', '94'] else 'inactive')
    chr22_shared_df['activeflag']  = chr22_shared_df['horclass'].apply(lambda x: 'active' if x in ['42', '94'] else 'inactive')


    chr13_shared_stat = chr13_shared_df.groupby(['chrom', 'hor', 'horstv', 'activeflag']).size().reset_index(name='count')
    chr21_shared_stat = chr21_shared_df.groupby(['chrom', 'hor', 'horstv', 'activeflag']).size().reset_index(name='count')
    chr14_shared_stat = chr14_shared_df.groupby(['chrom', 'hor', 'horstv', 'activeflag']).size().reset_index(name='count')
    chr22_shared_stat = chr22_shared_df.groupby(['chrom', 'hor', 'horstv', 'activeflag']).size().reset_index(name='count')
    print(chr13_shared_stat.head())
    print(chr21_shared_stat.head())
    print(chr14_shared_stat.head())
    print(chr22_shared_stat.head())

    merged_stat = pd.concat([chr13_shared_stat, chr21_shared_stat, chr14_shared_stat, chr22_shared_stat], ignore_index=True)
    merged_stat.loc[merged_stat['chrom'] == "chr13", 'chrom'] = "CHM13#chr13"
    merged_stat.loc[merged_stat['chrom'] == "chr21", 'chrom'] = "CHM13#chr21"
    merged_stat.loc[merged_stat['chrom'] == "chr14", 'chrom'] = "CHM13#chr14"
    merged_stat.loc[merged_stat['chrom'] == "chr22", 'chrom'] = "CHM13#chr22"
    merged_stat.columns = ['sample_hap_chrom', 'hor', 'horstv', 'activeflag', 'count']


    ##load cent info##
    df_cen = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", sep='\t', header=0)
    merged_stat = merged_stat.merge(df_cen, on="sample_hap_chrom", how='left')
    print(merged_stat.head())


    ##convert to wide##
    wide_stat = merged_stat.pivot_table(index=['sample_hap', 'hor', 'horstv', 'project', 'activeflag', 'filterflag'], columns='chrom', values='count', fill_value=0).reset_index()
    wide_stat.to_csv("chr13_chr21_chr14_chr22.HiCAT.shared.horstv.stat.xls", sep='\t', index=False, header=True)
    print(wide_stat.head())



def shared_hormon_stat():
    stv_cols = ['chrom', 'start', 'end', 'hor', 'horclass', 'horclass_color',  'horstv_color']
    chr13_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr13.graph.hordecomposition.final.xls"
    chr21_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr21.graph.hordecomposition.final.xls"
    chr14_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr14.graph.hordecomposition.final.xls"
    chr22_stv_file = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/chr22.graph.hordecomposition.final.xls"

    shared_horclass = ['8', '42', '49', '62', '65', '94', '133', '134', '149', '165',  '245', '254']

    chr13_df = pd.read_csv(chr13_stv_file, sep="\t", names=stv_cols)
    chr13_active_df = chr13_df[chr13_df['horclass'].isin(shared_horclass)].copy()
    print(chr13_active_df.head())


    chr21_df = pd.read_csv(chr21_stv_file, sep="\t", names=stv_cols)
    chr21_active_df = chr21_df[chr21_df['horclass'].isin(shared_horclass)].copy()
    print(chr21_active_df.head())

    chr14_df = pd.read_csv(chr14_stv_file, sep="\t", names=stv_cols)
    chr14_active_df = chr14_df[chr14_df['horclass'].isin(shared_horclass)].copy()
    print(chr14_active_df.head())

    chr22_df = pd.read_csv(chr22_stv_file, sep="\t", names=stv_cols)
    chr22_active_df = chr22_df[chr22_df['horclass'].isin(shared_horclass)].copy()
    print(chr22_active_df.head())


    chr13_active_horstvs = set(chr13_active_df['hor'].tolist())
    chr21_active_horstvs = set(chr21_active_df['hor'].tolist())
    chr14_active_horstvs = set(chr14_active_df['hor'].tolist())
    chr22_active_horstvs = set(chr22_active_df['hor'].tolist())

    sets = [chr13_active_horstvs, chr21_active_horstvs, chr14_active_horstvs, chr22_active_horstvs]
    element_counts = Counter()
    for s in sets:
        element_counts.update(s)
    shared_horstv = [element for element, count in element_counts.items() if count >= 2]

    print(f"total horstv in chr13: {len(chr13_active_horstvs)}")
    print(f"total horstv in chr21: {len(chr21_active_horstvs)}")
    print(f"total horstv in chr14: {len(chr14_active_horstvs)}")
    print(f"total horstv in chr22: {len(chr22_active_horstvs)}")
    print(f"shared horstv: {len(shared_horstv)}")


    chr13_shared_df  = chr13_active_df[chr13_active_df['hor'].isin(shared_horstv)].copy()
    chr21_shared_df  = chr21_active_df[chr21_active_df['hor'].isin(shared_horstv)].copy()
    chr14_shared_df  = chr14_active_df[chr14_active_df['hor'].isin(shared_horstv)].copy()
    chr22_shared_df  = chr22_active_df[chr22_active_df['hor'].isin(shared_horstv)].copy()

    chr13_shared_df['activeflag']  = chr13_shared_df['horclass'].apply(lambda x: 'active' if x in ['42', '94'] else 'inactive')
    chr21_shared_df['activeflag']  = chr21_shared_df['horclass'].apply(lambda x: 'active' if x in ['42', '94'] else 'inactive')
    chr14_shared_df['activeflag']  = chr14_shared_df['horclass'].apply(lambda x: 'active' if x in ['42', '94'] else 'inactive')
    chr22_shared_df['activeflag']  = chr22_shared_df['horclass'].apply(lambda x: 'active' if x in ['42', '94'] else 'inactive')

    chr13_shared_stat = chr13_shared_df.groupby(['chrom', 'hor', 'activeflag']).size().reset_index(name='count')
    chr21_shared_stat = chr21_shared_df.groupby(['chrom', 'hor', 'activeflag']).size().reset_index(name='count')
    chr14_shared_stat = chr14_shared_df.groupby(['chrom', 'hor', 'activeflag']).size().reset_index(name='count')
    chr22_shared_stat = chr22_shared_df.groupby(['chrom', 'hor', 'activeflag']).size().reset_index(name='count')

    print(chr13_shared_stat.head())
    print(chr21_shared_stat.head())

    merged_stat = pd.concat([chr13_shared_stat, chr21_shared_stat, chr14_shared_stat, chr22_shared_stat], ignore_index=True)
    merged_stat.loc[merged_stat['chrom'] == "chr13", 'chrom'] = "CHM13#chr13"
    merged_stat.loc[merged_stat['chrom'] == "chr21", 'chrom'] = "CHM13#chr21"
    merged_stat.loc[merged_stat['chrom'] == "chr14", 'chrom'] = "CHM13#chr14"
    merged_stat.loc[merged_stat['chrom'] == "chr22", 'chrom'] = "CHM13#chr22"
    merged_stat.columns = ['sample_hap_chrom', 'hor', 'activeflag', 'count']


    ##load cent info##
    df_cen = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", sep='\t', header=0)
    merged_stat = merged_stat.merge(df_cen, on="sample_hap_chrom", how='left')
    print(merged_stat.head())


    ##convert to wide##
    wide_stat = merged_stat.pivot_table(index=['sample_hap', 'hor', 'project', 'filterflag'], columns='chrom', values='count', fill_value=0).reset_index()
    wide_stat.to_csv("chr13_chr21_chr14_chr22.hormon.shared.horstv.stat.xls", sep='\t', index=False, header=True)
    print(wide_stat.head())




def main():
    shared_hicat_stat()
    shared_hormon_stat()


if __name__ == "__main__":
    main()

