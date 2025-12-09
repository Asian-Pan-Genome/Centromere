import pandas as pd

def stat_target_shared_stv():
    target_stvs = [
    "21663_29601",
    "21663_29601_21663_29759",
    "21391_29375_21533_29716_21663_29601_21663_29453",
    "21391_29453_21663_29601_21663_29716_21533_29375",
    "21377_29375_21533_29716_21663_29601_21663_29179",
    "21368_29375_21533_29716_21663_29601_21663_29179",
    "21368_29375_21533_29716_21663_29601_21663_29179_21377_29254_21663_29179",
    "21377_29254_21663_29179_21377_29375_21533_29716_21663_29601_21663_29179",
    "21377_29254_21663_29179",
    "21377_29375_21533_29601_21663_29205",
    "21368_29375_21533_29716_21663_29601_21663_29205",
    "21368_29375_21533_29601_21663_29205",
    "21363_29305_21368_29290_21878_29601_21663_29427_21448_28754",21663_29601
    "21363_29290_21878_29601_21663_29427_21448_28754"
    ]

    all_hicat_df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-10/statistics/all.HiCAT.hor.summary.final.xls", sep="\t", header=0)

    target_hicat_df = all_hicat_df[all_hicat_df['reorder_hor'].isin(target_stvs)].copy()

    print(target_hicat_df.head())

    target_hicat_df.to_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr1_5_16_19_mnphy/HiCAT_target_shared_stv_hicat_layered_info.xls", sep="\t", header=True, index=False)

    stat_repeat_count = target_hicat_df.groupby(['sample', 'hap', 'chromosome', 'reorder_hor'])['nrepeat'].sum().reset_index()
    
    print(stat_repeat_count)

    stat_repeat_count.to_csv("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr1_5_16_19_mnphy/HiCAT_target_shared_stv_hicat_layered_count.xls", sep="\t", header=True, index=False)

    
if __name__ == "__main__":
    stat_target_shared_stv()
