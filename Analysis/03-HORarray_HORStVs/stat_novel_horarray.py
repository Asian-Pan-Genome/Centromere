import pandas as pd

def stat_novel_horarray():

    ##get novel horclass id##
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/ancient_hor_array/Novel_horclass_index.xls", sep="\t", header=None, names=["HORclass"])
    novel_ids = df["HORclass"].tolist()
    novel_cols = [ str(id) + "(Novel)"  for id in novel_ids]
    novel_stat = {col : {} for col in novel_cols}
    print(novel_stat)

    ##get novel horclass shared assemblies and length ##
    chroms = [i for i in range(1,23)] + ["X", "Y"]

    for ichr in chroms:
        ##HiCAT##
        infile = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/horclass_length/chr{ichr}.HiCAT.horclass.wide.xls"
        df = pd.read_csv(infile, sep="\t", header=0)

        for col in novel_stat.keys():
            if col not in df.columns:
                continue
            tmpdf = df[df[col].notna()]
            sample_length = dict(zip(tmpdf['sample_hap_chrom'], tmpdf[col]))
            if novel_stat[col] == {}:
                novel_stat[col] = sample_length
            else:
                for sample, length in sample_length.items():
                    if sample not in novel_stat[col]:
                        novel_stat[col][sample]  = length 
                    elif novel_stat[col][sample] < length:
                        novel_stat[col][sample] = length
                    else:
                        pass
        ##HORmon##
        infile = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/horclass_length/chr{ichr}.HORmon.horclass.wide.xls"
        df = pd.read_csv(infile, sep="\t", header=0)

        for col in novel_stat.keys():
            if col not in df.columns:
                continue
            tmpdf = df[df[col].notna()]
            sample_length = dict(zip(tmpdf['sample_hap_chrom'], tmpdf[col]))
            if novel_stat[col] == {}:
                novel_stat[col] = sample_length
            else:
                for sample, length in sample_length.items():
                    if sample not in novel_stat[col]:
                        novel_stat[col][sample]  = length 
                    elif novel_stat[col][sample] < length:
                        novel_stat[col][sample] = length
                    else:
                        pass

    data = [(horclass, sample, length) for horclass, samples in novel_stat.items() for sample, length in samples.items()]
    outdf = pd.DataFrame(data, columns=["horclass", "sample_hap_chrom", "length"])

    def add_sample_hap(sample_hap_chrom, ichr):
        if sample_hap_chrom == ichr:
            return "CHM13"
        elif "#" in sample_hap_chrom:
            return "_".join(sample_hap_chrom.split('#')[:2])
        elif "_" in sample_hap_chrom:
            return "_".join(sample_hap_chrom.split('_')[:2])
        return sample_hap_chrom


    outdf['sample_hap'] = outdf['sample_hap_chrom'].apply(lambda x: add_sample_hap(x, ichr))
    outdf.to_csv("novelarray_sample_length.xls", sep="\t", header=True, index=False)
    print(outdf.head())

    result = outdf.groupby(['horclass'])['length'].mean().reset_index()
    print(result.shape)
    print(result.head())

    unique_sample_hap_count = outdf.groupby('horclass')['sample_hap'].nunique().reset_index()
    unique_sample_hap_count.rename(columns={'sample_hap': 'unique_sample_hap_count'}, inplace=True)


    merged_df = pd.merge(result, unique_sample_hap_count, on=['horclass'], how='inner')
    print(merged_df.shape)
    print(merged_df.head())
    merged_df.to_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/horclass_length/novelarray_sample_count_mean_length.xls", sep="\t", header=True, index=False)


if __name__ == "__main__":
    stat_novel_horarray()

    



                        

        


        

