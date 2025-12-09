import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

def get_trans(ichr):
    ##load cent_length file## 
    df_cen = pd.read_csv("cent_chrom.txt", sep="\t", header= 0)
    complete_df_cen = df_cen[df_cen['filterflag'] == 0][['sample_hap_chrom', 'start', 'end', 'len', 'project', 'filterflag']].copy()
    complete_df_cen.columns = ['chrom', 'cen_start', 'cen_end', 'cen_len', 'project', 'filterflag']

    ##load_cenanno##
    cenanno = pd.read_csv(ichr + ".merged.cenanno.bed", sep="\t", header=None)
    cenanno.columns = ["chrom", "start", "end", "satellite", "score", "strand", "s", "e", "color"]
    cenanno['color'] = cenanno['color'].apply(lambda x: "#891640" if x == "" else x)

    ##plot_cenanno_data_clean##
    plot_cenanno_df = pd.merge(cenanno, complete_df_cen, on="chrom", how='inner')
    plot_cenanno_df = plot_cenanno_df[( plot_cenanno_df['start'] >= plot_cenanno_df["cen_start"] ) & \
                                       (plot_cenanno_df['end'] <= plot_cenanno_df['cen_end'])].copy()
    plot_cenanno_df['modstart'] = plot_cenanno_df['start'] - plot_cenanno_df['cen_start'] + 1
    plot_cenanno_df['modend'] = plot_cenanno_df['end'] - plot_cenanno_df['cen_start'] + 1
    # print(plot_cenanno_df.head())

    ##cenanno summary##
    cenanno_summary = plot_cenanno_df[(plot_cenanno_df['satellite'] != 'rDNA') & \
                                       (plot_cenanno_df['end'] - plot_cenanno_df['start'] >= 50000)].copy()
    ##only APG samples##
    cenanno_summary = cenanno_summary[cenanno_summary['project'] == "APG"].copy()
    # print(cenanno_summary[cenanno_summary['chrom'].str.contains("C050-CHA-N10#Mat#chr1")])

    # print(cenanno_summary.head())
    transitions = []
    for chrom, subdf in cenanno_summary.groupby('chrom'):
        prev_satellite = None
        prev_end = None
        tmp = {'chrom' : chrom}
        for index, row in subdf.iterrows():            
            if prev_satellite is None:
                prev_satellite = row['satellite']
                prev_end = row['end']
            else:
                if (row['satellite'] != prev_satellite) and (row['start'] - prev_end <= 500000):
                    trans_type = f"{prev_satellite}_{row['satellite']}"
                    tmp[trans_type] = tmp.get(trans_type, 0) + 1
                    prev_satellite = row['satellite']
                    prev_end = row['end']
                else:
                    # if row['start'] - prev_end >= 100000:
                    #     prev_satellite = row['satellite']
                    #     prev_end = row['end']
                    prev_satellite = row['satellite']
                    prev_end = row['end']
                    # pass
            # print(tmp)
        # if chrom == "C050-CHA-N10#Mat#chr1":
        #     print(tmp)
        transitions.append(tmp)
    # print(transitions)
    transdf = pd.DataFrame(transitions)
    transdf.fillna(0, inplace=True)
    transdf['chr'] = ichr

    return transdf

def main():
    chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    dfs = []
    for ichr in chrs:
        transdf = get_trans(ichr)
        dfs.append(transdf)
    final_df = pd.concat(dfs, join='outer').fillna(0)
    tmp = final_df.copy()
    # print(final_df.head())
    # final_df.to_csv("satellite_transition.xls", sep="\t", index=True, header=True)


    final_df.drop(columns=['chrom'], inplace=True)
    grouped_df = final_df.groupby(list(final_df.columns)).size().reset_index(name='counts')
    grouped_df['chr'] = pd.Categorical(grouped_df['chr'], categories=chrs, ordered=True)
    grouped_df_desc = grouped_df.sort_values(by=['chr', 'counts'], ascending=[True, False])
    grouped_df_desc['index'] = range(1, len(grouped_df) + 1)
    grouped_df_desc['type'] = grouped_df_desc.apply(lambda row: f"{row['chr']}_{row['index']}_{row['counts']}", axis=1)
    merge_columns = [col for col in grouped_df_desc.columns if col not in ['counts', 'type', 'index']]
    merged_df = pd.merge(grouped_df_desc, tmp, on=merge_columns)
    merged_df.to_csv("satellite_transition.xls", sep="\t", index=True, header=True)
    print(merged_df.head())

    t = grouped_df.groupby('chr').size().reset_index(name='count_organization')
    
    print(t)
    # final_df.set_index('chr', inplace=True)
    
    # final_df = final_df.loc[:, ~final_df.columns.duplicated()]
    # final_df = final_df[~final_df.index.duplicated()]


    # sns.clustermap(final_df, cmap='viridis', row_cluster=True, col_cluster=True, annot=True, fmt='.2f')
    # plt.title('Heatmap with Clustering')
    # plt.savefig("translation_heatmap.png")
    # plt.show()
    # pca = PCA(n_components=2)
    # pca_result = pca.fit_transform(final_df)

    # pca_df = pd.DataFrame(data=pca_result, columns=['PCA1', 'PCA2'])
    # pca_df['chr'] = final_df.index

    # # Plot PCA results
    # plt.figure(figsize=(10, 8))
    # scatter = plt.scatter(pca_df['PCA1'], pca_df['PCA2'], c=pca_df['chr'].astype('category').cat.codes, cmap='viridis')
    # plt.colorbar(scatter, ticks=range(len(pca_df['chr'].unique())), label='Chromosome')
    # plt.xlabel('PCA1')
    # plt.ylabel('PCA2')
    # plt.title('PCA of Satellite Data')
    # plt.savefig("translation.png")
    # # plt.show()



if __name__ == "__main__":
    main()
