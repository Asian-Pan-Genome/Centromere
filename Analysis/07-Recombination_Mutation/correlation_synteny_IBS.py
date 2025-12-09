import pandas as  pd
from scipy.stats import pearsonr, linregress
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def change_id(x):
    if '-Mat' in x and '_2' not in x:
        sample_hap = x.split('_')[0]
    elif '-Pat' in x and '_2' not in x:
        sample_hap = x.split('_')[0]
    else:
        sample_hap = x
    return sample_hap


def load_ibs_data(ibs_file):
    df = pd.read_csv(ibs_file, sep='\s+', header=0)
    print(df.head())
    print(df.columns)

    df['1_IBS'] = 1 - df['DST']
    
    ibs_df = df[['IID1', 'IID2', '1_IBS']].copy()
    print(ibs_df.head())

    ibs_df.loc[ :, 'IID1_sample_hap'] = ibs_df['IID1'].apply(lambda x: change_id(x))
    ibs_df.loc[ :, 'IID2_sample_hap'] = ibs_df['IID2'].apply(lambda x: change_id(x))

    return ibs_df


def load_synteny_distance(distance_file):
    distance_df = pd.read_csv(distance_file, sep=',', header=0, index_col=0)
    print(distance_df.iloc[:10, :10])

    distance_df.reset_index(inplace=True)
    long_distance_df = distance_df.melt(id_vars=['source'], var_name='IID2_sample_hap', value_name='distance')
    long_distance_df.rename(columns={'source': 'IID1_sample_hap'}, inplace=True)

    print(long_distance_df.head())
    return long_distance_df


def normalize_data(df, columns):
    scaler = MinMaxScaler()
    df[columns] = scaler.fit_transform(df[columns])
    return df

def get_correlation(left_ibs_file, right_ibs_file, distance_file, out_png):
    left_ibs_df = load_ibs_data(left_ibs_file)
    left_ibs_df.rename(columns={"1_IBS": "left_IBS"}, inplace=True)
    right_ibs_df = load_ibs_data(right_ibs_file)
    right_ibs_df.rename(columns={"1_IBS": "right_IBS"}, inplace=True)

    merged_ibs_df = left_ibs_df.merge(right_ibs_df, left_on=['IID1_sample_hap', 'IID2_sample_hap'], right_on=['IID1_sample_hap', 'IID2_sample_hap'], how='inner')
    print(merged_ibs_df.shape)

    distance_df = load_synteny_distance(distance_file)
    final_df = merged_ibs_df.merge(distance_df, left_on=['IID1_sample_hap', 'IID2_sample_hap'], right_on=['IID1_sample_hap', 'IID2_sample_hap'], how='inner')
    final_df = final_df.dropna(subset=['left_IBS', 'right_IBS', 'distance'])
    print(len(final_df[['IID1_sample_hap', 'IID2_sample_hap']].drop_duplicates()))
    print(final_df.head())
    print(final_df.shape)

    # final_df = normalize_data(final_df, ['left_IBS', 'right_IBS', 'distance'])
    out = []
    
    # left vs. right correlation
    correlation_left_right, p_left_right = pearsonr(final_df['left_IBS'], final_df['right_IBS'])
    print(f'Pearson correlation and p-value between left_IBS and right_IBS: {correlation_left_right} {p_left_right}')
    out.append([correlation_left_right, p_left_right])

    # left vs. synteny correlation
    correlation_left_distance, p_left_distance = pearsonr(final_df['left_IBS'], final_df['distance'])
    print(f'Pearson correlation and p-value between left_IBS and distance: {correlation_left_distance} {p_left_distance}')
    out.append([correlation_left_distance, p_left_distance])

    # right vs. synteny correlation
    correlation_right_distance, p_right_distance = pearsonr(final_df['right_IBS'], final_df['distance'])
    print(f'Pearson correlation and p-value between right_IBS and distance: {correlation_right_distance} {p_right_distance}')
    out.append([correlation_right_distance, p_right_distance])

    #plot dotplot#
    fig, axs = plt.subplots(1, 3, figsize=(10, 3))
    sns.regplot(x='left_IBS', y='right_IBS', data=final_df, ax=axs[0], scatter_kws={'alpha':0.2})
    axs[0].set_title('left_IBS vs right_IBS')
    axs[0].text(0.95, 0.05, f'R={correlation_left_right:.2f}, p={p_left_right:.2e}', 
            ha='right', va='bottom', transform=axs[0].transAxes)
    
    sns.regplot(x='left_IBS', y='distance', data=final_df, ax=axs[1], scatter_kws={'alpha':0.2})
    axs[1].set_title('left_IBS vs distance')
    axs[1].text(0.95, 0.05, f'R={correlation_left_distance:.2f}, p={p_left_distance:.2e}', 
            ha='right', va='bottom', transform=axs[1].transAxes)

    sns.regplot(x='right_IBS', y='distance', data=final_df, ax=axs[2], scatter_kws={'alpha':0.2})
    axs[2].set_title('right_IBS vs distance')
    axs[2].text(0.95, 0.05, f'R={correlation_right_distance:.2f}, p={p_right_distance:.2e}', 
            ha='right', va='bottom', transform=axs[2].transAxes)
    
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    return out

def convert_to_matrix(outfile, cen_index):
    df = pd.read_csv(outfile, sep='\t', header=None)
    df.columns = ['rindex1', 'rindex2', 'r', 'p']
    print(df.head())

    all_indexs = set(list(df['rindex1']) + list(df['rindex2']))
    reindex = {index : i for i, index in enumerate(sorted(list(all_indexs)))}
    new_old_index = {i: ('CEN' if index == cen_index else index)  for index, i in reindex.items()}

    print(reindex)
    print(new_old_index)

    r_matrix = pd.DataFrame(0.0, index=range(len(reindex)), columns=range(len(reindex)))
    for _, row in df.iterrows():
        i = reindex[row['rindex1']]
        j = reindex[row['rindex2']]
        r_matrix.at[i, j] = row['r']
        r_matrix.at[j, i] =row['r']
    r_matrix.rename(index=new_old_index, columns=new_old_index, inplace=True)

    p_matrix = pd.DataFrame(0.0, index=range(len(reindex)), columns=range(len(reindex)))

    for _, row in df.iterrows():
        i = reindex[row['rindex1']]
        j = reindex[row['rindex2']]
        p_matrix.at[i, j] = row['p']
        p_matrix.at[j, i] =row['p']
    p_matrix.rename(index=new_old_index, columns=new_old_index, inplace=True)

    return r_matrix, p_matrix

def main(f_filtered_region, distance_file, cen_index, outfile, out_r, out_p):
    region_df = pd.read_csv(f_filtered_region, sep='\t', names = ['chrom', 'start', 'end', 'index', 'snp_num', 'sample_num'])
    indexs = list(region_df['index'])

    outf = open(outfile, 'w')

    for i in range(len(indexs)):
        for j in range(i+1, len(indexs)):
            print(indexs[i], indexs[j], " start compare")
            left_ibs_file = f"tmp/{indexs[i]}.matrix.genome"
            right_ibs_file = f"tmp/{indexs[j]}.matrix.genome"
            out_png = f"tmp/{indexs[i]}_vs_{indexs[j]}.dotplot.png"

            print(f"Compare {left_ibs_file} with {right_ibs_file}")

            out = get_correlation(left_ibs_file, right_ibs_file, distance_file, out_png)
            outf.write(f"{indexs[i]}\t{indexs[j]}\t{out[0][0]}\t{out[0][1]}\n")
        outf.write(f"{indexs[i]}\t{cen_index}\t{out[1][0]}\t{out[1][1]}\n")
    outf.close()
    
    r_matrix, p_matrix = convert_to_matrix(outfile, cen_index)
    r_matrix.to_csv(out_r, sep='\t', index=True, header=True)
    p_matrix.to_csv(out_p, sep='\t', index=True, header=True)

if __name__ == '__main__':
    f_filtered_region = sys.argv[1]
    distance_file = sys.argv[2]
    outfile = sys.argv[3]
    cen_index = float(sys.argv[4])
    out_r = sys.argv[5]
    out_p = sys.argv[6]
    main(f_filtered_region, distance_file, cen_index, outfile, out_r, out_p)
