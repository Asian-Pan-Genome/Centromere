import pandas as pd
import numpy as np

def unified_chromosome(row):
    if '#' in row['chrom']:
        ichr = row['chrom'].split('#')[-1]
    elif '_' in row['chrom']:
        ichr = row['chrom'].split('_')[-1]
    else:
        ichr = row['chrom']
    return ichr

def load_merge_config():
    replace_dict = {}
    with open('../../cluster_need_merge.xls','r') as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            replace_dict[tokens[1]] = tokens[0]
    return replace_dict

#replaced need-merged-cluster#
def replace_in_hor(ihor):
    replace_dict = load_merge_config()
    parts = ihor.split('_')
    parts = [replace_dict.get(part, part) for part in parts]
    return '_'.join(parts)

def load_original_HiCAT_summary(hicat_sum_file = "all.HiCAT.hor.summary.xls"):
    cols = ['sample', 'hap', 'blockindex',\
            'start', 'end', 'blocklen',\
            'mnnum', 'chrom', 'horstart', \
            'horend', 'hor_index_start', 'hor_index_end',
            'nrepeat', 'hor', 'layer', 'reorder_hor']
    df = pd.read_csv(hicat_sum_file, sep = "\t", names = cols)
    df['sample_hap_chrom'] = df['sample'] + '_' + df['hap'] + '_' + df['chrom']
    df['sample_hap'] = df['sample'] + '_' + df['hap']
    df['chromosome'] = df.apply(lambda row: unified_chromosome(row), axis=1)


    filterdf = df[(df['hor'] != "Na") &
                (df['chrom'].str.contains('chr')) &
                (~df['chrom'].str.contains('chrY-CEN'))].copy()
    filterdf['nrepeat'] = filterdf['nrepeat'].astype(int)
    filterdf['nmer'] = filterdf['reorder_hor'].apply(lambda x: len(x.split('_')))
    #print(filterdf['chromosome'].unique())
    filterdf.iloc[:, 13] = filterdf.iloc[:, 13].apply(replace_in_hor)
    filterdf.iloc[:, 15] = filterdf.iloc[:, 15].apply(replace_in_hor)
    #print(filterdf.head())
    #print(filterdf.iloc[8:15,8:15])

    ##define HOR and DR##
    groupdf = filterdf.groupby('reorder_hor')['nrepeat'].max().reset_index()
    #print(groupdf.head())
    groupdf.columns = ['reorder_hor', 'max_nrepeat']
    uniqDR = groupdf[groupdf['max_nrepeat'] < 3]
    dimer_repeat = list(uniqDR['reorder_hor'])
    uniqHOR = groupdf[groupdf['max_nrepeat'] >= 3]
    uniqHOR.to_csv("all.HiCAT.hor.summary.final.HOR_maxrepeat.xls", sep="\t", header=True, index=False)
    uniqDR.to_csv("all.HiCAT.hor.summary.final.DR_maxrepeat.xls", sep="\t", header=True, index=False)
    print(groupdf.shape)
    print(len(dimer_repeat))
    
    #filterdf['hor_flag'] = np.where(filterdf['reorder_hor'].isin(dimer_repeat), 'Dimer', 'HOR')
    filterdf['hor_flag'] = filterdf['reorder_hor'].apply(lambda x : 'Dimer' if x in dimer_repeat else 'HOR')
    print(filterdf.columns)
    print(filterdf.iloc[8:,8:])
    

    ##assign pan_HOR_type##
    dedup_df = filterdf[['reorder_hor', 'sample_hap', 'chromosome', 'hor_flag']].drop_duplicates()
    mainchr_df = dedup_df.groupby(['reorder_hor','chromosome']).size().reset_index(name='count')
    max_mainchr_df = mainchr_df.loc[mainchr_df.groupby('reorder_hor')['count'].idxmax()]
    max_mainchr_df.columns = ['reorder_hor', 'main_chromosome', 'count_main_chromosome_shared_by_haploid']
    mainchr_dedup_df = dedup_df.merge(max_mainchr_df, how='left', on=['reorder_hor'])
 
    haploid_shared_df = dedup_df.groupby('reorder_hor')['sample_hap'].nunique().reset_index(name='sample_hap_count')
    merged_df = mainchr_dedup_df.merge(haploid_shared_df, how='left', on=['reorder_hor'])
    merged_df['pan_hor_type'] = merged_df.apply(lambda row: assign_category(row), axis=1)   
    merged_df.to_csv("all.HiCAT.hor.summary.final.pan_hor_dr_type.xls", sep="\t", header=True, index=False)
     
    pan_hor_type_dict = dict(zip(merged_df['reorder_hor'], merged_df['pan_hor_type']))
    filterdf['pan_hor_type'] = df['reorder_hor'].map(pan_hor_type_dict)
    #finaldf = filterdf.merge(merged_df[['reorder_hor', 'pan_hor_type']], on='reorder_hor', how='left') 
    filterdf.to_csv("all.HiCAT.hor.summary.final.xls", sep="\t", header=True, index=False)

def assign_category(row):
    #chrX:312 assemblies; chrY:211 assemblies
    is_chrX = row['main_chromosome'] == 'chrX'
    is_chrY = row['main_chromosome'] == 'chrY'

    if row['sample_hap_count'] == 1:
        return 'Singleton'

    elif is_chrX:
        if 1 < row['sample_hap_count'] <= 3:
            return 'Rare'
        elif 3 < row['sample_hap_count'] <= 15:
            return 'Low frequency'
        elif 15 < row['sample_hap_count'] <= 280:
            return 'Dispensable'
        else:
            return 'Core'

    elif is_chrY:
        if 1 < row['sample_hap_count'] <= 2:
            return 'Rare'
        elif 2 < row['sample_hap_count'] <= 10:
            return 'Low frequency'
        elif 10 < row['sample_hap_count'] <= 190:
            return 'Dispensable'
        else:
            return 'Core'
 
    else:
        if 1 <  row['sample_hap_count'] <= 5 :
            return 'Rare'
        elif 5 < row['sample_hap_count'] <= 27 :
            return 'Low frequency'
        elif 27 < row['sample_hap_count'] <= 488:
            return 'Dispensable'
        else:
            return 'Core'
    return None    

if __name__ == "__main__":
    load_original_HiCAT_summary("all.HiCAT.hor.summary.dedup.xls")

#column_order = ['sample','hap','blockindex','chrom','start','end','blocklen','mnnum',
#                'horid','hor','n-mer','count','sample_hap_chrom','sample_hap','chromosome',
#                'sample_hap_count','sample_hap_chrom_count','total_count','max_count','pan_hor_type']
#out = out[column_order]
 
#out.to_csv("all.HiCAT.hors.final.updated.xls", sep="\t", header=True, index=False)

