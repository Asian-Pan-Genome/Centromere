import pandas as pd

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


df = pd.read_csv("all.graph.hordecomposition.summary.final.HOR.count.xls", sep="\t", names = ['reorder_hor', 'max_count'])
print(df.head())
print(df.shape)

sample_hap_df = pd.read_csv("all.graph.hordecomposition.summary.candidate.HOR.count.xls", sep="\t", names = ['sample', 'hap', 'chrom', 'reorder_hor', 'count'])

sample_hap_df = sample_hap_df[
    (sample_hap_df['chrom'].str.contains('chr')) & 
    (~sample_hap_df['chrom'].str.contains('chrY-CEN'))
].copy()


hordf = sample_hap_df.merge(df, on='reorder_hor', how='inner')
print(hordf.shape)
print(hordf.head())

mainchr_hordf = hordf.groupby(['reorder_hor','chrom']).size().reset_index(name='count')
print(mainchr_hordf.head())
print(mainchr_hordf.shape)

max_mainchr_df = mainchr_hordf.loc[mainchr_hordf.groupby('reorder_hor')['count'].idxmax()]
max_mainchr_df.columns = ['reorder_hor', 'main_chromosome', 'count_main_chromosome_shared_by_haploid']
print(max_mainchr_df.head())
print(max_mainchr_df.shape)

finaldf = hordf.merge(max_mainchr_df, how='left', on=['reorder_hor'])
finaldf['sample_hap'] = finaldf['sample'] + '_' + finaldf['hap']    
print(finaldf.head())
print(finaldf.shape)

haploid_shared_df = finaldf.groupby('reorder_hor')['sample_hap'].nunique().reset_index(name='sample_hap_count')
merged_df = finaldf.merge(haploid_shared_df, how='left', on=['reorder_hor'])
merged_df['pan_hor_type'] = merged_df.apply(lambda row: assign_category(row), axis=1)  
print(merged_df.head())
print(merged_df.shape)

merged_df.to_csv("all.graph.hordecomposition.summary.final.pan_hor_dr_type.xls", sep="\t", header=True, index=False)
