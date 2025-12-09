import pandas as pd
from pybedtools import BedTool
from collections import Counter
import os
import sys


def load_hor_class():
    df = pd.read_csv("../step-11/chrom_hor_phy/all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.color.xls", sep="\t", header=0)
    horclass_color = dict(zip(df['HOR_class'], df['color']))
    hor_class = dict(zip(df['HOR'], df['HOR_class']))

    return horclass_color, hor_class

def load_hor_stv():
    df = pd.read_csv("../step-11/chrom_hor_phy/all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.horcolor.xls", sep="\t", header=0)
    df['n-mer'] = df['HOR'].apply(lambda x : len(x.split('_')))
    df['index'] = range(1, len(df) + 1)
    df['HOR_class'] = df['HOR_class'].astype(str)
    df['index'] = df['index'].astype(str)
    df['n-mer'] = df['n-mer'].astype(str)
    df['horindex'] = "C" + df['HOR_class'] + "H" + df['index'] + "(" + df['n-mer'] + ")"
    print(df.head())

    #df.to_csv("../step-11/chrom_hor_phy/all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.horcolor.ext.xls", sep="\t", index=False, header=True)
    hor_horindex = dict(zip(df['HOR'], df['horindex']))
    hor_horcolor = dict(zip(df['HOR'], df['hor_color']))
    
    return hor_horindex, hor_horcolor


def load_project():
    cols = ['sample', 'pop', 'superpop', 'project']
    df = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/scripts/population.xls", sep="\t", names = cols)
    sample_project = dict(zip(df['sample'], df['project']))
    return sample_project


def load_hor_database(sample, hap, sample_project, outtmp):
    df = pd.read_csv("./statistics/all.HiCAT.hor.summary.final.xls", sep="\t", header=0)
    single_df = df[(df['sample'] == sample) & (df['hap'] == hap) & (df['layer'] == 'top')].copy()
        
    #generate HOR-tmpfile#
    project = sample_project[sample]
    if sample == "CHM13":
        outbed_df = single_df[['chrom', 'horstart', 'horend', 'reorder_hor', 'hor_flag', 'nrepeat', 'nmer']]
    elif (project == "APG") or (project == "REF"):
        #single_df['chr'] = single_df['sample'] + '#' + single_df['hap'] + '#' + single_df['chrom']
        outbed_df = single_df[['chrom', 'horstart', 'horend', 'reorder_hor', 'hor_flag', 'nrepeat', 'nmer']]
    else:
        #single_df['chr'] = single_df['sample'] + '_' + single_df['hap'] + '_' + single_df['chrom']
        outbed_df = single_df[['chrom', 'horstart', 'horend', 'reorder_hor','hor_flag', 'nrepeat', 'nmer']]

    outbed_df.to_csv(outtmp, sep="\t", header=False, index=False)

def merge_hor_mon(outtmp, f_filtermn, hor_class, horclass_color):
    a = BedTool(f_filtermn)
    b = BedTool(outtmp)

    mon_result = a.subtract(b)
    mon_df = mon_result.to_dataframe()
    mon_df['hor_flag'] = 'Mon'
    mon_df['strand_type'] = mon_df['strand']
    mon_df['nrepeat'] = 1
    mon_df['nmer'] = 1

    hor_result = a.intersect(b, wo=True)
    column_names = ['chrom', 'start', 'end', 'mn', 'score', 'strand', 
                    'chr', 'hor_start', 'hor_end', 'name', 
                    'hor_flag', 'nrepeat', 'nmer', 'overlap_len']
    hor_df = hor_result.to_dataframe(names=column_names)
    print(hor_df.head())

    #uniqhor_df = hor_df[['chr', 'hor_start', 'hor_end', 'hor', 'hor_flag']].drop_duplicates()
    #print(uniqhor_df.shape)
    
    def get_main_strand_and_type(strand_list):
        strand_count = Counter(strand_list)
    
        main_strand, count = strand_count.most_common(1)[0]
    
        if len(strand_count) > 1:
            strand_type = 'Mix'
            return main_strand, strand_type
        else:
            strand_type = main_strand
            return main_strand, strand_type
    df_grouped = hor_df.groupby(['chr', 'hor_start', 'hor_end', 'name', 'hor_flag', 'nrepeat', 'nmer'], as_index=False).agg({'strand': list})
    df_grouped[['main_strand', 'strand_type']] = df_grouped['strand'].apply(lambda x: pd.Series(get_main_strand_and_type(x)))
    
    df_grouped['score'] = 0
    uniqhor_df = df_grouped[['chr', 'hor_start', 'hor_end', 'name', 'score', 'main_strand', 'hor_flag', 'strand_type', 'nrepeat', 'nmer']]
    mon_df.columns = uniqhor_df.columns

    merged_df = pd.concat([mon_df, uniqhor_df], ignore_index=True)

    def get_chr_order(ichr):
        chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        if '#' in ichr:
            t = ichr.split('#')[-1]
            if t in chr_order:
                return chr_order.index(t)
            else:
                return 999
        elif '_' in ichr:
            t = ichr.split('_')[-1]
            if t in chr_order:
                return chr_order.index(t)
            else:
                return 999
        else:
            if ichr in chr_order:
                return chr_order.index(ichr)
            else:
                return 999
    merged_df['chr_key'] = merged_df['chr'].apply(lambda x: get_chr_order(x))
    sorted_df = merged_df.sort_values(by=['chr_key', 'hor_start']).drop(columns=['chr_key'])
    sorted_df.rename(columns={'chr': '#chr'}, inplace=True)
    sorted_df['HOR_class'] = sorted_df['name'].map(hor_class)
    sorted_df['color'] = sorted_df['HOR_class'].map(horclass_color)
    sorted_df['HOR_class'] = sorted_df['HOR_class'].fillna(999).astype(int)
    sorted_df['color'] = sorted_df['color'].fillna('NA').astype(str)
    sorted_df['name'] = sorted_df['name'].astype(str)
    print(sorted_df.head())
    return sorted_df


def add_color_class(df, outfile):

    cols = df.columns
    update_subdfs = []

    for chrom, subdf in df.groupby('#chr'):
        data = subdf.values.tolist()

        idx_name = 3
        idx_hor_class = -2
        idx_color = -1

        hor = []
        for i, line in enumerate(data):
            if line[6] == 'HOR':
                hor.append([i, line, set(line[idx_name].split('_'))])
                # print(i, line)
        
        if len(hor) == 0:
            for row in data:
                row[idx_color] = "#969696"
            updated_subdf = pd.DataFrame(data, columns=cols)
            update_subdfs.append(updated_subdf)
            continue

        pivot = 0
        for i, line in enumerate(data):
            if line[6] == 'Mon' or line[6] == 'Dimer':
                if pivot + 1 < len(hor):
                    while pivot + 1 < len(hor) and i > hor[pivot + 1][0]:
                        pivot += 1
                    if pivot + 1 < len(hor):
                        # print(pivot, i, len(hor))
                        if hor[pivot][0] < i < hor[pivot + 1][0]:
                            # print(hor[pivot][1][idx_hor_class], hor[pivot + 1][1][idx_hor_class])
                            if hor[pivot][1][idx_hor_class] == hor[pivot + 1][1][idx_hor_class]:
                                data[i][idx_hor_class] = hor[pivot][1][idx_hor_class]
                                data[i][idx_color] = hor[pivot][1][idx_color]
                                # print(hor[pivot][idx_color])
                                hor[pivot][2].update(data[i][idx_name].split('_'))

        pivot = 0
        for i, line in enumerate(data):
            if line[6] == 'Mon' or line[6] == 'Dimer':
                if pivot + 1 < len(hor):
                    while pivot + 1 < len(hor) and i > hor[pivot + 1][0]:
                        pivot += 1
                    if pivot + 1 < len(hor):
                        # print(pivot, i, len(hor))
                        if hor[pivot][0] < i < hor[pivot + 1][0]:
                            if hor[pivot][1][idx_hor_class] != hor[pivot + 1][1][idx_hor_class]:
                                if data[i][idx_name] in hor[pivot][2]:
                                    data[i][idx_hor_class] = hor[pivot][1][idx_hor_class]
                                    data[i][idx_color] = hor[pivot][1][idx_color]

                                    hor[pivot][2].update(data[i][idx_name].split('_'))

                                elif data[i][idx_name] in hor[pivot + 1][2]:
                                    data[i][idx_hor_class] = hor[pivot + 1][1][idx_hor_class]
                                    data[i][idx_color] = hor[pivot + 1][1][idx_color]

                                    hor[pivot + 1][2].update(data[i][idx_name].split('_'))

                                else:
                                    data[i][idx_hor_class] = '999'
                                    data[i][idx_color] = '#969696'

        for i in range(hor[0][0]):
            if data[i][idx_name] in hor[0][2]:
                data[i][idx_hor_class] = hor[0][1][idx_hor_class]
                data[i][idx_color] = hor[0][1][idx_color]

                hor[0][2].update(data[i][idx_name].split('_'))
            else:
                data[i][idx_hor_class] = '999'
                data[i][idx_color] = '#969696'

        for i in range(hor[-1][0] + 1, len(data)):
            if data[i][idx_name] in hor[-1][2]:
                data[i][idx_hor_class] = hor[-1][1][idx_hor_class]
                data[i][idx_color] = hor[-1][1][idx_color]

                hor[-1][2].update(data[i][idx_name].split('_'))
            else:
                data[i][idx_hor_class] = '999'
                data[i][idx_color] = '#969696'

        updated_subdf = pd.DataFrame(data, columns=cols)
        update_subdfs.append(updated_subdf)

    update_df = pd.concat(update_subdfs, ignore_index=True)
    update_df.to_csv(outfile, sep="\t", header=True, index=False)
    return update_df

def merge(intervals, distance):

    intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        if not merged or merged[-1][-1] + distance < interval[0]:
            ## if not merged or merged[-1][-1] - interval[0] < -1: (不连续区间合并)
            merged.append(interval)
        else:
            merged[-1][-1] = max(merged[-1][-1], interval[-1])
    return merged

def merge_same_horclass(df, distance, out_horclass):

    out_subdfs = []

    for chrom, subdf in df.groupby('#chr'):
        out_data = []
        data = subdf.values.tolist()

        idx_start = 1
        idx_end = 2
        idx_name = 3
        idex_strand = 5
        idx_hor_class = -2
        idx_color = -1

        tmp = []
        for i, line in enumerate(data):
            if i == 0:
                tmp.append(line)
                prev_hor_class = line[idx_hor_class]
                prev_color = line[idx_color]
                prev_strand = line[idex_strand]

            else:
                if (line[idx_hor_class] == prev_hor_class) and (line[idex_strand] == prev_strand):
                    tmp.append(line)
                else:
                    tmp_intervals = [[line[idx_start], line[idx_end]] for line in tmp]
                    tmp_merged_intervals = merge(tmp_intervals, distance)
                    for s, e in tmp_merged_intervals:
                        out_data.append([chrom, s, e, prev_hor_class, 0, prev_strand, prev_color])

                    tmp = [line]
                    prev_hor_class = line[idx_hor_class]
                    prev_color = line[idx_color]
                    prev_strand = line[idex_strand]

        tmp_intervals = [[line[idx_start], line[idx_end]] for line in tmp]
        tmp_merged_intervals = merge(tmp_intervals, distance)
        for s, e in tmp_merged_intervals:
            out_data.append([chrom, s, e, prev_hor_class, 0, prev_strand, prev_color])

        out_subdf = pd.DataFrame(out_data, columns=['#chr', 'start', 'end', 'hor_class', 'score', 'main_strand', 'color'])
        out_subdfs.append(out_subdf)
    out_df = pd.concat(out_subdfs, ignore_index=True)
    out_df.to_csv(out_horclass, sep="\t", header=False, index=False)

def add_stv_index_color(updated_df, hor_horindex, hor_horcolor, outfile):
    updated_df.loc[updated_df['hor_flag'] == 'HOR', 'horindex'] =  updated_df['name'].map(hor_horindex)
    updated_df.loc[updated_df['hor_flag'] == 'HOR', 'hor_color'] = updated_df['name'].map(hor_horcolor)
    
    updated_df['horindex'] = updated_df.apply(
        lambda row: row['hor_flag'] if pd.isna(row['horindex']) else row['horindex'], 
        axis=1
    )
    updated_df['hor_color'] = updated_df.apply(
        lambda row: (
            '#969696' if row['hor_flag'] == 'Dimer' else
            '#C3C3C3' if row['hor_flag'] == 'Mon' else
            row['hor_color']
        ),
        axis=1
    ) 

    print(updated_df.head())

    updated_df.to_csv(outfile, sep="\t", header=True, index=False)
    return updated_df



def main(hordir, sample, hap):
    sample_project = load_project()
    horclass_color, hor_class = load_hor_class()
    hor_horindex, hor_horcolor = load_hor_stv()

    outtmp = os.path.join(hordir, f"{sample}_{hap}.reorder.HOR.bed")
    f_filtermn = os.path.join(hordir, "filter.asat.bed")
    load_hor_database(sample, hap, sample_project, outtmp)
    sorted_df = merge_hor_mon(outtmp, f_filtermn, hor_class, horclass_color)

    outfile = os.path.join(hordir, f"{sample}_{hap}.HORmn.bed")
    outstv = os.path.join(hordir, f"{sample}_{hap}.HORmn.v2.bed")
    update_df = add_color_class(sorted_df, outfile)

    out_horclass = os.path.join(hordir, f"{sample}_{hap}.HORclass.bed")
    merge_same_horclass(update_df, 10000, out_horclass)

    update_stv_df = add_stv_index_color(update_df, hor_horindex, hor_horcolor, outstv) 
    out_merged_stv = os.path.join(hordir, f"{sample}_{hap}.HORstv.bed")
    merge_same_horclass(update_stv_df, 10000, out_merged_stv)

if __name__ == "__main__":
    hordir = sys.argv[1]
    sample = sys.argv[2]
    hap = sys.argv[3]
    main(hordir, sample, hap)
    #load_hor_stv()
