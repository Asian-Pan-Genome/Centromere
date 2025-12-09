import pandas as pd

def get_correlation_summary(ichr, \
                            cen_index, \
                            left_linked_start_index, \
                            left_linked_end_index, \
                            right_linked_start_index, \
                            right_linked_end_index, \
                            lr_start_index, \
                            lr_end_index, \
                            left_cent_start_index, \
                            left_cent_end_index, \
                            right_cent_start_index, \
                            right_cent_end_index):

    cols = ['chrom', 'start', 'end', 'index', 'snp_count']
    df = pd.read_csv(f"{ichr}/100kb_windows_with_snp_counts_filtered.bed", sep="\t", names = cols)
    correlation_df = pd.read_csv(f"{ichr}/100kb_windows_filtered_pearson_p.xls", sep="\t", names=['index1', 'index2', 'r', 'p'])

    def get_unilateral_linked_info(start_index, end_index):
        start = df.loc[df['index'] == start_index, 'start'].values[0]
        end = df.loc[df['index'] == end_index, 'end'].values[0]
        length = end - start
        target_indexs = list(df[(df['index'] >= start_index) & (df['index'] <= end_index)]['index'])
        linked_r = []
        for i in sorted(target_indexs):
            for j in sorted(target_indexs):
                if i < j:
                    r = correlation_df.loc[(correlation_df['index1'] == i) & (correlation_df['index2'] == j), 'r']
                    linked_r.append(0.0 if r.empty or r.values[0] < 0 else r.values[0])
        avg_r = sum(linked_r) / len(linked_r) if linked_r else 0.0
        return start, end, length, avg_r


    ##get left/right linked length and avg correlation##
    left_start, \
    left_end, \
    left_length, \
    left_avg_r = get_unilateral_linked_info(left_linked_start_index, left_linked_end_index)

    right_start, \
    right_end, \
    right_length, \
    right_avg_r = get_unilateral_linked_info(right_linked_start_index, right_linked_end_index)


    ##get left/right centromere length and avg correlation##
    def get_left_right_linked_info(lr_start_index, left_linked_end_index, right_linked_start_index, lr_end_index):

        left_start = df.loc[df['index'] == lr_start_index, 'start'].values[0]
        left_end = df.loc[df['index'] == left_linked_end_index, 'end'].values[0]
        left_length = left_end - left_start

        right_start = df.loc[df['index'] == right_linked_start_index, 'start'].values[0]
        right_end = df.loc[df['index'] == lr_end_index, 'end'].values[0]
        right_length = right_end - right_start

        left_indexs = list(df[(df['index'] >= lr_start_index) & (df['index'] <= left_linked_end_index)]['index'])
        right_indexs = list(df[(df['index'] >= right_linked_start_index) & (df['index'] <= lr_end_index)]['index'])
        linked_r = []
        for i in sorted(left_indexs):
            for j in sorted(right_indexs):
                r = correlation_df.loc[(correlation_df['index1'] == i) & (correlation_df['index2'] == j), 'r']
                linked_r.append(0.0 if r.empty or r.values[0] < 0 else r.values[0])
        avg_r = sum(linked_r) / len(linked_r) if linked_r else 0.0
        return left_start, left_end, right_start, right_end, left_length, right_length, avg_r
    lr_left_start, \
    lr_left_end, \
    lr_right_start, \
    lr_right_end, \
    lr_left_length, \
    lr_right_length, \
    lr_avg_r = get_left_right_linked_info(lr_start_index, left_linked_end_index, right_linked_start_index, lr_end_index)

    ##get left/right with cent length and avg correlation##
    def get_cent_unilateral_linked_info(cen_index, cent_start_index, cent_end_index):

        start = df.loc[df['index'] == cent_start_index, 'start'].values[0]
        end = df.loc[df['index'] == cent_end_index, 'end'].values[0]
        length = end - start

        target_indexs = list(df[(df['index'] >= cent_start_index) & (df['index'] <= cent_end_index)]['index'])
        # print(target_indexs)
        linked_r = []
        for i in sorted(target_indexs):
            r = correlation_df.loc[(correlation_df['index1'] == i) & (correlation_df['index2'] == cen_index), 'r']
            linked_r.append(0.0 if r.empty or r.values[0] < 0 else r.values[0])
        avg_r = sum(linked_r) / len(linked_r) if linked_r else 0.0
        print(linked_r)
        return start, end, length, avg_r

    left_cent_start, \
    left_cent_end, \
    left_cent_length, \
    left_cent_avg_r = get_cent_unilateral_linked_info(cen_index, left_cent_start_index, left_cent_end_index)

    right_cent_start, \
    right_cent_end, \
    right_cent_length, \
    right_cent_avg_r = get_cent_unilateral_linked_info(cen_index, right_cent_start_index, right_cent_end_index)
    print("-------------------------------")

    data = {
        "chrom" : [ichr] * 6,
        "start" : [left_start, right_start, lr_left_start, lr_right_start, left_cent_start, right_cent_start],
        "end" : [left_end, right_end, lr_left_end, lr_right_end, left_cent_end, right_cent_end],
        "length" : [left_length, right_length, lr_left_length, lr_right_length, left_cent_length, right_cent_length],
        "avg_r" : [left_avg_r, right_avg_r, lr_avg_r, lr_avg_r, left_cent_avg_r, right_cent_avg_r],
        "type" : ["left_linked", "right_linked", "lr_left_linked", "lr_right_linked", "left_cent", "right_cent"],
    }
    summary_df = pd.DataFrame(data)
    return summary_df

def get_lr_r(ichr, lr_start_index, left_linked_end_index, right_linked_start_index, lr_end_index):
    cols = ['chrom', 'start', 'end', 'index', 'snp_count']
    df = pd.read_csv(f"{ichr}/100kb_windows_with_snp_counts_filtered.bed", sep="\t", names = cols)
    correlation_df = pd.read_csv(f"{ichr}/100kb_windows_filtered_pearson_p.xls", sep="\t", names=['index1', 'index2', 'r', 'p'])

    left_indexs = list(df[(df['index'] >= lr_start_index) & (df['index'] <= left_linked_end_index)]['index'])
    right_indexs = list(df[(df['index'] >= right_linked_start_index) & (df['index'] <= lr_end_index)]['index'])

    out_correlation_df = correlation_df.loc[(correlation_df['index1'].isin(left_indexs)) & (correlation_df['index2'].isin(right_indexs))].copy()
    out_correlation_df['chrom'] = ichr
    return out_correlation_df



def main():
    input_index_df = pd.read_csv("chroms_linked_index_update.xls", sep="\t", header=0)
    dfs = []
    lr_dfs = []

    for index, row in input_index_df.iterrows():
        ichr = row['chrom']
        cen_index = row['cen_index']
        left_linked_start_index = row['left_linked_start_index']
        left_linked_end_index = row['left_linked_end_index']
        right_linked_start_index = row['right_linked_start_index']
        right_linked_end_index = row['right_linked_end_index']
        lr_start_index = row['lr_start_index']
        lr_end_index = row['lr_end_index']
        left_cent_start_index = row['left_cent_start_index']
        left_cent_end_index = row['left_cent_end_index']
        right_cent_start_index = row['right_cent_start_index']
        right_cent_end_index = row['right_cent_end_index']

        summary_df = get_correlation_summary(ichr, \
                                             cen_index, \
                                             left_linked_start_index, \
                                             left_linked_end_index, \
                                             right_linked_start_index, \
                                             right_linked_end_index, \
                                             lr_start_index, \
                                             lr_end_index, \
                                             left_cent_start_index, \
                                             left_cent_end_index, \
                                             right_cent_start_index, \
                                             right_cent_end_index)

        correlation_subdf = get_lr_r(ichr, \
                                     lr_start_index, \
                                     left_linked_end_index, \
                                     right_linked_start_index, \
                                     lr_end_index)
        
        dfs.append(summary_df)
        lr_dfs.append(correlation_subdf)
    final_df = pd.concat(dfs, ignore_index=True)
    final_df.to_csv("chroms_linked_correlation_summary_update.xls", sep="\t", index=False)

    final_lr_df = pd.concat(lr_dfs, ignore_index=True)
    final_lr_df.to_csv("chroms_lr_correlation_r_update.xls", sep="\t", index=False)
    

if __name__ == "__main__":
    main()