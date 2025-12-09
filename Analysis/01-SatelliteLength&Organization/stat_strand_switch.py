import pandas as pd
import os
import sys

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


def get_chr_order(ichr):
    chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    if '#' in ichr:
        t = ichr.split('#')[-1]
        if t in chr_order:
            return chr_order.index(t)
    elif '_' in ichr:
        t = ichr.split('_')[-1]
        if t in chr_order:
            return chr_order.index(t)
    elif ichr in chr_order:
        return chr_order.index(ichr)
    else:
        return 999


# def get_asat_inversion_breakpoint(cenanno, outfile,  outbreak, distance=10000,):
#     df = pd.read_csv(cenanno, sep="\t", names = ["chrom", "start", "end", "sat", "score", "strand", "hs", "he", "color"])

#     asatdf = df[df['sat'].str.startswith('S')].copy()

     
#     out_subdfs = []
#     out_breaks = []

#     for chrom, subdf in asatdf.groupby('chrom'):

#         tmp = []
#         out = []
#         breaks = [] 

#         for idx, (i, row) in enumerate(subdf.iterrows()):
#             if idx == 0:
#                 prev_strand = row['strand']
#                 prev_start = int(row['start'])
#                 prev_end = int(row['end'])
#                 tmp.append([prev_start, prev_end])
#             else:
#                 if row['strand'] == prev_strand:
#                     tmp.append([int(row['start']), int(row['end'])])
#                 else:
#                     merged = merge(tmp, distance)
#                     # print(chrom, i, merged)
#                     for s, e in merged:
#                         out.append([chrom, s, e, prev_strand])

#                     breaks.append([chrom, tmp[-1][-1], row['start'], row['start'] - merged[-1][-1], prev_strand+"/"+row['strand']])

#                     prev_strand = row['strand']
#                     prev_start = int(row['start'])
#                     prev_end = int(row['end'])
#                     tmp = [[prev_start, prev_end]]  

#         ##last##
#         merged = merge(tmp, distance)
#         for s, e in merged:
#             out.append([chrom, s, e, prev_strand])
    
#         out_subdf = pd.DataFrame(out, columns = ["chrom", "start", "end", "strand"])
#         out_subdfs.append(out_subdf)

#         out_break = pd.DataFrame(breaks, columns = ["chrom", "start", "end", "break_distance", "strand"])
#         out_breaks.append(out_break)

#     out_df = pd.concat(out_subdfs, ignore_index=True)
#     out_df['chr_key'] = out_df['chrom'].apply(lambda x: get_chr_order(x))
#     sorted_out_df = out_df.sort_values(by=['chr_key', 'start']).drop(columns=['chr_key'])
#     sorted_out_df.to_csv(outfile, sep="\t", header=False, index=False)

#     out_b = pd.concat(out_breaks, ignore_index=True)
#     out_b['chr_key'] = out_b['chrom'].apply(lambda x: get_chr_order(x))
#     sorted_out_b = out_b.sort_values(by=['chr_key', 'start']).drop(columns=['chr_key'])
#     sorted_out_b.to_csv(outbreak, sep="\t", header=False, index=False)


def get_sat_inversion_breakpoint(cenanno, outfile,  outbreak, distance=10000):
    df = pd.read_csv(cenanno, sep="\t", names = ["chrom", "start", "end", "sat", "score", "strand", "hs", "he", "color"])

    df["sat"] = df["sat"].apply(lambda x: "ASat" if x.startswith('S') else x)
    
    out_subdfs = []
    out_breaks = []

    for chrom, subdf in df.groupby('chrom'):

        tmp = []
        out = []
        breaks = [] 

        for idx, (i, row) in enumerate(subdf.iterrows()):
            if idx == 0:
                prev_strand = row['strand']
                prev_start = int(row['start'])
                prev_end = int(row['end'])
                prev_sat = row['sat']
                tmp.append([prev_start, prev_end])
            else:
                if (row['strand'] == prev_strand ) and (row['sat'] == prev_sat):
                    tmp.append([int(row['start']), int(row['end'])])
                else:
                    merged = merge(tmp, distance)
                    # print(chrom, i, merged)
                    for s, e in merged:
                        out.append([chrom, s, e, prev_strand, prev_sat])

                    breaks.append([chrom, tmp[-1][-1], row['start'], row['start'] - merged[-1][-1], prev_strand+"/"+row['strand'], prev_sat+"/"+row['sat']])

                    prev_strand = row['strand']
                    prev_start = int(row['start'])
                    prev_end = int(row['end'])
                    prev_sat = row['sat']
                    tmp = [[prev_start, prev_end]]  

        ##last##
        merged = merge(tmp, distance)
        for s, e in merged:
            out.append([chrom, s, e, prev_strand, prev_sat])
    
        out_subdf = pd.DataFrame(out, columns = ["chrom", "start", "end", "strand", "satellite"])
        out_subdfs.append(out_subdf)

        out_break = pd.DataFrame(breaks, columns = ["chrom", "start", "end", "break_distance", "strand", "satellite"])
        out_breaks.append(out_break)

    out_df = pd.concat(out_subdfs, ignore_index=True)
    out_df['chr_key'] = out_df['chrom'].apply(lambda x: get_chr_order(x))
    sorted_out_df = out_df.sort_values(by=['chr_key', 'start']).drop(columns=['chr_key'])
    sorted_out_df.to_csv(outfile, sep="\t", header=False, index=False)

    out_b = pd.concat(out_breaks, ignore_index=True)
    out_b['chr_key'] = out_b['chrom'].apply(lambda x: get_chr_order(x))
    sorted_out_b = out_b.sort_values(by=['chr_key', 'start']).drop(columns=['chr_key'])
    sorted_out_b.to_csv(outbreak, sep="\t", header=False, index=False)



def merge_horclass(infile, outfile):
    df = pd.read_csv(infile, sep="\t", names=["chrom", "start", "end", "horclass", "score", "strand", "color"])

    out_subdfs = []

    for chrom, subdf in df.groupby('chrom'):

        data = subdf.values.tolist()
        idx_hor_class = 3
        idx_color = 6
        idx_strand = 5
        idx_start = 1
        idx_end = 2 

        hor = []
        for i, line in enumerate(data):
            if line[idx_hor_class] != 999:
                hor.append([i, line])
        # print(hor[:10])
        
        pivot = 0
        for i, line in enumerate(data):
            if line[idx_hor_class] == 999:
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
        
        pivot = 0
        for i, line in enumerate(data):
            if line[idx_hor_class] == 999:
                if pivot + 1 < len(hor):
                    while pivot + 1 < len(hor) and i > hor[pivot + 1][0]:
                        pivot += 1
                    if pivot + 1 < len(hor):
                        # print(pivot, i, len(hor))
                        if hor[pivot][0] < i < hor[pivot + 1][0]:
                            if hor[pivot][1][idx_hor_class] != hor[pivot + 1][1][idx_hor_class]:
                                if (line[idx_strand] == hor[pivot][1][idx_strand]) and (line[idx_start] - hor[pivot][1][idx_end] < 100):
                                    data[i][idx_hor_class] = hor[pivot][1][idx_hor_class]
                                    data[i][idx_color] = hor[pivot][1][idx_color]
                                elif (line[idx_strand] == hor[pivot + 1][1][idx_strand]) and (hor[pivot + 1][1][idx_start] - line[idx_end] < 100):
                                    data[i][idx_hor_class] = hor[pivot + 1][1][idx_hor_class]
                                    data[i][idx_color] = hor[pivot + 1][1][idx_color]


        updated_subdf = pd.DataFrame(data, columns=["chrom", "start", "end", "horclass", "score", "strand", "color"])
        out_subdfs.append(updated_subdf)

    out_df = pd.concat(out_subdfs, ignore_index=True)
    out_df['chr_key'] = out_df['chrom'].apply(lambda x: get_chr_order(x))
    sorted_df = out_df.sort_values(by=['chr_key', 'start']).drop(columns=['chr_key'])
    sorted_df.to_csv(outfile, sep="\t", header=False, index=False) 


def get_horclass_inversion_breakpoint(outhorclass, outhorclass_strand_merge, out_horclass_break, distance=10000):
    df = pd.read_csv(outhorclass, sep="\t", names = ["chrom", "start", "end", "horclass", "score", "strand", "color"])
   
    out_subdfs = []
    out_breaks = []

    for chrom, subdf in df.groupby('chrom'):

        tmp = []
        out = []
        breaks = [] 

        for idx, (i, row) in enumerate(subdf.iterrows()):
            if idx == 0:
                prev_strand = row['strand']
                prev_start = int(row['start'])
                prev_end = int(row['end'])
                prev_horclass = row['horclass']
                tmp.append([prev_start, prev_end])
            else:
                if (row['strand'] == prev_strand)  and (row['horclass'] == prev_horclass):
                    tmp.append([int(row['start']), int(row['end'])])
                else:
                    merged = merge(tmp, distance)
                    # print(chrom, i, merged)
                    for s, e in merged:
                        out.append([chrom, s, e, prev_strand, prev_horclass])

                    breaks.append([chrom, tmp[-1][-1], row['start'], row['start'] - merged[-1][-1], prev_strand+"/"+row['strand'], str(prev_horclass)+"/"+str(row['horclass'])])

                    prev_strand = row['strand']
                    prev_start = int(row['start'])
                    prev_end = int(row['end'])
                    prev_horclass = row['horclass']
                    tmp = [[prev_start, prev_end]]  

        ##last##
        merged = merge(tmp, distance)
        for s, e in merged:
            out.append([chrom, s, e, prev_strand, prev_horclass])
    
        out_subdf = pd.DataFrame(out, columns = ["chrom", "start", "end", "strand", "horclass"])
        out_subdfs.append(out_subdf)

        out_break = pd.DataFrame(breaks, columns = ["chrom", "start", "end", "break_distance", "strand", "horclass"])
        out_breaks.append(out_break)

    out_df = pd.concat(out_subdfs, ignore_index=True)
    out_df['chr_key'] = out_df['chrom'].apply(lambda x: get_chr_order(x))
    sorted_out_df = out_df.sort_values(by=['chr_key', 'start']).drop(columns=['chr_key'])
    sorted_out_df.to_csv(outhorclass_strand_merge, sep="\t", header=False, index=False)

    out_b = pd.concat(out_breaks, ignore_index=True)
    out_b['chr_key'] = out_b['chrom'].apply(lambda x: get_chr_order(x))
    sorted_out_b = out_b.sort_values(by=['chr_key', 'start']).drop(columns=['chr_key'])
    sorted_out_b.to_csv(out_horclass_break, sep="\t", header=False, index=False)   

def main(sample, hap, cenanno):

    os.makedirs(f"{sample}/{hap}/inversion", exist_ok=True)

    #cenanno = f"{sample}/{hap}/{sample}_{hap}.v0.9.full.cenanno.updatedY.bed"
    #cenanno = f"{sample}/{hap}/{sample}_{hap}.full.cenanno.updatedY.bed"
    if sample != "CHM13":
        # outfile = f"{sample}/{hap}/inversion/{sample}_{hap}_asat_strand_merge.bed"
        # outbreak = f"{sample}/{hap}/inversion/{sample}_{hap}_asat_inversion_breakpoint.bed"
        outfile = f"{sample}/{hap}/inversion/{sample}_{hap}_sat_strand_merge.bed"
        outbreak = f"{sample}/{hap}/inversion/{sample}_{hap}_sat_inversion_breakpoint.bed"
    else:
        outfile = f"CHM13/inversion/CHM13_sat_strand_merge.bed"
        outbreak = f"CHM13/inversion/CHM13_sat_inversion_breakpoint.bed"

    if os.path.exists(cenanno):
        get_sat_inversion_breakpoint(cenanno, outfile, outbreak)
    else:
        #cenanno = f"{sample}/{hap}/{sample}_{hap}.v0.9.full.cenanno.bed"
        #cenanno = f"{sample}/{hap}/{sample}_{hap}.full.cenanno.bed"
        # get_sat_inversion_breakpoint(cenanno, outfile, outbreak)
        print(f"{cenanno} not exists, exit...")
        exit(1)

    # if sample != "CHM13":
    #     inhorclass = f"{sample}/{hap}/{sample}_{hap}.HORclass.bed"
    #     outhorclass = f"{sample}/{hap}/inversion/{sample}_{hap}.HORclass.merge.bed"
    # else:
    #     inhorclass = f"CHM13/CHM13.HORclass.bed"
    #     outhorclass = f"CHM13/inversion/CHM13.HORclass.merge.bed"
    # merge_horclass(inhorclass, outhorclass)

    # if sample != "CHM13":
    #     outhorclass_strand_merge = f"{sample}/{hap}/inversion/{sample}_{hap}.HORclass.strand.merge.bed"
    #     out_horclass_break = f"{sample}/{hap}/inversion/{sample}_{hap}.HORclass.inversion_breakpoint.bed"
    # else:
    #     outhorclass_strand_merge = f"CHM13/inversion/CHM13.HORclass.strand.merge.bed"
    #     out_horclass_break = f"CHM13/inversion/CHM13.HORclass.inversion_breakpoint.bed"
    # get_horclass_inversion_breakpoint(outhorclass, outhorclass_strand_merge, out_horclass_break)

if __name__ == "__main__":
    sample=sys.argv[1]
    hap=sys.argv[2]
    cenanno = sys.argv[3]
    main(sample, hap, cenanno)

