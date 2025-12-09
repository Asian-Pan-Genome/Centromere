import os
import subprocess
import pandas as pd
import sys
import pysam
import random
import re
from Bio.Seq import Seq

def run_bedtools_getfasta(genome, outfile, outfa):
    command = f"bedtools getfasta -fi {genome} -bed {outfile} -name -s -fo {outfa}"
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running bedtools: {e}")

def run_bedtools_intersection(mnbed, regionbed, outbed):
    command = f"bedtools intersect -a {mnbed} -b {regionbed} -wa -wb > {outbed}"
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running bedtools: {e}")

def run_bedtools_subtract(r1bed, dupbed, outbed):
    command = f"bedtools subtract -a {r1bed} -b {dupbed} > {outbed}"
    rmcommand = f"rm -f {r1bed}"
    try:
        subprocess.run(command, shell=True, check=True)
        subprocess.run(rmcommand, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running bedtools: {e}")



def extract_position(outprefix, outdimer, inverted_flag):

    if inverted_flag:
        sequence = ["21663", "29601"]
        strand_match = '-'
    else:
        sequence = ["29601", "21663"]
        strand_match = '+'

    current_index = 0
    first_start = None
    last_end = None
    index = 0

    with open(f"{outprefix}_mn.bed", 'r') as f, open(outdimer, 'w') as outf:
        for line in f:
            tokens = line.strip().split('\t')
            mn = tokens[3]
            strand = tokens[5]

            # Check if the current mn matches the expected sequence and strand is '-'
            if mn == sequence[current_index] and strand == strand_match:
                if mn == sequence[0]:
                    first_start = int(tokens[1])
                if mn == sequence[-1]:
                    last_end = int(tokens[2])

                current_index += 1  # Move to the next mn in the sequence

                if current_index == len(sequence):
                    outf.write(f"{tokens[0]}\t{first_start}\t{last_end}\t{outprefix}_{index}\t0\t{strand}\n")
                    current_index = 0
                    first_start = None
                    last_end = None
                    index += 1
            else:
                # Reset if the sequence is broken
                current_index = 0
                first_start = None
                last_end = None

            


def extract_dimer(regionbed, mnbed, genome, inverted_flag):

    r_basename = os.path.basename(regionbed)
    outprefix = os.path.splitext(r_basename)[0]
    print(regionbed, r_basename, outprefix)

    outbed = f"{outprefix}_mn.bed"
    outdimer = f"{outprefix}_dimer.bed"
    # tmpfa = f"{outprefix}_dimer.tmp.fa"
    outfa = f"{outprefix}_dimer.fa"

    run_bedtools_intersection(mnbed, regionbed, outbed)
    extract_position(outprefix, outdimer, inverted_flag)
    run_bedtools_getfasta(genome, outdimer, outfa)


def extract_target_dimer(horstv_bed, mnbed, genome, outprefix, chrom, inverted_flag):
    horstv = pd.read_csv(horstv_bed, sep="\t", names=["chrom", "start", "end", "hor", "horclass", "hccolor", "horstv", "horcolor"])
    horstv['horclass'] = horstv['horclass'].astype(str)

    r_chrom = outprefix

    if chrom == 'chr16':
        horstv = horstv[(horstv['chrom'] == r_chrom) & (horstv['horclass'].isin(['24', '25']))].copy()
    else:
        horstv = horstv[(horstv['chrom'] == r_chrom) & (horstv['horclass'] == '25')].copy()

    r0start = horstv['start'].min()
    r2end = horstv['end'].max()

    interval = 100000

    #r0#
    r0end = r0start + interval
    with open(f"{outprefix}_r0.bed", 'w') as out:
        out.write(f"{r_chrom}\t{r0start}\t{r0end}\n")
    r0bed = f"{outprefix}_r0.bed"
    print(r0bed)

    #r2#
    r2end = horstv['end'].max()
    r2start = r2end - interval
    with open(f"{outprefix}_r2.bed", 'w') as out:
        out.write(f"{r_chrom}\t{r2start}\t{r2end}\n")
    r2bed = f"{outprefix}_r2.bed"
    print(r2bed) 

    #r1#
    r1start = int((r0start + r2end) / 2 ) - int(interval / 2)
    r1end = int((r0start + r2end) / 2 ) + int(interval / 2)
    with open(f"{outprefix}_r1.bed", 'w') as out:
        out.write(f"{r_chrom}\t{r1start}\t{r1end}\n")
    r1bed = f"{outprefix}_r1.bed"
    print(r1bed)  
    

    if chrom in ['chr1', 'chr5', 'chr19']:
        extract_dimer(r0bed, mnbed, genome, inverted_flag)
        extract_dimer(r2bed, mnbed, genome, inverted_flag)
        extract_dimer(r1bed, mnbed, genome, inverted_flag)
    else:
        hicat_bed = f"/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-10/{sample}/{hap}/{sample}_{hap}.HORstv.bed"
        hicat_df = pd.read_csv(hicat_bed, sep="\t", names=["chrom", "start", "end", "hor", "score", "strand",  "horcolor", "layer"])
        chr16_dup = hicat_df[(hicat_df['chrom'] == outprefix) & (hicat_df['hor'] == 'C25H2403(2)')].copy()
        dupbed = f"{outprefix}_dup.bed"
        chr16_dup.to_csv(dupbed, sep="\t", header=False, index=False)

        r0_filtered_bed = f"{outprefix}_r0_filtered.bed"
        r2_filtered_bed = f"{outprefix}_r2_filtered.bed"
        r1_filtered_bed = f"{outprefix}_r1_filtered.bed"
        run_bedtools_subtract(r0bed, dupbed, r0_filtered_bed)
        run_bedtools_subtract(r2bed, dupbed, r2_filtered_bed)
        run_bedtools_subtract(r1bed, dupbed, r1_filtered_bed)

        extract_dimer(r0_filtered_bed, mnbed, genome, inverted_flag)
        extract_dimer(r2_filtered_bed, mnbed, genome, inverted_flag)
        extract_dimer(r1_filtered_bed, mnbed, genome, inverted_flag)
        extract_dimer(dupbed, mnbed, genome, inverted_flag)

 

def main(sample, hap, chrom):

    centdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", header=0, sep="\t")
    target_cent = centdf[(centdf['sample'] == sample) & (centdf['hap'] == hap) & (centdf['chrom'] == chrom)].copy()
    project = target_cent['project'].unique()[0]

    if project == "APG":
        genome = f"/share/home/zhanglab/user/sunyanqing/human/assembly/{sample}/{hap}/{sample}_{hap}.v0.9.fasta"
    else:
        genome = f"/share/home/zhanglab/user/sunyanqing/human/{project}/{sample}/{hap}/{sample}_{hap}_chrR.fasta"

    mnbed = f"/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-10/{sample}/{hap}/filter.asat.bed"
    outprefix = target_cent['sample_hap_chrom'].unique()[0]

    ##for inverted chr1##
    cols = ['sample_hap_chrom', 'start', 'end', 'strand', 'horclass', 'index', 'length']
    inverted_chr1_df = pd.read_csv("chr1_25_inversion_region_length.xls", names = cols, sep="\t")
    inverted_flag = True if outprefix in inverted_chr1_df['sample_hap_chrom'].values else False

    ##hicat_bed##
    horstv_bed = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HiCAT_bhlayer_update_split/{sample}_{hap}.HiCAT.horstv.bed"
    extract_target_dimer(horstv_bed, mnbed, genome, outprefix, chrom, inverted_flag)


def select_dimer_phy(n, output_file):

    dimer_fa = pysam.FastaFile(f"merged_seq/all_chr1-5-16-19.fa")

    merged_sequences = []
    seq_names = dimer_fa.references
    selected_names = random.sample(seq_names, min(n, len(seq_names)))
    for seq_name in selected_names:
        sequence = dimer_fa.fetch(seq_name)
        merged_sequences.append(f">{seq_name}\n{sequence}\n")


    with open(output_file, "w") as f:
        f.writelines(merged_sequences)



def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())



def inverted_genome():

    cols = ['sample_hap_chrom', 'start', 'end', 'strand', 'horclass', 'index', 'length']
    inverted_chr1_df = pd.read_csv("chr1_25_inversion_region_length.xls", names = cols, sep="\t")
    inverted_region_df = inverted_chr1_df[inverted_chr1_df['strand'] == '-'].copy()

    centdf = pd.read_csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/cent_chrom.txt", header=0, sep="\t")
    centdf = centdf[['sample_hap_chrom', 'sample_hap', 'sample', 'hap', 'chrom', 'project', 'filterflag']].copy()
    
    inverted_region_df = inverted_region_df.merge(centdf, on='sample_hap_chrom', how='left')
    complete_inverted_region_df = inverted_region_df[(inverted_region_df['filterflag'] == 0) & (inverted_region_df['project'] != 'Ref')].copy()
    print(len(complete_inverted_region_df['sample_hap_chrom'].unique()))

    print(complete_inverted_region_df.head())
    
    for name, group in complete_inverted_region_df.groupby(['sample_hap_chrom', 'sample', 'hap', 'project']):
        print(name)
        sample_hap_chrom = name[0]
        sample = name[1]
        hap = name[2]
        project = name[3]

        if project == "APG":
            genome = f"/share/home/zhanglab/user/sunyanqing/human/assembly/{sample}/{hap}/{sample}_{hap}.v0.9.fasta"
        else:
            genome = f"/share/home/zhanglab/user/sunyanqing/human/{project}/{sample}/{hap}/{sample}_{hap}_chrR.fasta"

        g = pysam.FastaFile(genome)
        chr1_seq = g.fetch(sample_hap_chrom)

        inverted_pos = []
        for index, row in group.iterrows():
            inverted_pos.append((row['start'], row['end']))
        inverted_pos.sort(key=lambda x: x[0])

        concat_seq = ""
        last_end = 0

        for start, end in inverted_pos:
            concat_seq += chr1_seq[last_end:start]
            
            inverted_seq = chr1_seq[start:end]
            modified_seq = reverse_complement(inverted_seq)
            
            concat_seq += modified_seq
            
            last_end = end
        concat_seq += chr1_seq[last_end:]

        
        with open(f"inverted_chr1_seq/{sample_hap_chrom}.mod.fa", 'w') as out:
            out.write(f">{sample_hap_chrom}\n{concat_seq}\n")
    





if __name__ == "__main__":
    # inverted_genome()

    #sample = sys.argv[1]
    #hap = sys.argv[2]
    #chrom = sys.argv[3]
    #main(sample, hap, chrom)

    n = int(sys.argv[1])
    output_file = sys.argv[2]
    select_dimer_phy(n, output_file)


    



