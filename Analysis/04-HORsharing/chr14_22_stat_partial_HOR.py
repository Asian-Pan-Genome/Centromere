import os
import subprocess
import pandas as pd
import sys
import pysam
import random
import re
import time
print("Current working directory:", os.getcwd())


# def run_bedtools_getfasta(genome, bed_file, output_fa):
#     command = f"bedtools getfasta -fi {genome} -bed {bed_file} -fo {output_fa}"
#     print(f"Running command: {command}")
#     result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
#     if result.returncode != 0:
#         print(f"Error running command: {result.stderr}")
#         raise RuntimeError(f"Bedtools command failed: {result.stderr}")
    
#     # Wait for the output file to be generated
#     for _ in range(10):  # Retry for up to 10 seconds
#         if os.path.exists(output_fa) and os.path.getsize(output_fa) > 0:
#             break
#         time.sleep(1)
#     else:
#         raise RuntimeError(f"Output file {output_fa} is empty or not generated.")

def run_bedtools_getfasta(genome, bedfile, outfa):
    genome = os.path.abspath(genome)
    bedfile = os.path.abspath(bedfile)
    outfa = os.path.abspath(outfa)

    print(f"[INFO] Genome fasta: {genome}")
    print(f"[INFO] BED file: {bedfile}")
    print(f"[INFO] Output fasta: {outfa}")

    if os.path.exists(bedfile):
        with open(bedfile) as f:
            for line in f:
                print("[INFO] First BED line:", line)
    else:
        print("[ERROR] BED file does not exist!")
        return

    try:
        result = subprocess.run(
            ["bedtools", "getfasta", "-fi", genome, "-bed", bedfile, "-fo", outfa],
            capture_output=True,
            text=True,
            check=True
        )
        print("[INFO] bedtools finished successfully.")
        if result.stdout:
            print("[STDOUT]\n", result.stdout)
        if result.stderr:
            print("[STDERR]\n", result.stderr)

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] bedtools failed with return code {e.returncode}")
        print("[STDOUT]\n", e.stdout)
        print("[STDERR]\n", e.stderr)

def run_bedtools_intersection(mnbed, regionbed, outbed):
    command = f"bedtools intersect -a {mnbed} -b {regionbed} -wa -wb -f 0.95 > {outbed}"
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running bedtools: {e}")

def concat_fasta_sequences(file_a, file_b):

    seq_a = None
    seq_b = None

    if not os.path.exists(file_a):
        print(f"File {file_a} not found.")
        seq_a = ""
    if not os.path.exists(file_b):
        print(f"File {file_b} not found.")
        seq_b = ""
    
    def get_seq(infile):
        with open(infile, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    continue
                else:
                    return line
    if seq_a is None:
        seq_a = get_seq(file_a)
    if seq_b is None:
        seq_b = get_seq(file_b)
    out_seq = seq_a + seq_b
    return out_seq
    

def extract_position(outprefix, mn_bed, ihor, genome, idx, ihor_flag):
    sequence = ihor.split('_')
    print("ihor:", sequence)

    # shared_sequences = ["3072", "24569", "3042", "25511", "3379", "2975"]

    original_positions = []
    mns = []
    strands = []
    start_index = None
    start_shared_index = None
    end_shared_index = None
    index = 0
    out_seq = None
    out_shared_seq = None
    oripos = pd.DataFrame()

    with open(mn_bed, 'r') as f:
        for line in f:
            tokens = line.strip().split('\t')
            mn = tokens[3]
            mns.append(tokens[3])
            strands.append(tokens[5])
            original_positions.append([int(tokens[1]), int(tokens[2])])

            if mn == sequence[0] and start_index is None:
                start_index = index

            index += 1

        print("original mns:", mns)
        print("start_index:", start_index)

        print("original positions:", original_positions[0][0], original_positions[-1][-1])


        if len(set(strands)) == 1:

            oripos_info = {
                "sample_hap_chrom": [outprefix],
                "start": [original_positions[0][0]],
                "end": [original_positions[-1][-1]],
                "name": f"{outprefix}_{ihor_flag}_{idx}"
            }
            oripos = pd.DataFrame(oripos_info)
            print(oripos)

            reorder_positions_s = original_positions[start_index:]
            reorder_positions_e = original_positions[:start_index]

            if start_index != 0:

                ##refine the boundary
                boundary = min(reorder_positions_s[0][0], reorder_positions_e[-1][1])

                with open(f"tmp/{outprefix}_tmp_s.bed", 'w') as out1:
                    out1.write(f"{outprefix}\t{boundary}\t{reorder_positions_s[-1][1]}\n")
                    out1.flush()
                    os.fsync(out1.fileno())
                    run_bedtools_getfasta(genome, f"tmp/{outprefix}_tmp_s.bed", f"tmp/{outprefix}_tmp_s.fa")

                with open(f"tmp/{outprefix}_tmp_e.bed", 'w') as out2:
                    out2.write(f"{outprefix}\t{reorder_positions_e[0][0]}\t{boundary}\n")
                    out2.flush()
                    os.fsync(out2.fileno())
                    run_bedtools_getfasta(genome, f"tmp/{outprefix}_tmp_e.bed", f"tmp/{outprefix}_tmp_e.fa")

            else:
                with open(f"tmp/{outprefix}_tmp_s.bed", 'w') as out1:
                    out1.write(f"{outprefix}\t{reorder_positions_s[0][0]}\t{reorder_positions_s[-1][1]}\n")
                    out1.flush()
                    os.fsync(out1.fileno())
                    run_bedtools_getfasta(genome, f"tmp/{outprefix}_tmp_s.bed", f"tmp/{outprefix}_tmp_s.fa")

            out_seq = concat_fasta_sequences(f"tmp/{outprefix}_tmp_s.fa", f"tmp/{outprefix}_tmp_e.fa")

            #delete tmp file
            if os.path.exists(f"tmp/{outprefix}_tmp_s.bed"):
                os.remove(f"tmp/{outprefix}_tmp_s.bed")
            if os.path.exists(f"tmp/{outprefix}_tmp_e.bed"):
                os.remove(f"tmp/{outprefix}_tmp_e.bed")
            if os.path.exists(f"tmp/{outprefix}_tmp_s.fa"):
                os.remove(f"tmp/{outprefix}_tmp_s.fa")
            if os.path.exists(f"tmp/{outprefix}_tmp_e.fa"):
                os.remove(f"tmp/{outprefix}_tmp_e.fa")
        else:
            print(f"Strands are not consistent, pass")
    return oripos, out_seq, out_shared_seq


def stat_target_hor(horstv_bed, mnbed, genome, outprefix):
    horstv = pd.read_csv(horstv_bed, sep="\t", names=["chrom", "start", "end", "hor", "horclass", "hccolor", "horstv", "horcolor"])
    horstv['horclass'] = horstv['horclass'].astype(str)

    r_chrom = outprefix

    horstv = horstv[(horstv['chrom'] == r_chrom) & (horstv['horclass'] == '94')].copy()


    r0start = horstv['start'].min()
    r2end = horstv['end'].max()

    total_length = r2end - r0start

    target_hor = {
    "1339_25692_2959_25361_1428_25301_1381_24708": "h1-8mer", 
    }
    # target_hor = {"2975_24578_3036_25494_3059_24579_3072_24569_3042_25511_3379" : "h2_11mer-24579"}

    horstv['horflag'] = horstv['hor'].map(target_hor).fillna("Others")
    horstv['length'] = horstv['end'] - horstv['start']
    print(horstv.head())

    ### extract sequences
    if ('CN1' not in horstv_bed) and ('HG002' not in horstv_bed) and ('YAO' not in horstv_bed):
        for ihor, ihor_flag in target_hor.items():
            outdf = horstv[horstv['horflag'] == ihor_flag].copy().reset_index(drop=True)
            print(f"{ihor_flag} dataframe: ", outdf.head())

            out_seqs = {}
            out_shared_seqs = {}
            oripos_dfs = []

            for idx, row in outdf.iterrows():
                region_bed = f"tmp/{outprefix}_{ihor_flag}_{idx}.bed"
                outbed = f"tmp/{outprefix}_{ihor_flag}_{idx}_mn.bed"
                out_seq = None
                with open(region_bed, 'w') as f:
                    f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\n")
                run_bedtools_intersection(mnbed, region_bed, outbed)
                oripos, out_seq, out_shared_seq = extract_position(outprefix, outbed, ihor, genome, idx, ihor_flag)
                if out_seq is not None:
                    out_seqs[f"{outprefix}_{ihor_flag}_{idx}"] = out_seq
                if out_shared_seq is not None:
                    out_shared_seqs[f"{outprefix}_{ihor_flag}_{idx}"] = out_shared_seq
                if not oripos.empty:
                    oripos_dfs.append(oripos)
                    
                os.remove(f"tmp/{outprefix}_{ihor_flag}_{idx}.bed")
                os.remove(f"tmp/{outprefix}_{ihor_flag}_{idx}_mn.bed")
            
            with open(f"tmp/{outprefix}_{ihor_flag}.fa", 'w') as f:
                for seq_name, sequence in out_seqs.items():
                    f.write(f">{seq_name}\n{sequence}\n")

            with open(f"tmp/{outprefix}_{ihor_flag}_shared.fa", 'w') as f:
                for seq_name, sequence in out_shared_seqs.items():
                    f.write(f">{seq_name}\n{sequence}\n")
            
            print("lenth of oripos_dfs:", len(oripos_dfs))
            if len(oripos_dfs) > 0:
                merged_oripos = pd.concat(oripos_dfs, ignore_index=True)
                merged_oripos.to_csv(f"tmp/{outprefix}_{ihor_flag}_oripos.bed", sep="\t", index=False)

    ##statistics
    stat_hor = horstv.groupby('horflag').agg(
    target_length=('length', 'sum'),
    target_count=('horflag', 'count')
).reset_index()
    stat_hor['length_percent'] = stat_hor['target_length'] / total_length * 100
    print(stat_hor)
    stat_hor.to_csv(f"tmp/{outprefix}_horstat.txt", sep="\t", index=False, header=False)
    return stat_hor


def main(sample, hap, chrom):

    mnbed = ""
    outprefix = ""
    outprefix = ""

    if sample == "CHM13":
        genome = "/share/home/zhanglab/user/sunyanqing/human/reference/CHM13/CHM13v2.fasta"
        mnbed = "/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-10/CHM13/filter.asat.bed"
        outprefix = chrom

    elif sample == "CN1":
        t="#"
        if hap == "Mat":
            genome = "/share/home/zhanglab/user/sunyanqing/human/reference/CN1v1.0/MF2_mat.v1.0.fa"
        else:
            genome = "/share/home/zhanglab/user/sunyanqing/human/reference/CN1v1.0/MF2_pat.v1.0.fa"
    elif sample == "HG002":
        t="#"
        if hap == "Mat":
            genome = "/share/home/zhanglab/user/sunyanqing/human/reference/HG002/hg002v1.0.1.mat.fasta"
        else:
            genome = "/share/home/zhanglab/user/sunyanqing/human/reference/HG002/hg002v1.0.1.pat.fasta"

    elif sample == "YAO":
        t="#"
        if hap == "Mat":
            genome = "/share/home/zhanglab/user/sunyanqing/human/reference/YAO/yao_mat.fasta"
        else:
            genome = "/share/home/zhanglab/user/sunyanqing/human/reference/YAO/yao_pat.fasta"
    elif sample.startswith('HG') or sample.startswith('NA'):
        t="_"
        if hap == "Mat" or hap == "Pat":
            genome = f"/share/home/zhanglab/user/sunyanqing/human/HPRC/{sample}/{hap}/{sample}_{hap}_chrR.fasta"
        else:
            genome = f"/share/home/zhanglab/user/sunyanqing/human/HGSVC/{sample}/{hap}/{sample}_{hap}_chrR.fasta"
    else:
        t="#"
        genome = f"/share/home/zhanglab/user/sunyanqing/human/assembly/{sample}/{hap}/{sample}_{hap}.v0.9.fasta"

    if mnbed == "":
        mnbed = f"/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-10/{sample}/{hap}/filter.asat.bed"
    if outprefix == "":
        outprefix = f"{sample}{t}{hap}{t}{chrom}"

    # horstv_bed = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HiCAT_bhlayer_update_split/{sample}_{hap}.HiCAT.horstv.bed"
    horstv_bed = f"/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HORmon_split/{sample}_{hap}.graph.hordecomposition.final.xls"
    stat_target_hor(horstv_bed, mnbed, genome, outprefix)


def select_phy(n, output_file):

    merged_fa = pysam.FastaFile(f"merged_seq/all_chr13-21_3mer.fa")

    merged_sequences = []
    
    seq_names = merged_fa.references
    selected_names = random.sample(seq_names, min(n, len(seq_names)))
    for seq_name in selected_names:
        sequence = merged_fa.fetch(seq_name)
        merged_sequences.append(f">{seq_name}\n{sequence}\n")


    with open(output_file, "w") as f:
        f.writelines(merged_sequences)

    

if __name__ == "__main__":
    sample = sys.argv[1]
    hap = sys.argv[2]
    chrom = sys.argv[3]
    main(sample, hap, chrom)

    # n = int(sys.argv[1])
    # output_file = sys.argv[2]
    # select_phy(n, output_file)


    



