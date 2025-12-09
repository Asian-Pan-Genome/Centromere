import pandas as pd
import random
import pysam
import time

def select_pairwise():
    ffai = "/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr13_21_mnphy/merged_hors/all_shared_filtered_dedup.aln.fa.fai"
    df = pd.read_csv(ffai, sep="\t", names = ['seqname', 'length', 'a', 'b', 'c'])

    h1_chr13 = []
    h2_chr13 = []
    h1_chr21 = []
    h2_chr21 = []
    h3_chr21 = []
    h3_chr13 = []
    h4_chr13 = []
    h4_chr21 = []

    '''
    pairwise compare

    #### intra divergence ### 
        h1_chr13 intra;
        h1_chr21 intra;
        h2_chr13 intra;
        h2_chr21 intra;
        h3_chr21 intra;
        h3_chr13 intra;
        h4_chr13 intra;
        h4_chr21 intra
    
    ###group one to confirm h3 origin###
        h1_chr13 vs. h3_chr21; 
        h2_chr13 vs. h3_chr21; 
        h2_chr13 vs. h3_chr21; 
        h2_chr21 vs. h3_chr21;


    ####group two to confirm h4 origin ###
        h1_chr13 vs h4_chr13
        h1_chr21 vs h4_chr13
        h2_chr13 vs h4_chr13
        h2_chr21 vs_h4_chr13
    '''

    for seq in list(df['seqname']):
        if "chr13_h1" in seq:
            h1_chr13.append(seq)
        elif "chr13_h2" in seq:
            h2_chr13.append(seq)
        elif "chr13_h3" in seq:
            h3_chr13.append(seq)
        elif "chr13_h4" in seq:
            h4_chr13.append(seq)
        elif "chr21_h1" in seq:
            h1_chr21.append(seq)
        elif "chr21_h2" in seq:
            h2_chr21.append(seq)
        elif "chr21_h3" in seq:
            h3_chr21.append(seq)
        elif "chr21_h4" in seq:
            h4_chr21.append(seq)
        else:
            pass


    print("total number of h1 on chr13: ",  len(h1_chr13))
    print("total number of h2 on chr13: ",  len(h2_chr13))
    print("total number of h3 on chr13: ",  len(h3_chr13))
    print("total number of h4 on chr13: ",  len(h4_chr13))
    print("total number of h1 on chr21: ",  len(h1_chr21))
    print("total number of h2 on chr21: ",  len(h2_chr21))
    print("total number of h3 on chr21: ",  len(h3_chr21))
    print("total number of h4 on chr21: ",  len(h4_chr21))

    target_list = [h1_chr13, h2_chr13, h4_chr13, h1_chr21, h2_chr21, h3_chr21, h4_chr21]
    random_list = []
    for ih in target_list:
        if len(ih) > 1000:
            selected_seqs = random.sample(ih, 1000)
            random_list.append(selected_seqs)
        else:
            random_list.append(ih)

    def get_divergent_snp(iseq, jseq):
        n = 0
        if len(iseq) == len(jseq):
            for i in range(len(iseq)):
                if iseq[i] != "-" and jseq[i] != "-" and iseq[i] != jseq[i]:
                    n += 1
        return n

    def intra_compare(ih):
        out = []
        for i in range(len(ih)):
            for j in range(i+1, len(ih)):
                iseq = fa.fetch(ih[i])
                jseq = fa.fetch(ih[j])
                div_n = get_divergent_snp(iseq, jseq)
                out.append((ih[i], ih[j], div_n))
        return out

    def inter_compare(ih, jh):
        out = []
        for i in range(len(ih)):
            for j in range(len(jh)):
                iseq = fa.fetch(ih[i])
                jseq = fa.fetch(jh[j])
                div_n = get_divergent_snp(iseq, jseq)
                out.append((ih[i], jh[j],div_n))
        return out

    fa = pysam.FastaFile("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr13_21_mnphy/merged_hors/all_shared_filtered_dedup.aln.fa")
     ##intra
    flag  = ['h1_chr13', 'h2_chr13', 'h4_chr13', 'h1_chr21', 'h2_chr21', 'h3_chr21', 'h4_chr21']
    for i in range(len(flag)):
        iflag = flag[i]
        ih = random_list[i]

        start_time = time.time()
        intra_div = intra_compare(ih)
        elapsed_time = time.time() - start_time
        print(f"complete for {iflag} in {elapsed_time} seconds")

        with open(f"{iflag}_intra.div", 'w') as outf:
            for idiv in intra_div:
                outf.write(f"{iflag}_intra\t{idiv[0]}\t{idiv[1]}\t{idiv[2]}\n")

    intercompare = [(0,5), (1,5), (3,5), (0,2), (1,2), (3,2)]

    for k in intercompare:
        i, j = k[0],k[1]
        iflag = flag[i]
        jflag = flag[j]

        ih = random_list[i]
        jh = random_list[j]

        start_time = time.time()
        inter_ih_jh = inter_compare(ih, jh)
        elapsed_time = time.time() - start_time
        print(f"complete for {iflag} vs. {jflag} in {elapsed_time} seconds")

        with open(f"{iflag}_vs_{jflag}_inter.div", 'w') as outf:
                for idiv in inter_ih_jh:
                    outf.write(f"{iflag}_vs_{jflag}\t{idiv[0]}\t{idiv[1]}\t{idiv[2]}\n")




if __name__ == "__main__":
    select_pairwise()

