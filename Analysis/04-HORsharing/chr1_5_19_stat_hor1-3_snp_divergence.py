import pandas as pd
import random
import pysam
import time

def select_random_sequence():
    ffai = "/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr1_5_16_19_mnphy/merged_hor1-3/all_hor1-3.uniq.selected.aln.fa.fai"
    fafile = "/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr1_5_16_19_mnphy/merged_hor1-3/all_hor1-3.uniq.selected.aln.fa"
    df = pd.read_csv(ffai, sep="\t", names=['seq', 'len', 'a', 'b', 'c'])
    fa = pysam.FastaFile(fafile)

    chr19_lst = list(df[df['seq'].str.contains('chr19_')]['seq'])
    chr5_lst = list(df[df['seq'].str.contains('chr5_')]['seq'])
    chr1_lst = list(df[df['seq'].str.contains('chr1_')]['seq'])

    print(len(chr19_lst))
    print(len(chr5_lst))
    print(len(chr1_lst))

    random_h1_lst = random.sample(chr19_lst, 1000) 
    random_h2_lst = random.sample(chr5_lst, 1000)
    random_h3_lst = random.sample(chr1_lst, 1000)

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

    flags = ['chr19', 'chr5', 'chr1']
    intra_lst = [random_h1_lst, random_h2_lst, random_h3_lst]

    for i, flag in enumerate(flags):
        flag = flags[i]
        ih = intra_lst[i]

        start_time = time.time()
        intra_div = intra_compare(ih)
        elapsed_time = time.time() - start_time
        print(f"complete for {flag} in {elapsed_time} seconds")

        with open(f"/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr1_5_16_19_mnphy/merged_hor1-3/{flag}_intra.div", 'w') as outf:
            for idiv in intra_div:
                outf.write(f"{flag}_intra\t{idiv[0]}\t{idiv[1]}\t{idiv[2]}\n")

    
    intercompare = [(0,1), (0,2), (1,2)]
    for k in intercompare:
        i, j = k[0],k[1]
        iflag = flags[i]
        jflag = flags[j]

        ih = intra_lst[i]
        jh = intra_lst[j]

        start_time = time.time()
        inter_ih_jh = inter_compare(ih, jh)
        elapsed_time = time.time() - start_time
        print(f"complete for {iflag} vs. {jflag} in {elapsed_time} seconds")

        with open(f"/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr1_5_16_19_mnphy/merged_hor1-3/{iflag}_vs_{jflag}_inter.div", 'w') as outf:
                for idiv in inter_ih_jh:
                    outf.write(f"{iflag}_vs_{jflag}\t{idiv[0]}\t{idiv[1]}\t{idiv[2]}\n")

if __name__ == "__main__":
    select_random_sequence()