import pandas as pd
import random
import pysam
import time

def select_random_sequence():
    ffai = "/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr8_mnphy/shared_sequence/all_shared_dedup.filter.fa.fai"
    fafile = "/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr8_mnphy/shared_sequence/all_shared_dedup.filter.aln.fa"
    df = pd.read_csv(ffai, sep="\t", names=['seq', 'len', 'a', 'b', 'c'])
    fa = pysam.FastaFile(fafile)

    h1_11mer_lst = list(df[df['seq'].str.contains('h1_11mer')]['seq'])
    h2_8mer_lst = list(df[df['seq'].str.contains('h2_8mer')]['seq'])
    h3_7mer_lst = list(df[df['seq'].str.contains('h3_7mer')]['seq'])

    print(len(h1_11mer_lst))
    print(len(h2_8mer_lst))
    print(len(h3_7mer_lst))

    random_h1_lst = random.sample(h1_11mer_lst, 1000) 
    random_h2_lst = random.sample(h2_8mer_lst, 1000)
    random_h3_lst = random.sample(h3_7mer_lst, 1000)

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

    flags = ['h1_11mer', 'h2_8mer', 'h3_7mer']
    intra_lst = [random_h1_lst, random_h2_lst, random_h3_lst]

    for i, flag in enumerate(flags):
        flag = flags[i]
        ih = intra_lst[i]

        start_time = time.time()
        intra_div = intra_compare(ih)
        elapsed_time = time.time() - start_time
        print(f"complete for {flag} in {elapsed_time} seconds")

        with open(f"{flag}_intra.div", 'w') as outf:
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

        with open(f"{iflag}_vs_{jflag}_inter.div", 'w') as outf:
                for idiv in inter_ih_jh:
                    outf.write(f"{iflag}_vs_{jflag}\t{idiv[0]}\t{idiv[1]}\t{idiv[2]}\n")

if __name__ == "__main__":
    select_random_sequence()