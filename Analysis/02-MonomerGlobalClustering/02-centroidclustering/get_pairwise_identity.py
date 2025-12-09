import sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from collections import defaultdict


def read_fa():
    sequences = list(SeqIO.parse("all.centroid.fasta", "fasta"))
    return sequences

def stat_cluster():
    cluster_seqs = {}
    with open("all.centroid.txt", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            cluster_seqs[tokens[1]] = tokens[8]
    #print(len(cluster_seqs))
    return cluster_seqs

def pairwise(iseq, jseq):
    aligner=PairwiseAligner()
    aligner.mode = 'global'

    alignments = aligner.align(iseq.seq, jseq.seq)
    alignment = alignments[0]
    matching = 0
    for (start_a, end_a), (start_b, end_b) in zip(*alignment.aligned):
        segment_a = iseq.seq[start_a:end_a]
        segment_b = jseq.seq[start_b:end_b]
        matching += sum(1 for a, b in zip(segment_a, segment_b) if a == b)
    shorter_length = min(len(iseq.seq), len(jseq.seq))
    #alignment_length = alignment.length
    identity = round(matching / shorter_length,2) if shorter_length > 0 else 0
    return identity

def one_vs_all(i, sequences, cluster_seqs, outfile):
    outf = open(outfile, "w")
     
    index = i+1
    iseq = sequences[i]

    for j in range(index, len(cluster_seqs)):
        jseq = sequences[j]
        print(i, j)
        identity = pairwise(iseq, jseq)
        outf.write(f"{i}\t{j}\t{identity}\n")
        
    outf.close()

if __name__== "__main__":
    i = int(sys.argv[1])
    sequences = read_fa()
    cluster_seqs = stat_cluster()
    outfile = f"tmp/centroid_{i}.pairwise.txt"
    one_vs_all(i, sequences, cluster_seqs, outfile)
