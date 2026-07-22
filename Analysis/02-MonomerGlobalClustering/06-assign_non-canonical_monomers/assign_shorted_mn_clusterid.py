import sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from collections import defaultdict


def read_fa(infile):
    sequences = list(SeqIO.parse(infile, "fasta"))
    return sequences


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
    identity = round(matching / shorter_length,3) if shorter_length > 0 else 0
    return identity


def load_nid_order():
    #order_nid = {}
    nid_order = {}
    with open("fasttree.order", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            nid = tokens[0]
            order = tokens[1]
            #order_nid[order] = nid
            nid_order[nid] = order
    print(list(nid_order.keys())[:10])
    return nid_order


def reassign_nonclassical_seq(fafile, nid_order):
    query_nid = {}
    query_sequences = read_fa(fafile)
    target_sequences = read_fa("marker_seq.fa")

    for query in query_sequences:
        tmp = []
        for target in target_sequences:
            identity = pairwise(query, target)
            tmp.append(identity)

        max_index = tmp.index(max(tmp))
        max_target = target_sequences[max_index]
        mapping_nid = max_target.id.split('@')[0]
        query_nid[query.id] = nid_order[mapping_nid]

    with open(f"{fafile}.index", "w") as outf:
        for query, order in query_nid.items():
            outf.write(f"{query}\t{order}\n")


if __name__ == "__main__":
    fafile = sys.argv[1]
    nid_order = load_nid_order()
    reassign_nonclassical_seq(fafile, nid_order)

