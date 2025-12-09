from get_continuous_block import *
from HiCAT_TD import *
from HORmon_graph import *
from graph_stat import *
from HiCAT_stat import *
from graph_decomposition_stat import *
import argparse
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Mining HOR using HORmon or HiCAT with classified clusters')
    parser.add_argument('-bed', dest="asatbed",help="bed-formatted decomposed alpha satellite",required=False)
    parser.add_argument('-c', dest="clstid", help="tsv-formatted cluster information", required=False)
    parser.add_argument('-O', dest="outdir", help="out directory", required=False)
    parser.add_argument('-d', dest="distance", help="maximum distance allowed for monomers to be merged.",
                            default=10000,required=False)
    parser.add_argument('-m', dest="method", help="method of HOR mining [HORmon or HiCAT]", required=False, default="stat")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cblocks = CountinuousBlock(args.asatbed, args.distance)
    seqid_clustid = MonomerCluster(args.clstid)
    outfile = os.path.join(args.outdir, "filter.asat.bed")
    OutFilterMn(cblocks, seqid_clustid, outfile)
    if args.method == "HiCAT":
        main_HiCAT(cblocks, seqid_clustid, args.outdir)
        hicat_hor_stat(args.outdir, args.asatbed, args.distance)
    elif args.method == "HORmon":
        #main_graph(cblocks, seqid_clustid, args.outdir)
        #graph_hor_stat(cblocks, args.outdir)
        graph_hordecomposition_stat(cblocks, seqid_clustid, args.outdir)
    else:
        pass


if __name__ == "__main__":
    main()
