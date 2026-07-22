#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80g



#python select_marker_seq.py
#mafft --auto marker_seq.fa > marker_seq.aln.fa
#sed "s/:/@/g"  marker_seq.aln.fa > marker_seq.aln.mod.fa
/share/home/zhanglab/user/sunyanqing/software/FastTree -nt marker_seq.aln.mod.fa > marker_seq.aln.mod.newick
/share/home/zhanglab/user/sunyanqing/software/iqtree-1.6.12-Linux/bin/iqtree -s marker_seq.aln.mod.fa -m MFP -bb 1000 -nt AUTO
 
