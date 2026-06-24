#!/bin/bash
grep H1L {sample}.cenSat_anno.bed | bedtools merge -i - -d 10000 | awk '$3-$2 > 50000' > asat.draft.bed 
sort -k1,1 -k2,2n {sample}.MethyFreq.freezed_dipasm.bed | awk 'BEGIN{OFS="\t"; print "chromosome","start","end","num_motifs_in_group","called_sites","called_sites_methylated","methylated_frequency","group_sequence"}{print $1,$2,$3,"1",$6,$7,$9,"split-group"}' > C002.tsv
awk 'NR > 3' {sample}_repeatmasker.out | sort -k5,5 -k6,6n |grep "ALR/Alpha" > {sample}.rm_clean.out 
snakemake
