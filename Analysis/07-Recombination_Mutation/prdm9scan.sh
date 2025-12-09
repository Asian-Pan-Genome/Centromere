#!/bin/bash
#SBATCH --job-name=PRDM9
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g

sample=$1
hap=$2

#if [[ $sample == HG* || $sample == NA* ]]; then
#    if [[ $hap == "hap1" || $hap == "hap2" ]]; then
#        echo "$sample $hap HGSVC"
#        fasta="/share/home/zhanglab/user/sunyanqing/human/HGSVC/$sample/$hap/${sample}_${hap}_chrR.fasta"
#    elif [[ $hap == "Mat" || $hap == "Pat" ]]; then
#        echo "$sample $hap HPRC"
#        fasta="/share/home/zhanglab/user/sunyanqing/human/HPRC/$sample/$hap/${sample}_${hap}_chrR.fasta"
#    else
#        echo "$sample $hap"
#    fi
#else
#    echo "$sample $hap APG"
#    fasta="/share/home/zhanglab/user/sunyanqing/human/assembly/${sample}/${hap}/${sample}_${hap}.v0.9.fasta"
#fi
#
#echo $fasta
#if [ ! -d ${sample}_${hap} ]; then mkdir ${sample}_${hap};fi
#sh fimo_scan.sh PRDM9_motifs.human.txt $fasta ${sample}_${hap}
#awk -F '\t' '{if ($9<0.3)print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$6"\t"$8"\t"$9"\t"$10}' ${sample}_${hap}/fimo.tsv | head -n -4 > ${sample}_${hap}/fimo.filtered.bed
#
#
#bedtools makewindows -g ${fasta}.fai -w 20000 > ${sample}_${hap}/genome.20kb.bed
#bedtools intersect -a ${sample}_${hap}/genome.20kb.bed -b ${sample}_${hap}/fimo.filtered.bed -wo  > ${sample}_${hap}/${sample}_${hap}.20kb.PRDM9.bed
#/share/home/zhanglab/user/sunyanqing/miniconda3/envs/R/bin/Rscript /share/home/zhanglab/user/sunyanqing/vscode_scripts/periphy/plot_PRDM9_density.R $sample $hap
/share/home/zhanglab/user/sunyanqing/miniconda3/envs/R/bin/Rscript /share/home/zhanglab/user/sunyanqing/vscode_scripts/periphy/plot_PRDM9_density_extcent.R $sample $hap
