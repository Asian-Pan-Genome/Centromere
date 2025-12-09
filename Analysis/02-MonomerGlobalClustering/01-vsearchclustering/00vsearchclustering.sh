#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100g


##filter non-canonical length
seqkit seq -m 120 -M 200 -g  all.mn.fasta -o all.mn.filtered.fasta
seqkit seq -M 119 -g all.mn.fasta -o all.mn.short.fasta
seqkit seq -m 201 -g all.mn.fasta -o all.mn.long.fasta

mnfasta="all.mn.filtered.fasta"

threshold=$1


##step-1 original clustering using vsearch

mkdir -p step-1 && cd step-1

vsearch --cluster_size $mnfasta --id 0.${threshold} --iddef 0 --centroids all.centroid.fasta  --uc all.uc --clusters original_seq/
samtools faidx all.centroid.fasta
clusternum=$(cat all.centroid.fasta.fai | wc -l)

# for i in {0..${clusternum}}
# do
#    ~/miniconda3/envs/viz/bin/samtools faidx original_seq/${i}
#    n=$(cat original_seq/${i}.fai | wc -l)
#    echo -e $i "\t" $n >> centroid.size
#    echo $i "done"
# done
