#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g

mkdir -p step-3 && cd step-3
ln -s Louvain_87-99_concatenated_clusters_sorted.txt .

cat col_threshold.txt | while read line
do
    arr=($line)
    col=${arr[0]}
    threshold=${arr[1]}

    tail -n +2 Louvain_87-99_concatenated_clusters_sorted.txt | sort -k${col},${col}n |  datamash -g $col count 1 collapse 1 > ${threshold}.sorted.txt
    python vs_chm13.py ${threshold}.sorted.txt ${threshold}.sorted.out
    sort -k2,2 -k3,3n ${threshold}.sorted.out | awk '{print $2"@"$3}' | datamash -g 1 count 1 | sed "s/@/\t/g" > ${threshold}.sorted.out.stat
done