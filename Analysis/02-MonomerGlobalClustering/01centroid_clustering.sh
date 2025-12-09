#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g

##step-2 centroid clustering 

mkdir -p step-2 && cd step-2

seq 0 ${clusternum} > clusterid.xls
split -l 2000 clusterid.xls cluster_part_

for file in cluster_part_*; do
    cat << EOF > ${file}.sh
#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH -J ${file}_centroid
#SBATCH -o ${file}.%j.out
#SBATCH -e ${file}.%j.err

date
echo "Processing file: ${file}"
while read i; do 
    python get_pairwise_identity.py "\$i" && echo "\$i finished"
done < "${file}"
date
EOF

    # sbatch ${file}.sh
done

python centroid_cluster.py $threshold $clusternum

##merge cluster results under different thresholds
cp scripts/merge_cluster_out.py cluster_out/
cd cluster_out
python merge_cluster_out.py