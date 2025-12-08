#!/bin/bash
#SBATCH --job-name=builddb
#SBATCH --partition=cpu64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g

makeblastdb -in MonomerT2T.consensus.fa  -dbtype nucl -out MonomerT2T.consensus
makembindex -iformat blastdb -input MonomerT2T.consensus
