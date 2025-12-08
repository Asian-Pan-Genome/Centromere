#!/bin/bash
#SBATCH --job-name=builddb
#SBATCH --partition=cpu64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g

makeblastdb -in MonomerGlobalDB.fasta  -dbtype nucl -out MonomerGlobalDB
makembindex -iformat blastdb -input MonomerGlobalDB
