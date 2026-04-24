#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g

cp    .
perl Assembly_HSat2and3_v3.pl C001-CHA-E01-01/Mat/C001-CHA-E01-01_Mat.v0.9.fasta .

