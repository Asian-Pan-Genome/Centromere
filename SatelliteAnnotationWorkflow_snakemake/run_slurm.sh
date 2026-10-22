#!/usr/bin/env bash
#SBATCH --job-name=satellite_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=satellite.%j.log
# Add your server's partition/account directives above.
set -euo pipefail

# Activate SatelliteAnnotationEnv BEFORE sbatch, or insert the appropriate
# conda.sh source and conda activate commands for your own server here.
# Submit from the workflow root. SLURM exports the activated environment by default.
cd "${SLURM_SUBMIT_DIR:?Submit this script from the workflow directory}"
command -v snakemake
command -v RepeatMasker
snakemake --cores "${SLURM_CPUS_PER_TASK}" --resources mem_mb=64000 --rerun-incomplete
