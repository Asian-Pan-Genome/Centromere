#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=150g

if [ ! -d "repeat/01_out_repeatmasker" ];then mkdir -p repeat/01_out_repeatmasker; fi

res=$(ls -A "repeat/01_out_repeatmasker")
if [ -z "$res" ];then
    RepeatMasker -species human -e rmblast  -s -pa 30 C001-CHA-E01-01/Mat/C001-CHA-E01-01_Mat.v0.9.fasta  -html -gff -dir repeat/01_out_repeatmasker
    echo "{log}: C001-CHA-E01-01 Mat repeatmasker done"
else
    echo "{log}: C001-CHA-E01-01 Mat repeatmasker have already been done."
conda deactivate
fi

### hsat1 ###
grep "Motif:SAR\"" repeat/01_out_repeatmasker/*fasta.out.gff | awk '{print$1"\t"$4-1"\t"$5"\t"$7"\tHSat1A"}' | sort -k1,1 -V -k2,2n > repeat/01_out_repeatmasker/hsat1a.raw.bed
grep "Motif:HSATI\"" repeat/01_out_repeatmasker/*fasta.out.gff | awk '{print$1"\t"$4-1"\t"$5"\t"$7"\tHSat1B"}' | sort -k1,1 -V -k2,2n  > repeat/01_out_repeatmasker/hsat1b.raw.bed
python  -i repeat/01_out_repeatmasker/hsat1a.raw.bed -o repeat/01_out_repeatmasker/hsat1a.bed
python  -i repeat/01_out_repeatmasker/hsat1b.raw.bed -o repeat/01_out_repeatmasker/hsat1b.bed
cat repeat/01_out_repeatmasker/hsat1a.bed repeat/01_out_repeatmasker/hsat1b.bed | sort -k1,1 -V -k2,2n > C001-CHA-E01-01.Mat.HSat1.bed

### bsat ###
egrep "BSAT|Beta|LSAU" repeat/01_out_repeatmasker/*fasta.out.gff | awk '{print$1"\t"$4-1"\t"$5"\tBSat\t0\t"$7"\t"$4-1"\t"$5"\t#EEA6B7"}' | sort -k1,1 -V -k2,2n > repeat/01_out_repeatmasker/bsat.bed
cp repeat/01_out_repeatmasker/bsat.bed C001-CHA-E01-01.Mat.bsat.bed

### gsat ###
egrep "GSAT|TAR1" repeat/01_out_repeatmasker/*fasta.out.gff | awk '{print$1"\t"$4-1"\t"$5"\tGSat\t0\t"$7"\t"$4-1"\t"$5"\t#5E665B"}' | sort -k1,1 -V -k2,2n  > repeat/01_out_repeatmasker/gsat.bed
cp repeat/01_out_repeatmasker/gsat.bed C001-CHA-E01-01.Mat.gsat.bed

