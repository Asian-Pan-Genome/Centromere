#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g

if [ $# -eq 0 ]; then
    echo "Error: globalflag parameter (True or False) is required"
    echo "Usage: $0 <globalflag>"
    exit 1
fi

globalflag=$1

# Validate parameter value
if [ "$globalflag" != "True" ] && [ "$globalflag" != "False" ]; then
    echo "Error: globalflag must be 'True' or 'False'"
    exit 1
fi


asat=false; hsat1=false; hsat23=false; bsat=false; gsat=false;

if [ "$globalflag" == "False" ]; then
    if [ ! -f "C001-CHA-E01-01_Mat.asat.bed" ]; then
        echo "C001-CHA-E01-01_Mat.asat.bed could not find!"
    else
        asat=true
        asatbed="C001-CHA-E01-01_Mat.asat.bed"
    fi
else
    echo "merge asat annotation based on MonomerGlobalDB!"
    if [ ! -f "C001-CHA-E01-01_Mat.asat.global.bed" ]; then
        echo "C001-CHA-E01-01_Mat.asat.global.bed could not find!"
    else
        asat=true
        asatbed="C001-CHA-E01-01_Mat.asat.global.bed"
    fi
fi



if [ ! -f "C001-CHA-E01-01.Mat.HSat1.bed" ]; then
    echo "C001-CHA-E01-01.Mat.HSat1.bed could not find!"
else
    hsat1=true
fi

if ls ${prefix}*.HSat2and3_Regions.bed 1> /dev/null 2>&1; then
    hsat23=true
else
    echo "${prefix}*.HSat2and3_Regions.bed could not find!"
fi

if [ ! -f "C001-CHA-E01-01.Mat.bsat.bed" ]; then
    echo "C001-CHA-E01-01.Mat.bsat.bed could not find!"
else
    bsat=true
fi

if [ ! -f "C001-CHA-E01-01.Mat.gsat.bed" ]; then
    echo "C001-CHA-E01-01.Mat.gsat.bed could not find!"
else
    gsat=true
fi

if [[ "$asat" = true && "$hsat1" = true && "$hsat23" = true && "$bsat" = true && "$gsat" = true ]] ; then
    awk '{if ($2 < $3) print}' C001-CHA-E01-01*.HSat2and3_Regions.bed > C001-CHA-E01-01.Mat.HSat2and3_Regions.tmp
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2"\t"$3"\t#C91D32"}' $asatbed > C001-CHA-E01-01_Mat.asat.tmp
    cat C001-CHA-E01-01_Mat.asat.tmp C001-CHA-E01-01.Mat.HSat1.bed C001-CHA-E01-01.Mat.HSat2and3_Regions.tmp C001-CHA-E01-01.Mat.bsat.bed C001-CHA-E01-01.Mat.gsat.bed| sort -k1,1 -V -k2,2n | sed "s/120,161,187/#00B0F0/g" | sed "s/51,51,102/#2D529F/g" > C001-CHA-E01-01.round0.cenanno
fi
rm C001-CHA-E01-01_Mat.asat.tmp
    
if [ ! -d "barplot/round1" ];then mkdir -p barplot/round1; fi

if [ -f "C001-CHA-E01-01.round0.cenanno" ] && [ -f "C001-CHA-E01-01.round0.cenpos" ];then
    python  -anno C001-CHA-E01-01.round0.cenanno \
                     -pos C001-CHA-E01-01.round0.cenpos \
                     -fai C001-CHA-E01-01_Mat.v0.9.fasta.fai \
                     -o C001-CHA-E01-01.round1.cenpos \
                     -p C001-CHA-E01-01.Mat \
                     -plotOutdir barplot/round1
    bedtools intersect -a C001-CHA-E01-01.round0.cenanno -b C001-CHA-E01-01.round1.cenpos -wa > C001-CHA-E01-01_Mat.cenanno.bed
else
    echo "could not find inputfiles!"
fi
python  C001-CHA-E01-01_Mat.cenanno.bed C001-CHA-E01-01_Mat.merged.cenanno.bed
bedtools intersect -a C001-CHA-E01-01_Mat.merged.cenanno.bed -b C001-CHA-E01-01.round1.cenpos -f 1.0 | sort -V -k1,1 -k2,2n | datamash -g 1 min 2 max 3 > C001-CHA-E01-01_Mat.cent.bed

