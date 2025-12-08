#!/bin/bash

source ./configure.yaml
set -xoe pipefail


prefix=$1
hap=$2 
if [ ! -d "$prefix/${hap}" ];then mkdir -p $prefix/${hap}; fi
cd $prefix/${hap}

## soft-link ##
genome="$assemblyDir/${prefix}/${hap}/${prefix}_${hap}.fasta"
if [ ! -f "${prefix}_${hap}.fasta.fai" ];then $samtools faidx $genome;fi
ffai="${prefix}_${hap}.fasta.fai"

## alpha satellite ##
echo  "#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40g

if [ ! -f \"${prefix}.AS\" ]; then
    date && blastn -db ${MonomerT2Tdb} \\
        -query $genome -out ${prefix}.AS \\
        -evalue 1e-10 -outfmt 6 -num_threads 10 -task megablast
        date && echo \"{log}: $prefix ${hap} blastn done\"
else
    echo \"{log}: $prefix ${hap} AS file have been exist\"
fi

infile=\"${prefix}.AS\"
tmpfile=\"${prefix}.as.tmp\"
echo \$infile

sort -V -k1,1 -k7,7n -k8,8n -k12,12nr \$infile | awk '
{
    key = \$1 \"\t\" \$7 \"\t\" \$8
    if (key != prev_key) {
        if (NR > 1) {
            if (pident >= 90 && alignlen >= 100) {
                if (qstart < qend && sstart < send) {
                    strand = \"+\"
                } else {
                    strand = \"-\"
                }
                print best_line \"\t\" strand
            }
        }
        prev_key = key
        bitscore = \$12
        qstart = \$7
        qend = \$8
        sstart = \$9
        send = \$10
        pident = \$3
        alignlen = \$4
        best_line = \$0
    } else {
        if (\$12 > bitscore) {
            bitscore = \$12
            qstart = \$7
            qend = \$8
            sstart = \$9
            send = \$10
            pident = \$3
            alignlen = \$4
            best_line = \$0
        }
    }
}
END {
    if (pident >= 90 && alignlen >= 100) {
        if (qstart < qend && sstart < send) {
            strand = \"+\"
        } else {
            strand = \"-\"
        }
        print best_line \"\t\" strand
    }
}
' | awk '{print \$1\"\t\"\$7\"\t\"\$8\"\t\"\$2\"\t\"\$12\"\t\"\$13}' > \$tmpfile 

intersect_tmpfile=\"${prefix}.as.intersect.tmp\"
unique_tmpfile=\"${prefix}.as.unique.tmp\"
overlap_id=\"${prefix}.as.dup.ids\"
overlap_tmpfile=\"${prefix}.as.dup.tmp\"
overlap_uniqfile=\"${prefix}.as.dedup.tmp\"
outfile=\"${prefix}_${hap}.asat.bed\"

bedtools intersect -a \$tmpfile -b \$tmpfile -wa -wb | awk '{print \$1\"@\"\$2\"@\"\$3\"@\"\$4\"@\"\$5\"@\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$9\"\t\"\$10\"\t\"\$11\"\t\"\$12}' > \$intersect_tmpfile
cat \$intersect_tmpfile | datamash -g 1 count 2 | awk '{if (\$2==1)print\$1}' | sed \"s/@/\t/g\" > \$unique_tmpfile
cat \$intersect_tmpfile | datamash -g 1 count 2 | awk '{if (\$2>1)print\$1}' > \$overlap_id
grep -F -f \$overlap_id \$intersect_tmpfile > \$overlap_tmpfile
awk '{
    if (!(\$1 in max) || \$6 > max[\$1]) {
        max[\$1] = \$6;
        line[\$1] = \$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7;
    }
} END {
    for (i in line) {
        print line[i];
    }
}' \$overlap_tmpfile | sort | uniq > \$overlap_uniqfile
cat \$unique_tmpfile \$overlap_uniqfile | sort -V -k1,1 -k2,2n -k3,3n -k12,12nr > \$outfile
rm \$intersect_tmpfile \$unique_tmpfile \$overlap_id \$overlap_tmpfile \$overlap_uniqfile

### choose gap:5 to estimate blurred peri-centromere region (round0) ###  
$bedtools merge -i \$outfile -d 5 > ${prefix}.g5.bed
datamash -g 1 min 2 max 3 < ${prefix}.g5.bed  > ${prefix}.g5.alpha.cenpos
$bedtools slop -i ${prefix}.g5.alpha.cenpos -g $ffai -b 5000000 > ${prefix}.g5.raw.cenpos
cut -f1 ${prefix}.g5.alpha.cenpos | grep \"chr\" | grep -v \"chrM\" | sort -t'#' -k3,3V > ${prefix}.g5.alpha.chrom
grep \"chr\" $ffai | grep -v \"chrM\" |cut -f1 |  sort -t'#' -k3,3V > ${prefix}.chrom
diff ${prefix}.g5.alpha.chrom ${prefix}.chrom > cenpos.check
if [ -s \"cenpos.check\" ];then
    cat cenpos.check | xargs echo
fi
mv ${prefix}.g5.raw.cenpos ${prefix}.round0.cenpos
echo \"blastn2bed done!\"
rm ${prefix}.g5.bed ${prefix}.g5.alpha.cenpos

date
" >0alpha.sh
#sbatch -J AS_${prefix}_${hap} -o asat.log -e asat.error 0alpha.sh 

#### Annotation_of_asat_based_on_MonomerGlobalDB ####
## alpha satellite ##
echo  "#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40g


if [ ! -f \"${prefix}.global.AS\" ]; then
    date && blastn -db ${MonomerGlobaldb} \\
        -query $genome -out ${prefix}.global.AS \\
        -evalue 1e-10 -outfmt 6 -num_threads 10 -task megablast
        date && echo \"{log}: $prefix ${hap} blastn done\"
else
    echo \"{log}: $prefix ${hap} AS file have been exist\"
fi

infile=\"${prefix}.global.AS\"
tmpfile=\"${prefix}.as.tmp\"
echo \$infile

sort -V -k1,1 -k7,7n -k8,8n -k12,12nr \$infile | awk '
{
    key = \$1 \"\t\" \$7 \"\t\" \$8
    if (key != prev_key) {
        if (NR > 1) {
            if (pident >= 90 && alignlen >= 100) {
                if (qstart < qend && sstart < send) {
                    strand = \"+\"
                } else {
                    strand = \"-\"
                }
                print best_line \"\t\" strand
            }
        }
        prev_key = key
        bitscore = \$12
        qstart = \$7
        qend = \$8
        sstart = \$9
        send = \$10
        pident = \$3
        alignlen = \$4
        best_line = \$0
    } else {
        if (\$12 > bitscore) {
            bitscore = \$12
            qstart = \$7
            qend = \$8
            sstart = \$9
            send = \$10
            pident = \$3
            alignlen = \$4
            best_line = \$0
        }
    }
}
END {
    if (pident >= 90 && alignlen >= 100) {
        if (qstart < qend && sstart < send) {
            strand = \"+\"
        } else {
            strand = \"-\"
        }
        print best_line \"\t\" strand
    }
}' | awk '{
    if (match(\$2, /^([0-9]+)_/, arr)) {
        extracted_num = arr[1]
    }else{
        extracted_num = \$2
    }
    print \$1\"\t\"\$7\"\t\"\$8\"\t\"extracted_num\"\t\"\$12\"\t\"\$13
}'> \$tmpfile

intersect_tmpfile=\"${prefix}.as.intersect.tmp\"
unique_tmpfile=\"${prefix}.as.unique.tmp\"
overlap_id=\"${prefix}.as.dup.ids\"
overlap_tmpfile=\"${prefix}.as.dup.tmp\"
overlap_uniqfile=\"${prefix}.as.dedup.tmp\"
outfile=\"${prefix}_${hap}.asat.global.bed\"

bedtools intersect -a \$tmpfile -b \$tmpfile -wa -wb | awk '{print \$1\"@\"\$2\"@\"\$3\"@\"\$4\"@\"\$5\"@\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$9\"\t\"\$10\"\t\"\$11\"\t\"\$12}' > \$intersect_tmpfile
cat \$intersect_tmpfile | datamash -g 1 count 2 | awk '{if (\$2==1)print\$1}' | sed \"s/@/\t/g\" > \$unique_tmpfile
cat \$intersect_tmpfile | datamash -g 1 count 2 | awk '{if (\$2>1)print\$1}' > \$overlap_id
grep -F -f \$overlap_id \$intersect_tmpfile > \$overlap_tmpfile
awk '{
    if (!(\$1 in max) || \$6 > max[\$1]) {
        max[\$1] = \$6;
        line[\$1] = \$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7;
    }
} END {
    for (i in line) {
        print line[i];
    }
}' \$overlap_tmpfile | sort | uniq > \$overlap_uniqfile
cat \$unique_tmpfile \$overlap_uniqfile | sort -V -k1,1 -k2,2n -k3,3n -k12,12nr > \$outfile
rm \$intersect_tmpfile \$unique_tmpfile \$overlap_id \$overlap_tmpfile \$overlap_uniqfile

### choose gap:5 to estimate blurred peri-centromere region (round0) ###  
$bedtools merge -i \$outfile -d 5 > ${prefix}.g5.global.bed
datamash -g 1 min 2 max 3 < ${prefix}.g5.global.bed  > ${prefix}.g5.global.alpha.cenpos
$bedtools slop -i ${prefix}.g5.global.alpha.cenpos -g $ffai -b 5000000 > ${prefix}.g5.global.raw.cenpos
cut -f1 ${prefix}.g5.global.alpha.cenpos | grep \"chr\" | grep -v \"chrM\" | sort -t'#' -k3,3V > ${prefix}.g5.global.alpha.chrom
grep \"chr\" $ffai | grep -v \"chrM\" |cut -f1 |  sort -t'#' -k3,3V > ${prefix}.chrom
diff ${prefix}.g5.global.alpha.chrom ${prefix}.chrom > cenpos.check
if [ -s \"cenpos.check\" ];then
    cat cenpos.check | xargs echo
fi
mv ${prefix}.g5.global.raw.cenpos ${prefix}.round0.global.cenpos
echo \"blastn2bed done!\"
rm ${prefix}.g5.bed ${prefix}.g5.alpha.cenpos
date
" >0alphaglobal.sh

echo "#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g

cp $hsat23pl $hsat2db $hsat3db .
perl Assembly_HSat2and3_v3.pl $genome .
" > 0hsat23.sh
#sbatch -J HSAT23_${prefix}_${hap} -o hsat23.log -e hsat23.error 0hsat23.sh

echo "#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=150g

#source /share/home/zhanglab/user/sunyanqing/miniconda3/bin/activate repeat
#if [ ! -d \"repeat/01_out_repeatmasker\" ];then mkdir -p repeat/01_out_repeatmasker; fi

#res=\$(ls -A \"repeat/01_out_repeatmasker\")
#if [ -z \"\$res\" ];then
#    $repeatmasker -species human -e rmblast  -s -pa 30 $genome  -html -gff -dir repeat/01_out_repeatmasker
#    echo \"{log}: $prefix ${hap} repeatmasker done\"
#else
#    echo \"{log}: $prefix ${hap} repeatmasker have already been done.\"
#conda deactivate
#fi

### hsat1 ###
#grep \"Motif:SAR\\\"\" repeat/01_out_repeatmasker/*fasta.out.gff | awk '{print\$1\"\t\"\$4-1\"\t\"\$5\"\t\"\$7\"\tHSat1A\"}' | sort -k1,1 -V -k2,2n > repeat/01_out_repeatmasker/hsat1a.raw.bed
#grep \"Motif:HSATI\\\"\" repeat/01_out_repeatmasker/*fasta.out.gff | awk '{print\$1\"\t\"\$4-1\"\t\"\$5\"\t\"\$7\"\tHSat1B\"}' | sort -k1,1 -V -k2,2n  > repeat/01_out_repeatmasker/hsat1b.raw.bed
#python $hsat1py -i repeat/01_out_repeatmasker/hsat1a.raw.bed -o repeat/01_out_repeatmasker/hsat1a.bed
#python $hsat1py -i repeat/01_out_repeatmasker/hsat1b.raw.bed -o repeat/01_out_repeatmasker/hsat1b.bed
#cat repeat/01_out_repeatmasker/hsat1a.bed repeat/01_out_repeatmasker/hsat1b.bed | sort -k1,1 -V -k2,2n > ${prefix}.${hap}.HSat1.bed

### bsat ###
egrep \"BSAT|Beta|LSAU\" repeat/01_out_repeatmasker/*fasta.out.gff | awk '{print\$1\"\t\"\$4-1\"\t\"\$5\"\tBSat\t0\t\"\$7\"\t\"\$4-1\"\t\"\$5\"\t#EEA6B7\"}' | sort -k1,1 -V -k2,2n > repeat/01_out_repeatmasker/bsat.bed
cp repeat/01_out_repeatmasker/bsat.bed ${prefix}.${hap}.bsat.bed

### gsat ###
egrep \"GSAT|TAR1\" repeat/01_out_repeatmasker/*fasta.out.gff | awk '{print\$1\"\t\"\$4-1\"\t\"\$5\"\tGSat\t0\t\"\$7\"\t\"\$4-1\"\t\"\$5\"\t#5E665B\"}' | sort -k1,1 -V -k2,2n  > repeat/01_out_repeatmasker/gsat.bed
cp repeat/01_out_repeatmasker/gsat.bed ${prefix}.${hap}.gsat.bed
" > 0repeat.sh
#sbatch -J RMA_${prefix}_${hap} -o repeat.log -e repeat.error 0repeat.sh
#sh 0repeat.sh

echo "#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g

asat=false; hsat1=false; hsat23=false; bsat=false; gsat=false

if [ ! -f \"${prefix}_${hap}.asat.bed\" ]; then
    echo \"${prefix}_${hap}.asat.bed could not find!\"
else
    asat=true
fi

if [ ! -f \"${prefix}.${hap}.HSat1.bed\" ]; then
    echo \"${prefix}.${hap}.HSat1.bed could not find!\"
else
    hsat1=true
fi

if ls \${prefix}*.HSat2and3_Regions.bed 1> /dev/null 2>&1; then
    hsat23=true
else
    echo \"\${prefix}*.HSat2and3_Regions.bed could not find!\"
fi

if [ ! -f \"${prefix}.${hap}.bsat.bed\" ]; then
    echo \"${prefix}.${hap}.bsat.bed could not find!\"
else
    bsat=true
fi

if [ ! -f \"${prefix}.${hap}.gsat.bed\" ]; then
    echo \"${prefix}.${hap}.gsat.bed could not find!\"
else
    gsat=true
fi

if [[ \"\$asat\" = true && \"\$hsat1\" = true && \"\$hsat23\" = true && \"\$bsat\" = true && \"\$gsat\" = true ]] ; then
    awk '{if (\$2 < \$3) print}' ${prefix}*.HSat2and3_Regions.bed > ${prefix}.${hap}.HSat2and3_Regions.tmp
    awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$2\"\t\"\$3\"\t#C91D32\"}' ${prefix}_${hap}.asat.bed > ${prefix}_${hap}.asat.tmp
    cat ${prefix}_${hap}.asat.tmp ${prefix}.${hap}.HSat1.bed ${prefix}.${hap}.HSat2and3_Regions.tmp ${prefix}.${hap}.bsat.bed ${prefix}.${hap}.gsat.bed| sort -k1,1 -V -k2,2n | sed \"s/120,161,187/#00B0F0/g\" | sed \"s/51,51,102/#2D529F/g\" > ${prefix}.round0.cenanno
fi
rm ${prefix}_${hap}.asat.tmp
    
if [ ! -d \"barplot/round1\" ];then mkdir -p barplot/round1; fi

if [ -f \"${prefix}.round0.cenanno\" ] && [ -f \"${prefix}.round0.cenpos\" ];then
    #python $barplotpy -anno ${prefix}.round0.cenanno \\
    #                 -pos ${prefix}.round0.cenpos \\
    #                 -fai $ffai \\
    #                 -o ${prefix}.round1.cenpos \\
    #                 -p ${prefix}.${hap} \\
    #                 -plotOutdir barplot/round1
    $bedtools intersect -a ${prefix}.round0.cenanno -b ${prefix}.round1.cenpos -wa > ${prefix}_${hap}.cenanno.bed
else
    echo \"could not find inputfiles!\"
fi
" > 0barplot.sh
#sbatch -J barplot_${prefix}_${hap} -o barplot.log -e barplot.error 0barplot.sh
#sh 0barplot.sh
