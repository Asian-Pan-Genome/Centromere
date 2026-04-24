#!/bin/bash
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40g

if [ ! -f "C001-CHA-E01-01.AS" ]; then
    date && blastn -db ./database/MonomerT2T.consensus \
        -query C001-CHA-E01-01/Mat/C001-CHA-E01-01_Mat.v0.9.fasta -out C001-CHA-E01-01.AS \
        -evalue 1e-10 -outfmt 6 -num_threads 10 -task megablast
        date && echo "{log}: C001-CHA-E01-01 Mat blastn done"
else
    echo "{log}: C001-CHA-E01-01 Mat AS file have been exist"
fi

infile="C001-CHA-E01-01.AS"
tmpfile="C001-CHA-E01-01.as.tmp"
echo $infile

sort -V -k1,1 -k7,7n -k8,8n -k12,12nr $infile | awk '
{
    key = $1 "\t" $7 "\t" $8
    if (key != prev_key) {
        if (NR > 1) {
            if (pident >= 90 && alignlen >= 100) {
                if (qstart < qend && sstart < send) {
                    strand = "+"
                } else {
                    strand = "-"
                }
                print best_line "\t" strand
            }
        }
        prev_key = key
        bitscore = $12
        qstart = $7
        qend = $8
        sstart = $9
        send = $10
        pident = $3
        alignlen = $4
        best_line = $0
    } else {
        if ($12 > bitscore) {
            bitscore = $12
            qstart = $7
            qend = $8
            sstart = $9
            send = $10
            pident = $3
            alignlen = $4
            best_line = $0
        }
    }
}
END {
    if (pident >= 90 && alignlen >= 100) {
        if (qstart < qend && sstart < send) {
            strand = "+"
        } else {
            strand = "-"
        }
        print best_line "\t" strand
    }
}
' | awk '{print $1"\t"$7"\t"$8"\t"$2"\t"$12"\t"$13}' > $tmpfile 

intersect_tmpfile="C001-CHA-E01-01.as.intersect.tmp"
unique_tmpfile="C001-CHA-E01-01.as.unique.tmp"
overlap_id="C001-CHA-E01-01.as.dup.ids"
overlap_tmpfile="C001-CHA-E01-01.as.dup.tmp"
overlap_uniqfile="C001-CHA-E01-01.as.dedup.tmp"
outfile="C001-CHA-E01-01_Mat.asat.bed"

bedtools intersect -a $tmpfile -b $tmpfile -wa -wb | awk '{print $1"@"$2"@"$3"@"$4"@"$5"@"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' > $intersect_tmpfile
cat $intersect_tmpfile | datamash -g 1 count 2 | awk '{if ($2==1)print$1}' | sed "s/@/\t/g" > $unique_tmpfile
cat $intersect_tmpfile | datamash -g 1 count 2 | awk '{if ($2>1)print$1}' > $overlap_id
grep -F -f $overlap_id $intersect_tmpfile > $overlap_tmpfile
awk '{
    if (!($1 in max) || $6 > max[$1]) {
        max[$1] = $6;
        line[$1] = $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;
    }
} END {
    for (i in line) {
        print line[i];
    }
}' $overlap_tmpfile | sort | uniq > $overlap_uniqfile
cat $unique_tmpfile $overlap_uniqfile | sort -V -k1,1 -k2,2n -k3,3n -k12,12nr > $outfile
rm $intersect_tmpfile $unique_tmpfile $overlap_id $overlap_tmpfile $overlap_uniqfile

### choose gap:5 to estimate blurred peri-centromere region (round0) ###  
bedtools merge -i $outfile -d 5 > C001-CHA-E01-01.g5.bed
datamash -g 1 min 2 max 3 < C001-CHA-E01-01.g5.bed  > C001-CHA-E01-01.g5.alpha.cenpos
bedtools slop -i C001-CHA-E01-01.g5.alpha.cenpos -g C001-CHA-E01-01_Mat.v0.9.fasta.fai -b 5000000 > C001-CHA-E01-01.g5.raw.cenpos
cut -f1 C001-CHA-E01-01.g5.alpha.cenpos | grep "chr" | grep -v "chrM" | sort -t'#' -k3,3V > C001-CHA-E01-01.g5.alpha.chrom
grep "chr" C001-CHA-E01-01_Mat.v0.9.fasta.fai | grep -v "chrM" |cut -f1 |  sort -t'#' -k3,3V > C001-CHA-E01-01.chrom
diff C001-CHA-E01-01.g5.alpha.chrom C001-CHA-E01-01.chrom > cenpos.check
if [ -s "cenpos.check" ];then
    cat cenpos.check | xargs echo
fi
mv C001-CHA-E01-01.g5.raw.cenpos C001-CHA-E01-01.round0.cenpos
echo "blastn2bed done!"
rm C001-CHA-E01-01.g5.bed C001-CHA-E01-01.g5.alpha.cenpos

date

