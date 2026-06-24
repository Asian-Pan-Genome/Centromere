suf=".fq.gz"
fastp -i ${igG_reads}"_1"${suf} -o igG_1.clean.fastq -I ${igG_reads}"_2"${suf} -O igG_2.clean.fastq -j igG.json -h igG.html -u 30 -q 20 -w 10
fastp -i ${chip_reads}"_1"${suf} -o chip_1.clean.fastq -I ${chip_reads}"_2"${suf} -O chip_2.clean.fastq -j chip.json -h chip.html -u 30 -q 20 -w 10
cat {sample}_Mat.fasta {sample}_Pat.fasta > {sample}_all.fasta  

bwa index {sample}_all.fasta
bwa mem -t 20 -k 50 -c 1000000 {sample}_all.fasta chip_1.clean.fastq chip_2.clean.fastq > chip.all.sam
samtools view -@20 -F 2308 -bh chip.all.sam > chip.all_test.bam
samtools sort chip.all.bam > chip.all.sorted.bam
samtools index chip.all.sorted.bam chip.all.sorted.bam.bai

bwa mem -t 20 -k 50 -c 1000000 {sample}_all.fasta igG_1.clean.fastq igG_2.clean.fastq > igG.all.sam 
samtools view -@20 -F 2308 -bh igG.all.sam > igG.all.bam
samtools sort -@20 igG.all.bam > igG.all.sorted.bam
samtools index igG.all.sorted.bam igG.all.sorted.bam.bai

bamCompare -b1 chip.all.sorted.bam -b2 igG.all.sorted.bam -o log2ratio.bedgraph --binSize 1000 --numberOfProcessors 2 --operation ratio --outFileFormat bedgraph 
