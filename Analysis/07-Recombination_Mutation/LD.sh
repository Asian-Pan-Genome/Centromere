#!/bin/bash
#SBATCH --job-name=SNP
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g


bcftools="/share/home/zhanglab/user/sunyanqing/miniconda3/envs/viz/bin/bcftools"
plink="/share/home/zhanglab/user/sunyanqing/software/plink/plink"

ichr="chr1"
#start_prefix=$1
#end_prefix=$2

####caculate LD for two-side pericentromere####

#$bcftools view -R ${ichr}_customized_region.bed ../CHM13-APGp1-HPRCp1-HGSVCp3_snp_diploid.vcf.gz -Oz -o ${ichr}_customized.vcf.gz
# tabix -p vcf ${ichr}_customized.vcf.gz

# $plink --vcf ${ichr}_customized.vcf.gz \
#        --double-id \
#        --set-missing-var-ids @:# \
#        --geno 0.1 \
#        --maf 0.02 \
#        --snps-only \
#        --make-bed \
#        --out ${ichr}_customized

#$plink --bfile ${ichr}_customized --ld-window-r2 0 --ld-window-kb 100000 --ld-window 99999999 --r2 --out ${ichr}_customized
#cat ${ichr}_customized.bim | cut -f2 | sort -R | tail -n 10000 > random.site
#$plink --bfile ${ichr}_customized --extract random.site --make-bed --out random
#$plink --bfile random --ld-window-r2 0 --ld-window-kb 10000  --ld-window 9999999 --r2 --out random
#cp ../chr8/02.plot.LD.R .
#~/miniconda3/envs/R/bin/Rscript 02.plot.LD.R random.ld ${ichr}_customized_random_abpos.png
# ~/miniconda3/envs/R/bin/Rscript 02.plot.LD.R ${ichr}_customized.ld ${ichr}_customized_abpos.png


####caculate LD for two-side pericentromere but exclude centromere####
# echo -e "chr1\t120519169\t121619172\nchr1\t142242033\t143242033" > chr1_customized_region_exclude_cent.bed
# $bcftools view -R ${ichr}_customized_region_exclude_cent.bed ../CHM13-APGp1-HPRCp1-HGSVCp3_snp_diploid.vcf.gz -Oz -o ${ichr}_customized_exclude_cent.vcf.gz
# tabix -p vcf ${ichr}_customized_exclude_cent.vcf.gz

# $plink --vcf ${ichr}_customized_exclude_cent.vcf.gz \
#       --double-id \
#       --set-missing-var-ids @:# \
#       --geno 0.1 \
#       --maf 0.02 \
#       --snps-only \
#       --make-bed \
#       --out ${ichr}_customized_exclude_cent

# $plink --bfile ${ichr}_customized_exclude_cent --ld-window-r2 0 --ld-window-kb 100000 --ld-window 99999999 --r2 --out ${ichr}_customized_exclude_cent
# cat ${ichr}_customized_exclude_cent.bim | cut -f2 | sort -R | tail -n 10000 > random_exclude_cent.site
# $plink --bfile ${ichr}_customized_exclude_cent --extract random_exclude_cent.site --make-bed --out random_exclude_cent
# $plink --bfile random_exclude_cent --ld-window-r2 0 --ld-window-kb 10000  --ld-window 9999999 --r2 --out random_exclude_cent
#cp ../chr5/02.plot.LD.relativepos.R .
#~/miniconda3/envs/R/bin/Rscript 02.plot.LD.R random_exclude_cent.ld ${ichr}_customized_random_abpos_exclude_cent.png
# ~/miniconda3/envs/R/bin/Rscript 02.plot.LD.relativepos.R  random_exclude_cent.ld ${ichr}_customized_random_abpos_exclude_cent.png ${ichr}_customized_random_abpos_index.csv 49115450
~/miniconda3/envs/R/bin/Rscript 02.plot.LD.relativepos.R  ${ichr}_customized_exclude_cent.ld ${ichr}_customized_exclude_cent.png ${ichr}_customized_exclude_cent.csv 121619172







