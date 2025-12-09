#!/bin/bash
#SBATCH --job-name=SNP
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g

input_region=$1
chr=$2
# bedtools makewindows -b $input_region -w 100000 > 100kb_windows.bed

# bcftools="/share/home/zhanglab/user/sunyanqing/miniconda3/envs/viz/bin/bcftools"
# plink="/share/home/zhanglab/user/sunyanqing/software/plink/plink"
# vcf_file="../../CHM13-APGp1-HPRCp1-HGSVCp3_snp.vcf.gz"

# tmp_dir="tmp"
# if [ ! -d "$tmp_dir" ]; then
#     mkdir "$tmp_dir"
# fi
# line_number=1

# output_bed="100kb_windows_with_snp_counts.bed"
# cp 100kb_windows.bed $output_bed

# # Read the input BED file line by line
# while IFS= read -r line; do
#     # Create a temporary BED file for the current line
#     echo "$line" > "$tmp_dir/tmp.bed"
    
#     # Extract variants using bcftools
#     $bcftools view -R "$tmp_dir/tmp.bed" "$vcf_file" -Oz -o "$tmp_dir/${line_number}.vcf.gz"
    
#     # Filter the VCF file using plink
#     $plink --vcf "$tmp_dir/${line_number}.vcf.gz" \
#           --set-missing-var-ids @:# \
#           --geno 0.2 \
#           --maf 0.02 \
#           --snps-only \
#           --recode vcf \
#           --out "$tmp_dir/${line_number}.filter"
    
#     # Count the number of SNPs in the filtered VCF file
#     if [ -f "$tmp_dir/${line_number}.filter.vcf" ]; then
#         snp_count=$(grep -v "^#" "$tmp_dir/${line_number}.filter.vcf" | wc -l)
#     else
#         snp_count=0
#     fi
    
#     # Append the SNP count to the corresponding line in the output BED file
#     # sed -i "${line_number}s/$/\t${snp_count}/" $output_bed
#     sed -i "${line_number}s/$/\t${line_number}\t${snp_count}/" $output_bed
    
#     # Increment the line number counter
#     line_number=$((line_number + 1))
# done < 100kb_windows.bed


# awk '{if ($5 > 50)print}' $output_bed > 100kb_windows_with_snp_counts_filtered.bed
# while IFS= read -r line; do
#     line_number=$(echo "$line" | awk '{print $4}')
    
#     filter_vcf="$tmp_dir/${line_number}.filter.vcf"
#     filter_prefix="$tmp_dir/${line_number}.filter"
#     matrix_output="$tmp_dir/${line_number}.matrix"
    
#     # Convert VCF to PLINK binary format
#     $plink --vcf $filter_vcf \
#             --double-id \
#            --snps-only \
#            --make-bed \
#            --out $filter_prefix
    
#     # Calculate pairwise IBS
#     $plink --bfile $filter_prefix \
#            --genome \
#            --distance square ibs flat-missing \
#            --out $matrix_output

# done < 100kb_windows_with_snp_counts_filtered.bed


# rm distance_matrix_${chr}.csv
# if [ ! -f "distance_matrix_${chr}.csv" ]; then
#     ln -s /share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/merge_distance/distance/distance_matrix_${chr}.csv .
# fi

#python ../correlation_synteny_IBS.py 100kb_windows_with_snp_counts_filtered.bed distance_matrix_${chr}.csv 100kb_windows_filtered_pearson_p.xls 11.5 out_pearson.xls out_pearson_p.xls > correlation.log


hor_distance="/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/distance/hor_distance_${chr}.csv"
syn_distance="/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/distance/merge_distance_${chr}.csv"

python ../correlation_synteny_IBS.py 100kb_windows_with_snp_counts_filtered.bed $hor_distance 100kb_windows_filtered_hor_pearson_p_hor.xls 11.5 out_pearson_hor.xls out_pearson_p_hor.xls
python ../correlation_synteny_IBS.py 100kb_windows_with_snp_counts_filtered.bed $syn_distance 100kb_windows_filtered_hor_pearson_p.xls 11.5 out_pearson.xls out_pearson_p.xls
