# Centromere Satellite Variation

Analysis of centromere size and satellite array length variation across human pan-genome assemblies (APGp1, HPRCy1, HGSVC3, Reference assemblies including T2T-CHM13, T2T-CN1, Q100-HG002 and T2T-YAO), including maternal–paternal haplotype comparisons, cross-chromosome distributions, and outlier detection.

## Directory Structure

    centromere-satellite-variation/
    ├── README.md
    ├── data/
    │   ├── all.sat.length.xls       # Satellite array lengths (all satellites, all samples)
    │   ├── cent_chrom.txt           # Centromere chromosome-level data with completion flag
    │   └── populaion.xls            # Sample-to-superpopulation mapping table
    ├── compute_mat_pat_comprehensive.R
    ├── plot_mat_vs_pat_satellite_length.R
    ├── plot_satellite_length_across_chroms.R
    ├── run_asat_signed_plots.R
    ├── run_hsat_signed_density.R
    ├── summarize_satlen.R
    ├── boxplot_satlen_superpop_downsampling.R
    └── batch_effect_APG_vs_nonAPG_EAS.R

## Scripts

### 1. `compute_mat_pat_comprehensive.R`

**Description:** Core data-processing script. Pairs maternal and paternal haplotypes (or haplotype1/haplotype2 for HGSVC3) across all projects, computes centromere size differences and satellite array length differences, performs per-chromosome outlier detection (mean ± 2 SD), and outputs a comprehensive summary table.

**Outputs:**

_Table (Supplementary Table 5):_

- `allsamples_mat_vs_pat_comprehensive.tsv` — Per-sample maternal–paternal centromere and satellite array comparison table

---

### 2. `plot_mat_vs_pat_satellite_length.R`

**Description:** Comprehensive visualization and statistical analysis of satellite array parental differences. Generates density plots (per chromosome × satellite), variance-share decomposition, cross-superpopulation comparisons, and maternal/paternal bias tests.

**Outputs:**

_Tables:_

- `apg_mat_vs_pat_satellite_length_summary.xls` — Per-sample APG satellite pair differences with outlier flags

- `apg_mat_vs_pat_satellite_combo_stats.xls` — Per (chromosome × satellite) statistics (n, mean, median, SD, outlier counts)

- `apg_mat_vs_pat_satellite_bias_test.xls` — One-sample t-test for maternal/paternal bias (BH-adjusted p-values)

- `apg_mat_vs_pat_satellite_variance_share.xls` — Variance-share decomposition: each satellite's contribution to total (Mat − Pat) centromere size variance per chromosome

- `allproj_mat_vs_pat_satellite_length_summary.xls` — All-project (APG + HPRC + HGSVC + Ref) satellite pair summary

_Cross-superpopulation tables:_

- `apg_mat_vs_pat_satellite_superpop_pair_counts.xls` — Pair counts per superpopulation × chromosome × satellite

- `apg_mat_vs_pat_satellite_superpop_pair_counts_wide.xls` — Wide-format version of the above

- `apg_mat_vs_pat_satellite_superpop_pair_summary.xls` — Superpopulation-level pair-completeness summary

- `apg_mat_vs_pat_satellite_superpop_downsample_wilcox.xls` — Downsampled Wilcoxon test: EAS vs other superpopulations

_Figures (Supplementary Fig. 7b):_

- `apg_mat_vs_pat_satellite_variance_share.pdf` — Stacked-bar chart: % of total Mat − Pat variance by satellite per chromosome

- `apg_mat_vs_pat_satellite_length_density.pdf` — Faceted density plots (chromosome × satellite) with BH-adjusted p-values and sample sizes annotated

- `apg_mat_vs_pat_satellite_superpop_density.pdf` — Superpopulation-stratified density plots with downsampled Wilcoxon p-values annotated

---

### 3. `plot_satellite_length_across_chroms.R`

**Description:** Visualize satellite array length distributions across chromosomes for the APG project. Generates violin plots, error-bar plots with outlier marking, and detailed outlier tables stratified by satellite, sample, and chromosome.

**Outputs:**

_Figures (Supplementary Fig. 2):_

- `apg_satellite_length_variation.pdf` — Violin plots with reference genome data points overlaid

- `apg_satellite_length_variation_errorbar.pdf` — Error-bar plots (mean ± SD) with outlier points (×) and reference genome points overlaid

---

### 4. `run_asat_signed_plots.R`

**Description:** Generate signed (Mat − Pat) density and histogram plots specifically for αSat (ASat). Includes mean ± 2 SD vertical reference lines with annotated values.

**Outputs:**

_Figures (Supplementary Fig. 7a):_

- `apg_mat_vs_pat_asat_length_signed_density.pdf` — αSat signed density plot

- `apg_mat_vs_pat_asat_length_signed_histogram.pdf` — αSat signed histogram

---

### 5. `run_hsat_signed_density.R`

**Description:** Generate signed (Mat − Pat) density and histogram plots for HSat1 (merged from HSat1A + HSat1B), HSat2, and HSat3. Each satellite uses appropriate x-axis limits based on its observed range (±8.5 Mb for HSat1, ±16 Mb for HSat2, ±14 Mb for HSat3).

**Outputs:**

_Figures (Supplementary Fig. 7b):_

- `apg_mat_vs_pat_hsat1_signed_density.pdf` — HSat1 signed density plot

- `apg_mat_vs_pat_hsat1_signed_histogram.pdf` — HSat1 signed histogram

- `apg_mat_vs_pat_hsat2_signed_density.pdf` — HSat2 signed density plot

- `apg_mat_vs_pat_hsat2_signed_histogram.pdf` — HSat2 signed histogram

- `apg_mat_vs_pat_hsat3_signed_density.pdf` — HSat3 signed density plot

- `apg_mat_vs_pat_hsat3_signed_histogram.pdf` — HSat3 signed histogram

---

### 6. `summarize_satlen.R`

**Description:** Generate descriptive statistics (mean, SD, Q1, median, Q3) for satellite array lengths across all projects, stratified by chromosome × satellite × project, plus a global (all-projects combined) summary.

**Outputs:**

_Table (Supplementary Table 3):_

- `satlen_summary_stats.tsv` — Descriptive statistics per chromosome × satellite × project and Global

---

### 7. `boxplot_satlen_superpop_downsampling.R`

**Description:** Compare satellite array lengths across superpopulations (AFR, AMR, EAS, EUR, SAS) for all five satellite types (ASat, BSat, HSat1, HSat2, HSat3). Uses downsampled pairwise Wilcoxon tests (K \= 10 iterations, min_n \= 5) to control for sample-size imbalance, and generates boxplots with significance brackets.

**Outputs:**

_Tables (Supplementary Table 4):_

- `downsampling_wilcox_ASat.tsv` — Pairwise downsampled Wilcoxon p-values for ASat

- `downsampling_wilcox_BSat.tsv` — Pairwise downsampled Wilcoxon p-values for BSat

- `downsampling_wilcox_HSat1.tsv` — Pairwise downsampled Wilcoxon p-values for HSat1

- `downsampling_wilcox_HSat2.tsv` — Pairwise downsampled Wilcoxon p-values for HSat2

- `downsampling_wilcox_HSat3.tsv` — Pairwise downsampled Wilcoxon p-values for HSat3

_Figures (Supplementary Fig. 4 and 5):_

- `boxplot_satlen_superpop_ASat.pdf` / `.png` — ASat boxplots across superpopulations (24-chromosome panel, 4 rows × 6 columns, with AFR-vs-others significance brackets)

- `boxplot_satlen_superpop_BSat_HSat.pdf` / `.png` — Combined boxplots for BSat, HSat1, HSat2, HSat3 with satellite-name side tags

---

### 8. `batch_effect_APG_vs_nonAPG_EAS.R`

**Description:** Detect potential batch effects within the EAS superpopulation by comparing satellite array lengths between APG project samples and non-APG (HPRC + HGSVC) samples. Uses downsampled Wilcoxon tests (K \= 10) per satellite × chromosome combination, with significance brackets on boxplots.

**Outputs:**

- `batch_effect_APG_vs_nonAPG_EAS_wilcox.tsv` — Per (satellite × chromosome) downsampled Wilcoxon results comparing APG vs non-APG EAS samples

_Figure (Supplementary Fig. 3):_

- `batch_effect_EAS_APG_vs_nonAPG_all_v6.pdf` — Multi-row faceted boxplot with significance brackets (ASat split into 2 rows, BSat/HSat1/HSat2/HSat3 each 1 row)

---

## Key Methodological Notes

- **HSat1 merging:** HSat1A and HSat1B arrays are summed to form a single HSat1 measurement in all scripts.

- **Filtering:** Centromeres with `filterflag != 0` are excluded. Satellite arrays on chromosomes with flagged centromeres are also excluded.

- **Reliable satellite × chromosome combinations:** Only combinations with mean satellite length > 500 kb are retained. Additionally, HSat1/2/3, BSat, and GSat on acrocentric chromosomes (chr13, chr14, chr15, chr21, chr22) are excluded due to known assembly difficulties in these regions.

- **HGSVC3 phasing:** HGSVC3 assemblies use hap1/hap2 (no parental phase information); hap1 − hap2 differences are treated analogously to Mat − Pat but with arbitrary sign.

- **Outlier detection:** Per-chromosome mean ± 2 SD of centromere size difference.
