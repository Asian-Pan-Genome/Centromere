# Strand Switch Analysis

Detection and statistical analysis of strand-switch (inversion) breakpoints within satellite DNA arrays across human pan-genome assemblies. Identifies where the strand orientation of satellite arrays flips along the chromosome, and characterizes the frequency, distribution, and correlation of these switches with satellite array length.

## Directory Structure

    strand_swicth/
    ├── README.md
    ├── stat_strand_switch.py              # Per-haplotype strand-switch breakpoint detection
    └── stat_switch_within_array.R         # Population-level statistical analysis and visualization

## Scripts

### 1. `stat_strand_switch.py`

**Description:** Detects strand-switch breakpoints within centromeric satellite arrays for a single haplotype.

**Input:**

- `{sample}/{hap}/{sample}_{hap}.full.cenanno.bed` — Per-haplotype centromere satellite annotation

**Outputs:**

- `{sample}/{hap}/inversion/{sample}_{hap}_sat_strand_merge.bed`

- `{sample}/{hap}/inversion/{sample}_{hap}_sat_inversion_breakpoint.bed` — Strand-switch breakpoints: `chrom, start, end, break_distance, strand(prev/curr), satellite(prev/curr)`

**Usage:**

```bash
python stat_strand_switch.py <sample> <hap> <cenanno_path>
```

---

### 2. `stat_switch_within_array.R`

**Description:** Population-level statistical analysis of satellite strand-switch breakpoints across multiple projects (APGp1, HPRCy1, HGSVC3). Computes general switching patterns per chromosome, tests the correlation between switch count and satellite array length (βSat on acrocentric chromosomes and HSat3 on chr9), classifies α-satellite (αSat) switches by HOR annotation type, and performs chr1 active-array-specific analysis.

**Input:**

- `apg_satellite_inversion_breakpoint_count.csv` — APG project strand-switch breakpoint counts

- `HPRC_HGSVC_satellite_inversion_breakpoint_count.csv` — HPRC/HGSVC project strand-switch breakpoint counts

- `cent_chrom.txt` — Centromere completeness metadata (used to filter incomplete centromeres)

- `all.sat.length.xls` — Satellite array lengths per haplotype

- `all.HORclass.inversion_breakpoint.bed` — α-satellite HOR-classified breakpoints

- `all.info.chr9.hsat3.breakpoint.stats.txt` — HSat3 breakpoint statistics on chr9

- `populaion.xls` — Sample-to-superpopulation mapping

- `chr1_25/chr1_25_inversion_region_length.xls` — chr1 active array (RID\=25) inversion region lengths

**Outputs:**

_Tables:_

- `all_satellite_strand_switch_general_pattern.xls` — Overall switching patterns (shared haploids, average count, SD) per chromosome and breakpoint preference (Supplementary Table 8)

- `apg_satellite_strand_switch_general_pattern.xls` — Same as above, restricted to APG project

- `global_bsat_strand_switch.xls` — βSat-specific switch data with switch frequency (switches per Mb) and length

- `global_hsat3_strand_switch.xls` — HSat3-specific switch data on chr9

- `global_asat_HORclass_switch_general_pattern.xls` — α-satellite HOR-class switch summary per chromosome (shared haploids, average count, SD)

_Figures:_

- `bsat_strand_switch_chrom.pdf` — Scatter plot of βSat switch count vs. βSat array length (Mb) for acrocentric chromosomes (chr1/13/14/15/21/22), with linear regression and Pearson correlation (Supplementary Fig 10)

- `hsat3_strand_switch_chrom.pdf` — Scatter plot of HSat3 switch count vs. HSat3 array length (Mb) on chr9, colored by project, with linear regression and Pearson correlation (Supplementary Fig 10)

- `chr1_25/rid_length_violin.pdf` — Violin plot of inversion region lengths by RID on chr1 (Supplementary Fig 11)
