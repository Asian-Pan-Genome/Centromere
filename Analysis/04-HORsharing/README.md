# HOR StV Sharing Across Non-Homologous Chromosomes

This directory contains the analysis scripts for investigating HOR structural variant (StV) sharing across non-homologous chromosomes in the pangenome. While most HOR arrays are chromosome-specific, we identified several HOR classes that are shared between two or more non-homologous chromosomes, suggesting inter-chromosomal exchange events. The analyses focus on three chromosome groups: **chr14 & chr22** , **chr13 & chr21** , and **chr1, chr5, chr19 & chr16**. For each group, we quantify sharing patterns, extract representative sequences, assess allelic divergence, and generate visualizations.

## Analysis Overview

The 17 scripts are organized into four functional categories:

| Category | Scripts | Description |
|----------|---------|-------------|
| [1] Global Sharing Statistics & Circos | 3 | Quantify HOR StV sharing across all chromosomes; generate Circos input data for genome-wide visualization |
| [2] chr14 & chr22 Sharing | 2 | Identify shared HOR StVs between chr14 and chr22; extract complete and core conserved sequences |
| [3] chr13 & chr21 Sharing | 5 | Characterize shared StVs across chr13, chr21 (and chr14, chr22); compare allelic divergence between chromosome-specific and shared variants |
| [4] chr1, chr5, chr19 & chr16 Sharing | 7 | Analyze complex multi-chromosome sharing patterns; handle chr1 inversions; extract dimer sequences and assess SNP divergence among three related HOR subtypes |

---

## Category 1: Global Sharing Statistics & Circos

#### `stat_sharing_stv_ratio_coverage.py`

| | |
|---|---|
| **Description** | Compute global statistics of HOR array and StV sharing across chromosomes for both HiCAT and HORmon methods: number of shared arrays, shared StV count and ratio per HOR class, and the fraction of each chromosome's array length occupied by shared StVs |
| **Input** | `all_hicat_horstv_summary.xls` / `all_graph_horstv_summary.xls` — genome-wide HOR StV summaries; `cent_chrom.txt` — centromere coordinates per haploid assembly; per-sample `.HiCAT.horstv.bed` / `.graph.hordecomposition.final.xls` — single-haploid StV annotation files |
| **Output** | `all_{method}_shared_array_haploids_num.xls` — per-class sharing haploid counts; `chr*-*_{method}_shared_horstv_ratio.xls` — per-class shared StV counts and ratios; `chr*-*_{method}_shared_stv_length.xls` — per-sample shared StV lengths and coverage ratios |

#### `stat_circis_config.R`

| | |
|---|---|
| **Description** | Tally the number of HOR StVs per HOR class per chromosome for both HiCAT and HORmon methods, generating input data for Circos karyotype and highlight tracks |
| **Input** | `hicat_hor_final_summary.xls` — HiCAT HOR summary with manual chromosome annotations; `graph_hor_final_summary.xls` — HORmon HOR summary |
| **Output** | `hicat_hor_num_each_chrom.xls` / `hormon_hor_num_each_chrom.xls` — four-column tables: chromosome, HOR class ID, StV count, and reported/novel annotation |

#### `generate_circos_conf_data.py`

| | |
|---|---|
| **Description** | Generate the full set of Circos configuration data: karyotype layout, HOR class highlight tracks (known vs. novel, suprachromosomal family coloring), and population-frequency scatter/heatmap tracks for superpopulation-specific analyses |
| **Input** | `hormon_hor_num_each_chrom.xls` — per-chromosome StV counts; `hormon_horclass_sf.xls` — SF family and color per HOR class; HORmon decomposition master table (`all.graph.hordecomposition.summary.v2.xls`); `populaion.xls` — sample/superpopulation metadata |
| **Output** | `hormon_karyotype.txt` — Circos karyotype definitions; `hormon_horclass_highlight.txt` — known HOR class color bands; `hormon_horclass_sf_highlight.txt` — SF family color bands; per-superpopulation frequency heatmap/scatter files; EAS-vs-AFR / EAS-vs-nonEAS differential scatter files |

---

## Category 2: chr14 & chr22 Sharing

#### `chr14_22_shared_horstv_stat.py`

| | |
|---|---|
| **Description** | Identify and quantify HOR StVs of class 94 that are shared between chr14 and chr22, for both HiCAT and HORmon methods; produce wide-format tables for downstream plotting |
| **Input** | Per-chromosome HiCAT `.HiCAT.horstv.bed` and HORmon `.graph.hordecomposition.final.xls` files for chr14 and chr22; `cent_chrom.txt` |
| **Output** | `chr14_chr22.HiCAT.shared.horstv.stat.xls` / `chr14_chr22.hormon.shared.horstv.stat.xls` — wide tables with sample, HOR string, StV index, and per-chromosome repeat counts |

#### `chr14_22_stat_partial_HOR.py`

| | |
|---|---|
| **Description** | Extract full monomer sequences and core conserved subsequences for HOR class 94 StVs on chr14 and chr22, handling circular permutation to reorder monomers starting from a defined reference position |
| **Input** | CLI arguments: `sample`, `hap`, `chrom`; per-sample `.HiCAT.horstv.bed` or `.graph.hordecomposition.final.xls`; `filter.asat.bed`; genome FASTA |
| **Output** | `tmp/{sample}_{horflag}.fa` — circularly permuted full HOR sequences; `tmp/{sample}_{horflag}_shared.fa` — core conserved subsequences; `tmp/{sample}_{horflag}_oripos.bed` — original genomic coordinates; `tmp/{sample}_horstat.txt` — per-subtype length and count statistics |

---

## Category 3: chr13 & chr21 Sharing

#### `chr13_21_shared_horstv_stat.py`

| | |
|---|---|
| **Description** | Identify HOR StVs shared among chr13, chr21, chr14, and chr22 (HOR class 42 and 94), compute sharing statistics with active/inactive status annotation |
| **Input** | Per-chromosome HiCAT and HORmon files for chr13, chr21, chr14, chr22; `cent_chrom.txt` |
| **Output** | `chr13_chr21_chr14_chr22.HiCAT.shared.horstv.stat.xls` / `.hormon.shared.horstv.stat.xls` — wide-format tables with per-chromosome counts and active/inactive flags |

#### `chr13_chr21_stat_partial_HOR.py`

| | |
|---|---|
| **Description** | Extract full and core conserved sequences for five class 42 HOR subtypes (h1–h5, 7–11mer) on chr13 and chr21; the core shared region comprises monomers 3072–24569–3042–25511–3379–2975 |
| **Input** | Same structure as `chr14_22_stat_partial_HOR.py` but targeting class 42 and five specific HOR subtypes |
| **Output** | Per-subtype full FASTA, shared-core FASTA, original position BED, and length statistics |

#### `chr13_21_divergent_allele_between_shared_stvs.py`

| | |
|---|---|
| **Description** | Quantify SNP divergence within and between chromosome-specific alleles of shared HOR subtypes (h1–h4) by sampling up to 1,000 sequences from a multiple sequence alignment and computing pairwise SNP counts |
| **Input** | `all_shared_filtered_dedup.aln.fa` — MSA of all shared HOR sequences; `.fai` index |
| **Output** | `{subtype}_intra.div` — intra-chromosomal pairwise SNP counts; `{subtypeA}_vs_{subtypeB}_inter.div` — inter-chromosomal pairwise SNP counts (three-column: seq1, seq2, SNP_count) |

#### `chr13_21_plot_shared_sqv.R`

| | |
|---|---|
| **Description** | Generate boxplots comparing per-sequence SNP variant (SqV) counts between chr13-resident and chr21-resident alleles of shared HOR subtypes h3 (10mer) and h4 (7mer), with Wilcoxon rank-sum tests |
| **Input** | `*_inter.div` and `*_intra.div` files from `chr13_21_divergent_allele_between_shared_stvs.py` |
| **Output** | `h3_snp_divergence.pdf` — boxplot: h1_chr13 vs h3_chr21, h2_chr13 vs h3_chr21, h1_chr21 vs h3_chr21 (with p-values); `h4_snp_divergence.pdf` — analogous comparison for the h4 subtype |

#### `plot_chr13_chr21_horstv_examples.R`

| | |
|---|---|
| **Description** | Produce a genome-browser-style visualization of HOR StV organization on chr21 across multiple haploid assemblies, using HORmon decomposition with recolored StV types |
| **Input** | `chr13.graph.hordecomposition.final.xls` / `chr21.graph.hordecomposition.final.xls`; `cent_chrom.txt`; `examples_chr13_chr21_id.xls` — list of sample assemblies to display; HOR StV color update tables |
| **Output** | `Figures/chr21_examples.cenhap.pdf` — horizontal stacked-rectangle plot: each row is one haploid assembly, each colored rectangle is one HOR StV, x-axis = genomic position (Mb) |

---

## Category 4: chr1, chr5, chr19 & chr16 Sharing

#### `chr1_5_19_16_select_dimer.py`

| | |
|---|---|
| **Description** | Extract dimer (2mer) sequences of HOR class 25 from three genomic intervals (r0/r1/r2: flank-proximal, center, flank-distal) per chromosome, handling chr1-specific inversion events and chr16-specific duplication masking |
| **Input** | Per-sample `.HiCAT.horstv.bed`; `filter.asat.bed`; genome FASTA; `chr1_25_inversion_region_length.xls` — coordinates of inverted regions on chr1; `cent_chrom.txt` |
| **Output** | `{sample}_r{0/1/2}_dimer.fa` — dimer sequences from each interval; for chr16, additionally `{sample}_dup_dimer.fa` (duplication region excluded) and filtered `r{0/1/2}_filtered_dimer.fa` |

#### `chr1_5_19_16_stat_target_shared_stv.py`

| | |
|---|---|
| **Description** | Extract the repeat count of 14 target shared HOR StVs (ranging from dimer to 12mer) across the four chromosomes from the HiCAT master summary table |
| **Input** | `all.HiCAT.hor.summary.final.xls` — HiCAT master table |
| **Output** | `HiCAT_target_shared_stv_hicat_layered_info.xls` — full HiCAT rows for the 14 target StVs; `HiCAT_target_shared_stv_hicat_layered_count.xls` — per-sample/haplotype/chromosome summed repeat counts |

#### `chr1_5_19_stat_hor1-3_snp_divergence.py`

| | |
|---|---|
| **Description** | Quantify SNP divergence within and between three related HOR subtypes (h1 = chr19 dimer, h2 = chr5 6mer, h3 = chr1 8mer) by sampling 1,000 aligned sequences each and computing pairwise SNP counts |
| **Input** | `all_hor1-3.uniq.selected.aln.fa` — MSA of the three subtypes |
| **Output** | `{chr}_intra.div` — intra-chromosomal SNP counts; `{chrA}_vs_{chrB}_inter.div` — inter-chromosomal SNP counts |

#### `chr1_5_19_16_upset_shared_horstv.R`

| | |
|---|---|
| **Description** | Generate an UpSet plot showing the intersection patterns of shared HOR StVs (class 25) among chr1, chr5, chr19, and chr16, and compute exact set intersections for each combination |
| **Input** | `hicat_hor_final_summary.xls` — filtered for class 25 entries with multi-chromosome annotations |
| **Output** | `chr1_5_16_19_shared_stv_hicat_upset_ordered.pdf` — UpSet plot |

#### `chr1_5_19_16_plot_shared_linked_boxplot.R`

| | |
|---|---|
| **Description** | For each of the 14 target shared StVs, produce linked boxplots showing repeat count distributions across the four chromosomes, with individual sample connecting lines |
| **Input** | `HiCAT_target_shared_stv_hicat_layered_count.xls` |
| **Output** | `HiCAT_shared_stv_figures/{stv}.shared_boxplot.pdf` — per-StV boxplot + jittered points + grey connecting lines across chromosomes |

#### `chr1_5_19_plot_hor1-3_hor_snp_div.R`

| | |
|---|---|
| **Description** | Visualize the mean ± SD of SqV counts for intra- and inter-chromosomal comparisons of the three HOR subtypes h1/h2/h3, with Wilcoxon tests |
| **Input** | `all_compare.div` — concatenated intra + inter divergence data |
| **Output** | `hor1-3_horstv_sqv.pdf` — mean ± SD error bar plot |

#### `plot_chr1_5_19_16_horstv_examples.R`

| | |
|---|---|
| **Description** | Produce a genome-browser-style visualization of HiCAT HOR StV organization on chr22 across selected haploid assemblies, focused on class 94 |
| **Input** | `chr22.HiCAT.horstv.bed`; `cent_chrom.txt`; `examples_chr22.xls` — sample list; color update tables |
| **Output** | `Figures/chr22_examples.cenhap.pdf` — horizontal stacked-rectangle plot |

---

## Data Flow

```
┌────────────────────────────────────────────────────────────┐
│         03-HORStVs&HORarray outputs                         │
│  (HiCAT / HORmon summaries, HOR class tables, BED files)   │
└────────────────────────┬───────────────────────────────────┘
                         │
         ┌───────────────┼───────────────┐
         ▼               ▼               ▼
   [1] Global       [2] chr14/22    [3] chr13/21   [4] chr1/5/19/16
   Statistics       Sharing         Sharing        Sharing
         │               │               │               │
         ▼               ▼               ▼               ▼
 stat_sharing_    chr14_22_      chr13_21_       chr1_5_19_16_
 stv_ratio_       shared_        shared_         select_dimer.py
 coverage.py      horstv_stat.py horstv_stat.py  chr1_5_19_16_
         │               │               │       stat_target_
         ▼               ▼               ▼       shared_stv.py
 stat_circis_    chr14_22_      chr13_chr21_    chr1_5_19_
 config.R        stat_partial_  stat_partial_   stat_hor1-3_
         │       HOR.py         HOR.py          snp_divergence.py
         ▼               │               │               │
 generate_                │               ▼               │
 circos_conf_             │      chr13_21_                │
 data.py                  │      divergent_               │
         │                │      allele_...py             │
         ▼                │               │               │
   Circos input           │               ▼               ▼
   files                  │      chr13_21_        R plotting:
         │                │      plot_shared_     upset, boxplot,
         ▼                │      sqv.R            snp_divergence,
   (Circos plot           │      plot_chr13_      examples
    generated             │      chr21_           
    externally)           │      horstv_          
                          │      examples.R       
                          │
                          ▼
                   R plotting:
                   (same structure)
```

---
