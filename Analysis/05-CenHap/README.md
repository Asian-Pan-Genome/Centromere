# CenHap Classification & Population Analysis

This directory contains scripts for centromeric haplotype (CenHap) classification and downstream population-genetic analyses. CenHaps are defined by clustering haplotype assemblies based on their centromeric HOR array composition (synteny distance matrix of HOR StV content), producing major haplotypes. The scripts cover SNP divergence analysis among chr8 HOR subtypes, superpopulation distribution and array length visualization, phylogeny-guided HOR StV organization plots, geographic mapping of CenHap frequencies, and synteny heatmap clustering.

## Analysis Overview

The 7 scripts are organized into five analysis groups:

| Group | Scripts | Description |
|-------|---------|-------------|
| chr8 HOR SNP Divergence | 2 | Quantify and visualize SqV divergence among three chr8 HOR subtypes |
| CenHap Population Distribution | 1 | Barplot of CenHap frequencies across superpopulations; α-satellite length boxplots per CenHap |
| Phylogeny + HOR StV Map | 1 | Tree-ordered HOR StV organization plots (left/right sides of centromere) with TE insertion marks and CenHap haplotype linking lines |
| Geographic Maps | 2 | World map scatter-pie charts showing CenHap subtype composition per geographic region |
| Synteny Heatmap | 1 | Distance-matrix heatmap clustering of haploid assemblies, with HOR StV tracks appended |

---

#### `stat_chr8_hor_snp_divergence.py`

| | |
|---|---|
| **Description** | Compute intra- and inter-subtype SNP divergence for three chr8 HOR subtypes (h1_11mer, h2_8mer, h3_7mer) by sampling 1,000 aligned sequences each |
| **Input** | `all_shared_dedup.filter.aln.fa` — MSA of shared chr8 HOR sequences |
| **Output** | `{subtype}_intra.div` — pairwise SNP counts within subtype; `{subtypeA}_vs_{subtypeB}_inter.div` — pairwise SNP counts between subtypes |

#### `plot_chr8_hor_snp_div.R`

| | |
|---|---|
| **Description** | Mean ± SD error-bar plot of SqV counts for intra- and inter-subtype comparisons of chr8 HORs |
| **Input** | `all_compare.div` — concatenated intra + inter divergence data |
| **Output** | `chr8_horstv_sqv.pdf` |

#### `plot_cenhap_superpopulation_length.R`

| | |
|---|---|
| **Description** | Multi-chromosome script (chr4, chr5, chr8, chr10, chr17): (1) stacked barplot of CenHap counts per superpopulation, (2) boxplot of α-satellite array length per CenHap, (3) boxplot of centromere length per CenHap, (4) boxplot of HSat1A insertion count per CenHap (chr4 only). chr17 additionally splits by insertion flag. Chromosome is selected by modifying the `chr` variable at the top of the script |
| **Input** | `{chr}_cenhap.xls` — CenHap classification; `cent_chrom.txt`; `populaion.xls`; `all.sat.length.xls`; `chr4_hsat1_insertion.xls` (chr4 only) |
| **Output** | `cenhap_classification/{chr}.superpopulation.barplot.pdf`; `{chr}.asat.boxplot.pdf`; `{chr}.cent.boxplot.pdf`; `{chr}.hsat1.boxplot.pdf` (chr4) |

#### `plot_phy_cenhap.R`

| | |
|---|---|
| **Description** | Generate a three-panel figure for a given chromosome: left side = HiCAT/HORmon HOR StV organization (ordered by left-side phylogeny), center = linking lines connecting left and right phylogeny tips per haplotype (colored by CenHap cluster), right side = HOR StV organization ordered by right-side phylogeny. TE insertion positions are overlaid as colored triangles |
| **Input** | `Phy/phynode_sample_hap.xls`; `Phy/{chr}/{chr}_left_ordered_phy.xls` / `_right_ordered_phy.xls`; `{chr}.HiCAT.horstv.bed` / `.graph.hordecomposition.final.xls`; `{chr}_cenhap.xls`; `TE/{chr}_te.bed`; `te_color.xls`; `populaion.xls`; `cent_chrom.txt` |
| **Output** | `Phy/{chr}/{chr}_left_ordered_out.xls` / `_right_ordered_out.xls`; `Phy/{chr}/{chr}_superpop.xls`; `Phy/{chr}/{chr}_cenhap_anno.xls`; `Phy/{chr}/{chr}_phy_HiCAT_cenhap.pdf`; `Phy/{chr}/{chr}_phy_HORmon_cenhap.pdf` |

#### `chr8_cenhap_map.R`

| | |
|---|---|
| **Description** | World map with scatter-pie charts at each geographic region showing the CenHap subtype composition (A1–A3, B1–B3, C1–C2) on chr8 |
| **Input** | `chr8_cenhap_map.txt` — sample, CenHap subtype, longitude, latitude |
| **Output** | `chr8_cenhap_map_v5.pdf` |

#### `chr4_cenhap_map.R`

| | |
|---|---|
| **Description** | Same as chr8 version but for chr4 CenHaps (A1–A3, B1–B4) |
| **Input** | `chr4_cenhap_map.txt` |
| **Output** | `chr4_cenhap_map_v3.pdf` |

#### `plot_synteny_matrix_heatmap.R`

| | |
|---|---|
| **Description** | Cluster haploid assemblies by centromeric HOR composition distance matrix (average-linkage hierarchical clustering), visualize as a heatmap with superpopulation color bar, and append HiCAT + HORmon HOR StV tracks aligned to the clustering order |
| **Input** | `distance/{hor|merge}_distance_{chr}.csv` — pairwise distance matrix; `{chr}.HiCAT.horstv.bed` / `.graph.hordecomposition.final.xls`; `cent_chrom.txt`; `populaion.xls` |
| **Output** | `Figures/{chr}.distance_{type}_matrix.heatmap.updated.pdf`; `Clustering/{chr}.{type}.distance_clustered.updated.csv`; `Clustering/{chr}.{type}.distance_clustered.result.updated.csv`; `Figures/{chr}.{type}.distance.updated.cenhap.pdf` (heatmap + StV tracks combo) |

---
