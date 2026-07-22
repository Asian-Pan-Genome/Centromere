# Recombination & Pericentromeric LD Analysis

This directory contains scripts for analyzing recombination, linkage disequilibrium (LD), and phylogenetic congruence in pericentromeric regions. The pipeline comprises: (1) scanning PRDM9 zinc-finger binding motifs and visualizing motif density around centromeres, (2) computing pairwise IBS from SNPs flanking centromeres and correlating genetic distances with HOR synteny distances, and (3) comparing left-arm vs. right-arm SNP phylogenies to test for asymmetric recombination across the centromere.

## Analysis Overview

The 18 scripts are organized into four groups:

| Group | Scripts | Description |
|-------|---------|-------------|
| PRDM9 Motif Scanning | 4 | Scan genomes for PRDM9 binding motifs (FIMO), tile into 20 kb windows, plot motif density with centromere/HOR StV annotation |
| LD & IBS Computation | 2 | Extract pericentromeric SNPs from a joint VCF, compute pairwise LD (r²) via PLINK; tile regions into 100 kb windows, compute per-window IBS distance matrices |
| Correlation & Visualization | 5 | Correlate pairwise IBS with HOR synteny distance (Pearson); generate LD heatmaps, correlation heatmaps, and summary bar/line plots across all chromosomes |
| Phylogenetic Tree Congruence | 7 | Build left/right pericentromeric SNP trees (IQ-TREE), compare topologies via MCI/Quartet/RF/Jaccard/ARI metrics, Mantel-test IBS congruence, and visualize with tanglegrams and genome-wide summary plots |

---

#### `fimo_scan.sh`

| | |
|---|---|
| **Description** | Run FIMO (MEME Suite) to scan a FASTA file for PRDM9 binding motifs at E-value threshold 1e-4 |
| **Input** | CLI args: `motif_file` (PRDM9 motif PWM), `fasta_file`, `output_dir` |
| **Output** | `{dir}/fimo.tsv` — genomic coordinates and statistics of motif hits |

#### `prdm9scan.sh`

| | |
|---|---|
| **Description** | Per-sample orchestration: run FIMO scanning, filter hits (q < 0.3), tile genome into 20 kb windows, intersect with PRDM9 hits, then call R to plot PRDM9 density with centromere and HOR StV annotation |
| **Input** | CLI args: `sample`, `hap`; genome FASTA; `PRDM9_motifs.human.txt`; centromere annotation BED; HiCAT HOR StV BED |
| **Output** | `{sample}_{hap}/fimo.filtered.bed` — filtered motif hits; `{sample}_{hap}/{sample}_{hap}.20kb.PRDM9.bed` — 20 kb window counts; per-chromosome `PRDM9_density_*.png` |

#### `plot_PRDM9_density.R`

| | |
|---|---|
| **Description** | Per-chromosome PRDM9 motif density line plot; for CHM13, overlaid with pericentromeric LD windows (black bars); for other assemblies, overlaid with centromere satellite annotation and HOR StV color tracks |
| **Input** | `{sample}_{hap}.20kb.PRDM9.bed`; centromere annotation BED; HiCAT HOR StV BED; `100kb_windows_with_snp_counts_filtered.bed` (CHM13 only) |
| **Output** | `{sample}_{hap}/PRDM9_density_{chrom}.png` |

#### `plot_PRDM9_density_extcent.R`

| | |
|---|---|
| **Description** | Same as above but zoomed to the extended centromere region (±1 Mb from centromere boundaries), with TE (α-satellite monomer) insertion marks overlaid as colored triangles |
| **Input** | Same as `plot_PRDM9_density.R` plus `TE/{chr}_te.bed` and `te_color.xls` |
| **Output** | `{sample}_{hap}/PRDM9_density_cent_{chrom}.pdf` |

---

#### `LD.sh`

| | |
|---|---|
| **Description** | Extract SNPs from a pan-genome joint VCF for customized pericentromeric regions, filter (MAF ≥ 0.02, genotyping rate ≥ 90%), and compute pairwise LD (r²) using PLINK. Supports both full-region and centromere-excluded modes |
| **Input** | `CHM13-APGp1-HPRCp1-HGSVCp3_snp_diploid.vcf.gz`; per-chromosome BED defining pericentromeric intervals |
| **Output** | `{chr}_customized.ld` — LD r² matrix; `random.ld` — LD for a random SNP subset |

#### `CorrelationHeatmap.sh`

| | |
|---|---|
| **Description** | Tile a pericentromeric BED into 100 kb windows, count SNPs per window (filter windows with >50 SNPs), compute pairwise IBS distance matrix per window via PLINK (`--distance square ibs`), then call `correlation_synteny_IBS.py` to correlate IBS with HOR synteny distance |
| **Input** | `{chr}_customized_region.bed`; pan-genome VCF; `distance_matrix_{chr}.csv` (synteny distance) or `hor_distance_{chr}.csv` / `merge_distance_{chr}.csv` |
| **Output** | `100kb_windows_with_snp_counts_filtered.bed`; per-window IBS matrices in `tmp/`; `out_pearson.xls` / `out_pearson_p.xls` (r and p-value matrices) |

---

#### `correlation_synteny_IBS.py`

| | |
|---|---|
| **Description** | Core computation: for each pair of 100 kb windows flanking the centromere, compute Pearson correlation between (1) left-vs-right IBS, (2) left-IBS vs HOR synteny distance, (3) right-IBS vs synteny distance. Also converts output to square r- and p-value matrices with "CEN" label for the centromeric row/column |
| **Input** | `100kb_windows_with_snp_counts_filtered.bed`; per-window `.matrix.genome` IBS files; synteny distance matrix CSV |
| **Output** | `100kb_windows_filtered_pearson_p.xls` — per-window-pair r and p; scatter plots `tmp/{i}_vs_{j}.dotplot.png`; `out_pearson.xls` / `out_pearson_p.xls` — square matrices |

#### `get_correlation_summary.py`

| | |
|---|---|
| **Description** | Read per-chromosome window coordinates and correlation results, compute summary statistics: left-side LD block length and average r, right-side LD block length and average r, left-to-right linkage average r, and left/right vs centromere flanking average r. Collates results across all chromosomes |
| **Input** | `chroms_linked_index_update.xls` — per-chromosome index definitions (cen_index, left/right linked start/end, lr start/end, cent start/end); `{chr}/100kb_windows_with_snp_counts_filtered.bed`; `{chr}/100kb_windows_filtered_pearson_p.xls` |
| **Output** | `chroms_linked_correlation_summary_update.xls` — 6 rows per chromosome (left_linked, right_linked, lr_left_linked, lr_right_linked, left_cent, right_cent) with start, end, length, avg_r; `chroms_lr_correlation_r_update.xls` — per-window-pair correlation values |

#### `02.plot.LD.relativepos.R`

| | |
|---|---|
| **Description** | Plot LD (r²) as a dot-plot heatmap: SNP index positions on x and y axes, colored by r² (Spectral palette). A black bar marks the centromere position |
| **Input** | `{chr}_customized_exclude_cent.ld` — PLINK LD output; `{chr}_customized_exclude_cent.csv` — SNP-to-index mapping; CLI arg: centromere start position |
| **Output** | `{chr}_customized_exclude_cent.png` — LD heatmap; `{chr}_customized_exclude_cent.csv` — BP-to-index mapping |

#### `plot_heatmap_correlation.R`

| | |
|---|---|
| **Description** | Triangular heatmap of Pearson r values between 100 kb window pairs (IBS vs. synteny correlations), using the Spectral color palette (white = NA, blue-to-red = 0→1). Only the lower triangle is shown |
| **Input** | `{chr}/out_pearson.xls` — r-value square matrix; `{chr}/out_pearson_p.xls` — p-value square matrix |
| **Output** | `{chr}/{chr}_correlation_heatmap_hor_update.pdf` or `_EAS.pdf` |

#### `plot_linkage_length_intensity.R`

| | |
|---|---|
| **Description** | Multi-panel summary figure across all chromosomes: (A) left-side LD block length (bar), (B) left-side average r (dot+line), (C) right-side LD block length, (D) right-side average r, (E) left-right inter-arm linkage average r (with dashed reference line at r=0.3), (F) centromere vs flank average r (left_cent and right_cent colored separately) |
| **Input** | `chroms_linked_correlation_summary_update.xls` — from `get_correlation_summary.py` |
| **Output** | `left_right_linked_pattern_update.pdf` — 4-panel left/right plot; `lr_linked_pattern_update.pdf` — inter-arm linkage plot; `cent_vs_lr_linked_pattern_update.pdf` — centromere-flank correlation plot |

---

### Phylogenetic Tree Congruence

#### `run_chrom_compare.sh`

| | |
|---|---|
| **Description** | Master orchestration for per-chromosome left-vs-right SNP tree comparison. Relabels IQ-TREE tips, drops low-quality samples, prunes to common tips, then runs: (1) `tree_congruence.py` (MCI + Quartet), (2) `compare_trees.R` (RF/Jaccard/ARI + tanglegram), (3) `ibs_congruence.py` (Mantel test). Two versions: `full` (drops GRCh38-mapped samples) and `complete` (additionally drops low-completeness assemblies). Finally generates tanglegrams for HiCAT and HORmon |
| **Input** | CLI arg: `chr`; `{chr}/{chr}_{side}.fasta.treefile` (IQ-TREE); `{chr}/{chr}.map.tsv`; `{chr}/{chr}.drop_full.txt` / `.drop_complete.txt`; side-specific FASTA files |
| **Output** | `{chr}/{chr}_{side}.{full,complete}.nwk`; `{chr}/compare/{chr}_treecong_{ver}.txt`; `{chr}/compare/{chr}_cmp_{ver}.report.txt` / `.tanglegram.pdf`; `{chr}/compare/{chr}_ibs_{ver}_scatter.pdf` / `_summary.txt`; `{chr}/figure/{chr}_tanglegram.{hicat,hormon}.pdf` |

#### `relabel_root_tree.R`

| | |
|---|---|
| **Description** | Relabel IQ-TREE tip names via an `old↦new` mapping table, optionally drop specified tips and root on an outgroup. Handles duplicate labels |
| **Input** | `--tree` (newick), `--map` (TSV: old↦new), `--out`; optional: `--drop-file`, `--root` |
| **Output** | Cleaned newick tree |

#### `tree_congruence.py`

| | |
|---|---|
| **Description** | Compare two trees using pure Python/NumPy: (1) normalized Mutual Clustering Information (MCI) — pairs splits via Hungarian algorithm, (2) Monte-Carlo Quartet similarity — random 4-taxon sampling. Supports bootstrap filtering (0/50/70/95) and permutation-based p-values |
| **Input** | `--tree-a` / `--tree-b` (newick); `--thresholds`; `--quartets`; `--perm`; `--out` |
| **Output** | Text table: threshold, n_splits_A/B, MCI normSim + p, Quartet resSim + strict + neutral% + p |

#### `ibs_congruence.py`

| | |
|---|---|
| **Description** | Compute pairwise IBS distances from left/right SNP FASTA alignments, test correlation via Spearman ρ with Mantel permutation (accounts for pairwise non-independence). Generates d_L vs d_R scatter plot |
| **Input** | `--left` / `--right` (FASTA); `--out-prefix`; `--perm`; `--seed` |
| **Output** | `{out_prefix}_scatter.pdf/.png`; `{out_prefix}_summary.txt` (n, pairs, ρ, r, Mantel p) |

#### `compare_trees.R`

| | |
|---|---|
| **Description** | Comprehensive tree comparison: Robinson-Foulds distance (raw/normalized, with bootstrap-filtered variants at 30/50/70), major-clade Jaccard (min-size 5/10/30/50), Adjusted Rand Index (K=2/3/4 via balanced K-cut). Outputs text report + 3-panel tanglegram |
| **Input** | `--tree-a` / `--tree-b`; `--chrom`; `--label-a` / `--label-b`; `--out-prefix`; optional drop/mismatch files |
| **Output** | `{out_prefix}.report.txt`; `{out_prefix}.tanglegram.pdf` |

#### `tanglegram.R`

| | |
|---|---|
| **Description** | Multi-panel tanglegram: left tree + bootstrap labels, superpopulation strip, project strip, cenhap left-aligned (HiCAT/HORmon HOR StV tracks, min_hor_start anchor), cross-tree links, cenhap right-end-aligned (max_hor_end anchor, x reversed), right tree. Supports great-ape outgroups (bonobo, chimp). Outputs cleaned newick files |
| **Input** | `chrom`; `--left-tree` / `--right-tree`; `--hor-source {hicat,hormon}`; `--cent-chrom`; `--pop`; `--hicat-dir` / `--hormon-dir`; `--color-map`; `--hor-index`; `--top-n-hor`; `--boot-thresh`; `--outdir` |
| **Output** | `{chrom}_{side}.cleaned.nwk`; `{chrom}_tanglegram.{source}.pdf` / `.preview.png`; `{chrom}_cenhap_missing.{source}.tsv` |

#### `plot_summary_v2.py`

| | |
|---|---|
| **Description** | Genome-wide summary bar chart of left-vs-right phylogeny congruence across 24 chromosomes (chr1–22, X, Y). Three metrics: IBS Spearman ρ (red), Quartet similarity ≥70 bootstrap (blue), MCI similarity (green). Two layouts: grouped bars and 3-panel faceted |
| **Input** | `compare_all_summary.tsv` — per-chromosome results; or auto-reads `{chr}/compare/{chr}_treecong_complete.txt` |
| **Output** | `compare_all_summary_grouped.pdf/.png`; `compare_all_summary_panels.pdf/.png` / `_mean.pdf/.png` |

---
