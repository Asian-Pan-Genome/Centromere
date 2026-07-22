# HOR StV Inference & HOR Array Characterization

This directory contains the complete downstream analysis pipeline for Higher-Order Repeat (HOR) structural variant (StV) inference and HOR array characterization across multi-ancestry pangenome assemblies. Two complementary HOR decomposition methods — **HiCAT** and **HORmon** — are used in parallel. Results from both are integrated, classified into StV types, annotated on genome coordinates, and analyzed for pan-genome frequency, diversity, and evolutionary conservation.

## Analysis Overview

The 14 scripts are organized into four functional categories:

| Category                                                       | Scripts | Description                                                                                                                                                                                                                      |
| -------------------------------------------------------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| \[1] HOR Database Construction                                 | 2       | Parse raw HiCAT outputs into unified summary tables; define HOR vs. Dimer; assign pan-genome frequency types                                                                                                                     |
| \[2] HOR StV Classification                                    | 3       | Group HOR StVs into HOR array by monomer composition similarity using Needleman–Wunsch alignment; unify HiCAT and HORmon classifications                                                                                         |
| \[3] HOR Array Annotation & Hierarchy                          | 2       | Annotate genome coordinates with HOR classes and colors; construct nested HOR hierarchy trees                                                                                                                                    |
| \[4] Pan-Genome Diversity, Novel Array & Evolutionary Analyses | 7       | Assign pan-genome types; compute cumulative discovery curves; characterize novel HOR arrays (statistics, prevalence vs. size, and synteny with great ape genomes); visualize chromosome distributions and repeat-number profiles |

---

## Category 1: HOR Database Construction

#### `get_hicat_hor_database.py`

|                 |                                                                                                                                                                                  |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Batch-merge per-haploid HiCAT output into a unified per-sample summary, iterating over all samples listed in a configuration file                                                |
| **Input**       | `input.xls` — three-column table: `sample`, `hap`, `hicat_output_directory`; per-haploid `{blockid}.all_layer.xls.reorder.xls` and `{blockid}.chain.xls` files in each directory |
| **Output**      | `HiCAT.hor.summary.xls` per haploid assembly — each row is a HOR call with sample, haplotype, block index, genomic coordinates, monomer count, HOR string, and layer annotation  |

#### `get_hicat_summary_table.py`

|                 |                                                                                                                                                                                                                                                                                                                                                                              |
| --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Load the raw HiCAT summary, merge redundant monomer cluster IDs, classify each HOR string as HOR (≥3 tandem repeats) or Dimer/DR (<3 repeats), and assign pan-genome frequency categories                                                                                                                                                                                    |
| **Input**       | `all.HiCAT.hor.summary.dedup.xls` — concatenated HiCAT output across all samples                                                                                                                                                                                                                                                                                             |
| **Output**      | `all.HiCAT.hor.summary.final.xls` — master HiCAT summary table; `all.HiCAT.hor.summary.final.HOR_maxrepeat.xls` — unique HOR structural variants (max nrepeat ≥ 3); `all.HiCAT.hor.summary.final.DR_maxrepeat.xls` — unique dimers (max nrepeat < 3); `all.HiCAT.hor.summary.final.pan_hor_dr_type.xls` — HOR/DR pan-genome type classifications                             |


---

## Category 2: HOR StV Classification

#### `alignment_v2.py`

|                 |                                                                                                                                                                                                         |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Custom implementation of the Needleman–Wunsch global alignment algorithm for comparing HOR monomer composition strings (e.g., `"5_12_5_12"`), used as the core similarity engine for HOR classification |
| **Input**       | Two sequences (lists of monomer cluster IDs)                                                                                                                                                            |
| **Output**      | Aligned sequences, edit distance, and match count. The caller uses `max_match / min_len` as the alignment identity ratio                                                                                |

#### `assign_hicat_horclass.py`

|                 |                                                                                                                                                                 |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Classify HiCAT-inferred HOR structural variants into groups (HOR classes) based on monomer composition similarity, using a hierarchical layer-by-layer approach |
| **Input**       | `all.HiCAT.hor.summary.final.HOR_maxrepeat.xls` — unique HORs with their maximum repeat counts and n-mer lengths                                                |
| **Output**      | `all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.xls` — mapping of each HOR string to its class group ID                                                   |

#### `assign_graph_horclass.py`

|                 |                                                                                                                                                                                                                                                             |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Apply the same HOR classification logic to HORmon (graph-based) HORs and integrate with HiCAT HOR class annotations for a unified classification                                                                                                            |
| **Input**       | `all.graph.hordecomposition.summary.final.pan_hor_dr_type.xls` — HORmon HOR summary; `all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.horindex.update.xls` — HiCAT HOR class assignments; `../../cluster.info.xls` — CHM13 monomer cluster annotations |
| **Output**      | `all.graph.hordecomposition.summary.final.horclass.05.xls` — HORmon HOR → group mapping; `*.hicat_horclass.xls` — graph HORs annotated with HiCAT HOR classes; `*.anno.xls` — each HOR annotated with its CHM13 monomer class composition string            |

---

## Category 3: HOR Array Annotation & Hierarchy

#### `merge_HOR_class.py`

|                 |                                                                                                                                                                                                                                                                                                                    |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Description** | Intersect monomer and HOR annotations with genomic coordinates, assign HOR class colors and structural variant indices, and produce BED files for genome browser visualization                                                                                                                                     |
| **Input**       | `all.HiCAT.hor.summary.final.xls` (HiCAT master table); per-sample `filter.asat.bed` (αSat monomer BED); HOR class color tables; HOR StV index and color tables                                                                                                                                                    |
| **Output**      | `{sample}_{hap}.reorder.HOR.bed` — raw HOR BED file; `{sample}_{hap}.HORmn.bed` — merged Monomer/Dimer/HOR annotation with HOR class and color; `{sample}_{hap}.HORclass.bed` — adjacent same-class regions merged (max gap 10 kb); `{sample}_{hap}.HORstv.bed` — same, with StV index and color for visualization |

#### `get_hor_hierarchy.py`

|                 |                                                                                                                                                                                                                                    |
| --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Build a tree-structured hierarchy of nested HORs for each chromosome, assign dot-separated hierarchy levels (e.g., `0.1.2`), and determine which HOR structural variant is dominant in each region for HiCAT HOR StV visualization |
| **Input**       | `all.HiCAT.hor.summary.final.xls` (HiCAT master table); HOR class/color/index tables; CLI arguments: `sample`, `hap`                                                                                                               |
| **Output**      | `HiCAT_bhlayer/{sample}_{hap}.HiCAT.hor.summary.final.hierarchy.bhlayer.xls` — original HiCAT rows augmented with `hierarchy` (nesting depth encoded as `top.sub.sub`) and `bhlayer` (1 \= dominant StV, 0 \= subordinate)         |

---

## Category 4: Pan-Genome Diversity, Novel Array & Evolutionary Analyses

#### `assign_graph_pantype.py`

|                 |                                                                                                                                                                                           |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Assign pan-genome frequency categories (Core / Dispensable / Low frequency / Rare / Singleton) to HORmon (graph-based) HOR structural variants                                            |
| **Input**       | `all.graph.hordecomposition.summary.final.HOR.count.xls` — total count per HOR; `all.graph.hordecomposition.summary.candidate.HOR.count.xls` — per-sample/haplotype/chromosome HOR counts |
| **Output**      | `all.graph.hordecomposition.summary.final.pan_hor_dr_type.xls` — HORmon HORs annotated with pan-genome type, main chromosome, and sample-haplotype count                                  |

#### `get_horstv_cumsum_curve_data.py`

|                 |                                                                                                                                                                                                                                                                                               |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Compute cumulative HOR StV discovery curves: as haploid assemblies are incrementally added (in various orderings), track the number of unique HOR StVs discovered (union) and the number shared by all assemblies examined so far (intersection). Also stratify by pan-genome type            |
| **Input**       | `cumsum.order.xls` — haploid assembly ordering (sorted, APG-first, HPRC-first); HiCAT and HORmon master summary tables; filtered HOR lists excluding novel single/dimeric monomers                                                                                                            |
| **Output**      | `hicat_*_cumsum_ordered_count.xls` / `graph_*_cumsum_ordered_count.xls` — index, sample_hap, union count, shared count; `*_pantype.xls` — same data stratified by pan-genome type; `*_project_uniq_hors.xls` — per-project (APG, HPRC, HGSVC, REF) unique HOR lists for Venn diagram analyses |

#### `plot_horstv_cumsum_curve.R`

|                 |                                                                                                                                                                                                                                                                                                                                                        |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Description** | Visualize cumulative HOR StV discovery as scatter plots and stacked bar charts                                                                                                                                                                                                                                                                         |
| **Input**       | `*_cumsum_ordered_count*.xls` and `*_pantype.xls` files generated by the Python script above; `cumsum.order.xls`                                                                                                                                                                                                                                       |
| **Output**      | `union_cumsum_curve_*_hicat_graph.pdf` — scatter plot: discovered HORs (y) vs. number of haploid assemblies (x), colored by project, shaped by method (HiCAT vs HORmon); `cumsum_hicat_pantype.pdf` / `cumsum_graph_pantype.pdf` — stacked bar charts of cumulative HOR discovery stratified by pan-genome type, with project-indicator bars at bottom |

#### `plot_pan_core_type.R`

|                 |                                                                                                                                                                                                                                                                                             |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Visualize the chromosome-level distribution and repeat-number profiles of HOR StVs stratified by pan-genome frequency type                                                                                                                                                                  |
| **Input**       | `*_hicat_hor_pan_core_type_across_chromosome_stat.xls` or `*_hormon_*` — per-chromosome counts of HOR StVs by pan-type; `*_pan_core_type_repeatnum_across_stvs_stat.xls` — per-StV mean and max repeat numbers by pan-type                                                                  |
| **Output**      | `*_pan_core_type_across_chromosome_stat.pdf` — stacked bar plot of HOR StV count per chromosome; `*_mean_repeatnum_across_stvs_stat.pdf` — ridge plot of log₁₀(mean repeat number) by pan-type; `*_max_repeatnum_across_stvs_stat.pdf` — ridge plot of log₁₀(max repeat number) by pan-type |

---

#### `stat_novel_horarray.py`

|                 |                                                                                                                                                                                            |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Description** | Compute per-sample array lengths and aggregate statistics (mean length, number of sharing haploid assemblies) for HOR arrays classified as novel (not reported in the T2T-CHM13 reference) |
| **Input**       | `Novel_horclass_index.xls` — list of novel HOR class IDs; per-chromosome HiCAT and HORmon `*.horclass.wide.xls` files (chr1–22, X, Y × 2 methods \= 48 files)                              |
| **Output**      | `novelarray_sample_length.xls` — per-novel-class, per-sample array length; `novelarray_sample_count_mean_length.xls` — aggregated mean length and unique sample count per novel class      |

#### `plot_novelarray_pattern.R`

|                 |                                                                                                                                                                                                                                                                   |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Description** | Scatter plot showing the relationship between novel HOR array prevalence and array size, with marginal density distributions                                                                                                                                      |
| **Input**       | `novelarray_sample_count_mean_length.xls` (from `stat_novel_horarray.py`); `Novel_horclass_info.xls` — metadata: suprachromosomal family (SF) assignment, maximum repeat number, CHM13-reported status                                                            |
| **Output**      | `Novelarray_pattern.svg` / `Novelarray_pattern_update.svg` — scatter plot: mean array length (Kb, y) vs. number of haploid assemblies sharing the array (x), with points colored by SF family and sized by max repeat number; marginal density plots on both axes |

#### `plot_novelarray_heatmap_synteny.R`

|                 |                                                                                                                                                                                |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Description** | Visualize pairwise monomer sequence identity matrices between human and great ape (gorilla, chimpanzee, bonobo) HOR arrays to assess evolutionary conservation of novel arrays |
| **Input**       | `*_vs_*_pairwise.xls` — tab-separated three-column files: `human_monomer_index`, `ape_monomer_index`, `sequence_identity`                                                      |
| **Output**      | `*.dotplot.pdf` — ComplexHeatmap with color gradient representing sequence identity                                                                                            |

---

## Data Flow

    ┌─────────────────────────────────────────────────────────────────┐
    │                   HiCAT raw output (per-block)                   │
    │    {blockid}.all_layer.xls.reorder.xls  │  {blockid}.chain.xls   │
    └────────────────────┬────────────────────────────────────────────┘
                         │
                         ▼
         ┌── get_hicat_hor_database.py ──┐
         │  get_hicat_summary_table.py    │
         └────────────────┬───────────────┘
                          │
             ┌────────────┴────────────┐
             ▼                         ▼
       HiCAT HORs               HORmon HORs
       (string-based)           (graph-based)
             │                         │
             ▼                         ▼
    assign_hicat_horclass.py   assign_graph_horclass.py
       (alignment_v2.py)         (alignment_v2.py)
             │                         │
             └─────────┬───────────────┘
                       │
                       ▼
              Unified HOR Classes
                       │
         ┌─────────────┼─────────────┐
         ▼             ▼             ▼
    merge_HOR      get_hor_      assign_graph
    _class.py      hierarchy.py  _pantype.py
         │             │             │
         ▼             ▼             ▼
    Genome BEDs   Hierarchy    Pan-genome
    (coordinates,  trees        types
    colors, StV    (bhlayer)    (Core/.../
    indices)                    Singleton)
                                    │
                        ┌───────────┴───────────┐
                        ▼                       ▼
              get_horstv_cumsum_        plot_pan_core_type.R
              curve_data.py             (chromosome distribution,
                        │                repeat-number ridges)
                        ▼
              plot_horstv_cumsum_curve.R
              (cumulative discovery plots)

     ┌─────────────────────────────────────────────────────────────┐
     │              Novel HOR Arrays & Evolutionary Analyses        │
     │                                                             │
     │  stat_novel_horarray.py  →  plot_novelarray_pattern.R       │
     │                           →  plot_novelarray_heatmap_synteny.R │
     └─────────────────────────────────────────────────────────────┘
