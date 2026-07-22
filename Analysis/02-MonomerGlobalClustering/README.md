# Global Monomer Clustering Pipeline

This directory contains the complete workflow for global clustering of α-satellite (αSat) monomer sequences extracted from multi-ancestry pangenome assemblies. The pipeline employs a hybrid clustering strategy combining VSEARCH-based greedy clustering with graph-based refinement to produce non-redundant, high-integrity monomer clusters.

## Overview

To mitigate potential bias introduced by limited sequence diversity in existing αSat databases, we developed a novel workflow for global monomer clustering. αSat monomer positions were identified through alignment to a consensus database generated for each monomeric class from the T2T-CHM13 assembly. Approximately 244 million αSat monomer sequences from APGp1, HPRCy1, HGSVC3, and reference-level complete assemblies (T2T-CHM13, T2T-CN1, Q100-HG002, and T2T-YAO) were processed through the following pipeline.

The pipeline consists of five main stages:

| Stage | Directory | Description |
|-------|-----------|-------------|
| 0 | `00-vsearchclustering/` | Length filtering and VSEARCH clustering |
| 1 | `01-centroidclustering/` | Graph-based centroid clustering with Louvain algorithm |
| 2 | `02-refineclustering/` | Hierarchical refinement and core-periphery assignment |
| 3 | `03-evaluation_with_T2T-CHM13/` | Evaluation against T2T-CHM13 reference annotations |
| 4 | `04-renameclusters/` | Phylogeny-guided cluster reordering and index generation |
| 5 | `05-assignnoncanicalmonomers/` | Assignment of non-canonical-length monomers to clusters |

---

## Stage 0: VSEARCH Clustering

**Script:** `00-vsearchclustering/00clustering_vsearch.sh`

### Description

Initial filtering and greedy clustering of αSat monomer sequences.

### Steps

1. **Length filtering:** Monomers shorter than 120 bp or longer than 200 bp are excluded using `seqkit seq -m 120 -M 200 -g`. Sequences failing this filter are retained separately as `all.mn.short.fasta` and `all.mn.long.fasta` for later assignment (see Stage 5).

2. **Deduplication and clustering:** Filtered sequences are clustered using VSEARCH (v2.30.0) at a specified identity threshold (87%–99%):

    ```
    vsearch --cluster_size <input> --id 0.${threshold} --iddef 0 \
            --centroids all.centroid.fasta --uc all.uc --clusters original_seq/
    ```

    - `--iddef 0`: Uses the standard identity definition (matches / shorter sequence length).
    - `--centroids`: Outputs one representative (centroid) sequence per cluster.
    - `--uc`: Produces a UCLUST-format cluster file.
    - `--clusters`: Writes each cluster's sequences to a separate file under `original_seq/`.

### Key Outputs

| Output | Description |
|--------|-------------|
| `all.mn.filtered.fasta` | Length-filtered monomers (120–200 bp) |
| `all.centroid.fasta` | Centroid sequences of initial VSEARCH clusters |
| `all.uc` | Cluster membership in UCLUST format |
| `original_seq/` | Per-cluster sequence files |

---

## Stage 1: Graph-Based Centroid Clustering

**Scripts:**
- `01-centroidclustering/01centroid_clustering.sh` — Main launcher
- `01-centroidclustering/get_pairwise_identity.py` — Pairwise identity computation
- `01-centroidclustering/centroid_cluster.py` — Similarity graph construction & Louvain clustering
- `01-centroidclustering/merge_cluster_out.py` — Multi-threshold result merging

### Description

Given the greedy nature of VSEARCH clustering, we implemented a graph-based refinement procedure. Centroid sequences from the initial VSEARCH clusters serve as nodes in a similarity graph, where edges connect node pairs whose sequence identity exceeds the initial clustering threshold (edge weight = similarity × 100). The Louvain algorithm is applied to detect community structure, regrouping centroids into more robust clusters.

### Steps

1. **Compute pairwise identities:** All-vs-all pairwise global alignments of centroid sequences are computed using BioPython's `PairwiseAligner`. For each centroid *i*, identities against centroids *j* > *i* are calculated and stored in `tmp/centroid_{i}.pairwise.txt`. The identity metric is defined as matches divided by the length of the shorter sequence, consistent with VSEARCH's `--iddef 0`.

2. **Build similarity graph and cluster:** `centroid_cluster.py` constructs a NetworkX graph from the pairwise identity matrix. Edges are retained when identity ≥ threshold. The Louvain community detection algorithm (`python-louvain`) partitions the graph into clusters. Results are written to `cluster_out/Louvain_{threshold}_cluster.out`.

3. **Merge across thresholds:** `merge_cluster_out.py` concatenates Louvain clustering results across all 13 identity thresholds (87%–99%). Clusters at each threshold are re-indexed by descending size, and the concatenated table is sorted by cluster indices across thresholds (lexicographic sort on `index_87`, `index_88`, …, `index_99`), producing a consistent ordering.

### Key Outputs

| Output | Description |
|--------|-------------|
| `tmp/centroid_{i}.pairwise.txt` | Pairwise identity scores for centroid *i* |
| `cluster_out/Louvain_{threshold}_cluster.out` | Louvain cluster assignments per threshold |
| `Louvain_87-99_concatenated_clusters_sorted.txt` | Merged cluster assignments across all thresholds |

---

## Stage 2: Hierarchical Refinement

**Scripts:**
- `02-refineclustering/inner_similarity_check.py` — Intra-community identity evaluation & splitting
- `02-refineclustering/intra_similarity_check.py` — Intra-cluster pairwise identity computation
- `02-refineclustering/direct_2nd_louvain_clustering.py` — Recursive Louvain sub-clustering
- `02-refineclustering/3nd_intra_subcluster_divergence.py` — Split quality visualization

### Description

To ensure cluster homogeneity and prevent over-inclusion of weakly linked sequences, a hierarchical refinement step evaluates cluster integrity by computing all pairwise identities between centroids within each community. Guided by the triangle inequality principle, a minimum pairwise identity threshold of **2 × threshold − 1** is defined for cluster coherence. Clusters failing this criterion are recursively split.

### Steps

1. **Inner similarity extraction:** `inner_similarity_check.py` extracts all pairwise identities among centroids assigned to the same merged cluster. The data are written to per-cluster files (`tmp/{cluster_id}.inner.txt`).

2. **Core-periphery strategy:** Within each merged cluster, VSEARCH centroids are ranked by size (number of monomer sequences). Clusters accounting for the top 90% of all monomers are designated as **core** centroids. The remaining (smaller) centroids are assigned to a core centroid if their minimum pairwise identity exceeds 95% (for the 87% threshold; equivalent to 2×threshold−1 at higher thresholds); otherwise they form new clusters.

3. **Recursive splitting:** If a cluster contains multiple large centroids with pairwise identity ≤ 93%, indicating inadequate cohesion, a similarity graph is re-constructed within that cluster using a stricter identity cutoff (95%). Louvain community detection (`direct_2nd_louvain_clustering.py`) identifies independent sub-clusters, which are then output as separate refined clusters.

4. **Split evaluation:** `3nd_intra_subcluster_divergence.py` generates boxplots comparing intra-cluster vs. inter-cluster identity distributions for selected clusters before and after splitting, enabling visual verification of improved cluster separation (see Supplementary Note Fig. S1 in the manuscript).

### Key Outputs

| Output | Description |
|--------|-------------|
| `tmp/{cluster_id}.inner.txt` | All pairwise centroid identities within a merged cluster |
| `final.linked.txt` | Refined cluster assignments linking merged cluster IDs to VSEARCH centroid IDs |
| `*.png` | Boxplots of intra- and inter-cluster identity distributions |

---

## Stage 3: Evaluation with T2T-CHM13

**Scripts:**
- `03-evaluation_with_T2T-CHM13/03vsCHM13.sh` — Evaluation launcher
- `03-evaluation_with_T2T-CHM13/vs_chm13.py` — Weighted Jaccard scoring
- `03-evaluation_with_T2T-CHM13/col_threshold.txt` — Column-to-threshold mapping

### Description

Clustering results at each identity threshold are evaluated against the T2T-CHM13 αSat monomer annotation, which provides a ground-truth monomer class label for each sequence.

### Steps

1. **Prepare input:** For each threshold (87%–99%), the concatenated cluster table is grouped by the corresponding column, and monomer sequences are collapsed per cluster.

2. **Compute weighted Jaccard score:** `vs_chm13.py` computes a weighted Jaccard similarity between the predicted clusters and CHM13 annotations. The weighting scheme accounts for both cluster purity (fraction of a cluster matching an annotation class) and annotation recovery (fraction of an annotation class captured by a cluster).

3. **Aggregate results:** Summary statistics are generated per threshold, reporting the mapping between clusters and known αSat monomer classes.

### Key Outputs

| Output | Description |
|--------|-------------|
| `{threshold}.sorted.out` | Cluster-to-annotation mapping per threshold |
| `{threshold}.sorted.out.stat` | Summary statistics of mapping concordance |

---

## Stage 4: Phylogeny-Guided Cluster Renaming

**Scripts:**
- `04-renameclusters/select_marker_seq.py` — Marker sequence selection
- `04-renameclusters/phy.sh` — Phylogenetic tree construction
- `04-renameclusters/reorder_linked.py` — Phylogeny-based reordering
- `04-renameclusters/update_mn_index.py` — Monomer index generation

### Description

Final clusters are assigned meaningful numeric identifiers consistent with their phylogenetic relationships. A representative (marker) sequence is chosen for each cluster, and a phylogenetic tree is constructed to guide cluster ordering.

### Steps

1. **Select marker sequences:** `select_marker_seq.py` selects the largest VSEARCH centroid within each refined cluster as its representative marker sequence. These are extracted and written to `marker_seq.fa` using `pysam`.

2. **Build phylogenetic tree:** `phy.sh` performs multiple sequence alignment of marker sequences with MAFFT, followed by phylogenetic tree construction using both FastTree (approximate maximum-likelihood) and IQ-TREE (maximum likelihood with 1,000 ultrafast bootstrap replicates, automatic model selection with ModelFinder Plus).

3. **Reorder clusters:** `reorder_linked.py` reorders the final clusters according to their position in the phylogenetic tree (FastTree order), producing `final.linked.ordered.v2.txt`.

4. **Generate monomer index:** `update_mn_index.py` updates the global monomer-to-cluster index. Each monomer sequence in every VSEARCH cluster is assigned the new phylogeny-ordered cluster ID. The output (`all.mn.updated.v2.index`) serves as the foundational annotation for all downstream HOR inference analyses.

### Key Outputs

| Output | Description |
|--------|-------------|
| `marker_seq.fa` | Representative sequences for each final cluster |
| `marker_seq.aln.mod.fa` | Multiple sequence alignment of marker sequences |
| `marker_seq.aln.mod.newick` | Phylogenetic tree (FastTree) |
| `final.linked.ordered.v2.txt` | Clusters reordered by phylogeny |
| `all.mn.updated.v2.index` | Final monomer-to-cluster ID mapping |

---

## Stage 5: Assignment of Non-Canonical Monomers

**Script:** `05-assignnoncanicalmonomers/assign_shorted_mn_clusterid.py`

### Description

Monomers that were excluded during length filtering (shorter than 120 bp or longer than 200 bp) are assigned to existing clusters based on sequence similarity to the marker sequences. Each non-canonical monomer is aligned against all marker sequences using global pairwise alignment (BioPython `PairwiseAligner`), and assigned to the cluster with the highest identity.

### Steps

1. Load marker sequences (`marker_seq.fa`) and the phylogenetic cluster order (`fasttree.order`).
2. For each non-canonical monomer, compute global pairwise identity against all markers.
3. Assign the monomer to the cluster of the best-matching marker.
4. Output a cluster index file for the non-canonical monomers.

### Key Outputs

| Output | Description |
|--------|-------------|
| `{input}.index` | Cluster ID assignments for non-canonical-length monomers |

---

## Dependencies

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| [VSEARCH](https://github.com/torognes/vsearch) | ≥2.30.0 | Initial greedy clustering & deduplication |
| [seqkit](https://bioinf.shenwei.me/seqkit/) | — | Sequence length filtering |
| [samtools](http://www.htslib.org/) | — | FASTA indexing |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | — | Multiple sequence alignment |
| [FastTree](http://www.microbesonline.org/fasttree/) | — | Approximate ML phylogeny |
| [IQ-TREE](http://www.iqtree.org/) | ≥1.6.12 | Maximum-likelihood phylogeny with bootstrap |
| [datamash](https://www.gnu.org/software/datamash/) | — | Command-line data grouping |
| [R](https://www.r-project.org/) | ≥4.1 | Statistical computing (optional) |

### Python Libraries

| Package | Purpose |
|---------|---------|
| `BioPython` | Sequence I/O, pairwise alignment |
| `networkx` | Graph construction and manipulation |
| `python-louvain` (`community`) | Louvain community detection |
| `igraph` + `leidenalg` | Alternative Leiden clustering (optional) |
| `numpy` | Numerical operations |
| `pandas` | Tabular data handling |
| `scikit-learn` | Silhouette score, F1 score |
| `pysam` | FASTA random access |
| `matplotlib` + `seaborn` | Visualization |

---
