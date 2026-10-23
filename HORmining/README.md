# HORmining

## Overview

HORmining is a bioinformatics pipeline that identifies Higher-Order Repeat (HOR) structures in alpha satellite DNA by integrating two complementary computational approaches: graph-based HORmon and hierarchical tandem repeat mining (HTRM)-based HiCAT.

---

## Installation

### Conda Environment Setup

The pipeline requires specific dependencies which can be installed via Conda:

```bash
# Create environment from the provided YAML file
conda env create -f environment.yml

# Activate the environment
conda activate HORmining
```

## Usage

```
##Identification by HiCAT
python HORmining.py -bed <input.bed> -c <clusters.tsv> -O <output_dir> -m HiCAT

##Identification by HORmon
python HORmining.py -bed <input.bed> -c <clusters.tsv> -O <output_dir> -m HORmon
```

## Minimal server smoke test

The `tests/data` directory contains a small artificial monomer array consisting of repeated `1-2-3` units (36 monomers; 12 repeats). Cluster identifiers must be integers because the summary code normalizes HORs numerically. This dataset is intended only to confirm installation and execution; it is not a biological benchmark.

Run both methods:

```bash
bash tests/run_smoke_test.sh
```

The script creates a timestamped directory under `tests/smoke_outputs`, runs HiCAT and HORmon, verifies their principal output files, and prints the output location. A successful run ends with:

```text
HORmining smoke tests passed
```

For a permanent record of the server test:

```bash
bash tests/run_smoke_test.sh 2>&1 | tee server_smoke_test.log
```

Please retain `server_smoke_test.log` and the generated output directory when reporting test results.


## Parameters

| Parameter | Description | Required | 
|-----------|-------------|----------|
| **BED file**  `-bed`, `--asatbed` | BED-formatted decomposed alpha satellite monomers | Yes |
| **Cluster file**  `-c`, `--clstid` | TSV-formatted monomer cluster information | Yes |
| **Output directory**  `-O`, `--outdir` | Directory for output files | Yes | 
| **Distance threshold**  `-d`, `--distance` | Maximum distance allowed for adjacent monomers to be merged (default: 10,000bp) | No |
| **Method selection**  `-m`, `--method` | HOR mining method: `HORmon` (graph-based) or `HiCAT` (HTRM-based) | Yes | 

## Input files
These files are typically generated using the `SatelliteAnnotationWorkflow`. However, any data structured in the required format is acceptable for the identification of HOR.

- **BED-formatted Decomposed Alpha Satellite Monomers**
   - **Format**: `chromosome\tstart\tend\tMonomer_seqname`
- **TSV-formatted Monomer Cluster Information**
   - **Format**: `Monomer_seqname\tMonomer_cluster_index`

To generate input files from data derived from the `SatelliteAnnotationWorkflow`, use the following command:

```bash
#<input.bed>
awk '{if ($4!~/Sat/)print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"("$6")"}' ${sample}/${hap}/${sample}_${hap}.cenanno.bed | sort -V -k1,1 -k2,2n > ${sample}_${hap}.asat.sort.bed
#<clusters.tsv>
awk '{if ($4!~/Sat/)print $1":"$2"-"$3"("$6")""\t"$4}' ${sample}/${hap}/${sample}_${hap}.cenanno.bed | sort -V -k1,1 -k2,2n > ${sample}_${hap}.mn.index
#run HORmining
python HORmining.py -bed  ${sample}_${hap}.asat.sort.bed -c ${sample}_${hap}.mn.index -O ${outdir} -m HORmon
```

## Output Files
**File:** `$outdir/graph/graph.HORdecomposition.stat.xls`

| Column Name | Description | 
| :--- | :--- | 
| **`blockID`** | Continuous block index. | 
| **`monomer_start_index`** | Start index of HOR in the block. |
| **`monomer_end_index`** | End index of HOR in the block. | 
| **`chromosome`** | Chromosome where the block is located. | 
| **`start`** | Genomic start position  |
| **`end`** | Genomic end position  |
| **`HOR(Mn)`** | The identified Higher-Order Repeat sequence using monomer IDs (e.g., 3_2_1). | 
| **`reordered_HOR(Mn)`** | The HOR sequence after reordering  (e.g., 1_2_3). | 

---
**File:** `$outdir/HiCAT/HiCAT.hor.summary.xls`

| Column Name | Description | 
| :--- | :--- |
| **`blockID`** | Continuous block index. |
| **`block_start`** | Block genomic start position. |
| **`block_end`** | Block genomic start position. |
| **`block_len`** | Total genomic length of the block (bp). |
| **`block_mn_num`** | Total number of monomers in the block. | 
| **`chromosome`** | Chromosome where the block is located. | 
| **`start`** | Genomic start position  | 
| **`end`** | Genomic end position | 
| **`monomer_start_index`** | Start index of HOR in the block. | 
| **`monomer_end_index`** | End index of HOR in the block. | 
| **`nrepeat`** | Number of times the HOR unit repeats. |
| **`HOR`** | The identified Higher-Order Repeat sequence using monomer IDs (e.g., 3_2_1). | 
| **`layer`** | The structural layer assigned to the HOR (top/cover). | 
| **`reordered_HOR`** | The HOR sequence after reordering  (e.g., 1_2_3). | 
