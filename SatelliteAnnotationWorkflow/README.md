# Workflow for Human Centromeric Satellite Annotation

## Overview

This workflow was developed for comprehensive annotation of centromeric satellites in human assemblies, including alpha (α), beta (β), and gamma (γ) satellites, as well as human satellites (HSat1, HSat2, and HSat3). The pipeline generates centromeric coordinates for each human genome assembly and visualizes the results.

---

## Annotation Strategy

- **αSat Annotation**: We provide two monomer databases for annotation:
  
  (1) a consensus database (**```MonomerT2TDB/MonomerT2TDB.fa```**) derived from the T2T-CHM13 assembly, representing monomeric classes from this reference genome;
  
  (2) a representative database (**```MonomerGlobalDB/MonomerGlobalDB.fa```**) constructed from global clustering of 543 human assemblies across diverse ancestries (including assemblies from APGp1, HPRCy1, HGSVC3, T2T-CHM13, T2T-CN1, Q100-HG002, and T2T-YAO).

  Annotation is performed using megablastn followed by stringent filtering. Although initial annotation of APGp1 used ```MonomerT2TDB```, we recommend **```MonomerGlobalDB```** for annotating future high-quality human genome assemblies. This database better captures the broader monomer diversity across populations, and the resulting annotation files serve as direct input for subsequent ```HORmining```, facilitating higher-order repeat identification.
  
- **HSat2 and HSat3 Annotation**: Follows the pipeline established for the T2T-CHM13 genome ([altemose/chm13_hsat](https://github.com/altemose/chm13_hsat)).
- **HSat1, βSat, and γSat Annotation**: Determined from RepeatMasker tracks using customized scripts based on the T2T-CHM13 method.

## Prerequisites

### Software Requirements
- **bedtools**
- **samtools**
- **RepeatMasker**

### Input Files Configuration

Edit `configure.yaml` to specify:

1. **Absolute path for human assemblies**: Assemblies must be formatted as:

```bash
${assemblyDir}/${sample}/${hap}/${sample}_${hap}.fasta
```

2. **Absolute paths for required software**: Specify locations for bedtools, samtools, and RepeatMasker.

## Usage

### Step 1: Build Index for Monomer Databases

```bash
cd MonomerT2TDB/
sh build_index.sh
cd ../MonomerGlobalDB/
sh build_index.sh
```
### Step 2: generate annotation scripts

For each phased assembly, generate annotation scripts using:
```bash
bash run.sh $sample $hap
```
This creates the following scripts in the ```${sample}/${hap}``` directory:

- **```0alpha.sh```** — αSat annotation using MonomerT2TDB

- **```0alphaglobal.sh```** — αSat annotation using MonomerGlobalDB

- **```0hsat23.sh```** — HSat2/HSat3 annotation (T2T-CHM13 pipeline)

- **```0repeat.sh```** — Annotation of HSat1, βSat, γSat

- **```0barplot.sh```** — Merge annotation tracks, generate centromere coordinates and visualization

### Step 3: Run Annotation Scripts

Run the scripts in the following order:
```
# Run annotation steps
bash 0alpha.sh   #MonomerT2TDB
bash 0alphaglobal.sh #MonomerGlobalDB
bash 0hsat23.sh
bash 0repeat.sh

# Final processing and visualization 
# using MonomerT2TDB αSat annotation
bash 0barplot.sh False

# using MonomerGlobalDB αSat annotation
bash 0barplot.sh True
```
## Output files
### Centromeric Coordinates
- **File**: `${sample}/${hap}/${sample}_${hap}.cent.bed`
- **Format**: `chromosome\tstart\tend`
- **Description**: BED file containing the genomic coordinates of the predicted centromeric regions for each haplotype.
  
### Satellite Annotation Tracks
   **Detailed Satellite Annotation**:
   - **File**: `${sample}/${hap}/${sample}_${hap}.cenanno.bed`
   - **Format**: `chromosome\tstart\tend\tsatellite\tscore\tstrand\tstart\tend\tcolor`
   - **Description**: Detailed, unmerged annotation track of satellite repeats within the centromeric regions.

   **Merged Satellite Annotation (for visualization)**:
   - **File**: `${sample}/${hap}/${sample}_${hap}.merged.cenanno.bed`
   - **Format**: `chromosome\tstart\tend\tsatellite\tscore\tstrand\tstart\tend\tcolor`
   - **Description**: Satellite annotation track where adjacent or closely spaced satellite units have been merged based on a distance threshold. This file is optimized for visualization in genome browsers, providing a clearer and more interpretable overview of satellite repeat organization.
  
### Visualization Plots for Manual Quality Control
- **Location**: `${sample}/${hap}/barplot/round1/`
- **Description**: Contains chromosome-level and genome-wide bar plots showing the distribution of satellite repeats within predicted centromeric regions. These plots are generated for manual inspection of the centromere coordinates.

## To Do
Integrate the entire workflow into a Snakemake pipeline for improved reproducibility, scalability, and automated execution.
