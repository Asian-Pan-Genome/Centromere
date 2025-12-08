# Workflow for Human Centromeric Satellite Annotation

## Overview

This workflow was developed for comprehensive annotation of centromeric satellites in human assemblies, including alpha (α), beta (β), and gamma (γ) satellites, as well as human satellites (HSat1, HSat2, and HSat3). The pipeline generates centromeric coordinates for each human genome assembly and visualizes the results.

---

## Annotation Strategy

- **αSat Annotation**: Uses a consensus database generated for each monomeric class from the T2T-CHM13 genome assembly. Identification is performed using megablastn with stringent filtering.
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

1. **MonomerT2TDB** (`MonomerT2TDB/MonomerT2TDB.fa`): Consensus sequences for each monomeric class annotated in T2T-CHM13. Used for initial αSat monomer annotation in APGp1.
To build the database index:

```bash
cd MonomerT2TDB/
sh build_index.sh
```
### Step 2: generate annotation scripts

For each phased assembly, generate annotation scripts using:
```bash
bash run.sh $sample $hap
```
This creates the following scripts in the ```${sample}/${hap}``` directory:

- **```0alpha.sh```** — αSat annotation using MonomerT2TDB

- **```0hsat23.sh```** — HSat2/HSat3 annotation (T2T-CHM13 pipeline)

- **```0repeat.sh```** — Annotation of HSat1, βSat, γSat

- **```0barplot.sh```** — Merge annotation tracks, generate centromere coordinates and visualization

### Step 3: Run Annotation Scripts

Run the scripts in the following order:
```
# Run annotation steps
bash 0alpha.sh   
bash 0hsat23.sh
bash 0repeat.sh

# Final processing and visualization
bash 0barplot.sh
```
## Output files
- **Centromeric Coordinates**:
  
  ```${sample}/${hap}/${sample}_${hap}.round1.cenpos```
      formated as ```chromosome\tstart\tend```
- **Satellite Annotation Track**:
  
  ```${sample}/${hap}/${sample}_${hap}.cenanno.bed```
  formated as ```chromosome\tstart\tend\tsatellite\tscore\tstrand\tstart\tend\tcolor```
  
- **Chromosome-level and genome-wide centromere visualization for manual QC**：
  located in ``` ${sample}/${hap}/barplot/round1/```

