# Structurally Complex Regions and Structural Haplotypes

This directory contains pipelines for epigenetic analysis for centromere. 

## Overview

### 1. CUT&Tag data processing

This part describes the pipeline used to generate normalized CUT&Tag enrichment profiles.
* Low-quality trimming and adapter removal were performed using (fastp)[https://github.com/opengene/fastp].
* Trimed reads were aligned to the corresponding diploid assembly with (BWA)[https://github.com/lh3/BWA] * Multi-mapped reads were discarded using (samtools)[https://github.com/samtools/samtools]
* Normalized CENP‑A enrichment was calculated with (bamCompare)[https://deeptools.readthedocs.io/en/develop/content/tools/bamCompare.html] 
The complete workflow is implemented in the script ```cut_tag.sh```.

### 2. CDR annotation

This part describes the workflow for annotating centromere dip regions (CDRs).
* HiFi or ONT reads are mapped to the diploid genome using (winnowmap2)[https://github.com/marbl/winnowmap];
* Methylation frequencies are called with (pb-CpG-tools)[https://github.com/PacificBiosciences/pb-CpG-tools](HiFi) or modkit[https://github.com/nanoporetech/modkit](ONT);
* CDRs are annotated from the resulting methylation frequency with (CDR-Finder)[https://github.com/arozanski97/CDR-Finder].
The complete pipeline is implemented in the script ```cdr_finder.sh```.

### 3. Centromere local identity

This part describes the pipelines for calculating local sequence identity of centromere region.
*  StainedGlass is used to compute pairwise identity across fixed‑size windows in centromere region.
* The custom script ```src/heat2bed.py``` processes the StainedGlass output to get local sequence identity values per window.
The complete pipeline is implemented in the script ```cdr_id.sh```.
