HOR Annotation Files

Two sets of higher-order repeat (HOR) annotation files are provided, generated using HiCAT and HORmon.

HiCAT (HOR_HiCAT)

For each haploid assembly, two output files are available:

1. ${sample}_${hap}_HiCAT.hor.summary.xls

This file contains the complete set of HOR structural variants (StVs), including both top-layer and cover-layer annotations.
Each row includes the following fields:

sample    hap    blockID    block_start    block_end    block_len    block_monomer_count
sample_hap_chromosome    HORstart    HORend    index_start    index_end
nrepeat    HOR/monomer_id    layer    reorder_HOR

2. ${sample}_${hap}.HiCAT.horstv.bed

This BED file provides decomposed HOR StVs, intended primarily for visualization.
Each row includes:

sample_hap_chromosome    start    end    HOR/monomer_id    HORindex    HOR_color

HORmon (HOR_HORmon)

For each haploid assembly, the HORmon output contains monomer-level annotations, with each row structured as follows:

sample_hap_chromosome    start    end    monomer_id    HORarray_index    HORarray_color    HOR_color
