## HOR Annotation Files

Two sets of higher-order repeat (HOR) annotation files are provided, generated using **HiCAT** and **HORmon**.

---

### **HiCAT (`HOR_HiCAT`)**

For each haploid assembly, two output files are included:

---

#### **1. `${sample}_${hap}_HiCAT.hor.summary.xls`**

Contains the full set of HOR structural variants (StVs), including both **top-layer** and **cover-layer** annotations.

**Columns:**

| **Column name**            | **Description** |
|----------------------------|------------------|
| **sample**                 | Sample ID |
| **hap**                    | Haploid (e.g., Mat/Pat) |
| **blockID**                | Block identifier |
| **block_start**            | Block start coordinate |
| **block_end**              | Block end coordinate |
| **block_length**              | Length of the block |
| **block_monomer_count**    | Number of monomers in the block |
| **sample_hap_chromosome**  | Chromosome identifier |
| **horstart**               | HOR start position |
| **horend**                 | HOR end position |
| **index_start**            | Start index of monomers |
| **index_end**              | End index of monomers |
| **nrepeat**                | Copy number of HOR unit |
| **HOR/monomer_id**         | HOR ID or monomer ID |
| **layer**                  | Annotation layer (top/cover) |
| **reorder_HOR**            | Reordered HOR (e. g. 2-3-4-1 ---> 1-2-3-4) |

---

#### **2. `${sample}_${hap}.HiCAT.horstv.bed`**

Provides **decomposed HOR structural variants** for visualization.

**Columns:**

| **Column name**            | **Description** |
|----------------------------|------------------|
| **sample_hap_chromosome**  | Chromosome identifier |
| **start**                  | Genomic start |
| **end**                    | Genomic end |
| **HOR/monomer_id**         | HOR ID or monomer ID |
| **HORindex**               | Index of the HOR/monomer |
| **HOR_color**              | Visualization color code |

---

### **HORmon (`HOR_HORmon`)**

For each haploid assembly, the HORmon output provides monomer-level annotations.

**Columns:**

| **Column name**            | **Description** |
|----------------------------|------------------|
| **sample_hap_chromosome**  | Chromosome identifier |
| **start**                  | Genomic start |
| **end**                    | Genomic end |
| **monomer_id**             | Monomer identifier |
| **HORrarray**               | HOR array assignment |
| **HORarray_color**         | Color code for HOR array |
| **HOR_color**              | Color code for monomer/HOR |
