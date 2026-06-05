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
| **Sample#hap#chromosome**  | Chromosome identifier |
| **Start**               | HOR start position |
| **End**                 | HOR end position |
| **Start_index**            | Start index of monomers |
| **End_index**              | End index of monomers |
| **Repeat_count**                | Copy number of HOR unit |
| **HOR(Dimer)**         | HOR ID or monomer ID |
| **Layer**                  | Annotation layer (top/cover) |
| **ReorderED_HOR(Dimer)**            | Reordered HOR (e. g. 2-3-4-1 ---> 1-2-3-4) |

---

#### **2. `${sample}_${hap}.HiCAT.horstv.bed`**

Provides **decomposed HOR structural variants** for visualization.

**Columns:**

| **Column name**            | **Description** |
|----------------------------|------------------|
| **Sample#hap#chromosome**  | Chromosome identifier |
| **Start**                  | Genomic start |
| **End**                    | Genomic end |
| **HORindex**               | Index of the HOR/monomer |
| **HOR_color**              | Visualization color code |

---

### **HORmon (`HOR_HORmon`)**

For each haploid assembly, the HORmon output provides monomer-level annotations.

**Columns:**

| **Column name**            | **Description** |
|----------------------------|------------------|
| **Sample#hap#chromosome**  | Chromosome identifier |
| **Start**                  | Genomic start |
| **End**                    | Genomic end |
| **HORindex**             | Monomer identifier |
| **HOR_color**              | Color code for monomer/HOR |
