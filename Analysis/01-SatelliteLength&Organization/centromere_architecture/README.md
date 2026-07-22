# Centromere Architecture (CenArch) Classification

Classification of centromeres based on the order of satellite arrays within each complete centromere across human pan-genome assemblies (APGp1, HPRCy1, HGSVC3, Reference assemblies including T2T-CHM13, T2T-CN1, Q100-HG002 and T2T-YAO). Defines the **CenArch (Centromere Architecture)** classification through satellite-array-order-based approaches, followed by manual curation into Final classes.

## Directory Structure

    Centromere_architecture/
    ├── README.md
    ├── satellite_seq_classification.py          # Satellite-array-order-based classification
    ├── make_class_barplots.py                   # Stacked barplots of CenArch class counts per chromosome
    ├── make_example_figures.py                  # Per-chromosome representative example figures
    └── satellite_seq_classification_all_final_fixed.tsv   # Final curated CenArch classification

## Scripts

### 1. `satellite_seq_classification.py`

**Description:** Classifies each centromere based on the order of its satellite array types along the chromosome (e.g., `ASat>HSat2>HSat3>BSat>HSat3>HSat2`).

**Input:**

- `data/*.full.merged.censat.bed.gz` — Satellite annotation within centromere coordinates (per-haplotype CenSat bed files, available at [CenSat Annotation](https://github.com/Asian-Pan-Genome/Centromere/tree/main/Annotation/CenSat))

- `data/cent_chrom.txt` — Centromere chromosome-level data with completion flag

- `data/populaion.xls` — Sample-to-superpopulation mapping table

**Outputs:**

_Tables:_

- `satellite_seq_classification_all_20kb.tsv` — Centromere architecture classification table of all 24 chromosomes

_Figures (per chromosome, 24 × 2 \= 48 PDFs; all PDFs available at _[Centromere Architecture Figures](https://github.com/Asian-Pan-Genome/Centromere/blob/main/Results/Centromere_Architecture_Figures.zip)):

- `figures/<chr>_seq_20kb_abs.pdf` — Satellite track plot in absolute coordinates (Mb), rows in class-frequency order

---

### 3. `make_class_barplots.py`

**Description:** Generate stacked barplots of complete centromere counts per chromosome, colored by the manually curated Final (CenArch) class.

**Outputs:**
_Figures (Fig. 1c and Supplementary Fig. 12):_

- `class_stacked_barplot_all_v2.pdf` — Stacked barplot across all samples

- `class_stacked_barplot_APG_v2.pdf` — Stacked barplot for APGp1 samples only

---

### 4. `make_example_figures.py`

**Description:** Generate representative centromeres for each Final (CenArch) class per chromosome .

**Outputs:**
_Figures (Supplementary Fig. 13):_

_24 PDFs:_

- `example_figure/<chr>_example.pdf`

---

## Key Result File

### `satellite_seq_classification_all_final_fixed.tsv`

This is the **final curated CenArch classification result file**. It contains 7,960 complete centromeres across 24 chromosomes and 4 projects. The `CenArch` column follows the format `{chr}_{rank}_{count}` (e.g., `chr1_1_173` \= the largest class on chr1, with 173 haplotypes).
