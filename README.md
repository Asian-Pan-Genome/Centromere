# Centromere and Satellite Annotation

This repository hosts a comprehensive data and workflow resource for **Centromere and Satellite Annotation** in the **APGp1** genome assembly.

---

## 1. Overview

This project provides precise characterization of centromeric regions, including:

* Complete Centromere Coordinates
* Detailed Satellite and Higher-Order Repeat (HOR) Annotation Tracks
* Robust Pipelines for HOR identification and Satellite Annotation.

---

## 2. Data Access

The `Annotation` directory contains the final, high-quality annotation data.

- Complete centromere coordinates in APGp1 : [APGp1_complete_centromere](https://github.com/Asian-Pan-Genome/Centromere/raw/refs/heads/main/Annotation/APGp1_complete_centromere_coordinates.bed)

- Satellite tracks for each phased assembly: [native_bed_format](https://github.com/Asian-Pan-Genome/Centromere/tree/main/Annotation/CenSat)

- HOR identified by HORmon: [HOR_HORmon](https://github.com/Asian-Pan-Genome/Centromere/raw/refs/heads/main/Annotation/HOR/HOR_HORmon.zip)

- HOR identified by HiCAT: [HOR_HiCAT](https://github.com/Asian-Pan-Genome/Centromere/blob/main/Annotation/HOR/HOR_HiCAT.zip)

## 3. Recommended Workflows

We provide and highly recommend the following two pipelines for future T2T human assemblies and related satellite studies.

### SatelliteAnnotationWorkflow

**Directory:** `SatelliteAnnotationWorkflow`

This workflow was developed for the comprehensive annotation of centromeric satellites (including $\alpha$, $\beta$, $\gamma$, HSat1, HSat2, and HSat3) in human assemblies. The pipeline outputs precise centromeric coordinates for each genome and provides scripts for result visualization.

* [**Access Workflow**](https://github.com/Asian-Pan-Genome/Centromere/tree/main/SatelliteAnnotationWorkflow)

### HORmining Pipeline

**Directory:** `HORmining`

HORmining is a bioinformatics pipeline designed for robust identification of Higher-Order Repeat (HOR) structures within alpha satellite DNA. It integrates two complementary computational approaches: the graph-based **HORmon** algorithm and the hierarchical tandem repeat mining (HTRM)-based **HiCAT** algorithm. 

* [**Access Workflow**](https://github.com/Asian-Pan-Genome/Centromere/tree/main/HORmining)

---

## 4. Integrated Analysis (Code Archive)

The `Analysis` directory contains the source code used for the final, integrated analysis of global centromere diversity, structural characterization, and evolutionary comparisons.
