# 🧬 Single-Cell RNA-seq Analysis of Breast Cancer Tumor Microenvironment (Seurat)

## Project Goal
Performed single-cell transcriptomic analysis to identify and characterize tumor and immune cell populations within the breast cancer tumor microenvironment (TME), capturing cellular heterogeneity and functional states.

---

## Overview
This project implements an end-to-end single-cell RNA-seq (scRNA-seq) analysis workflow using a public breast cancer dataset. The pipeline focuses on identifying distinct cellular populations and understanding tumor microenvironment composition.

---

## Pipeline Overview

1. Quality Control & Preprocessing  
2. Normalization & Feature Scaling  
3. Data Integration & Batch Correction  
4. Dimensionality Reduction (PCA, UMAP)  
5. Clustering & Cell Population Identification  
6. Cell-Type Annotation (cluster-based placeholder)  
7. Differential Expression Analysis  
8. Pathway Enrichment Analysis (GO)  

---

## Key Features

- scRNA-seq preprocessing, normalization, and batch correction  
- Integration of high-dimensional transcriptomic data  
- Identification of tumor and immune cell populations  
- Differential expression analysis across cell populations  
- Functional interpretation using pathway enrichment (clusterProfiler)  
- Reproducible and modular analysis pipeline  

---

## Project Structure

single-cell-breast-cancer-tme/  
│── scripts/  
│   ├── 1_qc_preprocessing.R  
│   ├── 2_integration_clustering_umap.R  
│   ├── 3_celltype_annotation_markers.R  
│   ├── 4_differential_expression.R  
│   ├── 5_pathway_enrichment.R  
│  
│── README.md  
│── environment.yml  

---

## Workflow Details

### 1️⃣ Quality Control
- Filter low-quality cells  
- Remove cells with high mitochondrial gene expression  

### 2️⃣ Normalization
- Normalize gene expression values  
- Scale features for downstream analysis  

### 3️⃣ Integration & Batch Correction
- Integrate datasets across batches (Seurat integration workflow)  
- Correct batch effects for consistent clustering  

### 4️⃣ Dimensionality Reduction
- Principal Component Analysis (PCA)  
- UMAP for visualization of cellular structure  

### 5️⃣ Clustering
- Graph-based clustering to identify cell populations  
- Detection of tumor and immune cell clusters  

### 6️⃣ Cell-Type Annotation
- Cluster-based placeholder annotation (to be refined with marker genes)  
- Supports identification of T cells, macrophages, tumor cells, etc.  

### 7️⃣ Differential Expression
- Identify marker genes between clusters  
- Compare transcriptional differences across cell populations  

### 8️⃣ Pathway Enrichment
- Gene Ontology (GO) enrichment using clusterProfiler  
- Functional interpretation of biological processes  

---

## 📈 Example Outputs

- UMAP visualization of cell clusters  
- Cluster marker gene tables  
- Differential expression results  
- GO pathway enrichment plots  

---

## ⚙️ Environment Setup

conda env create -f environment.yml  
conda activate scrna-seq-env  

---

## Tools & Technologies

- R: Seurat, dplyr, ggplot2  
- Bioinformatics: clusterProfiler, org.Hs.eg.db  
- Workflow: Reproducible pipeline with modular scripts  

---

## Skills Demonstrated

- Single-cell RNA-seq analysis  
- High-dimensional data integration and batch correction  
- Clustering and cell-type identification  
- Differential gene expression analysis  
- Pathway enrichment and biological interpretation  
- Reproducible computational workflows  

---

## Impact

This project reflects real-world bioinformatics workflows used in:

- Cancer genomics research  
- Tumor microenvironment analysis  
- Immuno-oncology  
- Precision medicine  

---

## Author

Divya Reddy  
MS Bioinformatics, Georgia Institute of Technology  
