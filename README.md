# Multi-Omics Atherosclerosis Benchmark Framework

**A Standardized Benchmarking Framework for Evaluating Computational Target Discovery Methods across DNA, RNA, and Protein Levels.**

## Project Description
Atherosclerosis remains a leading global cause of mortality, yet the discovery of clinically actionable therapeutic targets continues to be hindered by the fragmented use of single-omic computational approaches. Existing genomic-, transcriptomic-, and proteomic-level prediction methods often fail to capture cross-modal molecular consistency, resulting in targets with limited translational relevance.

To address this gap, we present a standardized benchmarking framework for evaluating computational target discovery methods. Using publicly accessible, data-matched multi-omics datasets and a curated **3-Tiered Benchmark Key** of validated atherosclerosis targets, we map all method outputs to HGNC gene symbols and assess their performance using the **Atherosclerosis Benchmark Score (ABS)**, dataset signal recovery, and cross-method consensus analysis.

This repository contains the implementation of both **Baseline Methods** (single-omic standard approaches) and **Complex Methods** (multi-omic integration and network analysis) to compare their efficacy in recovering biological signals.

---

## Repository Structure
The codebase is organized by method complexity and data type.

```text
.
├── Data/                       # Shared Patient-Matched Multi-Omics Input Data
│   ├── DNA/                    # Genotypes (PLINK) & GWAS Summary Statistics
│   ├── RNA/                    # RNA-seq Counts & TPM Matrices
│   └── Protein/                # Standardized Proteomics Data
│
├── Methods
│   ├── Baseline_Methods/           # Single-Omic Signal Recovery (Quality Control)
│       ├── DESEQ2/                 # Transcriptomic Baseline (Differential Expression)
│       ├── LIMMA/                  # Proteomic Baseline (Differential Abundance)
│       └── MAGMA/                  # Genomic Baseline (Gene-Level GWAS Association)
│
│   ├── Complex_Methods/            # Advanced Target Discovery Pipelines
│       ├── FUSION_TWAS/            # Transcriptome-Wide Association Study (DNA + RNA)
│       ├── DIABLO/                 # Supervised Multi-Omics Integration (DNA + RNA + Protein)
│       └── WGCNA/                  # Weighted Gene Co-expression Network Analysis (RNA)
│
└── README.md
```

# Pipeline

## Data Sources

All methods utilize a standardized, patient-matched multi-omics dataset derived from human carotid plaque samples (e.g., **PRJNA594843**).

### **Data/DNA/**
- **g1000_eur.bed/.bim/.fam**  
  Individual-level genotype data (1000 Genomes European subset) used for LD reference and genotype-based integration.
- **GWAS_Summary_Stats.tsv**  
  Public Coronary Artery Disease (CAD) summary statistics used for validation and TWAS.

### **Data/RNA/**
- **rna_seq_counts_matrix.csv.gz**  
  Raw RNA-seq counts used for **DESeq2** differential expression.
- **rna_seq_tpm_matrix.csv.gz**  
  Normalized TPM values used for **FUSION** and **WGCNA**.

### **Data/Protein/**
- **standardized_proteomics_matrix.csv**  
  Normalized protein abundance used for **LIMMA** and **DIABLO**.

---

## 🛠️ Prerequisites

### **Software**
- **R (v4.0+)**
  - DESeq2  
  - limma  
  - WGCNA  
  - mixOmics  
  - tidyverse  
  - data.table  

- **Python (v3.8+)**
  - pandas  
  - numpy  
  - scipy  

### **Bioinformatics Tools**
- **PLINK 1.9** — genotype QC and PCA  
- **FUSION** — TWAS analysis  
- **MAGMA** — gene-set analysis  

---

## Methods

## I. Baseline Methods
### **DESeq2 (RNA)**
- **Goal:** Identify differentially expressed genes between symptomatic vs. asymptomatic plaques  
- **Input:** `rna_seq_counts_matrix.csv`  
- **Script:** `DESEQ2_Analysis.R`

### **LIMMA (Protein)**
- **Goal:** Identify differentially abundant proteins using linear modeling  
- **Input:** `standardized_proteomics_matrix.csv`  
- **Script:** `Limma_Analysis.R`

### **MAGMA (DNA)**
- **Goal:** Aggregate SNP-level GWAS associations into gene-level Z-scores  
- **Input:** `GWAS_Summary_Stats.tsv`  
- **Script:** `MAGMA_analysis.sh`

---

## II. Complex Methods

These approaches integrate cross-omic signals and leverage network structure to uncover **novel or subtle molecular targets**.

### **FUSION / TWAS**
- **Goal:** Identify genetically regulated expression associated with CAD  
- **Process:** Integrates GWAS summary stats, eQTLs (GTEx v8 Artery–Coronary), and gene-expression prediction models  
- **Script:** `FUSION.assoc_test.R`

### **DIABLO (mixOmics)**
- **Goal:** Supervised multi-omics integration to identify correlated DNA/RNA/Protein signatures  
- **Script:** `diablo_analysis.R`

### **WGCNA**
- **Goal:** Build weighted co-expression networks to find **hub genes** and key modules correlated with plaque vulnerability  
- **Script:** `wgcna_consensus.R`

---

## Authors
- **Shuhan Lin**  
- **Vishank Raghavan**  
- **Vishnesh Jayanthi Ramanathan**  
- **Saksham Purbey**

**Course:** *CSE 6301 — Georgia Institute of Technology*

---

## References

- **FUSION:**  
  Gusev, A., et al. (2016). *Integrative approaches for large-scale transcriptome-wide association studies.* Nature Genetics.

- **DIABLO:**  
  Singh, A., et al. (2019). *DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays.* Bioinformatics.

- **WGCNA:**  
  Langfelder, P., & Horvath, S. (2008). *WGCNA: an R package for weighted correlation network analysis.* BMC Bioinformatics.

- **DESeq2:**  
  Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology.

- **MAGMA:**  
  de Leeuw, C. A., et al. (2015). *MAGMA: generalized gene-set analysis of GWAS data.* PLoS Computational Biology.

- **GTEx Consortium:**  
  GTEx Consortium. (2020). *The GTEx Consortium atlas of genetic regulatory effects across human tissues.* Science.

