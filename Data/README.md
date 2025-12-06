# Data Processing and Standardization

This directory contains the data files and scripts used to prepare the standardized multi-omic inputs for the benchmarking framework.

## Overview

The benchmark relies on a curated set of publicly available datasets where RNA-seq and Proteomics data were derived from the same patient cohort (matched samples), ensuring a fair evaluation of multi-omic integration methods without the confounding effects of "chimera" datasets.

## Datasets

### 1. Genomic Data (DNA)
* **Source:** CARDIoGRAMplusC4D Consortium (1000 Genomes-based meta-analysis).
* **Description:** GWAS summary statistics for Coronary Artery Disease (CAD) from >184,000 individuals.
* **Processing:**
    * Filtered for variants with valid P-values and effect sizes.
    * Converted to `sumstats` format (Z-scores) for FUSION/TWAS.
    * Mapped to NCBI37.3 gene locations for MAGMA.
* **Key Files:**
    * `sample_gwas_data.tsv`: The primary summary statistics input.

### 2. Transcriptomic Data (RNA)
* **Source:** NCBI SRA Accession **PRJNA594843**.
* **Description:** Bulk RNA-seq of human carotid atherosclerotic plaques and healthy arteries.
* **Processing:**
    * **Quality Control:** `fastp` (adapter trimming, quality filtering).
    * **Quantification:** `Salmon` (mapping to GENCODE v44).
    * **Aggregation:** Transcript-level counts aggregated to gene-level matrices.
* **Key Files:**
    * `rna_seq_counts_matrix.csv`: Raw integer counts (Input for DESeq2).
    * `rna_seq_tpm_matrix.csv`: TPM normalized counts (Input for DIABLO, WGCNA).
    * `SraRunTable.csv`: Metadata linking Run IDs to biological conditions.

### 3. Proteomic Data (Protein)
* **Source:** ProteomeXchange Accession **PXD056909**.
* **Description:** Data Independent Acquisition (DIA) mass spectrometry of the *same* carotid plaque cohort.
* **Processing:**
    * Filtered for high-confidence proteins (Q-value < 0.01).
    * Removed proteins with >50% missing values.
    * Imputed missing values (min-value imputation).
    * Log2-transformed abundances.
    * Mapped UniProt IDs to HGNC Gene Symbols.
* **Key Files:**
    * `standardized_proteomics_matrix.csv`: The processed protein expression matrix.

## Scripts

* **`sra_to_fastq.sh`**: Batch converts SRA files to FASTQ format.
* **`run_fastp.sh`**: Performs quality control and adapter trimming on raw reads.
* **`build_index.sh`**: Builds the Salmon index for the human transcriptome.
* **`run_salmon.sh`**: Quantifies gene expression from clean FASTQ files.
* **`aggregate_counts.py`**: Compiles individual Salmon outputs into a single gene expression matrix.
* **`process_proteomics.py`**: Cleans, filters, and pivots the raw Spectronaut proteomics report into a standardized matrix.
* **`format_fusion_from_pval.py`**: Converts GWAS summary statistics (Beta/P-value) into Z-scores for TWAS analysis.

## Usage

These scripts are typically run sequentially to build the input datasets from raw downloads.

**Example: Processing RNA-seq Data**
```bash
./sra_to_fastq.sh
./run_fastp.sh
./run_salmon.sh
python aggregate_counts.py
```
**Example: Processing Proteomics Data**
```bash
python process_proteomics.py
```
