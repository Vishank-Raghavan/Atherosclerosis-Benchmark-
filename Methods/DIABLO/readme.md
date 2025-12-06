# Multi-Omic Integration with DIABLO (mixOmics)

This directory contains the R script used to benchmark the **DIABLO** (Data Integration Analysis for Biomarker discovery using Latent cOmponents) method. DIABLO is a supervised multi-omics integration framework from the `mixOmics` package designed to identify highly correlated features (genes, proteins) across multiple datasets that discriminate between phenotypic groups.

## Overview

The script performs N-integration of transcriptomic (RNA-seq) and proteomic (Mass Spectrometry) data from matched patient samples to identify a multi-omic signature of atherosclerosis.

**Key Steps:**

1. **Data Loading:** Imports TPM-normalized RNA-seq data and log-transformed Proteomics data.

2. **Harmonization:** Matches samples between RNA (`SRR...`) and Protein (`A...`/`H...`) datasets using metadata from `SraRunTable.csv`.

3. **Preprocessing:** Filters out features with near-zero variance.

4. **Parameter Tuning:** Uses parallelized cross-validation to determine the optimal number of components (`ncomp`) and features to select (`keepX`).

5. **Model Execution:** Runs the final `block.splsda` model.

6. **Extraction:** Exports selected features to `diablo_significant_targets.csv`.

## Files

* **`run_diablo.R`**: The main R script for data processing, tuning, and model execution.

## Input Requirements

This script requires the following files to be present in the working directory (or paths adjusted):

* `rna_seq_tpm_matrix.csv`: Gene expression matrix (Rows: Genes, Cols: Samples).

* `standardized_proteomics_matrix.csv`: Protein abundance matrix (Rows: Proteins, Cols: Samples).

* `SraRunTable.csv`: Metadata linking SRA Run IDs to biological conditions and replicates.

## Usage

Run the script using R or RStudio.

**Command Line:**

```bash
Rscript run_diablo.R
```

**Dependencies:**

* `mixOmics` (Bioconductor)

* `tidyverse`

* `BiocParallel` (For parallel tuning)

## Output

* **`diablo_significant_targets.csv`**: A CSV file containing the selected features (Ensembl IDs for RNA, UniProt IDs for Protein) identified by the model as discriminative for Atherosclerosis vs. Healthy tissue.
