# Proteomic Analysis with LIMMA

This directory contains the R script used to benchmark the **LIMMA** (Linear Models for Microarray Data) method for identifying proteomic targets in atherosclerosis. Although originally developed for microarrays, LIMMA's empirical Bayes approach is highly effective for mass spectrometry-based proteomics data, especially for studies with limited sample sizes.

## Overview

The script performs differential abundance analysis on a standardized proteomics dataset to identify proteins that are significantly up- or down-regulated in atherosclerotic plaques compared to healthy tissue.

**Key Steps:**

1. **Data Loading:** Imports the log-transformed, standardized protein abundance matrix.

2. **Metadata Parsing:** Extracts biological conditions (Atheroma vs. Healthy) directly from the sample column names (e.g., `A1`, `H1`).

3. **Filtering:** Removes proteins with excessive missing values (>50%) to ensure robust statistical modeling.

4. **Modeling:** Fits a linear model to each protein using `lmFit`, accounting for the experimental design.

5. **Testing:** Computes moderated t-statistics and P-values using `eBayes` for the contrast `Atheroma - Healthy`.

6. **Extraction:** Exports two lists of targets:
   * **High-Confidence:** FDR adjusted P-value < 0.05.
   * **Suggestive:** Nominal P-value < 0.05.

## Files

* **`LIMMA Analysis.R`**: The main R script for data loading, linear modeling, and result extraction.

## Input Requirements

This script requires the following file to be present in the working directory:

* `standardized_proteomics_matrix.csv`: A matrix where rows are Protein IDs (UniProt) and columns are Samples. Values should be log-transformed abundances.
  * *Note:* Column names must start with 'A' (Atheroma) or 'H' (Healthy) for the automatic group parsing to work (e.g., `A1_DIA`, `H1_DIA`).

## Usage

Run the script using R or RStudio.

**Command Line:**

```bash
Rscript LIMMA Analysis.R
```
**Dependencies:**

* limma (Bioconductor)

* tidyverse

## Output
* `limma_significant_proteins_FDR05.csv`: A CSV file containing the list of proteins identified as significantly differentially abundant (FDR < 0.05), including their logFC and adjusted P-values.

* `limma_suggestive_proteins_P05.csv`: A broader list of proteins with nominal P-value < 0.05, useful for sensitivity analysis.
