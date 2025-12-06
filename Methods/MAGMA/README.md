# Genomic Target Identification with MAGMA

This directory contains the scripts and instructions used to benchmark **MAGMA** (Multi-marker Analysis of Genomic Annotation), a gene-based association analysis method. MAGMA aggregates SNP-level summary statistics from Genome-Wide Association Studies (GWAS) to identify significant gene-level targets associated with atherosclerosis.

## Overview

The workflow integrates large-scale GWAS summary statistics with a reference linkage disequilibrium (LD) panel to prioritize genes. It bridges the gap between variant-level genetic findings (Tier 2) and gene-level therapeutic targets.


## Files

* `MAGMA_analysis.sh`: main script to run MAGMA
* `annotation.genes.annot.gz` : gene annotation required for MAGMA
* `coronary_artery_disease.genes.out`: Unfiltered gene result from MAGMA
* `fdr_adjustment.py`: script to adjust p-vals and filter results
* `magma_fdr_signficant_genes.txt`: significant genes by FDR (q < 0.05)
* `signifcant_genes.txt`: signifcant genes by padj < 0.05

## Usage

```bash
./MAGMA_analysis.sh
python fdr_adjustment.py
```
