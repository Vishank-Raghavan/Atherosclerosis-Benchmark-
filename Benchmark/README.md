# Atherosclerosis Target Discovery Benchmark

This directory contains the code and data used to benchmark computational methods for identifying therapeutic targets in atherosclerosis. The framework evaluates methods across genomic (DNA), transcriptomic (RNA), and proteomic levels using a standardized scoring system.

## Directory Structure

### Scripts

* **`Atherosclerosis_Benchmark.py`**: The main driver script. It loads the mapped results from all methods, calculates the **Atherosclerosis Benchmark Score (ABS)**, performs dataset signal recovery (DTS) analysis, calculates overlap/consensus, and generates visualizations.

* **`map_magma_to_HGNC.py`**: Helper script to map Entrez Gene IDs from MAGMA output to HGNC Gene Symbols.

* **`convery_twas_IDS.py`**: Helper script to map Ensembl Gene IDs from TWAS/FUSION output to HGNC Gene Symbols.

* **`map_proteins_to_HGNC.py`**: Helper script to map UniProt IDs from Limma/Proteomics output to HGNC Gene Symbols.

* **`map_diablo.py`**: Helper script to map the mixed IDs (Ensembl transcripts + UniProt) from DIABLO output to HGNC Gene Symbols.

* **`map_rna_HGNC.py`**: Helper script to parse the complex ID strings from DESeq2 output and extract HGNC Gene Symbols.

### Input Data (Method Results)

These CSV files contain the raw or mapped significant targets identified by each method.

* `magma_significant_HGNC.csv`: Significant targets from MAGMA (Genomic).

* `TWAS_genes_mapped_filtered.csv`: Significant targets from TWAS/FUSION (Genomic/Transcriptomic).

* `deseq2_HGNC_genes.csv`: Significant differentially expressed genes from DESeq2 (Transcriptomic).

* `wgcna_significant_genes.csv`: Hub genes from WGCNA modules (Transcriptomic).

* `limma_suggestive_HGNC.csv`: Differentially abundant proteins from Limma (Proteomic).

* `diablo_targets_mapped.csv`: Selected features from DIABLO (Multi-omic Integration).

### Output Files

Running `run_benchmark.py` generates these results:

* **`final_benchmark_scores.tsv`**: A table ranking all methods by their ABS score, including details on Tier hits and Pathway coverage.

* **`dts_recovery_scores.tsv`**: A table showing the Dataset Signal Recovery (sensitivity) for each method.

* **`consensus_targets.txt`**: A text file listing specific genes found by multiple methods (Consensus).

* **Visualizations:**

  * `benchmark_abs_score.png`: Bar chart of the final ABS rankings.

  * `benchmark_tier_breakdown.png`: Stacked bar chart showing the number of Tier 1, 2, and 3 targets found by each method.

  * `benchmark_upset.png`: UpSet plot visualizing the intersection of targets across all methods.

## How to Run

### 1. Prerequisites

Ensure you have Python installed along with the following libraries:

```bash
pip install pandas matplotlib seaborn upsetplot mygene
```

### 2. Run the Benchmark

Execute the main script from this directory:

```bash
python Atherosclerosis_Benchmark.py
```

This will read the CSV files listed above, calculate all scores, and generate the output tables and plots.

## Benchmark Methodology

The benchmarking framework evaluates methods based on:

1. **Atherosclerosis Benchmark Score (ABS):** `ABS = TVS x PCM`

  * **Target Validation Score (TVS):** Weighted sum of targets found in our curated 3-Tier Key (Tier 1: Drug Targets, Tier 2: GWAS Loci, Tier 3: Mechanistic/Regulatory).

  * **Pathway Coverage Multiplier (PCM):** Multiplier (1x, 2x, 3x) based on coverage of the three core disease pathways (Lipid Metabolism, Inflammation, Vascular Wall Biology).

2. **Dataset Signal Recovery (DTS):** The percentage of "detectable" signals in the specific input dataset (defined by standard differential analysis) that the method successfully recovered.

3. **Consensus & Novelty:** Analysis of targets found by multiple independent methods versus those unique to a single approach.
