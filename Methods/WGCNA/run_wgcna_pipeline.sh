#!/usr/bin/env bash
# ==============================================================================
# WGCNA Pipeline - Master Execution Script
# ==============================================================================
# Purpose: Run the complete WGCNA analysis pipeline from start to finish
#
# Usage:
#   bash run_wgcna_pipeline.sh
#
# Runtime: Approximately 15-30 minutes
# ==============================================================================

set -e  # Exit on error

echo "==============================================================================="
echo "WGCNA NETWORK-BASED ANALYSIS PIPELINE"
echo "Atherosclerosis Benchmark Framework"
echo "==============================================================================="
echo ""

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R."
    exit 1
fi

# Check working directory
if [ ! -f "01_prepare_data.R" ]; then
    echo "ERROR: Please run this script from the Methods/WGCNA directory"
    exit 1
fi

# Check if input data exists
if [ ! -f "../../Data/RNA/rna_seq_counts_matrix.csv.gz" ]; then
    echo "ERROR: RNA-seq count matrix not found at ../../Data/RNA/rna_seq_counts_matrix.csv.gz"
    exit 1
fi

echo "✓ Prerequisites checked"
echo ""

# Step 1: Data Preparation
echo "==============================================================================="
echo "STEP 1/4: Data Preparation"
echo "==============================================================================="
echo "  - Aggregating transcripts to gene level"
echo "  - Filtering low-expression genes"
echo "  - Normalizing with VST"
echo ""

Rscript 01_prepare_data.R

if [ $? -ne 0 ]; then
    echo "ERROR: Step 1 failed"
    exit 1
fi

echo ""
echo "✓ Step 1 complete"
echo ""

# Step 2: Network Construction
echo "==============================================================================="
echo "STEP 2/4: Network Construction"
echo "==============================================================================="
echo "  - Choosing soft-thresholding power"
echo "  - Building co-expression network"
echo "  - Detecting gene modules"
echo ""
echo "  This step may take 10-15 minutes..."
echo ""

Rscript 02_wgcna_network_construction.R

if [ $? -ne 0 ]; then
    echo "ERROR: Step 2 failed"
    exit 1
fi

echo ""
echo "✓ Step 2 complete"
echo ""

# Step 3: Module-Trait Association
echo "==============================================================================="
echo "STEP 3/4: Module-Trait Association"
echo "==============================================================================="
echo "  - Correlating modules with atherosclerosis"
echo "  - Calculating gene significance"
echo "  - Calculating module membership"
echo ""

Rscript 03_module_trait_association.R

if [ $? -ne 0 ]; then
    echo "ERROR: Step 3 failed"
    exit 1
fi

echo ""
echo "✓ Step 3 complete"
echo ""

# Step 4: Hub Gene Identification
echo "==============================================================================="
echo "STEP 4/4: Hub Gene Identification"
echo "==============================================================================="
echo "  - Calculating gene connectivity"
echo "  - Identifying hub genes"
echo "  - Generating benchmark output"
echo ""

Rscript 04_hub_gene_identification.R

if [ $? -ne 0 ]; then
    echo "ERROR: Step 4 failed"
    exit 1
fi

echo ""
echo "✓ Step 4 complete"
echo ""

# Summary
echo "==============================================================================="
echo "PIPELINE COMPLETE!"
echo "==============================================================================="
echo ""
echo "Output files generated:"
echo ""
echo "Data files:"
echo "  ✓ gene_level_counts.csv"
echo "  ✓ normalized_expression.csv"
echo "  ✓ sample_metadata.csv"
echo ""
echo "Network results:"
echo "  ✓ module_assignments.csv"
echo "  ✓ module_eigengenes.csv"
echo "  ✓ wgcna_network.RData"
echo ""
echo "Module-trait analysis:"
echo "  ✓ module_trait_correlation.csv"
echo "  ✓ gene_significance.csv"
echo "  ✓ module_membership.csv"
echo ""
echo "Final results:"
echo "  ✓ wgcna_all_results.csv"
echo "  ★ wgcna_significant_genes.csv (BENCHMARK OUTPUT)"
echo "  ✓ hub_genes.csv"
echo ""
echo "Visualizations:"
echo "  ✓ sample_clustering.pdf"
echo "  ✓ network_topology_analysis.pdf"
echo "  ✓ module_dendrogram.pdf"
echo "  ✓ module_trait_heatmap.pdf"
echo "  ✓ hub_gene_scatterplots.pdf"
echo ""
echo "==============================================================================="
echo "Next steps:"
echo "  1. Review wgcna_significant_genes.csv for benchmark comparison"
echo "  2. Check PDF visualizations for quality control"
echo "  3. Compare with other methods (DESeq2, LIMMA, MAGMA)"
echo "==============================================================================="
