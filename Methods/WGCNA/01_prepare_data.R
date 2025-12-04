#!/usr/bin/env Rscript
# ==============================================================================
# WGCNA Pipeline - Step 1: Data Preparation
# ==============================================================================
# Purpose:
#   1. Load raw transcript-level RNA-seq counts
#   2. Aggregate transcripts to gene-level
#   3. Filter low-expression genes
#   4. Normalize using DESeq2's Variance Stabilizing Transformation
#   5. Prepare sample metadata
#
# Input:  ../../Data/RNA/rna_seq_counts_matrix.csv.gz
# Output: gene_level_counts.csv
#         normalized_expression.csv
#         sample_metadata.csv
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

# 2. Load Transcript-Level Count Data
# ------------------------------------------------------------------------------
cat("Loading transcript-level count matrix...\n")
counts_file <- "../../Data/RNA/rna_seq_counts_matrix.csv.gz"
transcript_counts <- read.csv(gzfile(counts_file), row.names = 1, check.names = FALSE)

cat(sprintf("Loaded %d transcripts across %d samples\n",
            nrow(transcript_counts), ncol(transcript_counts)))

# 3. Aggregate Transcripts to Gene Level
# ------------------------------------------------------------------------------
cat("\nAggregating transcripts to gene level...\n")

# Parse gene symbol from transcript ID
# Format: ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|1379|lncRNA|
# We want position 6 (DDX11L16)

parse_gene_symbol <- function(transcript_id) {
  # Split by pipe delimiter
  parts <- strsplit(transcript_id, "\\|")[[1]]

  # Extract gene symbol (6th element)
  if (length(parts) >= 6) {
    gene_symbol <- parts[6]
    # Handle cases where gene symbol might be empty
    if (gene_symbol == "" || gene_symbol == "-") {
      # Fall back to ENSG ID (2nd element)
      return(parts[2])
    }
    return(gene_symbol)
  } else {
    return(NA)
  }
}

# Extract gene symbols
gene_symbols <- sapply(rownames(transcript_counts), parse_gene_symbol)

# Remove transcripts with NA gene symbols
valid_idx <- !is.na(gene_symbols)
transcript_counts <- transcript_counts[valid_idx, ]
gene_symbols <- gene_symbols[valid_idx]

cat(sprintf("Removed %d transcripts with invalid gene symbols\n", sum(!valid_idx)))

# Aggregate to gene level by summing transcript counts
gene_counts <- transcript_counts %>%
  mutate(Gene = gene_symbols) %>%
  group_by(Gene) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("Gene")

cat(sprintf("Aggregated to %d genes\n", nrow(gene_counts)))

# Save gene-level counts
write.csv(gene_counts, "gene_level_counts.csv")
cat("Saved: gene_level_counts.csv\n")

# 4. Create Sample Metadata
# ------------------------------------------------------------------------------
cat("\nCreating sample metadata...\n")

# Define samples based on DESeq2 analysis
# Atheroma: SRR10687171-175 (5 samples)
# Healthy:  SRR10687176, 186, 197, 208, 209 (5 samples)

atheroma_samples <- paste0("SRR106871", c(71:75))
healthy_samples <- paste0("SRR106871", c(76, 86))
healthy_samples <- c(healthy_samples, paste0("SRR106872", c("08", "09")))
healthy_samples <- c(healthy_samples, "SRR10687197")

# Create metadata dataframe
all_samples <- c(atheroma_samples, healthy_samples)
sample_metadata <- data.frame(
  SampleID = all_samples,
  Condition = c(rep("atheroma", length(atheroma_samples)),
                rep("healthy", length(healthy_samples))),
  ConditionNumeric = c(rep(1, length(atheroma_samples)),
                       rep(0, length(healthy_samples))),
  row.names = all_samples
)

# Filter gene counts to keep only relevant samples
gene_counts_filtered <- gene_counts[, all_samples]

cat(sprintf("Filtered to %d relevant samples (%d atheroma, %d healthy)\n",
            ncol(gene_counts_filtered),
            sum(sample_metadata$Condition == "atheroma"),
            sum(sample_metadata$Condition == "healthy")))

# Save metadata
write.csv(sample_metadata, "sample_metadata.csv")
cat("Saved: sample_metadata.csv\n")

# 5. Filter Low-Expression Genes
# ------------------------------------------------------------------------------
cat("\nFiltering low-expression genes...\n")

# Keep genes with at least 10 counts in at least 5 samples (50% of samples)
min_count <- 10
min_samples <- 5

keep_genes <- rowSums(gene_counts_filtered >= min_count) >= min_samples
gene_counts_filtered <- gene_counts_filtered[keep_genes, ]

cat(sprintf("Kept %d genes passing filter (≥%d counts in ≥%d samples)\n",
            nrow(gene_counts_filtered), min_count, min_samples))
cat(sprintf("Removed %d low-expression genes\n", sum(!keep_genes)))

# 6. Normalize Using DESeq2 VST
# ------------------------------------------------------------------------------
cat("\nNormalizing data using Variance Stabilizing Transformation...\n")

# Round counts to integers (required by DESeq2)
gene_counts_int <- round(gene_counts_filtered)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = gene_counts_int,
  colData = sample_metadata,
  design = ~ Condition
)

# Apply Variance Stabilizing Transformation
# This produces approximately normal data suitable for WGCNA
vst_data <- vst(dds, blind = TRUE)
normalized_expr <- assay(vst_data)

# Transpose for WGCNA (rows = samples, columns = genes)
normalized_expr_t <- t(normalized_expr)

# Save normalized expression
write.csv(normalized_expr_t, "normalized_expression.csv")
cat("Saved: normalized_expression.csv\n")

# 7. Summary Statistics
# ------------------------------------------------------------------------------
cat("\n" ,rep("=", 70), "\n", sep = "")
cat("DATA PREPARATION COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("Final dataset: %d samples x %d genes\n",
            nrow(normalized_expr_t), ncol(normalized_expr_t)))
cat(sprintf("  Atheroma samples: %d\n", sum(sample_metadata$Condition == "atheroma")))
cat(sprintf("  Healthy samples:  %d\n", sum(sample_metadata$Condition == "healthy")))
cat("\nOutput files:\n")
cat("  - gene_level_counts.csv\n")
cat("  - sample_metadata.csv\n")
cat("  - normalized_expression.csv (ready for WGCNA)\n")
cat(rep("=", 70), "\n", sep = "")
