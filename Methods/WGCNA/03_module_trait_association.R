#!/usr/bin/env Rscript
# ==============================================================================
# WGCNA Pipeline - Step 3: Module-Trait Association
# ==============================================================================
# Purpose:
#   1. Load network data and module eigengenes
#   2. Correlate module eigengenes with atherosclerosis phenotype
#   3. Calculate gene-level trait significance
#   4. Calculate module membership for each gene
#   5. Visualize module-trait relationships
#
# Input:  wgcna_network.RData
#         sample_metadata.csv
#         module_eigengenes.csv
# Output: module_trait_correlation.csv
#         module_trait_heatmap.pdf
#         gene_significance.csv
#         module_membership.csv
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

# 2. Load Data
# ------------------------------------------------------------------------------
cat("Loading network data...\n")
load("wgcna_network.RData")

# Load sample metadata
sample_metadata <- read.csv("sample_metadata.csv", row.names = 1)

cat(sprintf("Loaded: %d samples, %d genes across %d modules\n",
            nrow(datExpr), ncol(datExpr), length(unique(moduleColors))))

# 3. Prepare Trait Data
# ------------------------------------------------------------------------------
cat("\nPreparing trait data...\n")

# Create trait matrix (atheroma = 1, healthy = 0)
trait <- as.data.frame(sample_metadata$ConditionNumeric)
rownames(trait) <- rownames(sample_metadata)
colnames(trait) <- "Atherosclerosis"

cat("Trait summary:\n")
print(table(sample_metadata$Condition))

# 4. Calculate Module-Trait Correlations
# ------------------------------------------------------------------------------
cat("\nCalculating module-trait correlations...\n")

# Get module eigengenes (drop SampleID column)
MEs_numeric <- MEs[, -1, drop = FALSE]

# Ensure column names are properly formatted
colnames(MEs_numeric) <- paste0("ME", colnames(MEs_numeric))

# Calculate correlations between MEs and traits
moduleTraitCor <- cor(MEs_numeric, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Save results
module_trait_results <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = moduleTraitCor[, 1],
  Pvalue = moduleTraitPvalue[, 1]
)
module_trait_results <- module_trait_results %>%
  arrange(Pvalue) %>%
  mutate(Significant = ifelse(Pvalue < 0.05, "Yes", "No"))

write.csv(module_trait_results, "module_trait_correlation.csv", row.names = FALSE)
cat("Saved: module_trait_correlation.csv\n")

# 5. Display Significant Modules
# ------------------------------------------------------------------------------
cat("\n", rep("=", 70), "\n", sep = "")
cat("MODULE-TRAIT CORRELATION RESULTS\n")
cat(rep("=", 70), "\n", sep = "")
cat("Modules significantly associated with Atherosclerosis (P < 0.05):\n\n")

sig_modules <- module_trait_results %>%
  filter(Pvalue < 0.05) %>%
  arrange(Pvalue)

if (nrow(sig_modules) > 0) {
  print(sig_modules, row.names = FALSE)
  cat(sprintf("\nTotal significant modules: %d\n", nrow(sig_modules)))
} else {
  cat("No modules reached significance threshold (P < 0.05)\n")
  cat("Showing top 5 modules by p-value:\n\n")
  print(head(module_trait_results, 5), row.names = FALSE)
}

# 6. Visualize Module-Trait Relationships
# ------------------------------------------------------------------------------
cat("\nGenerating module-trait heatmap...\n")

# Prepare text labels for heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Create heatmap
pdf("module_trait_heatmap.pdf", width = 6, height = max(8, length(unique(moduleColors)) * 0.3))
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(trait),
               yLabels = rownames(moduleTraitCor),
               ySymbols = gsub("ME", "", rownames(moduleTraitCor)),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))
dev.off()

cat("Saved: module_trait_heatmap.pdf\n")

# 7. Calculate Gene Significance (GS)
# ------------------------------------------------------------------------------
cat("\nCalculating gene-level trait significance...\n")

# Gene Significance = correlation of each gene with the trait
geneTraitSignificance <- as.data.frame(cor(datExpr, trait, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))

# Rename columns
names(geneTraitSignificance) <- "GS.Atherosclerosis"
names(GSPvalue) <- "GS.Pvalue"

# Create gene significance dataframe
gene_significance <- data.frame(
  Gene = colnames(datExpr),
  GS = geneTraitSignificance$GS.Atherosclerosis,
  GS_pvalue = GSPvalue$GS.Pvalue
)

write.csv(gene_significance, "gene_significance.csv", row.names = FALSE)
cat("Saved: gene_significance.csv\n")

# 8. Calculate Module Membership (MM)
# ------------------------------------------------------------------------------
cat("\nCalculating module membership...\n")

# Module Membership = correlation of each gene with its module eigengene
geneModuleMembership <- as.data.frame(cor(datExpr, MEs_numeric, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

# For each gene, extract its MM for its assigned module
module_membership_list <- list()

for (gene in colnames(datExpr)) {
  gene_module <- moduleColors[colnames(datExpr) == gene]
  me_name <- paste0("ME", gene_module)

  if (me_name %in% colnames(geneModuleMembership)) {
    mm_value <- geneModuleMembership[gene, me_name]
    mm_pvalue <- MMPvalue[gene, me_name]
  } else {
    mm_value <- NA
    mm_pvalue <- NA
  }

  module_membership_list[[gene]] <- data.frame(
    Gene = gene,
    Module = gene_module,
    MM = mm_value,
    MM_pvalue = mm_pvalue
  )
}

module_membership <- do.call(rbind, module_membership_list)
write.csv(module_membership, "module_membership.csv", row.names = FALSE)
cat("Saved: module_membership.csv\n")

# 9. Summary Statistics
# ------------------------------------------------------------------------------
cat("\n", rep("=", 70), "\n", sep = "")
cat("GENE SIGNIFICANCE SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("Total genes analyzed: %d\n", nrow(gene_significance)))
cat(sprintf("Genes with significant GS (P < 0.05): %d\n",
            sum(gene_significance$GS_pvalue < 0.05)))
cat(sprintf("Genes with high MM (MM > 0.7): %d\n",
            sum(module_membership$MM > 0.7, na.rm = TRUE)))
cat(sprintf("Genes with both GS P<0.05 AND MM>0.7: %d\n",
            sum(gene_significance$GS_pvalue < 0.05 &
                module_membership$MM > 0.7, na.rm = TRUE)))

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODULE-TRAIT ASSOCIATION COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat("Output files:\n")
cat("  - module_trait_correlation.csv\n")
cat("  - module_trait_heatmap.pdf\n")
cat("  - gene_significance.csv\n")
cat("  - module_membership.csv\n")
cat(rep("=", 70), "\n", sep = "")
