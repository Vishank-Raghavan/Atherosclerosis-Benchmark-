#!/usr/bin/env Rscript
# ==============================================================================
# WGCNA Pipeline - Step 4: Hub Gene Identification & Results Export
# ==============================================================================
# Purpose:
#   1. Integrate gene significance and module membership data
#   2. Identify hub genes in significant modules
#   3. Calculate gene connectivity measures
#   4. Generate final benchmark output (comparable to DESeq2, LIMMA, MAGMA)
#   5. Visualize hub genes
#
# Input:  wgcna_network.RData
#         gene_significance.csv
#         module_membership.csv
#         module_trait_correlation.csv
#         module_assignments.csv
# Output: wgcna_all_results.csv (all genes with metrics)
#         wgcna_significant_genes.csv (benchmark output)
#         hub_genes.csv (top hub genes per module)
#         hub_gene_scatterplots.pdf
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
cat("Loading all analysis results...\n")
load("wgcna_network.RData")

gene_significance <- read.csv("gene_significance.csv")
module_membership <- read.csv("module_membership.csv")
module_trait <- read.csv("module_trait_correlation.csv")
module_assignments <- read.csv("module_assignments.csv")

cat(sprintf("Loaded data for %d genes\n", nrow(gene_significance)))

# 3. Calculate Gene Connectivity
# ------------------------------------------------------------------------------
cat("\nCalculating gene connectivity measures...\n")

# Calculate adjacency matrix
adjacency <- adjacency(datExpr, power = softPower, type = "signed")

# Total connectivity (degree)
kTotal <- rowSums(adjacency) - 1  # Subtract 1 to exclude self-connection

# Within-module connectivity
kWithin <- rep(NA, ncol(datExpr))
names(kWithin) <- colnames(datExpr)

for (module in unique(moduleColors)) {
  module_genes <- colnames(datExpr)[moduleColors == module]
  if (length(module_genes) > 1) {
    kWithin[module_genes] <- rowSums(adjacency[module_genes, module_genes]) - 1
  }
}

# Create connectivity dataframe
connectivity <- data.frame(
  Gene = names(kTotal),
  kTotal = kTotal,
  kWithin = kWithin
)

cat("Connectivity calculated.\n")

# 4. Integrate All Gene-Level Data
# ------------------------------------------------------------------------------
cat("\nIntegrating all gene-level metrics...\n")

# Merge all dataframes
all_results <- gene_significance %>%
  left_join(module_membership, by = "Gene") %>%
  left_join(connectivity, by = "Gene") %>%
  left_join(module_trait %>% select(Module, Correlation, Pvalue) %>%
              rename(Module_Trait_Cor = Correlation,
                     Module_Trait_Pval = Pvalue),
            by = "Module")

# Sort by gene significance p-value
all_results <- all_results %>%
  arrange(GS_pvalue) %>%
  select(Gene, Module, GS, GS_pvalue, MM, MM_pvalue,
         kTotal, kWithin, Module_Trait_Cor, Module_Trait_Pval)

# Save all results
write.csv(all_results, "wgcna_all_results.csv", row.names = FALSE)
cat("Saved: wgcna_all_results.csv\n")

# 5. Identify Significant Genes (Benchmark Output)
# ------------------------------------------------------------------------------
cat("\nIdentifying significant genes for benchmark...\n")

# Criteria for significance (similar to other methods):
# 1. Gene Significance p-value < 0.05 (associated with atherosclerosis)
# 2. Module Membership > 0.7 (strong module member)
# 3. Exclude grey module (unassigned genes)

significant_genes <- all_results %>%
  filter(GS_pvalue < 0.05,
         MM > 0.7,
         Module != "grey") %>%
  arrange(GS_pvalue)

cat(sprintf("Identified %d significant genes\n", nrow(significant_genes)))

# Save significant genes (benchmark output)
write.csv(significant_genes, "wgcna_significant_genes.csv", row.names = FALSE)
cat("Saved: wgcna_significant_genes.csv (BENCHMARK OUTPUT)\n")

# 6. Identify Hub Genes per Module
# ------------------------------------------------------------------------------
cat("\nIdentifying hub genes in each module...\n")

# Get modules with significant trait association (P < 0.05)
sig_modules <- module_trait %>%
  filter(Pvalue < 0.05) %>%
  pull(Module) %>%
  gsub("ME", "", .)

if (length(sig_modules) == 0) {
  cat("No significant modules found. Using top 3 modules by correlation.\n")
  sig_modules <- module_trait %>%
    arrange(Pvalue) %>%
    head(3) %>%
    pull(Module) %>%
    gsub("ME", "", .)
}

cat(sprintf("Analyzing %d significant modules: %s\n",
            length(sig_modules), paste(sig_modules, collapse = ", ")))

# For each significant module, identify top hub genes
hub_genes_list <- list()

for (mod in sig_modules) {
  module_genes <- all_results %>%
    filter(Module == mod) %>%
    arrange(desc(MM), desc(kWithin)) %>%
    head(20)  # Top 20 hub genes per module

  module_genes$Rank <- 1:nrow(module_genes)

  hub_genes_list[[mod]] <- module_genes
}

hub_genes <- do.call(rbind, hub_genes_list)
write.csv(hub_genes, "hub_genes.csv", row.names = FALSE)
cat(sprintf("Saved: hub_genes.csv (%d genes)\n", nrow(hub_genes)))

# 7. Visualize Hub Genes (GS vs MM)
# ------------------------------------------------------------------------------
cat("\nGenerating hub gene scatterplots...\n")

pdf("hub_gene_scatterplots.pdf", width = 10, height = 8)

for (mod in sig_modules) {
  module_data <- all_results %>% filter(Module == mod)

  if (nrow(module_data) > 0) {
    # Get module-trait correlation for title
    mod_cor <- module_trait %>%
      filter(grepl(mod, Module)) %>%
      pull(Correlation)

    mod_pval <- module_trait %>%
      filter(grepl(mod, Module)) %>%
      pull(Pvalue)

    # Create scatterplot
    plot(module_data$MM, module_data$GS,
         xlab = paste("Module Membership in", mod, "module"),
         ylab = "Gene Significance for Atherosclerosis",
         main = sprintf("%s Module\n(Module-Trait Cor=%.3f, P=%.2e)",
                        mod, mod_cor, mod_pval),
         pch = 20, col = mod, cex = 1.2)

    # Add reference lines
    abline(h = 0, col = "grey", lty = 2)
    abline(v = 0.7, col = "red", lty = 2)
    abline(h = c(-0.3, 0.3), col = "blue", lty = 2)

    # Add correlation
    cor_test <- cor.test(module_data$MM, module_data$GS)
    text(min(module_data$MM, na.rm = TRUE),
         max(module_data$GS, na.rm = TRUE),
         sprintf("cor = %.3f\nP = %.2e", cor_test$estimate, cor_test$p.value),
         pos = 4, cex = 0.9)
  }
}

dev.off()
cat("Saved: hub_gene_scatterplots.pdf\n")

# 8. Summary Statistics
# ------------------------------------------------------------------------------
cat("\n", rep("=", 70), "\n", sep = "")
cat("HUB GENE IDENTIFICATION COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nSIGNIFICANT GENE SUMMARY:\n")
cat(sprintf("  Total genes analyzed: %d\n", nrow(all_results)))
cat(sprintf("  Genes with GS P-value < 0.05: %d\n",
            sum(all_results$GS_pvalue < 0.05, na.rm = TRUE)))
cat(sprintf("  Genes with MM > 0.7: %d\n",
            sum(all_results$MM > 0.7, na.rm = TRUE)))
cat(sprintf("  SIGNIFICANT GENES (GS P<0.05 & MM>0.7, not grey): %d\n",
            nrow(significant_genes)))

cat("\nMODULE BREAKDOWN OF SIGNIFICANT GENES:\n")
if (nrow(significant_genes) > 0) {
  module_breakdown <- significant_genes %>%
    group_by(Module) %>%
    summarise(Count = n(),
              Mean_GS = mean(abs(GS)),
              Mean_MM = mean(MM)) %>%
    arrange(desc(Count))
  print(module_breakdown, n = Inf)
} else {
  cat("  No significant genes found.\n")
  cat("\n  RELAXED CRITERIA (GS P<0.05 OR MM>0.8):\n")
  relaxed <- all_results %>%
    filter((GS_pvalue < 0.05 | MM > 0.8) & Module != "grey")
  cat(sprintf("  %d genes meet relaxed criteria\n", nrow(relaxed)))
}

cat("\nTOP 10 SIGNIFICANT GENES:\n")
if (nrow(significant_genes) > 0) {
  top_genes <- significant_genes %>%
    head(10) %>%
    select(Gene, Module, GS, GS_pvalue, MM, kTotal)
  print(top_genes, row.names = FALSE)
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("BENCHMARK OUTPUT READY\n")
cat(rep("=", 70), "\n", sep = "")
cat("Compare with other methods:\n")
cat("  - MAGMA:  Methods/MAGMA/significant_genes.txt\n")
cat("  - DESeq2: Methods/DESEQ2/deseq2_significant_genes.csv\n")
cat("  - LIMMA:  Methods/LIMMA/limma_significant_proteins_FDR05.csv\n")
cat("  - WGCNA:  Methods/WGCNA/wgcna_significant_genes.csv\n")
cat("\nOutput files:\n")
cat("  - wgcna_all_results.csv (all genes with metrics)\n")
cat("  - wgcna_significant_genes.csv (benchmark output)\n")
cat("  - hub_genes.csv (top hub genes per module)\n")
cat("  - hub_gene_scatterplots.pdf\n")
cat(rep("=", 70), "\n", sep = "")
