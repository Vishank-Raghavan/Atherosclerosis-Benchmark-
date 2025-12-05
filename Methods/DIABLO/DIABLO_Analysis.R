#!/usr/bin/env Rscript
# DIABLO Multi-Omics Integration
# Singh et al. (2019) Bioinformatics 35(17):3055-3062

# Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  options(repos = "https://cloud.r-project.org")
  install.packages("BiocManager")
}
if (!requireNamespace("mixOmics", quietly = TRUE)) {
  BiocManager::install("mixOmics")
}
library(mixOmics)

# Load data
protein <- read.csv("standardized_proteomics_matrix.csv", row.names = 1)
sra <- read.csv("SraRunTable.csv")

# Create phenotype mapping
protein_map <- data.frame(
  Sample = colnames(protein),
  Group = c(rep("Disease", 5), rep("Control", 10))
)

# Filter for significant proteins
#limma <- read.csv("../LIMMA/limma_significant_proteins_FDR05.csv", row.names = 1)
#protein_sig <- protein_sig[rownames(limma), ]

# Remove low variance features
vars <- apply(protein, 1, var, na.rm = TRUE)
protein <- protein[vars > 0.01 & !is.na(vars), ]

# Prepare data
X <- list(protein = t(protein))
Y <- factor(protein_map$Group, levels = c("Control", "Disease"))

# Design matrix
design <- matrix(0, 2, 2, dimnames = list(c("protein", "Y"), c("protein", "Y")))
design["protein", "Y"] <- design["Y", "protein"] <- 1

# Build model
model <- block.splsda(X = X, Y = Y, ncomp = 2, keepX = list(protein = c(10, 5)), design = design)

# Visualize
pdf("diablo_plots.pdf", width = 10, height = 8)
plotIndiv(model, legend = TRUE, title = "DIABLO Sample Plot")
plotVar(model, title = "Variable Plot")
plotLoadings(model, comp = 1, contrib = "max")
dev.off()

# Save selected features  
loadings <- model$loadings$protein[, 1]
selected_idx <- loadings != 0
selected <- data.frame(
  Protein = names(loadings)[selected_idx],
  Loading = loadings[selected_idx]
)
selected <- selected[order(abs(selected$Loading), decreasing = TRUE), ]
write.csv(selected, "diablo_selected_proteins.csv", row.names = FALSE)

cat("\nDIABLO analysis complete\n")
cat("Samples:", nrow(X$protein), "\n")
cat("Features:", ncol(X$protein), "\n")
cat("Selected:", nrow(selected), "proteins\n")
cat("\nOutputs:\n")
cat("  - diablo_plots.pdf\n")
cat("  - diablo_selected_proteins.csv\n")
