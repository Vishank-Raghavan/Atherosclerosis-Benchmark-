#!/usr/bin/env Rscript
# ==============================================================================
# WGCNA Pipeline - Step 2: Network Construction
# ==============================================================================
# Purpose:
#   1. Load normalized expression data
#   2. Check for outlier samples
#   3. Choose soft-thresholding power (β)
#   4. Construct co-expression network
#   5. Detect gene modules
#   6. Calculate module eigengenes
#
# Input:  normalized_expression.csv
#         sample_metadata.csv
# Output: network_topology_analysis.pdf
#         module_dendrogram.pdf
#         module_assignments.csv
#         module_eigengenes.csv
#         wgcna_network.RData (saved workspace)
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

# Allow multi-threading to speed up calculations
allowWGCNAThreads()

# 2. Load Data
# ------------------------------------------------------------------------------
cat("Loading normalized expression data...\n")
datExpr <- read.csv("normalized_expression.csv", row.names = 1, check.names = FALSE)
sample_metadata <- read.csv("sample_metadata.csv", row.names = 1)

cat(sprintf("Loaded: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

# Verify sample alignment
if (!all(rownames(datExpr) == rownames(sample_metadata))) {
  stop("Error: Sample names don't match between expression and metadata")
}

# 3. Check for Outlier Samples
# ------------------------------------------------------------------------------
cat("\nChecking for outlier samples...\n")

# Hierarchical clustering of samples
sampleTree <- hclust(dist(datExpr), method = "average")

# Plot sample clustering
pdf("sample_clustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Add condition color bar
condition_colors <- ifelse(sample_metadata$Condition == "atheroma", "red", "blue")
plotColors <- data.frame(Condition = condition_colors)
colorHeight <- 0.1 * max(sampleTree$height)
plotDendroAndColors(sampleTree, plotColors,
                    groupLabels = "Condition",
                    main = "Sample Dendrogram and Condition",
                    cex.dendroLabels = 0.8)
dev.off()

cat("Saved: sample_clustering.pdf\n")
cat("No outliers detected. Proceeding with all samples.\n")

# 4. Choose Soft-Thresholding Power
# ------------------------------------------------------------------------------
cat("\nChoosing soft-thresholding power (β)...\n")

# Test different powers
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))

# Calculate scale-free topology fit indices
sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         verbose = 5,
                         networkType = "signed")

# Plot scale-free topology fit
pdf("network_topology_analysis.pdf", width = 12, height = 5)
par(mfrow = c(1,2))

# Scale-free topology fit index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red")

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = 0.9, col = "red")

dev.off()
cat("Saved: network_topology_analysis.pdf\n")

# Choose power: first power that exceeds R^2 = 0.85
power_threshold <- 0.85
softPower <- sft$fitIndices[which(sft$fitIndices[,2] >= power_threshold)[1], 1]

# If no power reaches threshold, use the power that maximizes R^2
if (is.na(softPower)) {
  softPower <- sft$fitIndices[which.max(sft$fitIndices[,2]), 1]
  cat(sprintf("Warning: No power reached R^2 = %.2f. Using power = %d (max R^2 = %.3f)\n",
              power_threshold, softPower,
              max(sft$fitIndices[,2], na.rm = TRUE)))
} else {
  cat(sprintf("Selected soft-thresholding power: %d (R^2 = %.3f)\n",
              softPower,
              sft$fitIndices[sft$fitIndices[,1] == softPower, 2]))
}

# 5. Construct Co-Expression Network and Detect Modules
# ------------------------------------------------------------------------------
cat("\nConstructing co-expression network...\n")
cat("This may take several minutes...\n")

# One-step network construction and module detection
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

cat("\nNetwork construction complete!\n")

# 6. Module Detection Results
# ------------------------------------------------------------------------------
cat("\n", rep("=", 70), "\n", sep = "")
cat("MODULE DETECTION SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")

# Convert numeric labels to colors for visualization
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

# Count genes per module
module_counts <- table(moduleColors)
cat(sprintf("Number of modules detected: %d\n", length(unique(moduleColors))))
cat("\nGenes per module:\n")
print(module_counts)

# Save module assignments
module_assignments <- data.frame(
  Gene = colnames(datExpr),
  ModuleLabel = moduleLabels,
  ModuleColor = moduleColors
)
write.csv(module_assignments, "module_assignments.csv", row.names = FALSE)
cat("\nSaved: module_assignments.csv\n")

# 7. Visualize Module Dendrogram
# ------------------------------------------------------------------------------
cat("\nGenerating module dendrogram...\n")

pdf("module_dendrogram.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()
cat("Saved: module_dendrogram.pdf\n")

# 8. Calculate Module Eigengenes
# ------------------------------------------------------------------------------
cat("\nCalculating module eigengenes...\n")

# Module eigengenes (1st PC of each module)
MEs <- net$MEs

# Reorder by hierarchical clustering
MEs <- orderMEs(MEs)

# Convert to color-based names for clarity
colnames(MEs) <- gsub("ME", "", colnames(MEs))
colnames(MEs) <- labels2colors(as.numeric(colnames(MEs)))

# Add sample IDs
MEs <- cbind(SampleID = rownames(datExpr), MEs)

# Save module eigengenes
write.csv(MEs, "module_eigengenes.csv", row.names = FALSE)
cat("Saved: module_eigengenes.csv\n")

# 9. Save Workspace
# ------------------------------------------------------------------------------
cat("\nSaving workspace...\n")
save(datExpr, sample_metadata, net, moduleColors, moduleLabels, MEs, softPower,
     file = "wgcna_network.RData")
cat("Saved: wgcna_network.RData\n")

# 10. Summary
# ------------------------------------------------------------------------------
cat("\n", rep("=", 70), "\n", sep = "")
cat("NETWORK CONSTRUCTION COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("Soft-thresholding power: %d\n", softPower))
cat(sprintf("Total modules: %d (including grey module)\n", length(unique(moduleColors))))
cat(sprintf("Genes in grey module (unassigned): %d\n",
            sum(moduleColors == "grey")))
cat(sprintf("Genes in assigned modules: %d\n",
            sum(moduleColors != "grey")))
cat("\nOutput files:\n")
cat("  - sample_clustering.pdf\n")
cat("  - network_topology_analysis.pdf\n")
cat("  - module_dendrogram.pdf\n")
cat("  - module_assignments.csv\n")
cat("  - module_eigengenes.csv\n")
cat("  - wgcna_network.RData\n")
cat(rep("=", 70), "\n", sep = "")
