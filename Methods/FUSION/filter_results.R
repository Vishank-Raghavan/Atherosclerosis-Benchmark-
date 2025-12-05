
# Load the combined results
results <- read.table("All_Results.tsv", header=TRUE, stringsAsFactors=FALSE)

# 1. Remove the MHC Region (Standard QC)
# The Major Histocompatibility Complex (Chr 6: 25Mb - 35Mb) is very complex 
# and often produces false positives due to long-range LD. It is standard to remove it.
results <- results[ !(results$CHR == 6 & results$P0 > 25e6 & results$P1 < 35e6), ]

# 2. Calculate FDR (False Discovery Rate)
# This creates a new column "FDR" adjusting the raw TWAS.P value
results$FDR <- p.adjust(results$TWAS.P, method = "fdr")

# 3. Filter for Significance (Standard is FDR < 0.05)
# You can also use Bonferroni correction for a stricter list:
# threshold <- 0.05 / nrow(results)
sig_genes <- results[ results$FDR < 0.05, ]

# 4. Sort by Significance (most significant at the top)
sig_genes <- sig_genes[ order(sig_genes$FDR), ]

# 5. Save the Significant Targets
write.table(sig_genes, "Significant_TWAS_Targets_FDR05.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Print the number of significant genes found
cat("Number of significant targets found:", nrow(sig_genes), "\n")
