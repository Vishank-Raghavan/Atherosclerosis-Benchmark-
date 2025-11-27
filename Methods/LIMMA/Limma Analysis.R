# 1. Load Libraries
# If missing: BiocManager::install("limma")
library(limma)
library(tidyverse)

# 2. Load Data
# Read the log-transformed protein matrix created by your Python script
# format: Rows = Proteins, Cols = Samples
protein_data <- read.csv("standardized_proteomics_matrix.csv", row.names = 1)

# Check for missing values
# If >50% of samples are missing a protein, remove that protein
protein_data <- protein_data[rowMeans(is.na(protein_data)) < 0.5, ]

# 3. Create Metadata / Design Matrix
# We parse the condition from the column names (A1_DIA, A2_DIA, H1_DIA, H2_DIA, C1_DIA,...)
# The first letter indicates the group: A=Atheroma, H=Healthy, C=Complicated

sample_names <- colnames(protein_data)
# Extract the first character of each column name
group_codes <- substr(sample_names, 1, 1) 

# Create a mapping frame
meta_df <- data.frame(Sample = sample_names, Code = group_codes)

# Map single letters to full condition names
meta_df$Condition <- case_when(
  meta_df$Code == "A" ~ "Atheroma",
  meta_df$Code == "H" ~ "Healthy",
  meta_df$Code == "C" ~ "Complicated",
  TRUE ~ "Unknown"
)

# Filter to keep only Atheroma and Healthy for the benchmark comparison
# We remove "Complicated" (C) and any "Unknown" samples
keep_samples <- meta_df$Sample[meta_df$Condition %in% c("Atheroma", "Healthy")]
protein_data_clean <- protein_data[, keep_samples]
meta_clean <- meta_df[meta_df$Sample %in% keep_samples, ]

# Check if we have enough samples
if (nrow(meta_clean) < 4) {
  stop("Error: Too few samples remaining after filtering. Check column names.")
}

# Define groups for Limma
group_factor <- factor(meta_clean$Condition, levels = c("Healthy", "Atheroma"))

# 4. Create Design Matrix
# This tells Limma how to model the data (0 + group means we fit a mean for each group)
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

# 5. Fit Linear Model
# lmFit fits a linear model to every protein
fit <- lmFit(protein_data_clean, design)

# 6. Define Contrast (Comparison)
# We want to see what is different in Atheroma compared to Healthy
contrast_matrix <- makeContrasts(Diff = Atheroma - Healthy, levels = design)

# 7. Compute Differential Expression
# eBayes borrows information across proteins to stabilize variance (crucial for small N)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 8. Get Results
# topTable extracts the statistics
# number=Inf ensures we get ALL proteins, not just the top 10
results <- topTable(fit2, coef = "Diff", number = Inf)

# 9. Filter for Significance
# Standard threshold: Adjusted P-value < 0.05
sig_proteins <- subset(results, adj.P.Val < 0.05)

# Sort by significance (lowest P-value first)
sig_proteins <- sig_proteins[order(sig_proteins$adj.P.Val), ]

# 10. Save Results

# Set 1: Strict FDR (High Confidence)
write.csv(sig_proteins, "limma_significant_proteins_FDR05.csv")

# Set 2: Nominal P-value (Suggestive / Silver Standard)
# We accept P.Value < 0.05 (unadjusted)
suggestive_proteins <- subset(results, P.Value < 0.05)
suggestive_proteins <- suggestive_proteins[order(suggestive_proteins$P.Value), ]

write.csv(suggestive_proteins, "limma_suggestive_proteins_P05.csv")

print(paste("Total proteins tested:", nrow(results)))
print(paste("High-Confidence Targets (FDR < 0.05):", nrow(sig_proteins)))
print(paste("Suggestive Targets (P < 0.05):", nrow(suggestive_proteins)))
