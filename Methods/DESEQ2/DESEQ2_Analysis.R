# 1. Load Libraries
# If you don't have these, run: install.packages("BiocManager"); BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)

# 2. Load Data
# Read the counts matrix
counts_data <- read.csv("rna_seq_counts_matrix.csv", row.names = 1)

# Round counts to integers (Salmon outputs estimated floats)
counts_data <- round(counts_data)

# 3. Load Metadata
# Format: SampleID, Condition (e.g., Control, Disease)
# We use SraRunTable.csv which has a 'Run' column (SampleID) and 'treatment' column
raw_metadata <- read.csv("SraRunTable.csv")

# Clean the metadata
# We need to handle specific patterns like "atheroma replicate X" vs "8h treated replicate X"
col_data <- raw_metadata %>%
  select(Run, treatment) %>%
  # Create a clean 'Condition' column based on the text in 'treatment'
  mutate(Condition = case_when(
    grepl("atheroma", treatment) ~ "atheroma",
    grepl("healthy", treatment) ~ "healthy",
    grepl("complicated", treatment) ~ "complicated",
    TRUE ~ "exclude" # Label everything else (like "8h treated") as exclude
  )) %>%
  # Filter to keep only the relevant biological groups for the benchmark
  filter(Condition %in% c("atheroma", "healthy")) %>%
  column_to_rownames("Run") %>%
  select(Condition)

# Filter counts data to match the filtered metadata
# This automatically drops the "exclude" samples (8h, 16h, etc.) from the counts matrix
common_samples <- intersect(rownames(col_data), colnames(counts_data))
counts_data <- counts_data[, common_samples]
col_data <- col_data[common_samples, , drop = FALSE]

# Check alignment
print(paste("Number of samples after filtering:", ncol(counts_data)))
all(rownames(col_data) == colnames(counts_data)) # Should be TRUE

# 4. Create DESeqDataSet
# design = ~ Condition tells DESeq2 to test for differences based on the "Condition" column
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ Condition)

# 5. Run Analysis
# This single function runs normalization, dispersion estimation, and testing
dds <- DESeq(dds)

# 6. Get Results
# Contrast defines the comparison: "Condition", "Disease", "Control"
# We compare "atheroma" (Disease) vs "healthy" (Control)
# This means positive Log2FoldChange = Upregulated in Atheroma
res <- results(dds, contrast=c("Condition", "atheroma", "healthy"))

# 7. Filter and Save
# We use padj (FDR adjusted p-value) < 0.05
summary(res)
res_sig <- subset(res, padj < 0.05)

# Sort by most significant
res_sig <- res_sig[order(res_sig$padj), ]

# Save to CSV
write.csv(as.data.frame(res_sig), "deseq2_significant_genes.csv")

print(paste("Found", nrow(res_sig), "significant genes."))

