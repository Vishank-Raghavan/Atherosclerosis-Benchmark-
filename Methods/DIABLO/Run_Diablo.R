# ==============================================================================
# Multi-Omic Integration with DIABLO (mixOmics)
# ==============================================================================

# 1. Install/Load Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("mixOmics", quietly = TRUE))
  BiocManager::install("mixOmics")
if (!requireNamespace("BiocParallel", quietly = TRUE))
  BiocManager::install("BiocParallel")

library(mixOmics)
library(tidyverse)
library(BiocParallel)

# ==============================================================================
# 2. Load Data
# ==============================================================================

# A. Load RNA Data (TPM is preferred for integration over raw counts)
rna_data <- read.csv("rna_seq_tpm_matrix.csv", row.names = 1)
# Log-transform to normalize (log2(x+1))
rna_data <- log2(rna_data + 1)

# B. Load Protein Data (Already log-transformed in your python script)
prot_data <- read.csv("standardized_proteomics_matrix.csv", row.names = 1)

# C. Load Metadata
meta <- read.csv("SraRunTable.csv")

# ==============================================================================
# 3. Data Harmonization
# ==============================================================================

# --- Step 3a: Create Mapping from SRA Metadata ---
meta_mapping <- meta %>%
  select(Run, treatment) %>%
  mutate(
    # Extract condition type (A, H, C)
    Type = case_when(
      grepl("atheroma", treatment) ~ "A",
      grepl("healthy", treatment) ~ "H",
      grepl("complicated", treatment) ~ "C",
      TRUE ~ NA_character_
    ),
    # Extract replicate number
    Rep = str_extract(treatment, "[0-9]+"),
    # Construct Common ID (e.g., A5)
    CommonID = paste0(Type, Rep)
  ) %>%
  filter(!is.na(Type)) %>%
  column_to_rownames("Run")

# --- Step 3b: Update RNA Data with Common IDs ---
X_rna <- t(rna_data)
common_rna_srr <- intersect(rownames(X_rna), rownames(meta_mapping))
X_rna <- X_rna[common_rna_srr, ]
rownames(X_rna) <- meta_mapping[rownames(X_rna), "CommonID"]

# --- Step 3c: Update Protein Data with Common IDs ---
X_prot <- t(prot_data)
clean_prot_names <- gsub("_DIA", "", rownames(X_prot))
rownames(X_prot) <- clean_prot_names

# --- Step 3d: Final Match ---
common_samples <- intersect(rownames(X_rna), rownames(X_prot))

if (length(common_samples) < 2) {
  stop("Error: No matching sample IDs found.")
}

# Subset data to matched samples
X_rna <- X_rna[common_samples, ]
X_prot <- X_prot[common_samples, ]

# Define groups
Y_groups <- substr(common_samples, 1, 1)
Y_groups <- ifelse(Y_groups == "A", "Atheroma", 
                   ifelse(Y_groups == "H", "Healthy", "Other"))
Y <- as.factor(Y_groups)

# Filter for binary comparison
keep_indices <- Y %in% c("Atheroma", "Healthy")
X_rna <- X_rna[keep_indices, ]
X_prot <- X_prot[keep_indices, ]
Y <- droplevels(Y[keep_indices])

print(paste("Data harmonized. Samples:", length(Y)))

# ==============================================================================
# 4. NA Handling & Zero Variance
# ==============================================================================

# Helper function to clean a matrix
clean_matrix <- function(mat, name) {
  # 1. Check for Infinite values
  if (any(is.infinite(mat))) {
    print(paste("Fixing Infinite values in", name))
    mat[is.infinite(mat)] <- NA
  }
  
  # 2. Remove features with > 20% Missing Values
  na_prop <- colMeans(is.na(mat))
  keep_features <- na_prop < 0.2
  mat <- mat[, keep_features]
  print(paste("Removed", sum(!keep_features), "features with >20% NAs in", name))
  
  # 3. Impute remaining NAs
  if (any(is.na(mat))) {
    print(paste("Imputing remaining NAs in", name))
    mat <- apply(mat, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      return(x)
    })
  }
  
  return(mat)
}

# Clean both matrices
X_rna <- clean_matrix(X_rna, "RNA")
X_prot <- clean_matrix(X_prot, "Protein")

# Final Data List
data <- list(mRNA = X_rna, protein = X_prot)

# ==============================================================================
# 5. Parameter Tuning & DIABLO Run (PARALLELIZED FOR WINDOWS)
# ==============================================================================

# Setup Parallel Backend (Cross-Platform)
n_cores <- 24 
print(paste("Setting up parallel processing with", n_cores, "cores..."))

if (.Platform$OS.type == "windows") {
  # Windows uses SnowParam
  BPPARAM <- SnowParam(workers = n_cores)
} else {
  # Linux/Mac uses MulticoreParam
  BPPARAM <- MulticoreParam(workers = n_cores)
}

# Design matrix
design <- matrix(0.1, ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
diag(design) <- 0

# Check MixOmics version for function name
diablo_func <- if(exists("block.splsda")) block.splsda else block.spls

# Tune ncomp
print("Tuning ncomp (Parallel)...")
perf.diablo <- perf(diablo_func(X = data, Y = Y, ncomp = 5, design = design, near.zero.var = TRUE), 
                    validation = 'Mfold', folds = 5, nrepeat = 10, 
                    BPPARAM = BPPARAM) 

optimal_ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
if (is.na(optimal_ncomp)) optimal_ncomp <- 2 
print(paste("Optimal ncomp:", optimal_ncomp))

# Tune keepX
print("Tuning keepX (Parallel)...")
list.keepX <- list(mRNA = c(10, 25, 50, 100), protein = c(10, 25, 50, 100))

tune.diablo <- tune.block.splsda(X = data, Y = Y, ncomp = optimal_ncomp, 
                                 test.keepX = list.keepX, design = design,
                                 validation = 'Mfold', folds = 5, nrepeat = 10,
                                 near.zero.var = TRUE,
                                 BPPARAM = BPPARAM) 

optimal_keepX <- tune.diablo$choice.keepX

# Run Final Model
print("Running Final Model...")
final.diablo <- diablo_func(X = data, Y = Y, ncomp = optimal_ncomp, 
                            keepX = optimal_keepX, design = design,
                            near.zero.var = TRUE)

# ==============================================================================
# 6. Extract and Save
# ==============================================================================

selected_rna <- selectVar(final.diablo, block = 'mRNA', comp = 1)$mRNA$name
selected_prot <- selectVar(final.diablo, block = 'protein', comp = 1)$protein$name

if (optimal_ncomp > 1) {
  selected_rna <- c(selected_rna, selectVar(final.diablo, block = 'mRNA', comp = 2)$mRNA$name)
  selected_prot <- c(selected_prot, selectVar(final.diablo, block = 'protein', comp = 2)$protein$name)
}

results_df <- data.frame(
  ID = c(selected_rna, selected_prot),
  Type = c(rep("RNA", length(selected_rna)), rep("Protein", length(selected_prot))),
  Method = "DIABLO"
)

write.csv(results_df, "diablo_significant_targets.csv", row.names = FALSE)
print("DIABLO analysis complete.")
