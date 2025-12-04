# WGCNA Network-Based Analysis for Atherosclerosis

This directory contains a complete **Weighted Gene Co-expression Network Analysis (WGCNA)** pipeline for identifying atherosclerosis-associated genes through network-based methods.

## Overview

WGCNA identifies **modules** (clusters) of co-expressed genes and correlates them with disease phenotypes. Unlike traditional differential expression (DESeq2), WGCNA:
- Captures **biological context** - genes work in pathways, not isolation
- Identifies **hub genes** - key regulatory genes central to disease modules
- More **robust with small sample sizes** - aggregates signal across co-expressed genes
- Discovers **network rewiring** - changes in gene relationships, not just expression

---

## Pipeline Architecture

```
01_prepare_data.R
    ↓
    [gene_level_counts.csv, normalized_expression.csv, sample_metadata.csv]
    ↓
02_wgcna_network_construction.R
    ↓
    [module_assignments.csv, module_eigengenes.csv, wgcna_network.RData]
    ↓
03_module_trait_association.R
    ↓
    [module_trait_correlation.csv, gene_significance.csv, module_membership.csv]
    ↓
04_hub_gene_identification.R
    ↓
    [wgcna_significant_genes.csv ★ BENCHMARK OUTPUT ★]
```

---

## Input Data

**RNA-seq data:**
- `../../Data/RNA/rna_seq_counts_matrix.csv.gz`
  - 517,039 transcripts × 39 samples
  - Salmon quantification output

**Samples used:**
- **Atheroma (n=5):** SRR10687171-175
- **Healthy (n=5):** SRR10687176, 186, 197, 208, 209
- Tissue: Carotid artery from endarterectomy

---

## Pipeline Steps

### **Step 1: Data Preparation** (`01_prepare_data.R`)

**Purpose:** Aggregate transcripts to genes, filter, and normalize

**Operations:**
1. Load transcript-level counts (517K transcripts)
2. Parse gene symbols from transcript IDs
   - Format: `ENST...|ENSG...|...|GENE_SYMBOL|...|biotype|`
3. Aggregate to gene-level by summing transcript counts
4. Filter: Keep genes with ≥10 counts in ≥5 samples
5. Normalize: Variance Stabilizing Transformation (DESeq2::vst)
6. Create sample metadata (atheroma=1, healthy=0)

**Output:**
- `gene_level_counts.csv` - Raw gene-level counts
- `normalized_expression.csv` - VST-normalized, **ready for WGCNA**
- `sample_metadata.csv` - Sample phenotypes

**Expected:** ~15-20K genes after filtering

---

### **Step 2: Network Construction** (`02_wgcna_network_construction.R`)

**Purpose:** Build co-expression network and detect modules

**Operations:**
1. Check for sample outliers via hierarchical clustering
2. Choose soft-thresholding power (β)
   - Tests β=1-20, selects value achieving scale-free topology (R²>0.85)
3. Build signed co-expression network
   - Adjacency matrix: A_ij = |0.5 × (1 + cor(gene_i, gene_j))|^β
   - Topological Overlap Matrix (TOM): Measures network connectivity
4. Detect modules via hierarchical clustering
   - Min module size: 30 genes
   - Merge similar modules (cut height: 0.25)
5. Calculate module eigengenes (1st principal component)

**Output:**
- `sample_clustering.pdf` - Sample dendrogram with condition colors
- `network_topology_analysis.pdf` - Scale-free topology fit
- `module_dendrogram.pdf` - Gene clustering with module colors
- `module_assignments.csv` - Gene → Module mapping
- `module_eigengenes.csv` - Module expression profiles
- `wgcna_network.RData` - Complete workspace

**Expected:** 5-15 modules

---

### **Step 3: Module-Trait Association** (`03_module_trait_association.R`)

**Purpose:** Identify modules correlated with atherosclerosis

**Operations:**
1. Correlate module eigengenes with phenotype (atheroma vs healthy)
2. Calculate p-values for module-trait correlations
3. Calculate **Gene Significance (GS)**: correlation of each gene with trait
4. Calculate **Module Membership (MM)**: correlation of each gene with its module eigengene
5. Visualize module-trait relationships

**Output:**
- `module_trait_correlation.csv` - Modules ranked by disease association
- `module_trait_heatmap.pdf` - Correlation heatmap
- `gene_significance.csv` - Gene-level trait associations
- `module_membership.csv` - Gene-level module memberships

**Key Metrics:**
- **GS (Gene Significance):** How much a gene correlates with disease
- **GS p-value:** Statistical significance of GS
- **MM (Module Membership):** How central a gene is to its module
- **MM p-value:** Statistical significance of MM

---

### **Step 4: Hub Gene Identification** (`04_hub_gene_identification.R`)

**Purpose:** Extract significant genes for benchmark comparison

**Operations:**
1. Calculate gene connectivity (kTotal, kWithin)
2. Integrate all metrics: GS, MM, connectivity, module-trait correlation
3. **Filter significant genes:**
   - GS p-value < 0.05 (disease-associated)
   - MM > 0.7 (strong module member)
   - Exclude grey module (unassigned genes)
4. Identify hub genes in each significant module
5. Visualize GS vs MM scatterplots

**Output:**
- `wgcna_all_results.csv` - All genes with full metrics
- **`wgcna_significant_genes.csv`** ★ **BENCHMARK OUTPUT** ★
- `hub_genes.csv` - Top 20 hub genes per significant module
- `hub_gene_scatterplots.pdf` - GS vs MM plots per module

**Benchmark Format:**
```csv
Gene, Module, GS, GS_pvalue, MM, MM_pvalue, kTotal, kWithin, Module_Trait_Cor, Module_Trait_Pval
```

---

## Running the Pipeline

### **Prerequisites:**

```r
# Install required packages
install.packages("BiocManager")
BiocManager::install(c("WGCNA", "DESeq2"))
install.packages("tidyverse")
```

### **Execution:**

```bash
cd Methods/WGCNA

# Run full pipeline
Rscript 01_prepare_data.R
Rscript 02_wgcna_network_construction.R
Rscript 03_module_trait_association.R
Rscript 04_hub_gene_identification.R
```

**Or run all at once:**
```bash
Rscript 01_prepare_data.R && \
Rscript 02_wgcna_network_construction.R && \
Rscript 03_module_trait_association.R && \
Rscript 04_hub_gene_identification.R
```

**Estimated runtime:** 15-30 minutes (depending on CPU)

---

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Min counts | 10 | Minimum read count threshold |
| Min samples | 5 | Genes must pass threshold in ≥50% samples |
| Soft power (β) | Data-driven | Chosen to achieve scale-free topology |
| Network type | Signed | Preserves direction of correlation |
| Min module size | 30 | Minimum genes per module |
| Merge threshold | 0.25 | Similarity threshold for merging modules |
| GS significance | P < 0.05 | Gene-trait association threshold |
| MM threshold | 0.7 | Module membership threshold |

---

## Interpreting Results

### **Significant Genes:**
Genes in `wgcna_significant_genes.csv` meet **dual criteria:**
1. **Biologically associated** with atherosclerosis (GS p<0.05)
2. **Central to their module** (MM > 0.7)

This is **more stringent** than DESeq2 alone, as genes must also show coordinated behavior.

### **Module Colors:**
- Each module is assigned a color (turquoise, blue, brown, etc.)
- **Grey module** = genes not assigned to any module (excluded from results)
- Module size correlates with biological importance (larger ≠ better)

### **Hub Genes:**
- High MM (>0.8): Central regulators within modules
- High kWithin: Many connections within module
- **Therapeutic targets:** Hub genes often control entire pathways

### **Module-Trait Correlation:**
- **Positive correlation:** Module upregulated in atheroma
- **Negative correlation:** Module downregulated in atheroma
- **P < 0.05:** Statistically significant association

---

## Comparison with Other Methods

| Method | Data Type | Level | Output |
|--------|-----------|-------|--------|
| **MAGMA** | DNA (GWAS) | Gene sets | Population genetic associations |
| **DESeq2** | RNA | Individual genes | Differential expression |
| **LIMMA** | Protein | Individual proteins | Differential abundance |
| **WGCNA** | RNA | Gene modules | Network hubs & modules |

**WGCNA Advantages:**
- Identifies **regulatory hubs** (potential drug targets)
- Captures **pathway-level** effects
- More **robust** with small sample sizes
- Discovers genes DESeq2 might miss (co-expression without differential expression)

**WGCNA Limitations:**
- Requires ≥10 samples (we have 10 ✓)
- Computationally intensive
- Module interpretation requires biological knowledge

---

## Validation Strategy

### **Cross-Method Overlap:**
Compare WGCNA hub genes with:
1. DESeq2 significant genes (DEGs)
2. LIMMA significant proteins
3. MAGMA significant genes
4. Known atherosclerosis genes (literature)

### **Pathway Enrichment:**
Expected enrichments in hub genes:
- Inflammation (IL6, TNF, NFκB pathway)
- Lipid metabolism (APOE, LDLR, PCSK9)
- ECM remodeling (MMPs, collagens)
- Oxidative stress (SOD, GPX)

### **Hub Gene Validation:**
- Check if hub genes are:
  - Known atherosclerosis genes (GWAS catalog)
  - Transcription factors or signaling genes
  - Drugable targets

---

## Troubleshooting

### **Issue:** "No modules reached significance (P < 0.05)"
**Solution:** Check `module_trait_correlation.csv` - modules with P<0.1 may still be biologically relevant. Adjust threshold in Step 4.

### **Issue:** "Too few genes after filtering"
**Solution:** Lower filtering thresholds in `01_prepare_data.R`:
- `min_count <- 5` (instead of 10)
- `min_samples <- 3` (instead of 5)

### **Issue:** "Scale-free topology not achieved"
**Solution:** Common with small N. Accept lower R² (0.7-0.8) or use higher β.

### **Issue:** "Too many small modules"
**Solution:** Increase `minModuleSize` in `02_wgcna_network_construction.R` (e.g., 50 instead of 30)

---

## References

1. **WGCNA Tutorial:** https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
2. **Original Paper:** Langfelder & Horvath (2008). "WGCNA: an R package for weighted correlation network analysis." *BMC Bioinformatics*
3. **Signed Networks:** Langfelder et al. (2008). "Defining clusters from a hierarchical cluster tree." *Bioinformatics*

---

## Contact

For issues or questions about this pipeline, refer to the main project README or WGCNA documentation.

---

**Last Updated:** 2025-12-04
**Pipeline Version:** 1.0
