# FUSION Analysis
This is a Transcriptome-Wide Association Study (TWAS) using the FUSION framework to identify genes whose genetically regulated expression in coronary artery tissue is associated with Coronary Artery Disease.

## Run

```bash
FUSION.assoc_test.R
fusion.sh
```

## Outputs

- `Results.zip` folder - Raw results obtained from running fusion.sh
- `All_results.tsv` - Appended version of all the results from the Results folder
- `filter_results.R` - Filtered out the targets with a FDR value of less than 0.05
- `Significant_TWAS_Targets_FRD05.csv` - All the significant targets ranked from the most significant to the least significant 

## Method

This involved harmonizing GWAS summary statistics with GTEx v8 expression weights and the 1000 Genomes LD reference panel to overcome SNP mismatches, followed by executing the association test across all 22 chromosomes. Finally, the results were aggregated to filter for the most significant targets using False Discovery Rate (FDR) correction, creating a high-confidence gene list for benchmark comparison.

Reference: Gusev et al. “Integrative approaches for large-scale transcriptome-wide association studies” 2016 Nature Genetics

