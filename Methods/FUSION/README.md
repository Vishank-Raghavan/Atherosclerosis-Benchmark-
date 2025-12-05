# FUSION Analysis


## Run

```bash
fusion.sh
```

## Outputs

- `Results` folder - Raw results obtained from running fusion.sh
- `All_results.tsv` - Appended version of all the results from the Results folder
- `filter_results.R` - Filtered out the targets with a FDR value of less than 0.05
- `Significant_TWAS_Targets_FRD05.csv` - All the significant targets ranked from the most significant to the least significant 

## Method

FUSION builds predictive models of the genetic component of a functional/molecular phenotype and predicts and tests that component for association with disease using GWAS summary statistics.

Reference: Gusev et al. “Integrative approaches for large-scale transcriptome-wide association studies” 2016 Nature Genetics

