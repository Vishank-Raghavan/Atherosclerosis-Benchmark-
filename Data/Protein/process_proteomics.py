#!/usr/bin/env python3

import pandas as pd
import numpy as np

# 1. Load the dataset
# We use low_memory=False because these reports can be huge
file_path = 'DIA_Report_BGS_Factory_Report.tsv'
print("Loading proteomics report...")
df = pd.read_csv(file_path, sep = "\t", low_memory=False)

# 2. Filter for High-Confidence Hits
# Standard practice: Keep only proteins with a Q-value (False Discovery Rate) < 1%
# If PG.Qvalue is missing, we might skip this or use another metric, but it's usually there.
if 'PG.Qvalue' in df.columns:
    print("Filtering for Q-value < 0.01...")
    df_clean = df[df['PG.Qvalue'] < 0.01]
else:
    print("Warning: PG.Qvalue column not found. Skipping quality filtering.")
    df_clean = df

# 3. Select only the columns we need for the matrix
# We need: Protein ID, Sample ID, and Quantity
target_cols = ['PG.ProteinAccessions', 'R.FileName', 'PG.Quantity']
df_subset = df_clean[target_cols]

# 4. Pivot the Data (Long -> Wide)
# Index = Proteins
# Columns = Samples
# Values = Quantity
print("Pivoting data to create expression matrix...")
# We use mean() aggregation just in case there are duplicate entries for a protein in one run
protein_matrix = df_subset.pivot_table(
    index='PG.ProteinAccessions', 
    columns='R.FileName', 
    values='PG.Quantity', 
    aggfunc='mean'
)

# 5. Log-Transform (Standard for Proteomics)
# Proteomics data is usually log-normal. Most tools (DIABLO, Limma) expect log-transformed data.
# We add a small constant (+1) to avoid log(0) errors if there are zeros.
protein_matrix_log = np.log2(protein_matrix + 1)

# 6. Save the final matrix
output_file = 'standardized_proteomics_matrix.csv'
protein_matrix_log.to_csv(output_file)

print(f"Success! Processed matrix saved to {output_file}")
print(f"Matrix Dimensions: {protein_matrix_log.shape[0]} Proteins x {protein_matrix_log.shape[1]} Samples")
