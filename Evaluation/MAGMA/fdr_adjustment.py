#!/usr/bin/env python3

import pandas as pd
import numpy as np

# 1. Load the full MAGMA output
# Adjust filename if needed. MAGMA output is whitespace-separated.
print("Loading MAGMA results...")
df = pd.read_csv("coronary_artery_disease.genes.out", delim_whitespace=True)

# 2. Define the Benjamini-Hochberg FDR function
def benjamini_hochberg(p_values):
    n = len(p_values)
    # Sort p-values and keep original index
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    
    # Calculate rank (1 to n)
    ranks = np.arange(1, n + 1)
    
    # Calculate adjusted p-values: p * (n / rank)
    adjusted_p = sorted_p * n / ranks
    
    # Ensure monotonicity (the next p-value can't be lower than the previous)
    # We take the cumulative minimum from the back
    adjusted_p = np.minimum.accumulate(adjusted_p[::-1])[::-1]
    
    # Cap at 1.0
    adjusted_p[adjusted_p > 1] = 1
    
    # Restore original order
    original_order_q = np.zeros(n)
    original_order_q[sorted_indices] = adjusted_p
    
    return original_order_q

# 3. Calculate Q-values
# Ensure we are using the correct P-value column (usually 'P')
if 'P' in df.columns:
    print("Calculating FDR (Q-values)...")
    df['Q_value'] = benjamini_hochberg(df['P'].values)
else:
    print("Error: 'P' column not found in input file.")
    exit()

# 4. Filter for FDR < 0.05
significant_df = df[df['Q_value'] < 0.05]

# 5. Save results
output_file = "magma_fdr_significant_genes.txt"
significant_df.to_csv(output_file, sep='\t', index=False)

print(f"Analysis complete.")
print(f"Total genes tested: {len(df)}")
print(f"Genes with FDR < 0.05: {len(significant_df)}")
print(f"Results saved to: {output_file}")
