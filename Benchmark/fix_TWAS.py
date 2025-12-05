#!/usr/bin/env python3

import pandas as pd
import re
import os

# 1. Load the existing mapped results
input_file = "twas_targets_mapped.csv"

if not os.path.exists(input_file):
    print(f"Error: {input_file} not found. Please ensure the file is in the current directory.")
    exit()

print(f"Loading {input_file}...")
df = pd.read_csv(input_file)

# 2. Define the filter logic for clone-based IDs
def is_clone_id(symbol):
    """
    Returns True if the symbol looks like a clone ID (e.g., AC001234.1, AP000345)
    Standard HGNC symbols usually don't start with "AC/AP/AL" followed by digits
    or contain dots.
    """
    if not isinstance(symbol, str):
        return False
    
    # Pattern 1: Starts with AC, AP, or AL, followed by digits (e.g., AC001234)
    # Pattern 2: Contains a dot (e.g., RP11-123.4), often denoting versions/clones
    if re.match(r'^(AC|AP|AL)\d+', symbol) or '.' in symbol:
        return True
    return False

# 3. Apply Filtering
initial_count = len(df)

# Check if 'GeneSymbol' column exists
if 'GeneSymbol' not in df.columns:
    print("Error: 'GeneSymbol' column not found in the input file.")
    exit()

# Filter: Keep rows where GeneSymbol is NOT a clone ID
df_clean = df[~df['GeneSymbol'].apply(is_clone_id)]

filtered_count = len(df_clean)
removed_count = initial_count - filtered_count

print(f"Original row count: {initial_count}")
print(f"Rows removed (clone IDs): {removed_count}")
print(f"Final row count: {filtered_count}")

# 4. Save the filtered file
output_file = "twas_genes_mapped_filtered.csv"
df_clean.to_csv(output_file, index=False)

print(f"Filtered results saved to {output_file}")

if not df_clean.empty:
    print("Top 5 Mapped Genes (Filtered):")
    print(df_clean.head(5))