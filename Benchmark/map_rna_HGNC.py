#!/usr/bin/env python3

import pandas as pd
import sys

# 1. Load the RNA results
input_file = sys.argv[1]
output_file = sys.argv[2]
print(f"Loading {input_file}...")
df = pd.read_csv(input_file)

# The first column (index 0) contains the long ID string.
# Pandas might have named it "Unnamed: 0" or similar if it wasn't named in the CSV.
id_column = df.columns[0]

# 2. Function to extract Gene Symbol
def extract_symbol(long_id):
    try:
        # Split by pipe '|'
        parts = long_id.split('|')
        # The symbol is usually the 6th element (index 5)
        # Check if we have enough parts to avoid errors
        if len(parts) >= 6:
            return parts[5]
        else:
            return long_id # Fallback if format is weird
    except:
        return long_id

# 3. Apply extraction
df['GeneSymbol'] = df[id_column].apply(extract_symbol)

# 4. Clean up
# Select just the GeneSymbol and the statistics you care about
# We'll rename 'padj' to 'FDR' for clarity
clean_df = df[['GeneSymbol', 'log2FoldChange', 'padj']].copy()
clean_df.columns = ['Gene', 'Log2FC', 'FDR']

# remove any rows where Gene extraction failed (empty strings or NAs)
clean_df = clean_df[clean_df['Gene'] != ""]
clean_df = clean_df.dropna(subset=['Gene'])

# 5. Handle Duplicates
# Since DESeq2 analysis was transcript-level, you might have multiple transcripts mapping to one gene.
# For a gene-level benchmark, standard practice is to keep the most significant one (lowest FDR).
clean_df = clean_df.sort_values('FDR').drop_duplicates('Gene', keep='first')

# 6. Save
clean_df.to_csv(output_file, index=False)

print(f"Mapping complete. Found {len(clean_df)} unique significant genes.")
print(f"Saved to {output_file}")
print("Top 5 Genes found:")
print(clean_df.head(5))