#!/usr/bin/env python3

import pandas as pd
import mygene
import sys

# 1. Load the TWAS results
input_file = sys.argv[1]
print(f"Loading {input_file}...")
df = pd.read_csv(input_file)

# 2. Extract and Clean Ensembl IDs
# The 'ID' column contains values like 'ENSG00000112137.17'
# We need to remove the version suffix (.17) to get 'ENSG00000112137' for better mapping
def clean_ensembl(gene_id):
    if isinstance(gene_id, str):
        return gene_id.split('.')[0]
    return gene_id

df['Clean_ID'] = df['ID'].apply(clean_ensembl)
ensembl_ids = df['Clean_ID'].unique().tolist()

print(f"Found {len(ensembl_ids)} unique Ensembl IDs to map.")

# 3. Perform Mapping using MyGene.info
mg = mygene.MyGeneInfo()
print("Querying MyGene.info...")

# querymany is efficient for lists
# scopes='ensembl.gene': Tells it we are providing Ensembl Gene IDs
# fields='symbol': Tells it we want the HGNC Symbol back
results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')

# 4. Create a Mapping Dictionary
mapping = {}
for res in results:
    query = res.get('query')
    symbol = res.get('symbol')
    if query and symbol:
        mapping[query] = symbol

# 5. Apply Mapping to DataFrame
df['GeneSymbol'] = df['Clean_ID'].map(mapping)

# 6. Save Results
# We reorder columns to put GeneSymbol first for readability
# We also drop the temporary 'Clean_ID' column
cols = ['GeneSymbol'] + [c for c in df.columns if c not in ['GeneSymbol', 'Clean_ID']]
clean_df = df[cols]

# Optional: Drop rows that failed to map (if you only want mapped genes)
# clean_df = clean_df.dropna(subset=['GeneSymbol'])

output_file = sys.argv[2]
clean_df.to_csv(output_file, index=False)

print(f"Mapping complete. Mapped {clean_df['GeneSymbol'].notna().sum()} genes.")
print(f"Saved to {output_file}")
if not clean_df.empty:
    print("Top 5 Mapped Genes:")
    print(clean_df[['GeneSymbol', 'ID']].head(5))