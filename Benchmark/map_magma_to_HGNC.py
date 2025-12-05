#!/usr/bin/env python3

import pandas as pd
import mygene
import sys

# 1. Load the MAGMA results
# The file is tab-separated
input_file = sys.argv[1]
print(f"Loading {input_file}...")
df = pd.read_csv(input_file, sep='\t')

# 2. Extract Entrez IDs
# The 'GENE' column contains the Entrez IDs (integers)
# We convert them to strings for the query
entrez_ids = df['GENE'].astype(str).tolist()

print(f"Found {len(entrez_ids)} unique Entrez IDs to map.")

# 3. Perform Mapping using MyGene.info
mg = mygene.MyGeneInfo()
print("Querying MyGene.info...")

# querymany is efficient for lists
# scopes='entrezgene': Tells it we are providing Entrez IDs
# fields='symbol': Tells it we want the HGNC Symbol back
results = mg.querymany(entrez_ids, scopes='entrezgene', fields='symbol', species='human')

# 4. Create a Mapping Dictionary
mapping = {}
for res in results:
    query = res.get('query')
    symbol = res.get('symbol')
    if query and symbol:
        mapping[int(query)] = symbol # Store key as integer to match dataframe

# 5. Apply Mapping to DataFrame
df['GeneSymbol'] = df['GENE'].map(mapping)

# 6. Save Results
# We reorder columns to put GeneSymbol first for readability
cols = ['GeneSymbol'] + [c for c in df.columns if c != 'GeneSymbol']
df = df[cols]

output_file = sys.argv[2]
df.to_csv(output_file, index=False)

print(f"Mapping complete. Mapped {df['GeneSymbol'].notna().sum()} genes.")
print(f"Saved to {output_file}")
print("Top 5 Mapped Genes:")
print(df.head(5))