#!/usr/bin/env python3

import pandas as pd
import mygene
import os
import sys

# 1. Load the Limma results
input_file = sys.argv[1]
print(f"Loading {input_file}...")
df = pd.read_csv(input_file)

# The first column usually contains the IDs. Let's assume it's named "Unnamed: 0" if not specified.
id_col = df.columns[0]

# 2. Clean the IDs
# Take the first UniProt ID if there are multiple (e.g., "P01591;A0A..." -> "P01591")
def clean_uniprot(id_str):
    if pd.isna(id_str): return None
    return str(id_str).split(';')[0]

df['UniProt'] = df[id_col].apply(clean_uniprot)
unique_ids = df['UniProt'].dropna().unique().tolist()

print(f"Found {len(unique_ids)} unique UniProt IDs to map.")

# 3. Perform Mapping using MyGene.info
mg = mygene.MyGeneInfo()
print("Querying MyGene.info...")

# querymany is efficient for lists
results = mg.querymany(unique_ids, scopes='uniprot', fields='symbol', species='human')

# 4. Create a Mapping Dictionary
mapping = {}
failed_ids = []

for res in results:
    query = res.get('query')
    symbol = res.get('symbol')
    if query and symbol:
        mapping[query] = symbol
    else:
        if query not in mapping: # Avoid duplicates
             failed_ids.append(query)

print(f"Successfully mapped {len(mapping)} IDs via MyGene.info.")
print(f"Failed to map {len(failed_ids)} IDs.")

# 5. Handle Manual Mappings (Two-Step Process)

# Step A: Save failed IDs to a file for manual lookup
if failed_ids:
    with open(sys.argv[3], "w") as f:
        for uid in failed_ids:
            f.write(f"{uid}\n")
    print("\n[ACTION REQUIRED] Unmapped IDs saved to 'unmapped_ids.txt'.")
    print("1. Go to https://www.uniprot.org/id-mapping")
    print("2. Upload 'unmapped_ids.txt'")
    print("3. Map From 'UniProtKB_AC-ID' To 'Gene Name'")
    print("4. Download the result as 'manual_mapping.csv' (CSV format) and place it in this folder.")
    print("5. Run this script again to include them.\n")

# 5. Handle Manual Mappings (Using your provided TSV file)
manual_map_file = sys.argv[4]

if os.path.exists(manual_map_file):
    print(f"Found {manual_map_file}. Integrating manual mappings...")
    try:
        # Load the TSV file
        manual_df = pd.read_csv(manual_map_file, sep='\t')
        
        # Check for expected columns
        if 'From' in manual_df.columns and 'Gene Names' in manual_df.columns:
            mapped_count = 0
            for index, row in manual_df.iterrows():
                uid = row['From']
                gene_names = row['Gene Names']
                
                # If we have a gene name, use it
                if pd.notna(gene_names) and isinstance(gene_names, str):
                    # Take the first gene name if multiple are listed (separated by space)
                    symbol = gene_names.split(' ')[0]
                    
                    # Update mapping if not already present or if we want to overwrite
                    # (Usually good to keep MyGene.info as primary, but fill gaps with this)
                    if uid not in mapping:
                        mapping[uid] = symbol
                        mapped_count += 1
            print(f"Added {mapped_count} new mappings from manual file.")
        else:
            print(f"Warning: Columns 'From' and 'Gene Names' not found in {manual_map_file}.")
    except Exception as e:
        print(f"Error reading manual mapping file: {e}")
else:
    print(f"Manual mapping file '{manual_map_file}' not found. Skipping manual step.")

# 6. Apply Mapping to DataFrame
# This now uses both the automated MyGene results AND the loaded manual mappings
df['GeneSymbol'] = df['UniProt'].map(mapping)

# 7. Clean and Save
# We define the columns we want to keep
output_cols = ['GeneSymbol', 'UniProt', 'logFC', 'adj.P.Val', 'P.Value']

# Check if columns exist (limma output names might vary slightly, e.g. "B" stat)
available_cols = [c for c in output_cols if c in df.columns]

# If there are extra columns in the original file you want to keep, add them here
# For example, let's just keep everything from the original df plus GeneSymbol
clean_df = df.copy()

# Move GeneSymbol to the front for readability
cols = ['GeneSymbol'] + [c for c in clean_df.columns if c != 'GeneSymbol']
clean_df = clean_df[cols]

# Remove rows where mapping failed (even after manual attempt)
# If you want to keep unmapped proteins to investigate later, comment out the next line
clean_df = clean_df.dropna(subset=['GeneSymbol'])

# Save
output_file = "protein_genes_mapped.csv"
clean_df.to_csv(output_file, index=False)

print(f"Final mapping complete. Mapped {len(clean_df)} proteins to genes.")
print(f"Saved to {output_file}")
if not clean_df.empty:
    print("Top 5 Mapped Genes:")
    print(clean_df.head(5))