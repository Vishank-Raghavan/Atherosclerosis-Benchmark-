#!/usr/bin/env python3

import pandas as pd
import mygene
import os
import re

# 1. Load the DIABLO results
input_file = "diablo_significant_targets.csv"
print(f"Loading {input_file}...")
df = pd.read_csv(input_file)

# 2. Function to extract symbol from RNA ID string
# Example: "...|MIR3667HG|..." -> "MIR3667HG"
def extract_rna_symbol(long_id):
    try:
        parts = long_id.split('|')
        if len(parts) >= 6:
            return parts[5] # The symbol is usually at index 5
        return None
    except:
        return None

# 3. Function to clean Protein ID string
# Example: "P02790;A0A..." -> "P02790" (Take the first one)
def clean_protein_id(long_id):
    if pd.isna(long_id): return None
    return str(long_id).split(';')[0]

# 4. Process the DataFrame
mapped_symbols = []

# Collect Protein IDs for batch query
protein_ids_to_query = []
protein_indices = []

for index, row in df.iterrows():
    target_type = row['Type']
    target_id = row['ID']
    
    if target_type == 'RNA':
        # Extract directly from string
        symbol = extract_rna_symbol(target_id)
        mapped_symbols.append(symbol)
        
    elif target_type == 'Protein':
        # Clean ID and prepare for query
        clean_id = clean_protein_id(target_id)
        protein_ids_to_query.append(clean_id)
        protein_indices.append(index)
        mapped_symbols.append(None) # Placeholder, will fill later

# 5. Batch Query Protein IDs using MyGene.info
print(f"Querying MyGene.info for {len(protein_ids_to_query)} protein IDs...")
mg = mygene.MyGeneInfo()
results = mg.querymany(protein_ids_to_query, scopes='uniprot', fields='symbol', species='human')

# Create a mapping dict
prot_mapping = {}
failed_ids = []

for res in results:
    query = res.get('query')
    symbol = res.get('symbol')
    if query and symbol:
        prot_mapping[query] = symbol
    else:
        if query not in prot_mapping:
             failed_ids.append(query)

print(f"Successfully mapped {len(prot_mapping)} IDs via MyGene.info.")
print(f"Failed to map {len(failed_ids)} IDs via MyGene.info.")

# 6. Handle Manual Protein Mappings (Two-Step Process)

# Step A: Save failed IDs to a file for manual lookup
if failed_ids:
    with open("unmapped_protein_ids.txt", "w") as f:
        for uid in failed_ids:
            f.write(f"{uid}\n")
    print("\n[ACTION REQUIRED] Unmapped Protein IDs saved to 'unmapped_protein_ids.txt'.")
    print("If you have a manual mapping file (e.g., from UniProt ID Mapping tool), save it as 'manual_mapping.tsv' or 'manual_mapping.csv'.")

# Step B: Load manual mappings if the file exists
# We check for both TSV and CSV formats
manual_map_file_tsv = "manual_mapping.tsv"
manual_map_file_csv = "manual_mapping.csv"
manual_map_file = None

if os.path.exists(manual_map_file_tsv):
    manual_map_file = manual_map_file_tsv
    sep = '\t'
elif os.path.exists(manual_map_file_csv):
    manual_map_file = manual_map_file_csv
    sep = ','

if manual_map_file:
    print(f"Found {manual_map_file}. Integrating manual mappings...")
    try:
        manual_df = pd.read_csv(manual_map_file, sep=sep)
        
        # Check for expected columns (From/To logic)
        # UniProt TSV usually has 'From' and 'Gene Names'
        if 'From' in manual_df.columns:
            mapped_count = 0
            for index, row in manual_df.iterrows():
                uid = row['From']
                
                # Check for 'Gene Names' (TSV) or 'To' (CSV)
                symbol = None
                if 'Gene Names' in manual_df.columns:
                     gene_names = row['Gene Names']
                     if pd.notna(gene_names) and isinstance(gene_names, str):
                        symbol = gene_names.split(' ')[0] # Take first
                elif 'To' in manual_df.columns:
                    symbol = row['To']
                
                if symbol:
                    prot_mapping[uid] = symbol
                    mapped_count += 1
            print(f"Added {mapped_count} new mappings from manual file.")
        else:
            print(f"Warning: Columns not recognized in {manual_map_file}.")
    except Exception as e:
        print(f"Error reading manual mapping file: {e}")

# 7. Fill in Protein Symbols
for i in protein_indices:
    original_id = df.loc[i, 'ID']
    clean_id = clean_protein_id(original_id)
    
    # Try to find symbol in mapping
    if clean_id in prot_mapping:
        mapped_symbols[i] = prot_mapping[clean_id]
    else:
        # Fallback: keep clean ID if mapping fails
        mapped_symbols[i] = clean_id 

# Update the DataFrame with initial mapped symbols
df['GeneSymbol'] = mapped_symbols

# 8. Filter RNA IDs that look like Ensembl IDs
# Logic: If Type is RNA and GeneSymbol starts with "ENSG", set it to None (to be dropped or handled)
def is_ensembl_id(symbol):
    if isinstance(symbol, str) and symbol.startswith("ENSG"):
        return True
    return False

# Identify rows to drop (RNA rows with ENSG symbols)
rows_to_drop_rna = df[(df['Type'] == 'RNA') & (df['GeneSymbol'].apply(is_ensembl_id))].index

print(f"Filtering out {len(rows_to_drop_rna)} RNA targets that remained as Ensembl IDs (likely clone-based).")
df = df.drop(rows_to_drop_rna)

# 9. Filter Protein IDs that failed to map (still look like UniProt IDs)
# Logic: If Type is Protein and GeneSymbol matches the original cleaned ID (meaning it wasn't replaced by a symbol)
# Or if GeneSymbol matches typical UniProt patterns (e.g. starts with A0, P0, Q9 followed by alphanumeric)
def is_uniprot_id(symbol):
    # Basic check: UniProt IDs are alphanumeric and usually 6-10 chars.
    # HGNC symbols are usually shorter and don't start with numbers, but some do.
    # A robust check is simply: did the symbol change from the input ID?
    # But since we updated the column in place, we check for pattern.
    if isinstance(symbol, str) and re.match(r'^[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]', symbol):
         return True
    # Catch deleted/deprecated IDs like A0A...
    if isinstance(symbol, str) and (symbol.startswith("A0A") or symbol.startswith("deleted")):
        return True
    return False

rows_to_drop_prot = df[(df['Type'] == 'Protein') & (df['GeneSymbol'].apply(is_uniprot_id))].index
print(f"Filtering out {len(rows_to_drop_prot)} Protein targets that failed to map to a Gene Symbol.")
df = df.drop(rows_to_drop_prot)

# 9. Save
output_file = "diablo_targets_mapped.csv"
df.to_csv(output_file, index=False)

print(f"Mapping complete. Mapped {df['GeneSymbol'].notna().sum()} targets.")
print(f"Saved to {output_file}")
print(df.head())