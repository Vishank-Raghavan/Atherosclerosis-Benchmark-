#!/usr/bin/env python3

import pandas as pd
import os
import glob

# 1. Define setup
base_dir = 'salmon_quant'
output_counts = 'rna_seq_counts_matrix.csv'
output_tpm = 'rna_seq_tpm_matrix.csv'

# Find all quant.sf files
# Structure is salmon_quant/SRR.../quant.sf
files = glob.glob(os.path.join(base_dir, '*', 'quant.sf'))

print(f"Found {len(files)} quantification files.")

counts_list = []
tpm_list = []

for f in files:
    # Extract Sample ID from folder name (e.g., salmon_quant/SRR10687171/quant.sf -> SRR10687171)
    sample_id = os.path.basename(os.path.dirname(f))
    
    # Read the Salmon output file
    # Columns: Name (Transcript), Length, EffectiveLength, TPM, NumReads
    df = pd.read_csv(f, sep='\t', index_col='Name')
    
    # Store NumReads (Counts) and TPM
    counts_list.append(df['NumReads'].rename(sample_id))
    tpm_list.append(df['TPM'].rename(sample_id))

# 2. Concatenate into a Matrix (Rows=Transcripts, Cols=Samples)
print("Aggregating Counts...")
counts_matrix = pd.concat(counts_list, axis=1)

print("Aggregating TPM...")
tpm_matrix = pd.concat(tpm_list, axis=1)



# 3. Save
counts_matrix.to_csv(output_counts)
tpm_matrix.to_csv(output_tpm)

print(f"Saved Counts Matrix to {output_counts}")
print(f"Saved TPM Matrix to {output_tpm}")
print(f"Matrix Dimensions: {counts_matrix.shape}")
