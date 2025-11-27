#!/usr/bin/env bash

# 1. Create output directory
mkdir -p salmon_quant

# 2. Define range (matching your previous steps)
START=10687171
END=10687209
INDEX_DIR="human_index"

echo "Starting Salmon quantification..."

for (( i=$START; i<=$END; i++ ))
do
    ACC="SRR${i}"
    
    # Input file (Your cleaned single-end file)
    # Note: Check if your fastp output is .gz (the previous script made .gz)
    READS="${ACC}.clean.fastq.gz"
    
    # Output folder for this specific sample
    OUT_DIR="salmon_quant/${ACC}"

    if [ -f "$READS" ]; then
        echo "Processing ${ACC}..."
        
        # Run Salmon
        # quant: The command to quantify
        # -i: The index folder
        # -l A: Automatically detect library type (strandedness, etc.)
        # -r: Single-end reads (use -1 / -2 for paired)
        # --validateMappings: Improves accuracy
        # --gcBias: Corrects for GC content bias (important!)
        # -o: Output directory
        
        salmon quant \
            -i "$INDEX_DIR" \
            -l A \
            -r "$READS" \
            --validateMappings \
            --gcBias \
            -o "$OUT_DIR" \
            --threads 20

    else
        echo "Warning: Reads for ${ACC} not found."
    fi
done

echo "Quantification complete."
