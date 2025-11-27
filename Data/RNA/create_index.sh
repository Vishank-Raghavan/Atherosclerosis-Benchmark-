#!/usr/bin/env bash


# 1. Define your input transcriptome file
# (Make sure you have downloaded and unzipped it, or point to the .gz)
TRANSCRIPTOME="gentrome.fa.gz"

# 2. Build the index
# -t: transcript fasta
# -i: output index folder name
# -k: k-mer length (31 is standard for reads > 75bp)
echo "Building Salmon index..."
salmon index -t "$TRANSCRIPTOME" -d decoys.txt -p 20 -i salmon_index --gencode

echo "Index built in 'human_index/'."
