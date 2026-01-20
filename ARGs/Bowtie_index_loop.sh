#!/bin/bash

# Directories
FASTA_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/ARG_indexes_per_sample"
INDEX_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/bowtie_indexes_per_sample"

# Create output dir if it doesn't exist
mkdir -p "$INDEX_DIR"

# Loop over all FASTA files
for fasta in "$FASTA_DIR"/*.fasta; do
    # Extract base name (e.g., PC239_ARG_index)
    base=$(basename "$fasta" .fasta)

    echo "Building Bowtie2 index for $base..."
    
    # Run bowtie2-build
    bowtie2-build "$fasta" "$INDEX_DIR/$base"
done

echo "All Bowtie2 indices built."
