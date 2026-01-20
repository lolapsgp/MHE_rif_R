#!/bin/bash

#This script is used to find ARGs in the assembled contigs of Spanish cohort samples

INPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs"
OUTPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Resfinder_results"
DB_DIR="/fast/AG_Forslund/Lola/Sofware"

# Activate environment if needed
source ~/.bashrc
micromamba activate resfinder_py310

# Export required DB paths
export CGE_RESFINDER_RESGENE_DB="$DB_DIR/resfinder_db"
export CGE_RESFINDER_RESPOINT_DB="$DB_DIR/pointfinder_db"
export CGE_DISINFINDER_DB="$DB_DIR/disinfinder_db"
export CGE_BLASTN="/fast/home/l/lginerp/micromamba/envs/resfinder_py310/bin/blastn"

# Loop through each contig file
for fasta in "$INPUT_DIR"/*.fcontigs.fasta; do
    sample=$(basename "$fasta" .fcontigs.fasta)
    outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$outdir"

    echo "Running ResFinder on $sample"
    python -m resfinder \
        -ifa "$fasta" \
        -o "$outdir" \
        -s ecoli \
        --acquired \
        --point
done