#!/bin/bash
# Run ResFinder on Prevotella / Segatella MAGs

INPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/drep_output/dereplicated_genomes"
OUTPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/Resfinder/strains"
DB_DIR="/fast/AG_Forslund/Lola/Sofware"

# Activate environment
source ~/.bashrc
micromamba activate resfinder_py310

# Export DB paths
export CGE_RESFINDER_RESGENE_DB="$DB_DIR/resfinder_db"
export CGE_RESFINDER_RESPOINT_DB="$DB_DIR/pointfinder_db"
export CGE_DISINFINDER_DB="$DB_DIR/disinfinder_db"
export CGE_BLASTN="/fast/home/l/lginerp/micromamba/envs/resfinder_py310/bin/blastn"

# Loop through MAG FASTA files
for fasta in "$INPUT_DIR"/*.fa; do
    mag=$(basename "$fasta" .fa)
    outdir="$OUTPUT_DIR/$mag"
    mkdir -p "$outdir"

    echo "Running ResFinder on MAG: $mag"

    python -m resfinder \
        -ifa "$fasta" \
        -o "$outdir" \
        --acquired
done

echo "All MAGs processed."
