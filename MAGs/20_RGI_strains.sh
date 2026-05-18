#!/bin/bash
# =========================================
# Run RGI on MAGs (protein mode)
# =========================================

# Set input and output directories
INPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/prokka_dreplicated_genomes"
OUTPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/RGI_out/strains"

mkdir -p "$OUTPUT_DIR"

# Activate RGI environment
source ~/.bashrc
micromamba activate rgi  # Replace with your environment name if different

# Load CARD database (local mode)
rgi load --card_json /fast/AG_Forslund/shared/references/CARD/card.json --local

# Loop through all .faa files recursively
find "$INPUT_DIR" -type f -name "*.faa" | while read mag_file; do
    # Extract MAG name from filename (without path and .faa)
    mag_name=$(basename "$mag_file" .faa)

    echo "Running RGI for MAG: $mag_name"

    # Run RGI in protein mode
    rgi main \
      --input_sequence "$mag_file" \
      --output_file "$OUTPUT_DIR/${mag_name}_rgi" \
      --input_type protein \
      --local \
      --clean
done

echo "RGI protein mode processing complete!"
